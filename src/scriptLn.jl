#=
    scriptLn.jl — the census pipeline, top to bottom.

    One linear pass: calibrate the HMM active-site column, classify every aligned hit,
    dereplicate to the USEARCH centroids, resolve and attach taxonomy, attach habitat,
    then write the census and its manifest, and only then compute statistics.

    Two ordering rules this file exists to keep, both of which a linear script is easy to
    break:

      1. The full census CSV and its TOML manifest are written BEFORE anything is
         filtered. Rows without a TaxID, a kingdom or a habitat belong to the census; if
         they are dropped before the write, no reader can tell what was excluded
         downstream, and the rescued-taxonomy rows can never be audited.
      2. Nothing mutates the census table after that write. Every analysis subset is
         built with `filter`, so the table on disk and the table in memory stay equal.

    Run:  julia --project=. src/scriptLn.jl

    INPUTS ARE NOT IN THIS REPOSITORY. The alignment, the unaligned FASTA, the NCBI
    taxonomy dump and the USEARCH clustering are large third-party or derived files that
    the repository does not carry and this script does not regenerate. All four live under
    ONE data root — the directory named by LNFINDER_DATA, or dat/ in the clone when that is
    not set — in the layout docs/DATA.md documents:

      census/firstTry.fa            the unaligned sequences the hits came from
      census/clusters.uc            their USEARCH `-uc` clustering (see USEARCH_CMD)
      derived/aligned_to_model.afa  those sequences aligned to the profile, in model
                                    coordinates (\$LNFINDER_AFA names it somewhere else)
      derived/taxdump/              an unpacked NCBI taxonomy dump
                                    (\$LNFINDER_TAXDUMP names it somewhere else)

    A run stops on the first one that is not there, naming the input, where the layout
    expects it and which variable moves it. The two that answer to a variable of their own
    are the two a reader makes or fetches rather than unpacks, and they are often large
    enough to be kept outside a project tree.

      LNFINDER_CENSUS_OUT   directory this run writes the census and its manifest into

    Environment switches:
      RUN_HABITAT=false     skip the Omnicrobe lookups and the habitat statistics, for
                            reproduction runs that only need the census table.
      EXPECTED_CENTROIDS=N  abort unless the .uc file yields exactly N centroids. Set it
                            to the recorded centroid count when rebuilding an existing
                            census; left unset, the centroid count is only reported.
=#

#=
    Preamble & environment management
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

#=
    Packages
=#

using LnFinder
using CSV
using DataFrames
using Dates
using SHA
using Serialization
using TOML

#=
    Config
=#

const REPO = normpath(joinpath(@__DIR__, ".."))

# The one registry: every input this driver reads is a NAME there, resolved under the data
# root, and every refusal below is that registry's, naming the input rather than echoing a
# path off anybody's machine.
include(joinpath(REPO, "src", "ModuleDataInputs.jl"))

# A missing input is a sentence, not a sentence followed by several screens of Julia's own
# include machinery: when this file IS the program the refusal is printed and the process
# stops, and when it was included the exception is thrown for the caller to see.
abspath(PROGRAM_FILE) == (@__FILE__) &&
    DataInputs.require_or_exit("model_alignment", "hits_fasta", "taxdump",
                               "uc_clusters", "profile_hmm")

paths = (
    afa     = DataInputs.require_input("model_alignment"),
    fa      = DataInputs.require_input("hits_fasta"),
    taxdump = DataInputs.require_input("taxdump"),
    uc_file = DataInputs.require_input("uc_clusters"),
    hmm     = DataInputs.require_input("profile_hmm")
)

# Clustering identity of the .uc file above. USEARCH_CMD documents the clustering command;
# it is stored in the manifest and never run from here, so its -id argument and DEDUP_ID
# have to be kept in step by hand.
const DEDUP_ID = 0.90
const USEARCH_CMD = "usearch11 -sortbylength <in>.fa -fastaout <sorted>.fa ; " *
                    "usearch11 -cluster_fast <sorted>.fa -id 0.9 -centroids centroids.fa -uc clusters.uc"

const RUN_DATE = Dates.today()

# Where this run writes. Required rather than defaulted: a census written into a directory
# nobody chose is a census that overwrites, or hides beside, an earlier one.
const OUT_ROOT = get(ENV, "LNFINDER_CENSUS_OUT", "")
isempty(OUT_ROOT) && error("LNFINDER_CENSUS_OUT is not set. It must name the directory " *
                           "this run writes the census and its manifest into.")
const OUT_DIR = joinpath(OUT_ROOT, "census_rebuild_$(RUN_DATE)")
const CENSUS_CSV = joinpath(OUT_DIR, "final_deduplicated_data.csv")
const CENSUS_MANIFEST = joinpath(OUT_DIR, "final_deduplicated_data.manifest.toml")

const RUN_HABITAT = lowercase(get(ENV, "RUN_HABITAT", "true")) in ("1", "true", "yes")
const EXPECTED_CENTROIDS = haskey(ENV, "EXPECTED_CENTROIDS") ?
                           parse(Int, ENV["EXPECTED_CENTROIDS"]) : nothing

mkpath(OUT_DIR)
isfile(CENSUS_CSV) && error("$(CENSUS_CSV) already exists. A census is never overwritten: " *
                            "pick a fresh output directory, or move the old run aside.")

#=
    Manifest helpers
=#

"""
    file_sha256(path)

Hex sha256 of `path`, or a marker string when the file is absent, so a missing input
records itself in the manifest instead of aborting the run after the census is written.
"""
function file_sha256(path::AbstractString)
    isfile(path) || return "not found"
    return open(path) do io
        bytes2hex(sha256(io))
    end
end

"""
    direct_package_versions()

Version of every direct dependency of the activated project, for the manifest. Stdlibs
carry no version and are recorded as such.
"""
function direct_package_versions()
    deps = Pkg.dependencies()
    return Dict(name => (haskey(deps, uuid) && deps[uuid].version !== nothing ?
                         string(deps[uuid].version) : "stdlib")
                for (name, uuid) in Pkg.project().dependencies)
end

"""
    toml_ready(x)

Normalise `x` for `TOML.print`: dictionary keys become strings, `missing` and `nothing`
become empty strings, and anything else TOML cannot hold becomes its `string`. One absent
provenance field must not abort the manifest write.
"""
toml_ready(x::AbstractDict) = Dict(string(k) => toml_ready(v) for (k, v) in x)
toml_ready(x::AbstractVector) = [toml_ready(v) for v in x]
toml_ready(::Union{Missing,Nothing}) = ""
toml_ready(x::Union{AbstractString,Real}) = x
toml_ready(x) = string(x)

#=
    Active site calibration and parsing
=#

# Identifies the model index of the D-x-D motif, then classifies every aligned record.
idx = calibrate_active_site(paths.afa)
results_df = parse_hmmalign_results(paths.afa, idx)
rename!(results_df, :original_header => :id)

# `id` is the key against the .uc labels and the only handle on a row from here on.
n_dup_ids = nrow(results_df) - length(unique(results_df.id))
n_dup_ids == 0 || error("$(n_dup_ids) duplicated sequence ids in $(paths.afa); the .uc " *
                        "labels and the census id column have to be 1:1.")

#=
    Dereplication

    Dereplicate with USEARCH v11 cluster_fast (Edgar 2010) at 0.90 identity to reduce
    redundancy from near-identical sequences. This does not remove phylogenetic
    non-independence (handled by phylum standardisation and within-subfamily controls).
=#

if EXPECTED_CENTROIDS === nothing
    @warn "EXPECTED_CENTROIDS is unset, so the centroid count is not being checked. " *
          "Set it to the count of the run being reproduced to turn a silent change in " *
          "the clustering into an abort."
end
df_dedup = filter_centroids(results_df, paths.uc_file; expected = EXPECTED_CENTROIDS)
n_centroids = nrow(df_dedup)

#=
    Taxonomy

    Resolved on the centroids only: every accession dropped by dereplication would cost a
    UniProt round trip and never reach the census.
=#

# The three side files matter on a real run. accession_taxid_cache.csv is the primary
# layer: keyed on the bare accession, it makes a repeat run a local join costing no HTTP at
# all while the UniProt release is unchanged. uniparc_xref_cache.csv keeps the expensive
# per-entry cross-reference round trips resumable after an interruption. And the unresolved
# list is written before the rescue floor throws, so a collapse leaves behind the
# accessions to investigate rather than only a count.
#
# Do NOT seed the accession cache from this census's own TaxID column: a verification run
# that reads back the numbers it is meant to check has verified nothing.
tax = resolve_taxids(df_dedup.accession;
                     taxid_cache    = joinpath(OUT_DIR, "accession_taxid_cache.csv"),
                     xref_cache     = joinpath(OUT_DIR, "uniparc_xref_cache.csv"),
                     unresolved_csv = joinpath(OUT_DIR, "unresolved_accessions.csv"))
allunique(tax.accession) || error("resolve_taxids returned an accession more than once; " *
                                  "the join below would duplicate census rows.")

df_dedup = leftjoin(df_dedup, tax, on = :accession)
nrow(df_dedup) == n_centroids || error("the taxonomy join changed the row count " *
                                       "($(n_centroids) -> $(nrow(df_dedup))); the census " *
                                       "keeps exactly one row per centroid.")

# Lineages from the local NCBI dump. append_taxonomy! aborts when the dump's rank scheme
# leaves most rows that have a TaxID without a kingdom.
db = load_taxonomy_db(paths.taxdump)
append_taxonomy!(df_dedup, db; taxid_col = :TaxID, min_kingdom_coverage = 0.9)

#=
    Habitat
=#

if RUN_HABITAT
    fetch_omnicrobe_env!(df_dedup; taxid_col = :TaxID)
    # Heuristic rollup into ecological pillars. Taxa Omnicrobe never mapped keep an
    # explicit label instead of a missing, so the column stays a plain Vector{String}.
    df_dedup.report_group = [ismissing(h) ? "_Not Mapped_" : agg_habitat(h)
                             for h in df_dedup.habitat_specific]
else
    @info "RUN_HABITAT is off: skipping the Omnicrobe lookups and the habitat statistics."
    df_dedup.report_group = fill("_Not Mapped_", nrow(df_dedup))
end

#=
    Motif columns over the whole census
=#

# Second motif position, the residue that separates the classes.
df_dedup.pos2_aa = [length(s) >= 2 ? string(s[2]) : "-" for s in df_dedup.active_site]

# Residue frequencies, computed before any subsetting.
stats_df = count_residues_per_position(df_dedup.active_site)

#=
    CENSUS OUTPUT — written before anything is filtered

    Every centroid is written here, rows without a TaxID, kingdom or habitat included.
    The manifest beside it carries the provenance that tells two runs apart. Everything
    below this point only reads df_dedup.
=#

CSV.write(CENSUS_CSV, df_dedup)

uniprot_releases = unique(string.(skipmissing(tax.uniprot_release)))
manifest = Dict(
    "run" => Dict(
        "date" => string(RUN_DATE),
        "script" => "src/scriptLn.jl",
        "julia_version" => string(VERSION),
        "run_habitat" => RUN_HABITAT,
        "census_rows" => nrow(df_dedup)
    ),
    "inputs" => Dict(
        "afa" => paths.afa,
        "afa_sha256" => file_sha256(paths.afa),
        "fa" => paths.fa,
        "fa_sha256" => file_sha256(paths.fa),
        "uc" => paths.uc_file,
        "uc_sha256" => file_sha256(paths.uc_file),
        "hmm" => paths.hmm,
        "hmm_sha256" => file_sha256(paths.hmm),
        "active_site_model_index" => idx
    ),
    "dereplication" => Dict(
        "identity" => DEDUP_ID,
        "usearch_cmd" => USEARCH_CMD,
        "expected_centroid_count" => EXPECTED_CENTROIDS,
        "census_rows_after_dedup" => n_centroids
    ),
    "taxonomy" => Dict(
        "uniprot_release" => join(uniprot_releases, ";"),
        "taxdump_dir" => paths.taxdump,
        "taxdump" => taxdump_provenance(paths.taxdump)
    ),
    "packages" => direct_package_versions()
)

open(CENSUS_MANIFEST, "w") do io
    TOML.print(io, toml_ready(manifest))
end

println("Census written:   $(CENSUS_CSV) ($(nrow(df_dedup)) rows)")
println("Manifest written: $(CENSUS_MANIFEST)")

#=
    Analysis subsets

    `filter` returns a copy, so the census written above is never modified. Rows without a
    kingdom are excluded here, and only here.
=#

ln_users = filter(r -> r.classification == :Lanthanide && !ismissing(r.kingdom), df_dedup)

# FASTA for the downstream tree, headers enriched from the full census.
export_tree_fasta(df_dedup, paths.fa, joinpath(OUT_DIR, "tree_input_deduped.fasta"))

#=
    Statistics: domain and ecology against the second motif residue
=#

arch_vs_bac_df, p_val_arch = compare_archaea_bacteria(ln_users)

CSV.write(joinpath(OUT_DIR, "archaea_comparison.csv"), arch_vs_bac_df)
CSV.write(joinpath(OUT_DIR, "motif_conservation_stats.csv"), stats_df)

# Serialize the raw statistical objects for later visualization in Julia/R.
if RUN_HABITAT
    matrix_df, obs, z_scores, pval = calculate_enrichment_stats(ln_users)

    # Convert the contingency table to a z-score table with the same layout.
    z_df = copy(matrix_df)
    for (i, colname) in enumerate(names(z_df)[2:end])
        z_df[!, colname] = round.(z_scores[:, i], digits = 2)
    end

    CSV.write(joinpath(OUT_DIR, "niche_enrichment.csv"), matrix_df)
    CSV.write(joinpath(OUT_DIR, "niche_z_scores.csv"), z_df)

    serialize(joinpath(OUT_DIR, "statistical_results.jls"),
        Dict(:p_value_arch => p_val_arch, :p_value_niche => pval,
             :obs => obs, :total => sum(obs)))
else
    serialize(joinpath(OUT_DIR, "statistical_results.jls"),
        Dict(:p_value_arch => p_val_arch))
end
