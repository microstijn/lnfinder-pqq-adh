#=
    prokaryote_background.jl — the denominator: where do prokaryotes live IN GENERAL?

    "Where do lanthanide binders live?" cannot be answered by counting their habitats.
    OmniCrobe is text-mined from publications, and a model organism accumulates habitat
    relations by orders of magnitude faster than an unstudied one, so a raw habitat tally
    measures research effort, not ecology. The distribution only becomes interpretable
    against a background: if 35% of all studied prokaryotes are aquatic, then 35% of Ln
    hits being aquatic is not a finding.

    This builds that background by running a generic sample of prokaryotes through the
    IDENTICAL pipeline the census used (fetch_relations → aggregate_habitats), which is
    what makes the two directly comparable.

    TAXON POOL — UniProt reference proteomes for Bacteria (taxonomy_id:2) and Archaea
    (2157). Chosen deliberately: the census hits are themselves UniProt proteins from
    sequenced organisms, so the pool shares the census's strongest bias (only cultured,
    sequenced, named organisms are represented) and it therefore cancels in the contrast.
    An "all of NCBI taxonomy" pool would NOT share it and would manufacture enrichment.

    The `lineage` field is fetched alongside so a PHYLUM-MATCHED background can be built
    downstream. That control matters: census taxa are those carrying a PQQ-ADH, so an
    apparent habitat enrichment could just be the protein family's phylogenetic range.

    ⚠ UniProt gzips responses even when identity encoding is requested — see
    `ValidateClassifier._http_fetch` in src/validate_classifier.jl. Without the defensive
    gunzip every parse fails silently and the whole background comes back empty.

    Resumable: taxa already present in the output CSV are skipped, so an interrupted run
    costs nothing. The OBT ontology cache is shared with the census run.

    Run:  julia --project=. scripts/prokaryote_background.jl   (it takes no options)
=#

# ---------------------------------------------------------------------------
# Command line. This step takes no options, but it still reads the command line, before
# any package is loaded and before any file or URL is touched. Every other step in the
# pipeline answers --help, so a reader probes this one the same way; without this block
# `--help` falls through to the body and begins a multi-thousand-request UniProt fetch and
# a write into the reader's working tree. An unrecognised argument is refused rather than
# ignored, for the same reason: a mistyped flag must not silently become a full run.
# ---------------------------------------------------------------------------
const USAGE = """
Usage: julia --project=. scripts/prokaryote_background.jl [-h|--help]

Habitat step 0 — the generic-prokaryote background the later steps are contrasted
against. Takes no options.

Reads:  UniProt reference proteomes over the network, then the habitat relations of each
        sampled taxon (network, cached on disk).
Writes, under the data root (\$LNFINDER_DATA, else dat/ in the clone; docs/DATA.md):
        habitat/prokaryote_taxa.csv        the sampled taxon pool
        habitat/prokaryote_background.csv  one habitat row per taxon
        habitat/obt_path_cache.csv         the shared ontology-path cache

Resumable: a rerun keeps the cached pool and skips taxa already present in the output, so
an interrupted run costs nothing. Delete the two CSVs to start over.
"""
if abspath(PROGRAM_FILE) == (@__FILE__) && !isempty(ARGS)
    if any(a -> a in ("-h", "--help"), ARGS)
        print(USAGE)
        exit(0)
    end
    println(stderr, "prokaryote_background.jl: this step takes no options; got ",
            join(ARGS, " "), "\n")
    print(stderr, USAGE)
    exit(1)
end

const REPO = normpath(joinpath(@__DIR__, ".."))
using Pkg; Pkg.activate(REPO)
using DataFrames, CSV, HTTP, JSON3, Printf, Random, CodecZlib, ProgressMeter

include(joinpath(REPO, "src", "ModuleEnvironmentTools.jl")); using .EnvironmentTools

# Names in the one registry, resolved under the data root (src/ModuleDataInputs.jl). This
# step WRITES all three, and step 1 reads them back from exactly these names.
include(joinpath(REPO, "src", "ModuleDataInputs.jl"))

const OUT           = DataInputs.input_path("pool_background")
const TAXA_CACHE    = DataInputs.input_path("pool_taxa")
const OBT_CACHE_CSV = DataInputs.input_path("obt_path_cache")

# This step is the first thing that writes under the data root, so it is the step that has
# to be able to create the directory the layout puts its outputs in. Done here rather than
# at load, so that loading this file still writes nothing.
make_output_dir() = mkpath(dirname(OUT))
const MAX_TAXA      = 8000      # take everything below this; subsample above it
const SUBSAMPLE     = 6000
const SEED          = 20260810  # recorded so the sample is reproducible

# ---- OBT path cache, shared with every habitat run in this repository ------------
# The loader and the writer live in ModuleEnvironmentTools so that every habitat run shares
# one on-disk cache and one ordering; local copies drift, and a copy that writes the CSV in
# Dict order reintroduces the non-determinism the habitat work exists to remove.
prime_cache!() = EnvironmentTools.load_obt_cache!(OBT_CACHE_CSV)
save_cache!() = EnvironmentTools.save_obt_cache(OBT_CACHE_CSV)

# ---- HTTP with the UniProt gunzip guard ------------------------------------------
function fetch_text(url; timeout = 180)
    try
        r = HTTP.get(url, ["Accept-Encoding" => "identity"];
                     status_exception = false, readtimeout = timeout)
        r.status == 200 || (@warn("HTTP $(r.status)", url); return "")
        b = r.body
        (length(b) >= 2 && b[1] == 0x1f && b[2] == 0x8b) &&
            return String(transcode(GzipDecompressor, b))
        return String(b)
    catch e
        @warn "fetch failed" url exception = e
        return ""
    end
end

# ---- taxon pool -------------------------------------------------------------------
"""Reference-proteome taxa for a domain: (taxid, organism, phylum)."""
function reference_taxa(taxid_root::Int, label::AbstractString)
    url = "https://rest.uniprot.org/proteomes/stream?query=" *
          HTTP.escapeuri("(taxonomy_id:$taxid_root) AND (reference:true)") *
          "&format=tsv&fields=organism_id,organism,lineage"
    body = fetch_text(url)
    rows = NamedTuple[]
    for (i, line) in enumerate(split(body, '\n'))
        i == 1 && continue                      # header
        isempty(strip(line)) && continue
        p = split(line, '\t')
        length(p) >= 2 || continue
        tid = tryparse(Int, strip(p[1])); tid === nothing && continue
        lineage = length(p) >= 3 ? String(p[3]) : ""
        # UniProt lineage is comma-separated, ranks unlabelled; phylum is the token
        # ending in -ota / -etes / -bacteria under the domain. Take the first plausible.
        phy = ""
        for t in strip.(split(lineage, ','))
            (endswith(t, "ota") || endswith(t, "etes") || endswith(t, "bacteria") ||
             endswith(t, "archaeota")) || continue
            t in ("Bacteria", "Archaea") && continue
            phy = String(t); break
        end
        push!(rows, (taxid = tid, organism = String(strip(p[2])), phylum = phy, domain = label))
    end
    @printf("  %-8s %5d reference-proteome taxa\n", label, length(rows))
    return rows
end

if isfile(TAXA_CACHE)
    taxa = CSV.read(TAXA_CACHE, DataFrame)
    @printf("Loaded cached taxon pool: %d taxa\n", nrow(taxa))
else
    println("Fetching UniProt reference-proteome taxa …")
    pool = vcat(reference_taxa(2, "Bacteria"), reference_taxa(2157, "Archaea"))
    isempty(pool) && error("no taxa retrieved — UniProt unreachable, or the gunzip guard is missing")
    taxa = DataFrame(pool)
    unique!(taxa, :taxid)
    if nrow(taxa) > MAX_TAXA
        @printf("Pool of %d exceeds %d — subsampling %d (seed %d)\n",
                nrow(taxa), MAX_TAXA, SUBSAMPLE, SEED)
        taxa = taxa[shuffle(MersenneTwister(SEED), 1:nrow(taxa))[1:SUBSAMPLE], :]
    end
    make_output_dir()
    CSV.write(TAXA_CACHE, taxa)
    @printf("Taxon pool: %d taxa → %s\n", nrow(taxa), TAXA_CACHE)
end

# ---- resume: skip taxa already annotated -----------------------------------------
done = Set{Int}()
existing = DataFrame()
if isfile(OUT)
    existing = CSV.read(OUT, DataFrame)
    done = Set(Int.(existing.taxid))
    @printf("Resuming: %d taxa already annotated\n", length(done))
end
todo = taxa[.!in.(taxa.taxid, Ref(done)), :]
@printf("To fetch: %d taxa (~%.0f min at 0.3 s each)\n", nrow(todo), 0.3nrow(todo) / 60)

n = prime_cache!(); n > 0 && @printf("Primed OBT path cache with %d entries\n", n)

# ---- annotate ---------------------------------------------------------------------
rows = NamedTuple[]
@showprogress desc="Habitat-mapping background taxa " for r in eachrow(todo)
    # `nothing` means the fetch failed, which is not the same as "no habitat"; this
    # script tolerates a transient failure and records the taxon as unannotated.
    rels = something(EnvironmentTools.fetch_relations(r.taxid), EnvironmentTools.Relation[])
    agg = aggregate_habitats(rels, EnvironmentTools.get_obt_path)
    push!(rows, (taxid = r.taxid, organism = r.organism, phylum = r.phylum, domain = r.domain,
                 habitat_type = ismissing(agg.habitat_type) ? "" : String(agg.habitat_type),
                 habitat_category = ismissing(agg.category) ? "" : String(agg.category),
                 habitat_chemistry = agg.chemistry,
                 habitat_specific = ismissing(agg.name) ? "" : String(agg.name),
                 confidence = round(agg.confidence, digits = 3),
                 n_relations = agg.n_relations, n_distinct = agg.n_distinct))
end

out = isempty(rows) ? existing : (isempty(existing) ? DataFrame(rows) : vcat(existing, DataFrame(rows); cols = :union))
make_output_dir()
CSV.write(OUT, out)
@printf("Saved OBT path cache (%d entries)\n", save_cache!())

# ---- summary ----------------------------------------------------------------------
mapped = out[.!ismissing.(out.habitat_type) .& (string.(coalesce.(out.habitat_type, "")) .!= ""), :]
@printf("\n%s\n", "="^72)
@printf("BACKGROUND: %d taxa queried, %d with a habitat (%.1f%%)\n",
        nrow(out), nrow(mapped), 100nrow(mapped) / max(nrow(out), 1))
@printf("%s\n", "="^72)

println("\n=== Habitat type distribution (background) ===")
tt = sort(collect(Dict(String(v) => count(==(v), mapped.habitat_type) for v in unique(mapped.habitat_type))),
          by = kv -> -kv[2])
for (k, v) in tt
    @printf("  %-18s %5d  (%5.1f%%)\n", k, v, 100v / nrow(mapped))
end

println("\n=== Habitat chemistry (background) ===")
ch = Dict{String,Int}()
for c in mapped.habitat_chemistry
    s = ismissing(c) ? "" : String(c); isempty(s) && continue
    for f in split(s, '+'); ch[String(f)] = get(ch, String(f), 0) + 1; end
end
for (k, v) in sort(collect(ch), by = kv -> -kv[2])
    @printf("  %-18s %5d  (%5.1f%%)\n", k, v, 100v / nrow(mapped))
end

println("\n=== Domain / top phyla ===")
for d in unique(mapped.domain)
    @printf("  %-10s %5d\n", d, count(==(d), mapped.domain))
end
ph = sort(collect(Dict(String(v) => count(==(v), mapped.phylum) for v in unique(mapped.phylum) if !ismissing(v) && String(v) != "")),
          by = kv -> -kv[2])
println("  top phyla: ", join(["$k=$v" for (k, v) in first(ph, 8)], "  "))

println("\n✅ wrote $OUT")
