module TaxonomyTools

export load_taxonomy_db, get_lineage, append_taxonomy!, merge_taxonomy_with_hits, filter_centroids
export taxdump_provenance

using DataFrames
using ProgressMeter
using SHA

#=
    NCBI reorganised the top of the taxonomy in March 2025: "superkingdom" became "domain",
    a real "kingdom" rank was inserted below it (Pseudomonadati, Metazoa, …) and viruses
    were moved under the "acellular root". Code that keys on "superkingdom" alone silently
    returns an empty kingdom for every row of a post-March-2025 dump, and code that keys on
    rank == "kingdom" picks up Pseudomonadati instead of Bacteria. Both schemes are handled
    here, and append_taxonomy! refuses to return a census whose kingdom column is mostly
    empty rather than letting the gap travel downstream.

    Dumps before March 2025 use rank "superkingdom" for Bacteria/Archaea/Eukaryota; later
    dumps use "domain" (Cox et al. 2025, doi:10.1093/nar/gkae967; NCBI Insights 2025-02-27).
=#

# Ranks copied into the lineage. "superkingdom" (pre-2025) and "domain" (post-2025) both
# land in :kingdom; NCBI's own "kingdom" rank is deliberately absent, see KINGDOM_RANKS.
const TARGET_RANKS = ["superkingdom", "domain", "phylum", "class", "order", "family", "genus", "species"]

# The two spellings of the top-level cellular division, whichever dump vintage is loaded.
const KINGDOM_RANKS = ("superkingdom", "domain")

# Post-2025 dumps hang viruses off an "acellular root" node instead of giving them a
# superkingdom. Recording it as :kingdom keeps virus rows comparable across dump vintages.
const ACELLULAR_ROOT_RANK = "acellular root"
const VIRUSES_NAME = "Viruses"

"""
    TaxonomyDB
    Holds the parent-child relationships, names, merged ids and deleted ids in memory.
"""
struct TaxonomyDB
    parents::Dict{Int, Int}
    ranks::Dict{Int, String}
    names::Dict{Int, String}
    merged::Dict{Int, Int}
    deleted::Set{Int}
end

"""
    TaxonomyDB(parents, ranks, names)
    Back-compat constructor for callers built before merged.dmp/delnodes.dmp were parsed:
    no id is then known to be merged or deleted.
"""
TaxonomyDB(parents::Dict{Int, Int}, ranks::Dict{Int, String}, names::Dict{Int, String}) =
    TaxonomyDB(parents, ranks, names, Dict{Int, Int}(), Set{Int}())

"""
    load_taxonomy_db(dump_path::String)
    Parses NCBI nodes.dmp, names.dmp and, when present, merged.dmp and delnodes.dmp.
"""
function load_taxonomy_db(dump_path::String)
    nodes_file = joinpath(dump_path, "nodes.dmp")
    names_file = joinpath(dump_path, "names.dmp")
    merged_file = joinpath(dump_path, "merged.dmp")
    delnodes_file = joinpath(dump_path, "delnodes.dmp")

    if !isfile(nodes_file) || !isfile(names_file)
        error("NCBI dump files not found in $dump_path. Ensure nodes.dmp and names.dmp are there.")
    end

    parents = Dict{Int, Int}()
    ranks = Dict{Int, String}()
    names_map = Dict{Int, String}()
    merged = Dict{Int, Int}()
    deleted = Set{Int}()

    println("📂 Loading NCBI Taxonomy Nodes...")
    # nodes.dmp format: tax_id | parent_id | rank | ...
    open(nodes_file, "r") do io
        for line in eachline(io)
            isempty(strip(line)) && continue
            parts = split(line, "|")
            length(parts) < 3 && continue
            tid = parse(Int, strip(parts[1]))
            pid = parse(Int, strip(parts[2]))
            rnk = strip(parts[3])
            parents[tid] = pid
            ranks[tid] = rnk
        end
    end

    println("📂 Loading NCBI Taxonomy Names...")
    # names.dmp format: tax_id | name_txt | unique_name | name_class |
    open(names_file, "r") do io
        for line in eachline(io)
            isempty(strip(line)) && continue
            parts = split(line, "|")
            length(parts) < 4 && continue
            # Test the name CLASS field, not the whole line: synonyms such as
            # "no scientific name known" contain the phrase in their text and, matched by
            # occursin, overwrite the scientific name of whatever id they belong to.
            strip(parts[4]) == "scientific name" || continue
            tid = parse(Int, strip(parts[1]))
            nam = strip(parts[2])
            names_map[tid] = nam
        end
    end

    # merged.dmp format: old_tax_id | new_tax_id |
    if isfile(merged_file)
        open(merged_file, "r") do io
            for line in eachline(io)
                isempty(strip(line)) && continue
                parts = split(line, "|")
                length(parts) < 2 && continue
                merged[parse(Int, strip(parts[1]))] = parse(Int, strip(parts[2]))
            end
        end
    else
        @warn "merged.dmp not found in $dump_path; retired TaxIDs cannot be followed to their current id."
    end

    # delnodes.dmp format: tax_id |
    if isfile(delnodes_file)
        open(delnodes_file, "r") do io
            for line in eachline(io)
                isempty(strip(line)) && continue
                parts = split(line, "|")
                push!(deleted, parse(Int, strip(parts[1])))
            end
        end
    else
        @warn "delnodes.dmp not found in $dump_path; deleted TaxIDs cannot be told apart from unknown ones."
    end

    println("✅ Taxonomy DB loaded ($(length(parents)) nodes, $(length(merged)) merged, $(length(deleted)) deleted).")
    return TaxonomyDB(parents, ranks, names_map, merged, deleted)
end

"""
    resolve_taxid(taxid::Int, db::TaxonomyDB)
    Follows merged.dmp to the id that is current in this dump.
    Returns (resolved_id, status, merged_to) with status ∈ {"ok", "merged", "deleted"};
    merged_to is `missing` unless the id was actually redirected.
"""
function resolve_taxid(taxid::Int, db::TaxonomyDB)
    taxid in db.deleted && return (taxid, "deleted", missing)

    current_id = taxid
    merged_to = missing
    # merged.dmp is normally one hop, but a chain costs nothing to follow; the cap is there
    # so a cyclic dump cannot hang the census.
    hops = 0
    while haskey(db.merged, current_id) && hops < 10
        current_id = db.merged[current_id]
        merged_to = current_id
        hops += 1
    end

    ismissing(merged_to) && return (current_id, "ok", missing)
    current_id in db.deleted && return (current_id, "deleted", merged_to)
    return (current_id, "merged", merged_to)
end

"""
    lineage_of(taxid::Int, db::TaxonomyDB)
    Walks parent links from an already-resolved id to the root, collecting TARGET_RANKS.
    Returns `missing` when no target rank is met.
"""
function lineage_of(taxid::Int, db::TaxonomyDB)
    lineage = Dict{Symbol, String}()

    current_id = taxid
    # Safety counter to prevent infinite loops on corrupted dmp files
    max_depth = 50

    while current_id != 1 && max_depth > 0
        rank = get(db.ranks, current_id, "no rank")
        name = get(db.names, current_id, "Unknown")

        if rank in TARGET_RANKS
            # Both dump vintages collapse onto :kingdom. NCBI's own "kingdom" rank is not
            # in TARGET_RANKS, so Pseudomonadati/Metazoa can never land there.
            key = rank in KINGDOM_RANKS ? :kingdom : Symbol(rank)
            lineage[key] = name
        elseif rank == ACELLULAR_ROOT_RANK && name == VIRUSES_NAME
            lineage[:kingdom] = VIRUSES_NAME
        end

        parent_id = get(db.parents, current_id, 1)
        if parent_id == current_id; break; end # Reached root
        current_id = parent_id
        max_depth -= 1
    end

    return isempty(lineage) ? missing : lineage
end

"""
    lineage_with_status(taxid::Int, db::TaxonomyDB)
    get_lineage plus the provenance the census needs per row:
    (lineage, status, merged_to) with status ∈ {"ok", "merged", "deleted", "no_domain"}.
    "no_domain" flags a live id whose walk produced no :kingdom — the signature of a dump
    whose rank scheme this code does not handle. A merged id keeps the "merged" status so
    the redirect stays visible.
"""
function lineage_with_status(taxid::Int, db::TaxonomyDB)
    resolved_id, status, merged_to = resolve_taxid(taxid, db)
    status == "deleted" && return (missing, "deleted", merged_to)

    lineage = lineage_of(resolved_id, db)
    if status == "ok" && (ismissing(lineage) || !haskey(lineage, :kingdom))
        status = "no_domain"
    end
    return (lineage, status, merged_to)
end

"""
    get_lineage(taxid::Int, db::TaxonomyDB)
    Reconstructs the full lineage (Kingdom to Species) for a given TaxID, following
    merged.dmp first. Deleted ids return `missing`.
"""
function get_lineage(taxid::Int, db::TaxonomyDB)
    lineage, _, _ = lineage_with_status(taxid, db)
    return lineage
end

"""
    append_taxonomy!(df::DataFrame, db::TaxonomyDB; taxid_col=:TaxID, min_kingdom_coverage=0.9)
    Modifies the DataFrame in-place to add taxonomy columns, plus :taxid_lineage_status and
    :taxid_merged_to. Errors when fewer than `min_kingdom_coverage` of the rows that have a
    TaxID come back with a kingdom.
"""
function append_taxonomy!(df::DataFrame, db::TaxonomyDB; taxid_col=:TaxID,
                          min_kingdom_coverage::Float64=0.9)
    println("🧬 Appending Taxonomic Lineages...")

    # Initialize columns as Union{String, Missing}
    ranks = [:kingdom, :phylum, :class, :order, :family, :genus, :species]
    for r in ranks
        df[!, r] = Vector{Union{String, Missing}}(missing, nrow(df))
    end
    df[!, :taxid_lineage_status] = Vector{String}(undef, nrow(df))
    df[!, :taxid_merged_to] = Vector{Union{Int, Missing}}(missing, nrow(df))

    n_with_taxid = 0
    n_with_kingdom = 0

    # Progress bar for long datasets
    @showprogress for i in 1:nrow(df)
        tid = df[i, taxid_col]

        # skipmissing handles our 'missing' TaxIDs from earlier steps
        if ismissing(tid)
            df[i, :taxid_lineage_status] = "no_taxid"
        else
            n_with_taxid += 1
            lin, status, merged_to = lineage_with_status(Int(tid), db)
            df[i, :taxid_lineage_status] = status
            df[i, :taxid_merged_to] = merged_to

            if !ismissing(lin)
                for r in ranks
                    if haskey(lin, r)
                        df[i, r] = lin[r]
                    end
                end
                haskey(lin, :kingdom) && (n_with_kingdom += 1)
            end
        end
    end

    # Fail loudly here rather than let a kingdom-less census reach the figures: a
    # kingdom-less census is exactly what the March 2025 rank change produces, and nothing
    # downstream notices it.
    if n_with_taxid == 0
        @warn "No row carries a TaxID; kingdom coverage cannot be checked."
    else
        coverage = n_with_kingdom / n_with_taxid
        if coverage < min_kingdom_coverage
            error("Kingdom coverage $(round(100 * coverage, digits=1))% " *
                  "($n_with_kingdom/$n_with_taxid rows with a TaxID) is below the required " *
                  "$(round(100 * min_kingdom_coverage, digits=1))%. Likely cause: a rank-scheme " *
                  "mismatch (the dump uses neither 'superkingdom' nor 'domain' at the top) or a " *
                  "taxdump that does not match the TaxIDs in this table (check merged.dmp and " *
                  "delnodes.dmp, and the :taxid_lineage_status column for 'no_domain'/'deleted').")
        end
        println("✅ Kingdom coverage: $n_with_kingdom / $n_with_taxid rows with a TaxID.")
    end

    return df
end

"""
    taxdump_provenance(dir::AbstractString)
    sha256 and mtime of the taxdump files a census run depends on, for the run manifest.
    Absent files are recorded as present=false rather than skipped, so a manifest shows
    which dumps the run did not have.
"""
function taxdump_provenance(dir::AbstractString)
    provenance = Dict{String, Any}()
    for fname in ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp")
        path = joinpath(dir, fname)
        if isfile(path)
            digest = open(sha256, path)
            provenance[fname] = Dict{String, Any}(
                "present" => true,
                "sha256" => bytes2hex(digest),
                "mtime" => mtime(path),
                "bytes" => filesize(path),
            )
        else
            provenance[fname] = Dict{String, Any}("present" => false)
        end
    end
    return provenance
end

"""
    merge_taxonomy_with_hits(hits_df, tax_df; on_col=:TaxID)

Joins the HMM results with the Omnicrobe/NCBI Taxonomy metadata.
We use a left join to ensure no protein hits are lost during the merge.
"""
function merge_taxonomy_with_hits(hits_df::DataFrame, tax_df::DataFrame; on_col=:TaxID)
    # 1. Ensure the join column is the same type (Int64)
    # This prevents 'No matching method' errors during join
    hits_df[!, on_col] = Int64.(hits_df[!, on_col])
    tax_df[!, on_col] = Int64.(tax_df[!, on_col])

    println("🔗 Merging $(nrow(hits_df)) hits with taxonomic metadata...")

    # 2. Perform Left Join
    # hits_df is the 'Left' table (the primary data)
    # tax_df is the 'Right' table (the metadata)
    merged_df = leftjoin(hits_df, tax_df, on=on_col)

    # 3. Post-merge Cleanup
    # The habitat summary below indexes :habitat_specific, so a metadata frame without it
    # would throw two lines after the warning that exists to soften exactly that case.
    if !all(c -> c in names(merged_df), ("habitat_category", "habitat_specific"))
        @warn "Metadata columns not found in merged result. Check tax_df headers."
        return merged_df
    end

    total_mapped = count(.!ismissing.(merged_df.habitat_specific))
    println("✅ Merge complete. $total_mapped / $(nrow(merged_df)) hits now have habitat data.")

    return merged_df
end

"""
    filter_centroids(df, uc_path; expected=nothing)
Filters the main DataFrame to keep only the cluster representatives from USEARCH.
Errors if the .uc holds no seed records, or if `expected` is given and the seed count
differs from it.
"""
function filter_centroids(df, uc_path; expected::Union{Int, Nothing}=nothing)
    # A Set, not a Vector: `row.id in centroids` over a Vector is a linear scan per row, so
    # a full census pays a quadratic number of string comparisons for the lookup.
    centroids = Set{String}()
    for line in readlines(uc_path)
        parts = split(line, '\t')
        length(parts) < 9 && continue
        if parts[1] == "S"  # S indicates a Centroid (Seed)
            # Column 9 is the sequence label in UC format
            push!(centroids, String(parts[9]))
        end
    end

    # An empty set silently filters the census down to zero rows, which is how a wrong or
    # truncated .uc path would otherwise pass unnoticed.
    if isempty(centroids)
        error("No 'S' (seed) records in $uc_path. Either the file is not USEARCH .uc output, " *
              "or the clustering run wrote no centroids.")
    end
    if expected !== nothing && length(centroids) != expected
        error("$uc_path holds $(length(centroids)) centroids, expected $expected. " *
              "The clustering run does not match the one this census was built from.")
    end

    # Filter the original dataframe
    dedup_df = filter(row -> row.id in centroids, df)

    println("✂️ Deduplication Summary:")
    println("  Total sequences: $(nrow(df))")
    println("  Cluster Centroids: $(nrow(dedup_df))")
    println("  Reduction: $(round(100 * (1 - nrow(dedup_df)/nrow(df)), digits=1))%")

    return dedup_df
end

end # module
