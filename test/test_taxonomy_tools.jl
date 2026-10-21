#=
    Offline self-test for the NCBI taxonomy layer. No network and no real dump: both
    fixtures are written into an mktempdir() in .dmp syntax, with a handful of synthetic
    nodes.

    Three failures this suite exists to prevent:

    1. The March-2025 rank change. NCBI renamed "superkingdom" to "domain", inserted a real
       "kingdom" rank (Pseudomonadati, Metazoa, …) below it and moved viruses under an
       "acellular root". A lineage walk matching "superkingdom" only brings every row of a
       current dump back without a kingdom — and append_taxonomy! reports nothing, having
       no idea what coverage to expect. The obvious repair (accept rank "kingdom") is
       worse: it fills the column with Pseudomonadati instead of Bacteria.

    2. Retired TaxIDs. Ids that NCBI merged into another node, or deleted outright, looked
       up as-is: the merged ones fall out of the lineage walk silently, and the deleted
       ones are indistinguishable from ids that simply carry no rank.

    3. names.dmp matched with occursin("scientific name", line), which tests the whole line
       rather than the name-class field. A synonym whose TEXT contains that phrase then
       overwrites the scientific name of its node.

    The OLD-scheme block is a regression fence: its expectations are the lineages the
    pre-2025 rank scheme yields on this same fixture, so a dump of that vintage must still
    yield them byte for byte.

    Run:  julia --project=. test/test_taxonomy_tools.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
using DataFrames
using SHA
include(joinpath(@__DIR__, "..", "src", "ModuleTaxonomyTools.jl"))
using .TaxonomyTools

# ---------------------------------------------------------------------------
# Fixture writers. Rows are (tax_id, parent_id, rank, scientific name); the field layout
# is NCBI's own: values separated by "\t|\t", every line closed with "\t|".
# ---------------------------------------------------------------------------
function write_dump(dir, rows; merged = nothing, deleted = nothing, extra_name_lines = String[])
    open(joinpath(dir, "nodes.dmp"), "w") do io
        for (tid, pid, rank, _) in rows
            println(io, "$tid\t|\t$pid\t|\t$rank\t|\t\t|\t0\t|\t1\t|")
        end
    end
    open(joinpath(dir, "names.dmp"), "w") do io
        for (tid, _, _, name) in rows
            println(io, "$tid\t|\t$name\t|\t\t|\tscientific name\t|")
        end
        for line in extra_name_lines
            println(io, line)
        end
    end
    if merged !== nothing
        open(joinpath(dir, "merged.dmp"), "w") do io
            for (old_id, new_id) in merged
                println(io, "$old_id\t|\t$new_id\t|")
            end
        end
    end
    if deleted !== nothing
        open(joinpath(dir, "delnodes.dmp"), "w") do io
            for tid in deleted
                println(io, "$tid\t|")
            end
        end
    end
    return dir
end

# Pre-March-2025 vintage: one superkingdom rank for cellular life and for viruses.
const OLD_ROWS = [
    (1,      1,      "no rank",      "root"),
    (131567, 1,      "no rank",      "cellular organisms"),
    (2,      131567, "superkingdom", "Bacteria"),
    (1224,   2,      "phylum",       "Pseudomonadota"),
    (28211,  1224,   "class",        "Alphaproteobacteria"),
    (356,    28211,  "order",        "Hyphomicrobiales"),
    (119045, 356,    "family",       "Methylobacteriaceae"),
    (407,    119045, "genus",        "Methylobacterium"),
    (408,    407,    "species",      "Methylobacterium extorquens"),
    (2157,   131567, "superkingdom", "Archaea"),
    (10239,  1,      "superkingdom", "Viruses"),
]

# Current vintage: "domain" instead of "superkingdom", an NCBI "kingdom" node wedged
# between domain and phylum, and viruses under the acellular root.
const NEW_ROWS = [
    (1,       1,       "no rank",        "root"),
    (131567,  1,       "cellular root",  "cellular organisms"),
    (2,       131567,  "domain",         "Bacteria"),
    (3379134, 2,       "kingdom",        "Pseudomonadati"),
    (1224,    3379134, "phylum",         "Pseudomonadota"),
    (28211,   1224,    "class",          "Alphaproteobacteria"),
    (356,     28211,   "order",          "Hyphomicrobiales"),
    (119045,  356,     "family",         "Methylobacteriaceae"),
    (407,     119045,  "genus",          "Methylobacterium"),
    (408,     407,     "species",        "Methylobacterium extorquens"),
    (2157,    131567,  "domain",         "Archaea"),
    (10239,   1,       "acellular root", "Viruses"),
]

# Neither top-level spelling: the shape of a dump this code does not understand.
const HEADLESS_ROWS = [
    (1,      1,      "no rank", "root"),
    (131567, 1,      "no rank", "cellular organisms"),
    (2,      131567, "no rank", "Bacteria"),
    (1224,   2,      "phylum",  "Pseudomonadota"),
    (408,    1224,   "species", "Methylobacterium extorquens"),
]

# A synonym whose TEXT contains "scientific name" — the line that wins wherever occursin
# is matched anywhere on the line instead of on the name-class field.
const DECOY_SYNONYM = "2\t|\tno scientific name assigned\t|\t\t|\tsynonym\t|"

# ---------------------------------------------------------------------------
# Frozen lineages for OLD_ROWS under the PRE-2025 rank scheme, as that scheme yields them
# on this fixture; any drift here means a pre-2025 dump is no longer annotated the way this
# baseline records.
# ---------------------------------------------------------------------------
const OLD_BASELINE = Dict{Int, Any}(
    1      => missing,
    131567 => missing,
    2      => Dict(:kingdom => "Bacteria"),
    1224   => Dict(:kingdom => "Bacteria", :phylum => "Pseudomonadota"),
    28211  => Dict(:kingdom => "Bacteria", :phylum => "Pseudomonadota",
                   :class => "Alphaproteobacteria"),
    356    => Dict(:kingdom => "Bacteria", :phylum => "Pseudomonadota",
                   :class => "Alphaproteobacteria", :order => "Hyphomicrobiales"),
    119045 => Dict(:kingdom => "Bacteria", :phylum => "Pseudomonadota",
                   :class => "Alphaproteobacteria", :order => "Hyphomicrobiales",
                   :family => "Methylobacteriaceae"),
    407    => Dict(:kingdom => "Bacteria", :phylum => "Pseudomonadota",
                   :class => "Alphaproteobacteria", :order => "Hyphomicrobiales",
                   :family => "Methylobacteriaceae", :genus => "Methylobacterium"),
    408    => Dict(:kingdom => "Bacteria", :phylum => "Pseudomonadota",
                   :class => "Alphaproteobacteria", :order => "Hyphomicrobiales",
                   :family => "Methylobacteriaceae", :genus => "Methylobacterium",
                   :species => "Methylobacterium extorquens"),
    2157   => Dict(:kingdom => "Archaea"),
    10239  => Dict(:kingdom => "Viruses"),
)

const TMP = mktempdir()
const OLD_DIR = write_dump(mkpath(joinpath(TMP, "old")), OLD_ROWS)
const NEW_DIR = write_dump(mkpath(joinpath(TMP, "new")), NEW_ROWS;
                           merged = [(999999, 408)], deleted = [888888],
                           extra_name_lines = [DECOY_SYNONYM])
const HEADLESS_DIR = write_dump(mkpath(joinpath(TMP, "headless")), HEADLESS_ROWS;
                                merged = Tuple{Int, Int}[], deleted = Int[])

# The old fixture has no merged.dmp/delnodes.dmp, which is also how the loader announces a
# dump it cannot use to follow retired ids.
db_old = @test_logs (:warn,) (:warn,) match_mode = :any load_taxonomy_db(OLD_DIR)
db_new = load_taxonomy_db(NEW_DIR)
db_headless = load_taxonomy_db(HEADLESS_DIR)

@testset "TaxonomyTools — dump vintages, retired ids, centroids" begin

    @testset "OLD scheme: regression against the pre-2025 lineages" begin
        for (tid, expected) in OLD_BASELINE
            @test isequal(get_lineage(tid, db_old), expected)
        end
        # The two anchors stated separately, so a failure names the rank that broke.
        @test get_lineage(408, db_old)[:kingdom] == "Bacteria"
        @test get_lineage(408, db_old)[:phylum] == "Pseudomonadota"
        @test get_lineage(10239, db_old)[:kingdom] == "Viruses"
        @test db_old.merged == Dict{Int, Int}()
        @test isempty(db_old.deleted)

        # Callers written before merged.dmp/delnodes.dmp were parsed still build a db with
        # three arguments, and get the same lineages out of it.
        db_3arg = TaxonomyTools.TaxonomyDB(db_old.parents, db_old.ranks, db_old.names)
        @test isempty(db_3arg.merged) && isempty(db_3arg.deleted)
        @test isequal(get_lineage(408, db_3arg), OLD_BASELINE[408])
    end

    @testset "NEW scheme: 'domain' is the kingdom, NCBI's 'kingdom' is not" begin
        lin = get_lineage(408, db_new)
        @test lin[:kingdom] == "Bacteria"
        @test lin[:phylum] == "Pseudomonadota"
        # Pseudomonadati sits between domain and phylum and must land in no rank at all —
        # accepting rank "kingdom" would have put it in :kingdom.
        @test !("Pseudomonadati" in values(lin))
        @test get_lineage(2157, db_new)[:kingdom] == "Archaea"
        # Viruses have no domain in this vintage; the acellular root stands in for it, so
        # virus rows stay comparable with the old dump's.
        @test get_lineage(10239, db_new)[:kingdom] == "Viruses"
        # Same species, same lineage, either vintage.
        @test isequal(get_lineage(408, db_new), get_lineage(408, db_old))
    end

    @testset "merged.dmp and delnodes.dmp" begin
        @test isequal(get_lineage(999999, db_new), get_lineage(408, db_new))
        lin, status, merged_to = TaxonomyTools.lineage_with_status(999999, db_new)
        @test status == "merged"
        @test merged_to == 408
        @test lin[:species] == "Methylobacterium extorquens"

        @test ismissing(get_lineage(888888, db_new))
        _, del_status, del_merged = TaxonomyTools.lineage_with_status(888888, db_new)
        @test del_status == "deleted"
        @test ismissing(del_merged)

        # A live id is neither, and says so.
        @test TaxonomyTools.lineage_with_status(408, db_new)[2] == "ok"
    end

    @testset "names.dmp: the name CLASS decides, not the line text" begin
        @test db_new.names[2] == "Bacteria"
        @test get_lineage(408, db_new)[:kingdom] == "Bacteria"
    end

    @testset "append_taxonomy!: columns, status and the coverage guard" begin
        make_df() = DataFrame(id = ["c1", "c2", "c3", "c4", "c5"],
                              TaxID = Union{Int, Missing}[408, 999999, 888888, missing, 10239])

        df = make_df()
        # Three of the four rows that carry a TaxID resolve to a kingdom (the deleted one
        # cannot), so the run needs a threshold below that to get past the guard.
        append_taxonomy!(df, db_new; min_kingdom_coverage = 0.5)
        @test isequal(df.kingdom, ["Bacteria", "Bacteria", missing, missing, "Viruses"])
        @test df.species[1] == "Methylobacterium extorquens"
        @test df.taxid_lineage_status == ["ok", "merged", "deleted", "no_taxid", "ok"]
        @test isequal(df.taxid_merged_to, Union{Int, Missing}[missing, 408, missing, missing, missing])
        # :TaxID itself is never rewritten: the merged target is reported alongside it.
        @test isequal(df.TaxID, Union{Int, Missing}[408, 999999, 888888, missing, 10239])

        # The same frame under the default threshold: the guard has to bite.
        @test_throws ErrorException append_taxonomy!(make_df(), db_new)

        # A live id that reaches no kingdom is reported as "no_domain" in the column, which
        # is what the guard's error message tells the reader to look for. Nine rows in ten
        # resolve here, so the default threshold is met exactly.
        df_edge = DataFrame(TaxID = Union{Int, Missing}[fill(408, 9); 131567])
        append_taxonomy!(df_edge, db_new)
        @test df_edge.taxid_lineage_status[10] == "no_domain"
        @test ismissing(df_edge.kingdom[10])

        # A dump whose top rank this code does not know: every lineage lacks a kingdom.
        df_headless = DataFrame(id = ["c1", "c2"], TaxID = Union{Int, Missing}[408, 1224])
        @test_throws ErrorException append_taxonomy!(df_headless, db_headless;
                                                     min_kingdom_coverage = 0.5)
        @test TaxonomyTools.lineage_with_status(408, db_headless)[2] == "no_domain"
    end

    @testset "filter_centroids: seeds only, and a .uc that is not one" begin
        uc_path = joinpath(TMP, "clusters.uc")
        # USEARCH .uc: S = centroid/seed, H = member hit, C = cluster summary. Only the S
        # labels are representatives; C repeats them and would double the set.
        open(uc_path, "w") do io
            println(io, "S\t0\t250\t*\t*\t*\t*\t*\tseq_A\t*")
            println(io, "H\t0\t250\t99.2\t+\t0\t0\t250M\tseq_B\tseq_A")
            println(io, "H\t0\t250\t97.5\t+\t0\t0\t250M\tseq_C\tseq_A")
            println(io, "S\t1\t180\t*\t*\t*\t*\t*\tseq_D\t*")
            println(io, "H\t1\t180\t95.0\t+\t0\t0\t180M\tseq_E\tseq_D")
            println(io, "C\t0\t3\t*\t*\t*\t*\t*\tseq_A\t*")
            println(io, "C\t1\t2\t*\t*\t*\t*\t*\tseq_D\t*")
        end
        df = DataFrame(id = ["seq_A", "seq_B", "seq_C", "seq_D", "seq_E"], score = 1:5)

        kept = filter_centroids(df, uc_path)
        @test kept.id == ["seq_A", "seq_D"]
        @test filter_centroids(df, uc_path; expected = 2).id == ["seq_A", "seq_D"]
        @test_throws ErrorException filter_centroids(df, uc_path; expected = 3)

        # Hit records only: with no seed record to key on, the census filters down to zero
        # rows and says nothing.
        hits_only = joinpath(TMP, "no_seeds.uc")
        open(hits_only, "w") do io
            println(io, "H\t0\t250\t99.2\t+\t0\t0\t250M\tseq_B\tseq_A")
            println(io, "C\t0\t3\t*\t*\t*\t*\t*\tseq_A\t*")
        end
        @test_throws ErrorException filter_centroids(df, hits_only)
    end

    @testset "taxdump_provenance" begin
        prov = taxdump_provenance(NEW_DIR)
        @test Set(keys(prov)) == Set(["nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"])
        @test all(prov[f]["present"] for f in keys(prov))
        @test prov["nodes.dmp"]["sha256"] ==
              bytes2hex(sha256(read(joinpath(NEW_DIR, "nodes.dmp"))))
        @test prov["nodes.dmp"]["bytes"] == filesize(joinpath(NEW_DIR, "nodes.dmp"))
        @test prov["nodes.dmp"]["mtime"] > 0
        # The old fixture has no merged.dmp: a manifest must record the absence rather than
        # leave the reader guessing which files the run had.
        prov_old = taxdump_provenance(OLD_DIR)
        @test prov_old["merged.dmp"]["present"] == false
        @test prov_old["names.dmp"]["present"] == true
    end

    @testset "merge_taxonomy_with_hits joins instead of raising UndefVarError" begin
        hits = DataFrame(id = ["c1", "c2"], TaxID = [408, 2157])
        meta = DataFrame(TaxID = [408], habitat_category = ["Terrestrial"],
                         habitat_specific = ["soil"])
        merged = merge_taxonomy_with_hits(hits, meta)
        @test nrow(merged) == 2
        @test count(.!ismissing.(merged.habitat_specific)) == 1
    end
end
