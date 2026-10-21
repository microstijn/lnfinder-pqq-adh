#=
    Offline self-test for scripts/habitat_enrichment.jl — the consumer that turns the
    background carrier scan, its three persisted records and the habitat calls into the
    codes a downstream text step reads. No network, no HMMER and no repository data file:
    every input below is built in a temp directory at the size of the property under test.

    The script is included behind its own `PROGRAM_FILE` guard, so including it defines the
    functions and runs nothing. One testset starts a copy as a subprocess, because what the
    driver does when the file IS the program can only be measured from outside.

    EVERY FIXTURE IS ABSTRACT. Habitats and phyla are role names such as "cat-a" and
    "phy-1", taxids are small arbitrary integers, and bars and cutoffs are arbitrary values
    chosen so the property under test is decided. Nothing here is a quantity measured on the
    data this script analyses: a suite that encodes a run's own numbers passes only until the
    next run, and it publishes what the run found.

    The failures this suite exists to prevent:

    1. A LICENCE WIRED BY CODE NO TEST CAN REACH. Guards written as top-level script code
       can be loaded by no suite: the evidence verdict replaced by `true`, the file-state
       gate dropped, or the record's own pass field forwarded in place of the re-derived
       verdict would each license the strong code on evidence that fails, with every suite
       still green. Every one of those inputs is passed through one function, and each
       mutation of it has to fail a case below.

    2. A TABLE HASHED FROM A SECOND READ. The digest is what lets a record be refused for
       describing another table, so it has to come from the very bytes that were parsed. With
       two reads, a scan run promoting in between leaves the run analysing one table while its
       records vouch for another, and every digest agrees. The reader is injectable here for
       exactly that reason: the case below answers differently on a second read, so an
       implementation that reads twice cannot pass it.

    3. A RECORD READ TWICE, HALF OF IT CHECKED. The threshold and recall records are read two
       ways — typed, for the fields the licence rests on, and as written, for the provenance
       keys the summary echoes. Two reads are two moments, and only one of them was ever
       digest-checked.

    4. A PERFECT RATE ON A TINY SAMPLE. A sample that recovered every carrier it flagged has a
       rate of one at any size. The verdict is taken on the exact lower confidence bound of
       the same two counts, so the same counts must fall back under a bar the bound cannot
       reach and license under one it does. Both cases are below, and neither names a number
       measured on any data.

    5. A RECORD'S OTHER NUMBERS TAKEN AS WRITTEN. With every digest intact, a rate that
       disagrees with its own counts, a gated acceptance row that was cut out, a flagged count
       the table read does not have, and a verdict word that disagrees with the numbers beside
       it each licensed the strong code. Each is spoiled on its own below.

    6. AN EXCLUSION IMPLEMENTED TWICE. The any-hit exclusion is the producer's own tested
       predicate; a second copy of it here would drift, and the copy analysis ran would not be
       the one the producer's suite tests.

    7. A LEVEL THAT MOVES NOBODY, OR MOVES THE SAME TAXON TWICE. The move-in level is the one
       specification that ADDS taxa to the foreground, so it is the one that can break a
       result rather than merely move it — and it has two silent failures. Its taxa carry
       their domain in `domain` and leave the census's `kingdom` column empty, so a
       prokaryote test read off `kingdom` drops every one of them and reports the unmoved
       specification under this level's name; and a moved taxon left in the background would
       put an arm against a background holding part of itself. The cases below pin both, and
       pin that no other mode goes anywhere near the new path.

    8. A SEED NOBODY READ. The background pool is a seeded subsample, so its seed belongs in
       this step's provenance. A seed typed into this script, or a default filled in when the
       generator cannot be read, would state the seed of a subsample the script did not take.
       The seed is read out of the generator's source, only from a single unambiguous line,
       and every other outcome is a named state with the value left missing.

    Run:  julia --project=. test/test_habitat_enrichment.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Test, DataFrames, CSV, SHA, Printf

# The script goes into a module of its own, the way it includes the producer: the aggregate
# runner includes every suite into one process, and a script's constants must not collide
# with another suite's. `HabitatConsumer.ScanRecords` is the producer it brought in with it.
module HabitatConsumer end
Base.include(HabitatConsumer, joinpath(@__DIR__, "..", "scripts", "habitat_enrichment.jl"))
const HE = HabitatConsumer
const SR = HabitatConsumer.ScanRecords

@testset "habitat_enrichment consumer" begin

# ---------------------------------------------------------------- fixtures
# A frame of six taxa: two carry hits and a family signature, one carries neither but was
# searched, and three were never searched. Roles, not a sample of anything.
fixture_scan(; taxa = [101, 102, 103, 104, 105, 106], flagged = [101, 106],
               family = [101, 106]) =
    DataFrame(taxid = taxa,
              n_candidates = [t in flagged ? 2 : (t == 102 ? 1 : 0) for t in taxa],
              n_hits = [t in flagged ? 2 : 0 for t in taxa],
              hit_classes = [t in flagged ? "Lanthanide" : "" for t in taxa],
              in_ipr017512_ge500 = [t in family for t in taxa],
              n_hits_ln_or_ca = [t in flagged ? 1 : 0 for t in taxa])

# A recall sample in which the full-proteome scan flagged `n` taxa and the prefilter
# recovered `recovered` of them. The bar is an argument because the two cases this round
# turns on — a clean sample whose bound falls short of the bar, and the same counts under a
# bar the bound reaches — differ in nothing else.
fixture_recall(; n_flagged = 9, recovered = 9, bar = 0.85, n_sampled = 30) =
    SR.recall_summary(n_sampled = n_sampled, full_flagged = collect(1:n_flagged),
                      pre_flagged = collect(1:recovered), n_never_searched = 4,
                      bar = bar, seed = 12345)

FIXTURE_BITS    = 12.5
FIXTURE_RELEASE = "REL-A"

fixture_threshold_pairs(; threshold_bits = FIXTURE_BITS) =
    ["rule" => "fixture", "threshold_bits" => @sprintf("%.4f", threshold_bits),
     "threshold_source" => "fixture", "active_site_column" => "17",
     "hmm" => "fixture.hmm", "hmm_sha256" => repeat("a", 64), "recovery" => "1.000000"]

# The acceptance table a passing scan run over `scan` would write: every check the consumer
# requires, built with the producer's own `accept!` so each verdict is derived from the row's
# numbers, the two recall rows from the producer's `recall_acceptance!`, and one disclosure.
# `flagged_value` plants a flagged count the table does not have; `drop` removes rows by
# check text.
function fixture_acceptance(scan, recall; flagged_value = nothing, drop = String[])
    acc = SR.acceptance_table()
    n_family = count(scan.in_ipr017512_ge500)
    n_flagged = count(>(0), scan.n_hits)
    SR.accept!(acc, HE.FAMILY_CHECK; value = n_family, denominator = n_family,
               bar = n_family, note = "bar computed by the fixture")
    SR.accept!(acc, "calibration recovery"; value = 1.0, denominator = 8,
               bar = SR.MIN_RECOVERY, note = "fixed bar")
    SR.recall_acceptance!(acc, recall)
    SR.accept!(acc, HE.FLAGGED_CHECK; value = something(flagged_value, n_flagged),
               denominator = nrow(scan), bar = n_family, note = "bar computed by the fixture")
    SR.accept!(acc, "hit sequences classed at the active-site column"; value = 4,
               denominator = 4, bar = 4, note = "bar computed by the fixture")
    for (t, name) in SR.POSITIVE_CONTROLS
        SR.accept!(acc, "positive control $name ($t)"; value = 1, bar = 1, note = "fixed bar")
    end
    SR.accept!(acc, "a quantity with no bar"; value = 1, denominator = nrow(scan),
               note = "no bar, disclosed")
    isempty(drop) && return acc
    return acc[[!(c in drop) for c in acc.check], :]
end

# One published set — the table and its three records — written by the producer's own
# `publish_scan`, so every digest and identity stamp is the producer's and not this suite's.
function publish_fixture(dir; scan = fixture_scan(), recall = fixture_recall(),
                         acceptance = nothing, threshold_bits = FIXTURE_BITS,
                         release = FIXTURE_RELEASE, run_id = "RUN-FIXTURE")
    mkpath(dir)
    paths = (scan = joinpath(dir, "background_pqq_scan.csv"),
             threshold = joinpath(dir, "threshold.txt"),
             recall = joinpath(dir, "recall.txt"),
             acceptance = joinpath(dir, "acceptance.tsv"),
             release_record = joinpath(dir, "uniprot_release.txt"))
    # The release record the scan writes beside its search-space download; it is not part of
    # the published set, so a promotion neither writes nor displaces it.
    SR.write_keyvals(paths.release_record, ["X-UniProt-Release" => release])
    pub = SR.publish_scan(; out_csv = paths.scan, csv = SR.csv_bytes(scan),
                          threshold_path = paths.threshold,
                          threshold_pairs = fixture_threshold_pairs(; threshold_bits),
                          recall_path = paths.recall,
                          recall_pairs = SR.recall_keyvals(recall; threshold_bits = threshold_bits,
                                                           uniprot_release = release),
                          acceptance_path = paths.acceptance,
                          acceptance = something(acceptance,
                                                 fixture_acceptance(scan, recall)),
                          run_id = run_id)
    return merge(paths, (promoted = pub.promoted, sha256 = pub.scan_sha256, release = release))
end

# What the licence reads off an `analyse` result: the configuration's taxa, the frame AS THE
# SCAN SAW IT, and the analysed background arm. A tuple of that shape keeps each case at the
# size of the property under test; `analyse` itself is exercised on a call table further
# down, and the two are tied together by the leftover and frame cases.
#
# `frame_scanned`, not `frame`: the move-in level hands part of the frame to the foreground,
# and the scan's coverage has to be judged against what the scan was asked about, or every
# run of that level would report its own definition as a hole in the scan.
fixture_arms(; frame = [101, 102, 103, 104, 105, 106], arm = [102, 103, 104, 105]) =
    (df = DataFrame(taxid = frame), frame_scanned = trues(length(frame)),
     frame = trues(length(frame)), bgm = DataFrame(taxid = arm))

# The chain the driver runs, from one published set to the licence, in the driver's own order.
# `pre_registered_bar` defaults to the fixture's own recorded bar, because a fixture sample is
# small and its bar is chosen per case; the production default is the scan's constant, and the
# case below plants a disagreement between the two to show the check is live.
function licence_of(paths; arms = fixture_arms(), scan_path = paths.scan,
                    pre_registered_bar = SR.read_recall(paths.recall).bar)
    state = HE.scan_table_state(scan_path)
    table = HE.read_scan_table(scan_path)
    records = HE.read_scan_records(paths.threshold, paths.recall, paths.acceptance,
                                   paths.release_record)
    return HE.scan_evidence_and_licence(arms; table_read = table, table_state = state,
                                        records = records,
                                        pre_registered_bar = pre_registered_bar)
end

# One key/value record line rewritten in place, order kept, so a planted disagreement leaves
# every other field and every digest exactly as the producer wrote them.
function edit_record!(path, key, value)
    lines = readlines(path)
    hit = [startswith(l, key * "\t") for l in lines]
    any(hit) || error("no $key line in $(basename(path))")
    write(path, join([hit[i] ? "$key\t$value" : lines[i] for i in eachindex(lines)], "\n") * "\n")
    return path
end

reason_for(res, key) = any(r -> startswith(r, key * ":"), res.reasons)

# ---------------------------------------------------------------- the category selection
# One selection file, written with abstract role names, so nothing here is a category of the
# data this script analyses and nothing here can be mistaken for the shipped choice — there is
# no shipped choice any more, which is the property these tests exist to hold.
SEL_TIER    = "cat-t"
SEL_GENUS   = "cat-g"
SEL_SHARE   = "cat-s"
SEL_GRID    = ["cat-b", "cat-q"]
SEL_PHYLUM  = "phy-2"

"""Write a category-selection criteria file; `rows` replaces the default row list wholesale."""
function selection_csv(dir::AbstractString; rows = nothing, name = "category_selection.csv")
    r = rows === nothing ?
        vcat([("expect", "selections", string(1 + 1 + 1 + length(SEL_GRID) + 1), "size note")],
             [("select", "depletion_tier", SEL_TIER, "note"),
              ("select", "genus_clause", SEL_GENUS, "note"),
              ("select", "sharing", SEL_SHARE, "note")],
             [("select", "grid_count", c, "note") for c in SEL_GRID],
             [("select", "standardisation_phylum", SEL_PHYLUM, "note")]) : rows
    path = joinpath(dir, name)
    CSV.write(path, DataFrame(kind = [x[1] for x in r], role = [x[2] for x in r],
                              value = [x[3] for x in r], note = [x[4] for x in r]))
    return path
end

SEL_PATH = selection_csv(mktempdir())
SELECTED = HE.read_category_selection(SEL_PATH)

@testset "the category selection is supplied, never defaulted, and its keys derive from it" begin
    # What the file says is what the run is about, and the key prefixes are DERIVED.
    @test SELECTED.depletion_tier == SEL_TIER
    @test SELECTED.genus_clause == SEL_GENUS
    @test SELECTED.sharing == SEL_SHARE
    @test SELECTED.grid_count == Tuple(SEL_GRID)
    @test SELECTED.standardisation_phylum == SEL_PHYLUM
    # The governed set is exactly the categories the supplied rules are about, in role order,
    # and is not a set this script may hold an opinion about.
    @test SELECTED.governed == (SEL_TIER, SEL_GENUS, SEL_SHARE)
    @test SELECTED.key.depletion_tier == "catt" && SELECTED.key.sharing == "cats"
    @test SELECTED.key.standardisation_phylum == "phy2"
    @test SELECTED.digest == bytes2hex(sha256(read(SEL_PATH)))
    # The prefix rule, on its own.
    @test HE.key_slug("Host-plant") == "hostplant"
    @test HE.key_slug("  Mixed Case 9 ") == "mixedcase9"
    @test_throws ErrorException HE.key_slug("---")
    @test_throws ErrorException HE.key_slug("9lives")
    # Every declared role must be filled, in both directions, and the file states its own size.
    dir = mktempdir()
    full = [("expect", "selections", "5", "n"),
            ("select", "depletion_tier", SEL_TIER, "n"),
            ("select", "genus_clause", SEL_GENUS, "n"),
            ("select", "sharing", SEL_SHARE, "n"),
            ("select", "grid_count", "cat-b", "n"),
            ("select", "standardisation_phylum", SEL_PHYLUM, "n")]
    @test HE.read_category_selection(selection_csv(dir; rows = full, name = "a.csv")).grid_count ==
          ("cat-b",)
    bad(rows, name) = HE.read_category_selection(selection_csv(dir; rows = rows, name = name))
    # a role nobody fills
    @test_throws ErrorException bad([r for r in full if r[2] != "sharing"], "b.csv")
    # a role this script does not declare
    @test_throws ErrorException bad(vcat(full, [("select", "colour", "x", "n")]), "c.csv")
    # a single-valued role filled twice
    @test_throws ErrorException bad(vcat(full, [("select", "sharing", "cat-x", "n")]), "d.csv")
    # a row with no note: nothing would say where the selection came from
    @test_throws ErrorException bad([r[2] == "sharing" ? (r[1], r[2], r[3], "") : r for r in full],
                                   "e.csv")
    # a row with no value: there is no default to fall back to
    @test_throws ErrorException bad([r[2] == "sharing" ? (r[1], r[2], "", r[4]) : r for r in full],
                                   "f.csv")
    # a count that does not match the rows, and a file that states no count at all
    @test_throws ErrorException bad([r[2] == "selections" ? (r[1], r[2], "4", r[4]) : r
                                     for r in full], "g.csv")
    @test_throws ErrorException bad([r for r in full if r[1] != "expect"], "h.csv")
    # two selections that reduce to one key prefix would file both under one family of keys
    @test_throws ErrorException bad([r[2] == "genus_clause" ? (r[1], r[2], "Cat T", r[4]) : r
                                     for r in full], "i.csv")
    # an unknown row kind
    @test_throws ErrorException bad(vcat(full, [("assume", "sharing", "x", "n")]), "j.csv")
    @test_throws ErrorException HE.read_category_selection(joinpath(dir, "absent.csv"))
end

# ---------------------------------------------------------------- the command line
@testset "the command line is read once, into typed values" begin
    # THE ONE OPTION WITH NO DEFAULT AT ALL. The categories each recorded rule speaks for were
    # selected from this data set's own estimates, so a run that names no criteria file is
    # refused rather than given a set this software picked.
    @test_throws ErrorException HE.options_from(Dict{String,String}())
    base() = Dict("--categories" => SEL_PATH)
    o = HE.options_from(base())
    @test o.categories == SEL_PATH
    @test o.config == "primary"
    @test o.bg_exclusion == HE.PRIMARY_MODE
    @test (o.bg_exclusion_source, o.bg_exclusion_role) == (:default, :primary)
    @test o.chem_col === :chemistry_flags
    @test o.plan_bg_size_below === nothing            # a plan bound is never defaulted
    @test o.plan_bg_phylum_pct_max === nothing
    @test endswith(o.anchor_out, "_anchor.csv")       # beside --out, not beside the repo's
    # A mode that was typed is recorded as typed, and a non-primary mode as a sensitivity run.
    s = HE.options_from(merge(base(), Dict("--bg-exclusion" => "lnca", "--out" => "/tmp/x.csv")))
    @test (s.bg_exclusion, s.bg_exclusion_source, s.bg_exclusion_role) ==
          ("lnca", :explicit, :sensitivity)
    @test s.anchor_out == "/tmp/x_anchor.csv"
    @test HE.options_from(merge(base(), Dict("--config" => "legacy"))).chem_col === :chemistry
    @test HE.options_from(merge(base(), Dict("--boot" => "0"))).nboot == 0
    # The move-in pool belongs to the one level that moves it. Handed to any other mode it
    # would be read by nothing, so the run would quietly be a different one from the command
    # that was typed — which is the failure this parser exists to stop.
    @test isempty(HE.options_from(base()).movein_pool)
    m = HE.options_from(merge(base(), Dict("--bg-exclusion" => HE.MOVEIN_MODE,
                                           "--movein-pool" => "p.csv")))
    @test m.movein_pool == "p.csv" && m.bg_exclusion_role === :sensitivity
    @test_throws ErrorException HE.options_from(merge(base(), Dict("--movein-pool" => "p.csv")))
    @test_throws ErrorException HE.options_from(merge(base(),
                                                     Dict("--bg-exclusion" => "lntype",
                                                          "--movein-pool" => "p.csv")))
    @test_throws ErrorException HE.options_from(merge(base(),
                                                     Dict("--bg-exclusion" => "hmmsearch_movein!")))
    # A bad value is an error, never a fall-back to the default under the name that was typed.
    @test_throws ErrorException HE.options_from(merge(base(), Dict("--bg-exclusion" => "HMMSEARCH")))
    @test_throws ErrorException HE.options_from(merge(base(), Dict("--chem-col" => "bogus")))
    @test_throws ErrorException HE.options_from(merge(base(), Dict("--boot" => "-1")))
    @test_throws ErrorException HE.options_from(merge(base(), Dict("--boot" => "2.5")))
    @test_throws ErrorException HE.options_from(merge(base(), Dict("--seed" => "abc")))
    @test_throws ErrorException HE.options_from(merge(base(), Dict("--plan-bg-size-below" => "0")))
    @test_throws ErrorException HE.options_from(merge(base(),
                                                     Dict("--plan-bg-phylum-pct-max" => "101")))
    @test_throws ErrorException HE.options_from(merge(base(), Dict("--config" => "")))
    @test_throws ErrorException HE.options_from(Dict("--categories" => ""))
    # Every option the parser defines is one the usage text names, and an unknown token stops.
    @test all(o -> occursin(o, HE.USAGE), keys(HE.OPTION_SPEC))
    @test_throws ErrorException HE.parse_options(["--bg-exlusion", "lnca"], HE.OPTION_SPEC)
    @test_throws ErrorException HE.parse_options(["--scan", "--boot", "0"], HE.OPTION_SPEC)
    @test_throws ErrorException HE.parse_options(["primary"], HE.OPTION_SPEC)
end

# ---------------------------------------------------------------- panels and category sets
@testset "the row panels a category set is read off" begin
    # Bounds by role: one interval that excludes zero, one that includes it, one that was
    # never computed, and one row outside the units a category's standing is read off.
    row(contrast, unit, habitat, lo, hi; corrected = false, zero_cell = 0, a = 3.0, b = 1.0) =
        (contrast = contrast, unit = unit, habitat = habitat,
         ci_lo_std = lo, ci_hi_std = hi, ci_lo = lo, ci_hi = hi,
         continuity_corrected = corrected, n_zero_stratum_cells = zero_cell, a = a, b = b)
    rows = [row(HE.LN_VS_BG, "taxa", "cat-a", 0.2, 0.9),
            row(HE.LN_VS_BG, "taxa", "cat-b", -0.3, 0.4),
            row(HE.LN_VS_BG, "taxa (chemistry)", "cat-c", NaN, NaN),
            row(HE.LN_VS_BG, "genera", "cat-d", 0.5, 1.5),
            row(HE.CA_VS_BG, "taxa", "cat-e", 0.1, 0.2; corrected = true, zero_cell = 2, b = 0.0)]
    sets = HE.interval_sets(rows, :std)
    @test sets.excludes == ["cat-a"]                    # the genera row is not a primary unit
    @test sets.undecidable == ["cat-c"]
    @test sets.labels == ["cat-a", "cat-b", "cat-c"]
    @test sets.n_undecidable_rows == 1
    @test HE.interval_sets(rows, :crude).excludes == ["cat-a"]
    # THE REGRESSION. A category the anchor measured and the primary panel has no row for is
    # undecidable in the primary run, not a category that stopped excluding zero.
    @test HE.undecidable_categories(sets, ["cat-a", "cat-b", "cat-c", "cat-z"]) ==
          ["cat-c", "cat-z"]
    @test HE.undecidable_categories(sets, String[]) == ["cat-c"]
    # A disclosed arm's row counts in the plain counters and never in the registered ones.
    @test length(HE.registered_rows(rows)) == 4
    @test count(r -> r.continuity_corrected, rows) == 1
    @test count(r -> r.continuity_corrected, HE.registered_rows(rows)) == 0
    @test count(r -> r.n_zero_stratum_cells > 0, HE.registered_rows(rows)) == 0
    @test count(r -> r.a == 0 || r.b == 0, HE.registered_rows(rows)) == 0
    @test HE.pick(rows, HE.LN_VS_BG, "taxa", "cat-b").ci_lo_std == -0.3
    @test HE.pick(rows, HE.LN_VS_BG, "taxa", "cat-z") === nothing
    @test HE.row_bounds(nothing, :std) === (NaN, NaN)
    @test HE.zero_std(nothing) === :undecidable          # no row is no measurement
end

# ---------------------------------------------------------------- the arm builder
@testset "the arms, and the arm the leftovers are counted over" begin
    # Three Ln carriers, two other census carriers, five frame taxa, one of them without a
    # usable term. The taxids are the scan fixture's, so the exclusion set means something.
    spec = [(101, true, false, 2, 2, "cat-a"), (102, true, false, 1, 1, "cat-a"),
            (103, true, false, 3, 1, "cat-b"), (104, true, false, 2, 0, "cat-a"),
            (105, true, false, 1, 0, "cat-b"), (201, false, true, 0, 0, "cat-a"),
            (202, false, true, 0, 0, "cat-b"), (203, false, true, 0, 0, "cat-a"),
            (204, false, true, 0, 0, ""), (205, false, true, 0, 0, "cat-b")]
    calls = DataFrame(config = fill("primary", length(spec)),
                      taxid = [s[1] for s in spec], in_census = [s[2] for s in spec],
                      in_pool = [s[3] for s in spec], n_census_seqs = [s[4] for s in spec],
                      n_ln_seqs = [s[5] for s in spec],
                      habitat_type = [s[6] for s in spec],
                      type_weights = [isempty(s[6]) ? "" : s[6] * ":1.0" for s in spec],
                      kingdom = fill("Bacteria", length(spec)),
                      phylum = [isodd(s[1]) ? "phy-1" : "phy-2" for s in spec],
                      genus = ["gen-$(s[1] % 3)" for s in spec],
                      organism = ["Org $(s[1])" for s in spec],
                      chemistry_flags = fill("chem-x", length(spec)))
    A = HE.analyse(calls, "primary", :chemistry_flags, Set([201]); ca_class_taxids = nothing)
    @test A.n1 == 3                                   # the Ln carriers with a usable term
    @test A.n2 == 3                                   # the frame, minus 201 and the blank term
    @test sort(A.bgm.taxid) == [202, 203, 205]
    @test sort(A.df.taxid[A.frame]) == [201, 202, 203, 204, 205]
    @test !isempty(A.rows)
    @test isempty(A.ca) && isempty(A.caw)             # no census: the disclosed arm is left out
    @test_throws ErrorException HE.analyse(calls, "nosuch", :chemistry_flags, Set{Int}();
                                           ca_class_taxids = nothing)
    # The leftovers are counted over the ANALYSED arm, not over the frame or the scan table.
    scan = fixture_scan(taxa = [201, 202, 203, 204, 205], flagged = [201, 202],
                        family = [203])
    left = HE.leftover_counts(scan, A.bgm.taxid)
    @test left.hmm_hit == 1                           # 202 stayed in the arm
    @test left.ipr017512_ge500 == 1                   # 203, which no hit flagged
    @test left.union == 2
    @test HE.leftover_counts(scan, A.df.taxid[A.frame]).hmm_hit == 2   # the frame, for contrast
    # Nothing was moved, so the frame the scan saw IS the frame that was analysed, and this
    # call table has no domain column at all — the shape an older table has.
    @test A.n_moved == 0 && A.frame == A.frame_scanned
    @test A.n1 == count(A.ln_census .& A.usable)
end

@testset "a taxon with no phylum is in no arm" begin
    # Phylum is the stratum of the standardised contrast, so a taxon with no phylum — NCBI's
    # catch-all placeholders are the case that arises — belongs to no stratum and enters no
    # arm. One of each role the arm builder assembles, each with a usable term: a census Ln
    # carrier (102), a census non-Ln carrier (104), a frame taxon (202), and a frame taxon the
    # move-in pool names (204).
    spec = [(101, true,  false, 2, 2, "cat-a", "phy-1"),
            (102, true,  false, 1, 1, "cat-b", ""),
            (103, true,  false, 2, 2, "cat-b", "phy-2"),
            (104, true,  false, 1, 0, "cat-a", ""),
            (105, true,  false, 2, 0, "cat-b", "phy-1"),
            (201, false, true,  0, 0, "cat-a", "phy-1"),
            (202, false, true,  0, 0, "cat-b", ""),
            (203, false, true,  0, 0, "cat-b", "phy-2"),
            (204, false, true,  0, 0, "cat-a", ""),
            (205, false, true,  0, 0, "cat-a", "phy-1")]
    calls = DataFrame(config = fill("primary", length(spec)),
                      taxid = [s[1] for s in spec], in_census = [s[2] for s in spec],
                      in_pool = [s[3] for s in spec], n_census_seqs = [s[4] for s in spec],
                      n_ln_seqs = [s[5] for s in spec],
                      habitat_type = [s[6] for s in spec],
                      type_weights = [s[6] * ":1.0" for s in spec],
                      kingdom = fill("Bacteria", length(spec)),
                      phylum = [s[7] for s in spec],
                      genus = ["gen-$(s[1] % 3)" for s in spec],
                      organism = ["Org $(s[1])" for s in spec],
                      chemistry_flags = fill("chem-x", length(spec)))
    A = HE.analyse(calls, "primary", :chemistry_flags, Set{Int}(); ca_class_taxids = nothing)
    @test sort(A.fg.taxid)  == [101, 103]
    @test sort(A.ac.taxid)  == [101, 103, 105]
    @test sort(A.nl.taxid)  == [105]
    @test sort(A.bgm.taxid) == [201, 203, 205]
    @test !any(isempty, A.fg_phy) && !any(isempty, A.bg_phy)
    # The frame the scan saw still holds the phylum-less taxa: they leave the arms, not the
    # record of what was scanned.
    @test 202 in A.df.taxid[A.frame] && 204 in A.df.taxid[A.frame]
    # A pool that names a phylum-less frame taxon does not move it: it stays out of both sides.
    M = HE.analyse(calls, "primary", :chemistry_flags, Set{Int}(); ca_class_taxids = nothing,
                   moved_in = Set([204, 205]))
    @test M.n_moved == 1
    @test sort(M.fg.taxid)  == [101, 103, 205]
    @test sort(M.bgm.taxid) == [201, 203]
end

# ---------------------------------------------------------------- the move-in level
# A call table shaped like the real one in the one respect this level turns on: `kingdom` is
# filled from the CENSUS lineage, so every taxon the census never registered — which is every
# frame taxon — has an empty kingdom cell and carries its domain in `domain`. A prokaryote
# test read off kingdom alone therefore drops each of them at the mask.
function movein_calls()
    spec = [(11,  true,  false, 2, 2, "cat-a", "Bacteria", ""),
            (12,  true,  false, 1, 1, "cat-b", "Bacteria", ""),
            (13,  true,  false, 1, 1, "cat-a", "Archaea",  ""),
            (101, false, true,  0, 0, "cat-a", "", "Bacteria"),
            (102, false, true,  0, 0, "cat-b", "", "Bacteria"),
            (103, false, true,  0, 0, "cat-a", "", "Bacteria"),
            (104, false, true,  0, 0, "",      "", "Bacteria"),
            (105, false, true,  0, 0, "cat-b", "", "Eukaryota")]
    n = length(spec)
    return DataFrame(config = fill("primary", n), taxid = [s[1] for s in spec],
                     in_census = [s[2] for s in spec], in_pool = [s[3] for s in spec],
                     n_census_seqs = [s[4] for s in spec], n_ln_seqs = [s[5] for s in spec],
                     habitat_type = [s[6] for s in spec],
                     type_weights = [isempty(s[6]) ? "" : s[6] * ":1.0" for s in spec],
                     kingdom = [s[7] for s in spec], domain = [s[8] for s in spec],
                     phylum = fill("phy-1", n), genus = ["gen-$(s[1])" for s in spec],
                     organism = ["Org $(s[1])" for s in spec],
                     chemistry_flags = fill("chem-x", n))
end

@testset "the moved taxa are in the foreground and not in the background" begin
    calls = movein_calls()
    arms(pool; excluded = Set{Int}()) =
        HE.analyse(calls, "primary", :chemistry_flags, excluded;
                   ca_class_taxids = nothing, moved_in = pool)
    base = arms(Set{Int}())
    @test (base.n1, base.n2) == (3, 4)          # three census carriers; the frame, less a blank
    # 101 is annotated and prokaryotic, 104 is prokaryotic with no usable call, 105 is a frame
    # taxon that is not a prokaryote, and 999 is not a frame taxon at all.
    moved = arms(Set([101, 104, 105, 999]))
    @test moved.n_moved == 2 && moved.n_moved_usable == 1
    @test 101 in moved.fg.taxid                        # it joined the foreground ...
    @test !(101 in moved.bgm.taxid)                    # ... and left the analysed background
    @test !(101 in moved.df.taxid[moved.frame])        # ... and the frame it was drawn from
    # THE INVARIANT THE ARM BUILDER ASSERTS, held here as well: no taxon on both sides, and
    # none entering the foreground twice. A level that left a moved taxon in the background
    # would compare an arm with a background containing part of itself.
    @test isempty(intersect(Set(moved.fg.taxid), Set(moved.bgm.taxid)))
    @test isempty(intersect(Set(moved.fg.taxid), Set(moved.df.taxid[moved.frame])))
    @test allunique(moved.fg.taxid)
    # n1 grows by the moved taxa that carry a usable call, and by nothing else; n2 loses the
    # one the exclusion had left in.
    @test moved.n1 == base.n1 + moved.n_moved_usable == 4
    @test moved.n2 == base.n2 - 1 == 3
    @test sort(moved.fg.taxid) == [11, 12, 13, 101]
    @test sort(moved.bgm.taxid) == [102, 103, 105]
    # The frame the scan saw is untouched: the level reassigns part of it, it does not change
    # which taxa were scanned.
    @test moved.df.taxid[moved.frame_scanned] == base.df.taxid[base.frame_scanned]
    @test sort(moved.df.taxid[moved.frame_scanned]) == [101, 102, 103, 104, 105]
    # THE PER-SEQUENCE UNIT DOES WHAT THE COMMENT SAYS. A moved taxon has no census sequence —
    # that is what made it background — so it contributes no per-sequence unit and the series
    # is the census foreground against the post-move background.
    @test moved.n_moved_census_sequences == 0
    @test length(moved.seq_idx) == length(base.seq_idx) == 4
    @test HE.pick(moved.rows, HE.LN_VS_BG, "sequences", "cat-a").n1 == 4
    @test HE.pick(moved.rows, HE.LN_VS_BG, "taxa", "cat-a").n1 == 4        # the taxa unit has it
    @test HE.pick(moved.rows, HE.LN_VS_BG, "sequences", "cat-a").n2 == 3   # the moved background
    # The genus unit does take it, because a genus label is not a census column.
    @test length(moved.gfg.genus) == length(base.gfg.genus) + 1
    # A pool that names no frame taxon moves nobody and changes nothing — a census member
    # named there is already in the foreground, and moving it would count it twice.
    outside = arms(Set([11, 999]))
    @test (outside.n1, outside.n2, outside.n_moved) == (base.n1, base.n2, 0)
    @test count(outside.named_in_pool) == 0
    # A LEVEL THAT MOVES NOBODY IS REFUSED. The pool names a frame taxon and none is moved:
    # that is what a prokaryote test taken off the empty `kingdom` column does to every moved
    # taxon, and it is also what a pool of non-prokaryotes does. Either way the run would
    # report the unmoved specification under this level's name.
    @test_throws ErrorException arms(Set([105]))
    @test_throws ErrorException arms(Set([105, 999]))
end

@testset "the move-in pool: read by rule, counted at run time" begin
    dir = mktempdir()
    hits = joinpath(dir, "hits.csv")
    CSV.write(hits, DataFrame(id = ["s1", "s2", "s3", "s4"],
                              classification = ["Lanthanide", "Calcium", "Lanthanide",
                                                "Lanthanide"],
                              TaxID = [101, 102, 103, missing]))
    p = HE.read_movein_pool(hits)
    @test p.taxids == Set([101, 103])                  # the class decides, not the row count
    @test (p.rows, p.class_rows, p.rows_without_a_taxon) == (4, 3, 1)
    @test p.sha256 == bytes2hex(sha256(read(hits)))    # the digest of the bytes that were read
    @test "classification" in p.columns
    # The repository's own spelling of the taxon column is accepted beside the census's.
    lower = joinpath(dir, "lower.csv")
    CSV.write(lower, DataFrame(taxid = [7, 8], classification = ["Lanthanide", "Gap"]))
    @test HE.read_movein_pool(lower).taxids == Set([7])
    # A file that is not a hit table carries a number and no definition, so it is refused
    # rather than read as a pool of nobody.
    plain = joinpath(dir, "list.csv"); CSV.write(plain, DataFrame(taxid = [1, 2]))
    @test_throws ErrorException HE.read_movein_pool(plain)
    @test_throws ErrorException HE.read_movein_pool(joinpath(dir, "nothing-here.csv"))
    # The rule is named in the result, never assumed by the caller.
    scan = fixture_scan()
    byfile = HE.movein_pool(; scan = scan, path = hits)
    @test byfile.rule == HE.REGISTERED_MOVEIN_POOL == HE.MOVEIN_POOL_RULES[1]
    @test byfile.rule_source === :explicit && byfile.taxids == Set([101, 103])
    @test byfile.source == "hits.csv"
    byscan = HE.movein_pool(; scan = scan, path = "")
    @test byscan.rule == HE.MOVEIN_POOL_RULES[2] && byscan.rule_source === :default
    # The scan rule is the lntype set, from the same function that mode uses: a second copy
    # of the predicate would drift away from the one the mode ran.
    @test byscan.taxids == HE.exclusion_taxids(scan, "lntype")
    @test byscan.sha256 === nothing
    # No table, no rule. `nothing` is not an empty pool: an empty pool is a level that moved
    # nobody, which is the one outcome this variant must never report by accident.
    @test HE.movein_pool(; scan = nothing, path = "").taxids === nothing
    @test_throws ErrorException HE.movein_pool_note("bogus")
    @test occursin("census", HE.movein_pool_note(HE.MOVEIN_POOL_RULES[1]))
    @test occursin("DIFFERENT", HE.movein_pool_note(HE.MOVEIN_POOL_RULES[2]))
end

# ---------------------------------------------------------------- the exclusion
@testset "the any-hit exclusion is the producer's own predicate" begin
    scan = fixture_scan()
    @test HE.exclusion_taxids(scan, "hmmsearch") == Set([101, 106])
    # The producer's tested function, called on the same table, and not a copy of its rule.
    @test HE.exclusion_taxids(scan, "hmmsearch") ==
          Set(SR.exclude_carriers(DataFrame(taxid = scan.taxid), scan).removed.taxid)
    @test HE.exclusion_taxids(scan, "none") == Set{Int}()
    @test HE.exclusion_taxids(scan, "lntype") == Set([101, 106])
    @test HE.exclusion_taxids(scan, "lnca") == Set([101, 106])
    @test HE.exclusion_taxids(scan, "ipr017512") == Set([101, 106])
    @test HE.exclusion_taxids(nothing, "hmmsearch") === nothing
    @test HE.exclusion_taxids(nothing, "none") == Set{Int}()   # a mode needing no table
    # A table without the column a mode reads cannot answer it, and an empty set would look
    # like a mode that excluded nobody.
    @test HE.exclusion_taxids(select(scan, Not(:n_hits_ln_or_ca)), "lnca") === nothing
    @test HE.exclusion_taxids(select(scan, Not(:in_ipr017512_ge500)), "ipr017512") === nothing
    @test_throws ErrorException HE.exclusion_taxids(scan, "bogus")
    @test_throws ErrorException HE.exclusion_note("bogus")
    @test occursin("hmmsearch", HE.exclusion_note("hmmsearch"))
    # The move-in level REMOVES what the registered primary removes — it is that exclusion
    # plus a reassignment — so it shares the predicate instead of carrying a second copy.
    @test HE.exclusion_taxids(scan, HE.MOVEIN_MODE) == HE.exclusion_taxids(scan, "hmmsearch")
    @test occursin("foreground", HE.exclusion_note(HE.MOVEIN_MODE))
    # And it is not one of the modes a run tallies beside the mode in force: those counts are
    # background sizes over the frame of the run that is executing, and a level which hands
    # taxa to the other arm has neither that n2 nor that n1.
    @test !(HE.MOVEIN_MODE in HE.SENSITIVITY_MODES)
    @test HE.SENSITIVITY_MODES == Tuple(m for m in HE.EXCLUSION_MODES if m != HE.MOVEIN_MODE)
    @test HE.PRIMARY_MODE != HE.MOVEIN_MODE
end

# ---------------------------------------------------------------- one read of the table
@testset "the table's digest is of the bytes that were parsed" begin
    dir = mktempdir()
    p = publish_fixture(dir)
    parsed = read(p.scan)
    @test p.promoted
    # A promotion lands after the read: the file now holds another table, and a reader that
    # answers with the bytes of the first read is what a single read looks like from outside.
    write(p.scan, SR.csv_bytes(fixture_scan(taxa = [101, 102, 103, 104, 105, 106, 107])))
    reads = Ref(0)
    one_read(path) = (reads[] += 1; reads[] == 1 ? parsed : read(path))
    got = HE.read_scan_table(p.scan; read_bytes = one_read)
    @test reads[] == 1                                       # one read, so one moment
    @test got.sha256 == bytes2hex(sha256(parsed))            # the digest of what was parsed
    @test got.table.taxid == fixture_scan().taxid            # the table those bytes hold
    @test got.sha256 != bytes2hex(sha256(read(p.scan)))      # not the table now on disk
    @test got.sha256 == p.sha256                             # which is what the records name
    # The plain reader agrees with the producer's own digest of the bytes it published.
    q = publish_fixture(mktempdir())
    plain = HE.read_scan_table(q.scan)
    @test plain.sha256 == q.sha256 == bytes2hex(sha256(read(q.scan)))
    @test nrow(plain.table) == nrow(fixture_scan())
    # No file is no table and no digest — never an empty table, which reads as a clean frame.
    absent = HE.read_scan_table(joinpath(dir, "nothing-here.csv"))
    @test absent.table === nothing && absent.sha256 === nothing
end

# ---------------------------------------------------------------- one read of each record
@testset "each record is read once, and both views come from that read" begin
    dir = mktempdir()
    p = publish_fixture(dir)
    pinned = Dict(f => read(f) for f in (p.threshold, p.recall, p.release_record))
    # The records on disk are replaced, as a promotion would replace them, before the read.
    other = publish_fixture(mktempdir(); recall = fixture_recall(bar = 0.65),
                            run_id = "RUN-OTHER", release = "REL-OTHER")
    for (mine, theirs) in ((p.threshold, other.threshold), (p.recall, other.recall),
                           (p.release_record, other.release_record))
        write(mine, read(theirs))
    end
    reads = Dict{String,Int}()
    once(path) = (reads[path] = get(reads, path, 0) + 1;
                  reads[path] == 1 && haskey(pinned, path) ? pinned[path] : read(path))
    r = HE.read_scan_records(p.threshold, p.recall, p.acceptance, p.release_record;
                             read_bytes = once)
    @test reads[p.threshold] == 1 && reads[p.recall] == 1 && reads[p.release_record] == 1
    @test r.release_kv["X-UniProt-Release"] == FIXTURE_RELEASE && r.release_present
    # The typed fields and the record as written are the SAME read: both are the first one.
    @test r.recall.run_id == "RUN-FIXTURE" == r.recall_kv["run_id"]
    @test r.threshold.run_id == "RUN-FIXTURE" == r.threshold_kv["run_id"]
    @test r.recall.bar == 0.85 && r.recall_kv["bar"] == "0.850000"
    @test r.recall.scan_sha256 == p.sha256 == r.recall_kv["scan_sha256"]
    @test r.recall.lower_bound_95 === SR.clopper_pearson_lower(9, 9)
    @test r.threshold.threshold_bits == FIXTURE_BITS
    @test r.recall_present && r.threshold_present && r.acceptance_present
    @test isempty(r.acceptance_error)
    @test nrow(r.acceptance) == nrow(fixture_acceptance(fixture_scan(), fixture_recall()))
    @test unique(r.acceptance.scan_sha256) == [p.sha256]
    # A record of an older shape is a record of something else: it reads as unknown, with the
    # producer's own message, and does not end a run that can still report its numbers.
    old = joinpath(dir, "old.tsv")
    write(old, "check\tvalue\tverdict\tnote\n" * "a floor\t1\tPASS\tolder shape\n")
    o = HE.read_scan_records(p.threshold, p.recall, old, p.release_record)
    @test o.acceptance === nothing
    @test occursin("acceptance columns", o.acceptance_error)
    @test o.acceptance_present
    # Absent records are unknown in every field, never a default.
    a = HE.read_scan_records(joinpath(dir, "no.txt"), joinpath(dir, "no2.txt"),
                             joinpath(dir, "no3.tsv"), joinpath(dir, "no4.txt"))
    @test !a.threshold_present && !a.recall_present && !a.acceptance_present
    @test !a.release_present && isempty(a.release_kv)
    @test a.recall.pass === nothing && a.recall.scan_sha256 === nothing
    @test a.recall.lower_bound_95 === nothing && isempty(a.recall_kv)
    @test nrow(a.acceptance) == 0 && isempty(a.acceptance_error)
end

# ---------------------------------------------------------------- what the file name says
@testset "what the file --scan names IS" begin
    dir = mktempdir()
    p = publish_fixture(dir)
    published = HE.scan_table_state(p.scan)
    @test published.state === :published
    @test isempty(published.reasons)
    @test !published.pending_beside_it && !published.displaced_copy_on_disk
    # A displaced copy is internally consistent — it was a passing run once — so only its
    # name says that a later run refused the frame.
    stale = p.scan * SR.STALE_EXT
    cp(p.scan, stale)
    displaced = HE.scan_table_state(stale)
    @test displaced.state === :displaced
    @test length(displaced.reasons) == 1 && occursin("displaced", only(displaced.reasons))
    @test HE.scan_table_state(p.scan).displaced_copy_on_disk
    @test isempty(HE.scan_table_state(p.scan).reasons)   # a copy beside a newer run is no fault
    pending = p.scan * SR.PENDING_EXT
    cp(p.scan, pending)
    @test HE.scan_table_state(pending).state === :unpromoted
    @test occursin("never promoted", only(HE.scan_table_state(pending).reasons))
    beside = HE.scan_table_state(p.scan)
    @test beside.state === :published && beside.pending_beside_it
    @test occursin("may be superseded", only(beside.reasons))
    gone = HE.scan_table_state(joinpath(dir, "nothing-here.csv"))
    @test gone.state === :absent && occursin("does not exist", only(gone.reasons))
end

# ---------------------------------------------------------------- the table's own counts
@testset "the counts the acceptance table is held to come from the parsed table" begin
    scan = fixture_scan()
    c = HE.scan_table_counts(scan)
    @test c.n_flagged == 2 && c.n_family == 2
    @test HE.scan_table_counts(fixture_scan(flagged = Int[], family = Int[])) ==
          (n_flagged = 0, n_family = 0)
    @test HE.scan_table_counts(nothing) == (n_flagged = nothing, n_family = nothing)
    # A column that is absent, or that holds anything but whole counts or flags, is unknown
    # rather than zero: zero would read as a frame in which nothing was found.
    @test HE.scan_table_counts(select(scan, Not(:in_ipr017512_ge500))).n_family === nothing
    @test HE.scan_table_counts(select(scan, Not(:n_hits))).n_flagged === nothing
    text = copy(scan); text.n_hits = string.(scan.n_hits)
    @test HE.scan_table_counts(text).n_flagged === nothing
    gap = copy(scan); gap.n_hits = [missing; scan.n_hits[2:end]]
    @test HE.scan_table_counts(gap).n_flagged === nothing
end

# ---------------------------------------------------------------- the required gated rows
@testset "the required gated checks are the scan's own rows" begin
    src = read(joinpath(@__DIR__, "..", "scripts", "background_pqq_scan.jl"), String)
    # The family row and the controls are built from the producer's constants, so its source
    # is checked for their templates; the other three are spelled out in the consumer, and it
    # is the CONSUMER'S spelling that is looked for in the producer's rows. The fixtures below
    # build their rows from that same spelling, so a reworded check would stay consistent with
    # them and only this comparison with the producer can notice it.
    @test occursin("\"\$IPR_FAMILY frame taxa flagged " *
                   "(>= \$MIN_FAMILY_LEN aa, whole-frame membership)\"", src)
    @test occursin("\"positive control \$name (\$t)\"", src)
    spelled = [c for c in HE.REQUIRED_GATED_CHECKS
               if c != HE.FAMILY_CHECK && c != HE.RECALL_CHECK && !startswith(c, "positive control ")]
    @test length(spelled) == 3 && HE.FLAGGED_CHECK in spelled
    for check in spelled
        @test occursin("accept!(accept, \"$check\"", src)
    end
    @test HE.REQUIRED_GATED_CHECKS[end-length(SR.POSITIVE_CONTROLS)+1:end] ==
          ["positive control $name ($t)" for (t, name) in SR.POSITIVE_CONTROLS]
    # The recall row's text comes from a function the producer exposes, so it is checked by
    # calling it: the gate is on the rate row, and the bound row beside it carries no bar.
    acc = SR.recall_acceptance!(SR.acceptance_table(), fixture_recall())
    @test acc.check[1] == HE.RECALL_CHECK
    @test !isempty(acc.bar[1]) && isempty(acc.bar[2])
    @test HE.FAMILY_CHECK == "$(SR.IPR_FAMILY) frame taxa flagged " *
                             "(>= $(SR.MIN_FAMILY_LEN) aa, whole-frame membership)"
    @test length(HE.REQUIRED_GATED_CHECKS) == 5 + length(SR.POSITIVE_CONTROLS)
    @test HE.FLAGGED_CHECK in HE.REQUIRED_GATED_CHECKS
    @test issubset(HE.REQUIRED_GATED_CHECKS,
                   fixture_acceptance(fixture_scan(), fixture_recall()).check)
    # The two rows the producer gates only when there is something to gate are NOT required:
    # a first scan run, with no earlier sample and no recorded download, writes them as
    # disclosures, and requiring them would withhold every licence on a table nobody faults.
    for conditional in ("recall sample extends the earlier sample",
                        "recall proteome downloads recorded from the search-space release")
        @test occursin("accept!(accept, \"$conditional\"", src)
        @test !(conditional in HE.REQUIRED_GATED_CHECKS)
    end
end

# ---------------------------------------------------------------- an acceptance row's verdict
@testset "an acceptance row's verdict has to follow its own numbers" begin
    scan = fixture_scan(); rec = fixture_recall()
    # The table as a consumer reads it back, with numeric value, denominator and bar cells.
    acc = fixture_acceptance(scan, rec)
    published = SR.read_acceptance(SR.write_acceptance(joinpath(mktempdir(), "a.tsv"), acc))
    @test isempty(HE.acceptance_row_reasons(published, (n_full_flagged = 9, n_recovered = 9)))
    # Unknown recall counts leave the recall row nothing to be held against; the verdicts are
    # still re-derived, and the module reports the unknown counts on its own.
    @test isempty(HE.acceptance_row_reasons(published, SR.read_recall(tempname())))
    # A gated row that says PASS beside a value below its bar, and a gated row printed as a
    # disclosure, are each refused: the module is handed verdict words, not the numbers.
    for (col, value) in ((:verdict, "PASS"), (:verdict, "DISCLOSED"))
        spoilt = copy(published)
        spoilt.value[1] = spoilt.bar[1] - 1
        spoilt[1, col] = value
        rs = HE.acceptance_row_reasons(spoilt, (n_full_flagged = 9, n_recovered = 9))
        @test length(rs) == 1 && occursin(HE.FAMILY_CHECK, only(rs))
    end
    disclosure = copy(published)
    disclosure.bar[end] = 1.0                       # a bar on a row printed as a disclosure
    @test length(HE.acceptance_row_reasons(disclosure, (n_full_flagged = 9, n_recovered = 9))) == 1
    # The recall row and the recall record describe one check, so their counts must agree.
    row = findfirst(==(HE.RECALL_CHECK), published.check)
    wrong_n = copy(published); wrong_n.denominator[row] = 8
    @test any(r -> occursin("counts 8 flagged taxa", r),
              HE.acceptance_row_reasons(wrong_n, (n_full_flagged = 9, n_recovered = 9)))
    wrong_rate = copy(published); wrong_rate.value[row] = 0.5
    rs = HE.acceptance_row_reasons(wrong_rate, (n_full_flagged = 9, n_recovered = 9))
    @test any(r -> occursin("gives a rate of 0.5", r), rs)
    @test isempty(HE.acceptance_row_reasons(nothing, (n_full_flagged = 9, n_recovered = 9)))
    @test HE.acceptance_cell(published, HE.FLAGGED_CHECK, :value) == 2
    @test HE.acceptance_cell(published, "positive control", :value) === nothing   # no such row
    @test HE.acceptance_cell(published, HE.FLAGGED_CHECK, :denominator) == nrow(scan)
    @test HE.acceptance_cell(nothing, HE.FLAGGED_CHECK, :value) === nothing
    # A doubled row has no single verdict, so no cell of it may be read.
    doubled = vcat(published, published[[row], :])
    @test HE.acceptance_cell(doubled, HE.RECALL_CHECK, :bar) === nothing
end

# ---------------------------------------------------------------- the frame comparison
@testset "the scanned set against the analysed frame" begin
    scan = fixture_scan()
    same = HE.frame_comparison([101, 102, 103, 104, 105, 106], scan)
    @test same.matches && isempty(same.missing_from_scan) && isempty(same.extra_in_scan)
    # A carrier whose row is missing from the table reads as a non-carrier to every count
    # taken from it, so a frame the scan does not cover withholds the licence.
    short = HE.frame_comparison([101, 102, 103, 104, 105, 106, 107], scan)
    @test !short.matches && short.missing_from_scan == [107] && isempty(short.extra_in_scan)
    wide = HE.frame_comparison([101, 102], scan)
    @test !wide.matches && wide.extra_in_scan == [103, 104, 105, 106]
    none = HE.frame_comparison([101], nothing)
    @test none.matches === nothing                    # unknown, not "matches nothing"
end

# ---------------------------------------------------------------- the licence
@testset "the licence: a set that vouches for its table, and every way it fails to" begin
    # THE CASE THIS ROUND TURNS ON. Every flagged taxon was recovered, so the rate is one and
    # the scan's own gate passes; the exact lower bound of those counts is below the bar, so
    # the stronger background statement is not licensed.
    lb = HE.recall_lower_bound(9, 9)
    @test lb < 0.85 <= 1.0                            # the premise, from the module's own bound
    small = publish_fixture(mktempdir(); recall = fixture_recall(bar = 0.85))
    @test small.promoted                              # the scan itself publishes this set
    @test SR.read_recall(small.recall).pass === true  # and its own gate, on the rate, passed
    below = licence_of(small)
    @test below.code === :fallback
    @test below.evidence.recall_pass_derived === false
    @test below.evidence.recall_rate == 1.0
    @test below.evidence.recall_lower_bound == lb
    @test length(below.reasons) == 1                  # the bound, and nothing else, withholds it
    @test startswith(only(below.reasons), "recall_n_recovered:") &&
          occursin("lower bound", only(below.reasons))
    @test !below.ok
    # A bound that does not reach the bar is an honest measurement, not a broken record: the
    # records still vouch for the table while the licence is withheld.
    @test below.records_ok
    # The licence's own reasons carry the re-derived verdict, not the record's pass field:
    # forwarding the field would leave this reason out while the record still says pass.
    @test any(r -> startswith(r, "recall_pass:"), below.removal.reasons)

    # THE SAME COUNTS under a bar the bound reaches: every input clean, and the strong code.
    @test lb >= 0.65
    above = licence_of(publish_fixture(mktempdir(); recall = fixture_recall(bar = 0.65)))
    @test above.code === :strong
    @test above.ok && isempty(above.reasons) && above.removal.complete
    @test above.evidence.recall_pass_derived === true
    @test above.frame.matches && above.leftovers.hmm_hit == 0 &&
          above.leftovers.ipr017512_ge500 == 0
    @test above.table_counts == (n_flagged = 2, n_family = 2)
    clean() = publish_fixture(mktempdir(); recall = fixture_recall(bar = 0.65))

    # A RECORD FROM ANOTHER TABLE. The recall record of a set published over a different
    # table, beside this one: its digest names that table, and nothing else is wrong with it.
    p = clean()
    q = publish_fixture(mktempdir(); recall = fixture_recall(bar = 0.65),
                        scan = fixture_scan(taxa = [101, 102, 103, 104, 105, 106, 107]))
    @test p.sha256 != q.sha256
    cp(q.recall, p.recall; force = true)
    foreign = licence_of(p)
    @test foreign.code === :fallback
    @test reason_for(foreign, "record_sha256s") &&
          any(r -> occursin("recall record describes scan table", r), foreign.reasons)

    # A RATE THAT DISAGREES WITH ITS OWN COUNTS, every digest intact.
    p = clean()
    edit_record!(p.recall, "rate", "0.500000")
    wrong_rate = licence_of(p)
    @test wrong_rate.code === :fallback
    @test reason_for(wrong_rate, "recall_rate_field")
    @test wrong_rate.evidence.recall_pass_derived === true   # the counts still pass; the field is refused

    # THE BAR THE BOUND IS HELD TO IS ITSELF EVIDENCE. Every other input is checked against a
    # second copy; the bar's second copy is the scan's own constant, so a record set that agrees
    # with itself on a lowered bar is still refused. This is the one number a two-cell edit of
    # the records would otherwise move, and the stronger statement turns on it.
    @test SR.MIN_RECALL == 0.98              # the recorded bar, pinned where a change must show
    p = clean()                               # a self-consistent set at a bar its bound reaches
    @test licence_of(p).code === :strong
    planted = licence_of(p; pre_registered_bar = SR.MIN_RECALL)
    @test planted.code === :fallback
    @test reason_for(planted, "bar_in_other_records") &&
          any(r -> occursin("recorded recall bar", r), planted.reasons)
    @test !planted.records_ok                 # a bar that disagrees is a record failure
    # A BOUND FIELD THAT DISAGREES WITH ITS OWN COUNTS, every digest intact: the module derives
    # the bound from the counts and refuses a field that does not match it.
    p = clean()
    edit_record!(p.recall, "lower_bound_95", "0.990000")
    wrong_bound = licence_of(p)
    @test wrong_bound.code === :fallback
    @test reason_for(wrong_bound, "recall_lower_bound_field")

    # A MISSING GATED ACCEPTANCE ROW: the table is internally consistent, published by the
    # producer, and simply lacks the row. Every required row gets its own case, because a list
    # that silently loses an entry is exactly how a dropped check stops being checked.
    for cut in HE.REQUIRED_GATED_CHECKS
        acc = fixture_acceptance(fixture_scan(), fixture_recall(bar = 0.65); drop = [cut])
        missing_row = licence_of(publish_fixture(mktempdir(); recall = fixture_recall(bar = 0.65),
                                                 acceptance = acc))
        @test missing_row.code === :fallback
        @test any(r -> occursin("gated check \"$cut\" is absent", r), missing_row.reasons)
    end

    # A FLAGGED-ROW VALUE THE PARSED TABLE DOES NOT HAVE. One more flagged taxon than the table
    # holds, still at or above its bar, so the verdict re-derives PASS and only the count can
    # refuse it.
    acc = fixture_acceptance(fixture_scan(), fixture_recall(bar = 0.65); flagged_value = 3)
    @test acc.verdict[findfirst(==(HE.FLAGGED_CHECK), acc.check)] == "PASS"
    flagged = licence_of(publish_fixture(mktempdir(); recall = fixture_recall(bar = 0.65),
                                         acceptance = acc))
    @test flagged.code === :fallback
    @test reason_for(flagged, "flagged_row_value")
    @test isempty(flagged.row_reasons)

    # THE CROSS-CHECK IS AGAINST THE TABLE THAT WAS PARSED: a family member added to the table
    # the records were published over, with the records re-stamped to its digest, is refused
    # by count and by leftover — the digest alone cannot see it.
    fam = fixture_scan(family = [101, 103, 106])
    famset = publish_fixture(mktempdir(); recall = fixture_recall(bar = 0.65), scan = fam,
                             acceptance = fixture_acceptance(fixture_scan(), fixture_recall(bar = 0.65)))
    refused = licence_of(famset)
    @test refused.code === :fallback
    @test reason_for(refused, "family_row_denominator")
    @test refused.leftovers.ipr017512_ge500 == 1
    @test any(r -> startswith(r, "left_ipr_ge500:"), refused.removal.reasons)

    # A VERDICT WORD THAT DISAGREES WITH ITS NUMBERS: a gated row printed PASS below its bar.
    p = clean()
    tbl = SR.read_acceptance(p.acceptance)
    i = findfirst(==("hit sequences classed at the active-site column"), tbl.check)
    tbl.value[i] = tbl.bar[i] - 1
    lines = readlines(p.acceptance)
    cells = split(lines[i + 1], '\t'); cells[2] = SR.num_cell(tbl.value[i])
    lines[i + 1] = join(cells, '\t'); write(p.acceptance, join(lines, "\n") * "\n")
    words = licence_of(p)
    @test words.code === :fallback
    @test words.evidence.ok                           # the module sees only the verdict words
    @test length(words.row_reasons) == 1 && words.reasons == words.row_reasons

    # A MISS, A DISPLACED SET, AN UNPROMOTED SIBLING: each refused on its own. A set with a MISS
    # is never promoted, so the one way it reaches plain names is a hand copy of its .pending
    # files — which no file name can tell from a published set, and only its verdicts refuse.
    promote_by_hand!(set) = (for f in (set.scan, set.threshold, set.recall, set.acceptance)
                                 mv(f * SR.PENDING_EXT, f)
                             end; set)
    acc = fixture_acceptance(fixture_scan(), fixture_recall(bar = 0.65))
    SR.accept!(acc, "a floor that failed"; value = 0, denominator = 3, bar = 1, note = "fixture")
    missed_set = publish_fixture(mktempdir(); recall = fixture_recall(bar = 0.65), acceptance = acc)
    @test !missed_set.promoted
    missed = licence_of(promote_by_hand!(missed_set))
    @test missed.code === :fallback
    @test reason_for(missed, "acceptance_verdicts") && isempty(missed.row_reasons)
    @test HE.scan_table_state(missed_set.scan).state === :published
    p = clean()
    for f in (p.scan, p.threshold, p.recall, p.acceptance)
        mv(f, f * SR.STALE_EXT)
    end
    stale = (scan = p.scan * SR.STALE_EXT, threshold = p.threshold * SR.STALE_EXT,
             recall = p.recall * SR.STALE_EXT, acceptance = p.acceptance * SR.STALE_EXT,
             release_record = p.release_record)
    displaced = licence_of(stale)
    @test displaced.code === :fallback
    @test displaced.evidence.ok                     # the set agrees with itself: only its name
    @test !displaced.ok && any(r -> occursin("displaced", r), displaced.reasons)
    p = clean()
    cp(p.scan, p.scan * SR.PENDING_EXT)
    sibling = licence_of(p)
    @test sibling.code === :fallback && sibling.evidence.ok && !sibling.ok

    # A RECALL THAT MEASURED NOTHING, and a record written before the bound existed.
    empty_set = publish_fixture(mktempdir();
                                recall = fixture_recall(n_flagged = 0, recovered = 0, bar = 0.65))
    @test !empty_set.promoted                         # the scan refuses it as a MISS itself
    nothing_flagged = licence_of(promote_by_hand!(empty_set))
    @test nothing_flagged.code === :fallback
    @test nothing_flagged.evidence.recall_pass_derived === false
    @test reason_for(nothing_flagged, "recall_n_full_flagged")
    @test isnan(nothing_flagged.evidence.recall_lower_bound)
    p = clean()
    write(p.recall, join(filter(l -> !startswith(l, "lower_bound_95\t"), readlines(p.recall)), "\n") * "\n")
    older = licence_of(p)
    @test older.code === :fallback && reason_for(older, "recall_lower_bound_field")

    # THE OTHER RECORDS HOLD THE RECALL RECORD'S BAR, CUTOFF AND RELEASE.
    p = clean()
    edit_record!(p.recall, "bar", "0.400000")
    lowered = licence_of(p)
    @test lowered.code === :fallback && reason_for(lowered, "bar_in_other_records")
    p = clean()
    edit_record!(p.threshold, "threshold_bits", "13.5000")
    @test reason_for(licence_of(p), "threshold_bits_in_other_records")
    p = clean()
    edit_record!(p.release_record, "X-UniProt-Release", "REL-B")
    @test reason_for(licence_of(p), "uniprot_release_in_other_records")
    rm(p.release_record)
    @test reason_for(licence_of(p), "uniprot_release_in_other_records")

    # AN ACCEPTANCE RECORD OF ANOTHER SHAPE ends nothing: the evidence is unknown.
    p = clean()
    write(p.acceptance, "check\tvalue\tverdict\n" * "a floor\t1\tPASS\n")
    oldshape = licence_of(p)
    @test oldshape.code === :fallback
    @test any(r -> startswith(r, "acceptance record:"), oldshape.reasons)

    # THE FRAME AND THE LEFTOVERS ARE READ OFF THE ARMS PASSED IN.
    p = clean()
    @test licence_of(p; arms = fixture_arms(arm = [101, 102, 103])).code === :fallback
    @test licence_of(p; arms = fixture_arms(frame = [101, 102, 103, 104, 105])).code === :fallback
    @test licence_of(p).code === :strong             # and the clean set is still clean
    # NO TABLE AT ALL.
    gone = licence_of(p; scan_path = joinpath(dirname(p.scan), "nothing-here.csv"))
    @test gone.code === :fallback && gone.frame.matches === nothing
end

# ---------------------------------------------------------------- the anchor
@testset "the anchor-status ladder, rung by rung" begin
    ok = (frames_match = true, n_undecidable_rows = 0, reference_supplied = true,
          reference_exists = true, reference_readable = true,
          anchor_categories = ["cat-a", "cat-b"], reference_categories = ["cat-a", "cat-b"])
    status(; kw...) = HE.anchor_status_of(; merge(ok, (; kw...))...)
    @test status() === :reproduces_reference
    @test status(reference_categories = ["cat-a"]) === :differs_from_reference
    @test status(reference_readable = false) === :undecidable_reference_unreadable
    @test status(reference_exists = false) === :undecidable_reference_missing
    @test status(reference_supplied = false, reference_exists = false) === :computed_unverified
    @test status(n_undecidable_rows = 1) === :undecidable_anchor_intervals
    @test status(frames_match = false) === :undecidable_frame_mismatch
    # Each earlier rung wins over every later one, so a failing frame is never reported as a
    # reference problem, and an unreadable reference never as agreement.
    @test status(frames_match = false, n_undecidable_rows = 2, reference_readable = false) ===
          :undecidable_frame_mismatch
    @test status(n_undecidable_rows = 2, reference_exists = false) ===
          :undecidable_anchor_intervals
    @test status(reference_exists = false, reference_readable = false,
                 reference_categories = String[]) === :undecidable_reference_missing
    @test status(reference_readable = false, reference_categories = ["cat-a", "cat-b"]) ===
          :undecidable_reference_unreadable
    # Only an anchor that stands, verified or merely computed, may feed the partition.
    @test HE.anchor_partition_available(:reproduces_reference)
    @test HE.anchor_partition_available(:computed_unverified)
    for s in (:differs_from_reference, :undecidable_reference_unreadable,
              :undecidable_reference_missing, :undecidable_anchor_intervals,
              :undecidable_frame_mismatch, :undecidable_config_absent)
        @test !HE.anchor_partition_available(s)
    end
end

# ---------------------------------------------------------------- the background pool
# The pool the background frame was drawn from is a SEEDED subsample, so its seed is provenance
# of this step. It is read out of the generator's source — never typed into this script — and
# every way that read can fail has to be a recorded state rather than a guessed value.
@testset "the background pool's seed is read from its generator, or recorded as absent" begin
    recorded(; kw...) = (out = Pair{String,Any}[];
                         HE.background_pool_provenance!((k, v) -> push!(out, k => v); kw...);
                         out)
    as_dict(out) = Dict(out)
    mktempdir() do dir
        gen = joinpath(dir, "generator.jl")
        pool = joinpath(dir, "pool.csv")
        frame = joinpath(dir, "frame.csv")
        # A generator with each constant on a line of its own, among lines that must not
        # match: a comment, a longer name that starts with the same letters, and a constant
        # that is not a plain integer.
        write(gen, """
                   # const SEED = 1
                   const SEEDS     = 2
                   const MAX_TAXA  = 30   # a cap
                   const SUBSAMPLE = 4
                   const SEED      = 57
                   const OTHER     = 0x10
                   """)
        @test HE.pool_constant(gen, "SEED") == (value = 57, line = 5)
        @test HE.pool_constant(gen, "SUBSAMPLE") == (value = 4, line = 4)
        @test HE.pool_constant(gen, "MAX_TAXA") == (value = 30, line = 3)
        @test HE.pool_constant(gen, "OTHER") === nothing
        @test HE.pool_constant(gen, "ABSENT") === nothing
        @test HE.pool_constant(joinpath(dir, "no_such_generator.jl"), "SEED") === nothing
        # UNIQUENESS: two definitions of one constant are no provenance at all.
        twice = joinpath(dir, "twice.jl")
        write(twice, "const SEED = 5\nconst SUBSAMPLE = 4\nconst SEED = 6\n")
        @test HE.pool_constant(twice, "SEED") === nothing
        @test HE.pool_constant(twice, "SUBSAMPLE") == (value = 4, line = 2)

        # EVERYTHING PRESENT: fourteen keys, in the order a run writes them.
        write(pool, "taxid,phylum\n1,phy-1\n2,phy-1\n3,phy-2\n4,phy-2\n")
        write(frame, "taxid,habitat_type\n1,cat-a\n1,cat-a\n2,cat-b\n")
        full = recorded(generator = gen, pool = pool, frame = frame)
        @test first.(full) == ["background_pool_provenance_state", "background_pool_seed",
                               "background_pool_seed_origin", "background_pool_subsample_size",
                               "background_pool_subsample_above", "background_pool_generator",
                               "background_pool_generator_sha256", "background_pool_file",
                               "background_pool_sha256", "background_pool_taxa",
                               "background_frame_source_file", "background_frame_source_sha256",
                               "background_frame_source_taxa",
                               "background_pool_subsample_applied"]
        f = as_dict(full)
        @test f["background_pool_provenance_state"] == "generator_source"
        @test f["background_pool_seed"] == 57
        @test f["background_pool_seed_origin"] == "scripts/prokaryote_background.jl:5 (const SEED)"
        @test f["background_pool_subsample_size"] == 4
        @test f["background_pool_subsample_above"] == 30
        @test f["background_pool_generator_sha256"] == bytes2hex(sha256(read(gen)))
        @test f["background_pool_file"] == "pool.csv"
        @test f["background_pool_sha256"] == bytes2hex(sha256(read(pool)))
        @test f["background_pool_taxa"] == 4
        # the frame source counts distinct taxa, not rows
        @test f["background_frame_source_taxa"] == 2
        @test f["background_frame_source_sha256"] == bytes2hex(sha256(read(frame)))
        # four pool rows against SUBSAMPLE = 4: the subsample was taken
        @test f["background_pool_subsample_applied"] === true
        write(pool, "taxid,phylum\n1,phy-1\n2,phy-1\n3,phy-2\n")
        @test as_dict(recorded(generator = gen, pool = pool,
                               frame = frame))["background_pool_subsample_applied"] === false

        # GENERATOR ABSENT: no seed, no origin, no generator keys, and nothing the subsample
        # check could be made against.
        absent = recorded(generator = joinpath(dir, "no_such_generator.jl"), pool = pool,
                          frame = frame)
        a = as_dict(absent)
        @test a["background_pool_provenance_state"] == "generator_absent"
        @test ismissing(a["background_pool_seed"])
        @test a["background_pool_seed_origin"] == ""
        @test ismissing(a["background_pool_subsample_size"])
        @test !haskey(a, "background_pool_generator")
        @test !haskey(a, "background_pool_generator_sha256")
        @test ismissing(a["background_pool_subsample_applied"])
        # the tables are still described: their provenance does not depend on the generator's
        @test a["background_pool_sha256"] == bytes2hex(sha256(read(pool)))

        # SEED NOT FOUND: a generator that defines the seed twice has no seed to cite.
        nf = as_dict(recorded(generator = twice, pool = pool, frame = frame))
        @test nf["background_pool_provenance_state"] == "seed_not_found"
        @test ismissing(nf["background_pool_seed"])
        @test nf["background_pool_seed_origin"] == ""
        @test haskey(nf, "background_pool_generator_sha256")

        # SUBSAMPLE_APPLIED IS MISSING when either side of the check is unknown: no pool file,
        # or no single SUBSAMPLE line.
        nopool = as_dict(recorded(generator = gen, pool = joinpath(dir, "no_pool.csv"),
                                  frame = frame))
        @test ismissing(nopool["background_pool_subsample_applied"])
        @test !haskey(nopool, "background_pool_sha256")
        nosub = joinpath(dir, "nosub.jl")
        write(nosub, "const SEED = 57\n")
        ns = as_dict(recorded(generator = nosub, pool = pool, frame = frame))
        @test ns["background_pool_provenance_state"] == "generator_source"
        @test ismissing(ns["background_pool_subsample_size"])
        @test ismissing(ns["background_pool_subsample_applied"])

        # A table with no taxid column is described by digest and counts no taxa.
        write(frame, "id,habitat_type\n1,cat-a\n")
        @test ismissing(as_dict(recorded(generator = gen, pool = pool,
                                         frame = frame))["background_frame_source_taxa"])
    end

    # THIS REPOSITORY'S GENERATOR defines each constant exactly once, so a run on it cites the
    # seed by line; a second definition added there would turn every run's state into
    # seed_not_found. The tables it writes are the registry's.
    for name in ("SEED", "SUBSAMPLE", "MAX_TAXA")
        @test HE.pool_constant(HE.POOL_GEN, name) !== nothing
    end
    @test HE.POOL_FILE == HE.DataInputs.input_path("pool_taxa")
    @test HE.FRAME_FILE == HE.DataInputs.input_path("pool_background")
end

# ---------------------------------------------------------------- the driver, from outside
# The script copied into a tree of its own with the project files it activates, the three
# modules and the producer it includes, the pool generator whose seed it reads, a call table
# and one published scan set. Its REPO is that tree, so every default path it reads is a file
# of the fixture and nothing of this repository. Two runs over one table, differing only in
# the recall record's bar: the licence a run WRITES has to follow the bound, and has to come
# out strong when the bound allows it.
@testset "the driver writes the licence its inputs decide" begin
    repo = joinpath(@__DIR__, "..")
    tree = mktempdir()
    for f in ("Project.toml", "Manifest.toml",
              joinpath("scripts", "habitat_enrichment.jl"),
              joinpath("scripts", "background_pqq_scan.jl"),
              joinpath("scripts", "prokaryote_background.jl"),
              joinpath("src", "ModuleEnvironmentTools.jl"),
              joinpath("src", "ModuleHabitatStats.jl"), joinpath("src", "ModuleOutcomeRules.jl"),
              joinpath("src", "ModuleDataInputs.jl"))
        mkpath(dirname(joinpath(tree, f)))
        cp(joinpath(repo, f), joinpath(tree, f))
    end
    # The data root of this tree is dat/ under it, because LNFINDER_DATA is not set and the
    # registry falls back to the clone: the record below therefore goes where the LAYOUT puts
    # it, not where a repository-relative path would put it.
    mkpath(joinpath(tree, "dat", "habitat", "background_scan"))
    SR.write_keyvals(joinpath(tree, "dat", "habitat", "background_scan", "uniprot_release.txt"),
                     ["X-UniProt-Release" => FIXTURE_RELEASE])
    # The background pool's two tables, where the layout puts them, and the generator copied
    # above: the run describes all three in its summary and reads nothing else from them.
    pool_taxa = joinpath(tree, "dat", "habitat", "prokaryote_taxa.csv")
    write(pool_taxa, "taxid,phylum,domain,organism\n101,phy-1,Bacteria,Org 101\n" *
                     "102,phy-1,Bacteria,Org 102\n103,phy-1,Bacteria,Org 103\n")
    write(joinpath(tree, "dat", "habitat", "prokaryote_background.csv"),
          "taxid,habitat_type,habitat_chemistry\n101,cat-a,\n102,cat-a,\n")
    generator = joinpath(tree, "scripts", "prokaryote_background.jl")
    # Ln carriers and other census carriers by role, and the scan fixture's frame.
    ln = [(11, "cat-a"), (12, "cat-a"), (13, "cat-b"), (14, "cat-a"), (15, "cat-b"), (16, "cat-a")]
    other = [(21, "cat-b"), (22, "cat-a"), (23, "cat-b")]
    frame = [(101, "cat-a"), (102, "cat-a"), (103, "cat-b"), (104, "cat-a"),
             (105, "cat-b"), (106, "cat-b")]
    spec = vcat([(t, h, true, false, 1, 1) for (t, h) in ln],
                [(t, h, true, false, 1, 0) for (t, h) in other],
                [(t, h, false, true, 0, 0) for (t, h) in frame])
    n = length(spec)
    # kingdom is a CENSUS lineage column, so the frame taxa carry their domain instead — the
    # shape the move-in run below depends on, and the one an implementation that reads
    # kingdom alone turns into a level that moves nobody.
    calls = DataFrame(config = fill("primary", n), taxid = [s[1] for s in spec],
                      in_census = [s[3] for s in spec], in_pool = [s[4] for s in spec],
                      n_census_seqs = [s[5] for s in spec], n_ln_seqs = [s[6] for s in spec],
                      kingdom = [s[3] ? "Bacteria" : "" for s in spec],
                      domain = [s[3] ? "" : "Bacteria" for s in spec],
                      phylum = fill("phy-1", n),
                      genus = ["gen-$(s[1] % 4)" for s in spec],
                      organism = ["Org $(s[1])" for s in spec],
                      habitat_type = [s[2] for s in spec],
                      type_weights = [s[2] * ":1.0" for s in spec],
                      chemistry = fill("chem-x", n), chemistry_flags = fill("chem-x", n),
                      n_relations = fill(5, n), n_distinct = fill(3, n),
                      n_truncated = fill(0, n), n_tied = fill(1, n),
                      resolved_by = fill("none", n), n_src_record = fill(2, n),
                      n_src_text = fill(3, n), n_isolation_terms = fill(0, n),
                      n_generic_terms = fill(0, n), isolation_only = fill(false, n),
                      generic_only = fill(false, n))
    input = joinpath(tree, "input"); mkpath(input)
    CSV.write(joinpath(input, "calls.csv"), calls)
    # A census whose classes put ONE other-census carrier in the disclosed calcium arm. It is
    # here so that a run has rows of a disclosed contrast: the row counters a re-check reads
    # have to count the registered contrasts alone, and with no disclosed row the two forms
    # of every counter would be equal whatever counted them.
    census = DataFrame(id = ["seq-$(t)" for t in vcat([l[1] for l in ln], 21)],
                       classification = vcat(fill("Lanthanide", length(ln)), "Calcium"),
                       kingdom = fill("Bacteria", length(ln) + 1),
                       TaxID = vcat([l[1] for l in ln], 21))
    CSV.write(joinpath(input, "census.csv"), census)
    # The selection this driver run is held to. The script ships none, so the run has to be
    # handed one; these are the fixture's own abstract categories, and the chemistry role is
    # filled by a chemistry flag because the tier rule reads the chemistry unit.
    catsfile = selection_csv(input;
                             rows = [("expect", "selections", "5", "fixture"),
                                     ("select", "depletion_tier", "chem-x", "fixture"),
                                     ("select", "genus_clause", "cat-a", "fixture"),
                                     ("select", "sharing", "cat-b", "fixture"),
                                     ("select", "grid_count", "cat-a", "fixture"),
                                     ("select", "standardisation_phylum", "phy-1", "fixture")])
    # Both sets carry the recorded bar, because that is what the driver holds the bound to: no
    # CLI option can lower it. They differ only in how many taxa the full search flagged, which
    # is what the interval is made of — nine recovered of nine cannot show the bar is met, and
    # the smallest all-recovered sample that can needs 149 of them.
    below = publish_fixture(joinpath(input, "below"); recall = fixture_recall(bar = SR.MIN_RECALL))
    above = publish_fixture(joinpath(input, "above");
                            recall = fixture_recall(bar = SR.MIN_RECALL, n_flagged = 149,
                                                    recovered = 149, n_sampled = 200))
    @test read(below.scan) == read(above.scan)       # one table, two record sets
    # The pool for the third run: a census hit table naming two frame taxa lanthanide-type —
    # one the scan flagged (101) and one it did not (103) — and one calcium-class row that
    # must not be moved. The counts below are of THIS file, read by the run.
    poolfile = joinpath(input, "pool.csv")
    CSV.write(poolfile, DataFrame(id = ["s-101", "s-102", "s-103"],
                                  classification = ["Lanthanide", "Calcium", "Lanthanide"],
                                  TaxID = [101, 102, 103]))
    script = joinpath(tree, "scripts", "habitat_enrichment.jl")
    # The runs start together and are read back once all have ended: each pays the package
    # load on its own, and there is no reason to pay them one after the other.
    function start_driver(set, label; extra = String[])
        out = joinpath(tree, "out-$label"); mkpath(out)
        args = ["--categories", catsfile,
                "--calls", joinpath(input, "calls.csv"), "--scan", set.scan,
                "--threshold-record", set.threshold, "--recall", set.recall,
                "--acceptance", set.acceptance, "--census", joinpath(input, "census.csv"),
                "--boot", "0", "--out", joinpath(out, "rows.csv"),
                "--summary", joinpath(out, "summary.csv"), extra...]
        log = joinpath(out, "run.log")
        p = run(pipeline(ignorestatus(`$(Base.julia_cmd()) --startup-file=no $script $args`);
                         stdout = log, stderr = log); wait = false)
        return (process = p, out = out, log = log)
    end
    function finish_driver(started)
        wait(started.process)
        sumfile = joinpath(started.out, "summary.csv")
        summary = isfile(sumfile) ?
            Dict(r.key => (ismissing(r.value) ? "" : string(r.value))
                 for r in eachrow(CSV.read(sumfile, DataFrame; types = String))) :
            Dict{String,String}()
        return (code = started.process.exitcode, summary = summary,
                log = isfile(started.log) ? read(started.log, String) : "", out = started.out)
    end
    started = (start_driver(below, "below"), start_driver(above, "above"),
               start_driver(above, "movein";
                            extra = ["--bg-exclusion", HE.MOVEIN_MODE,
                                     "--movein-pool", poolfile]))
    b, a, c = finish_driver.(started)
    @test b.code == 0
    @test get(b.summary, "background_licence_code", "") == "fallback"
    @test get(b.summary, "background_recall_pass_derived", "") == "false"
    @test get(b.summary, "background_prefilter_recall_pass", "") == "true"   # the scan's own gate
    @test occursin("lower bound", get(b.summary, "background_scan_evidence_reasons", ""))
    @test parse(Float64, get(b.summary, "background_prefilter_recall_lower_bound_derived", "NaN")) ==
          HE.recall_lower_bound(9, 9)
    @test parse(Float64, get(b.summary, "background_prefilter_recall_lower_bound", "NaN")) ==
          HE.recall_lower_bound(9, 9)
    @test get(b.summary, "background_prefilter_recall_rate_derived", "") == "1.0"
    @test get(b.summary, "background_scan_sha256", "") == below.sha256
    @test get(b.summary, "background_scan_table_taxa_flagged", "") == "2"
    @test get(b.summary, "uniprot_release", "") == FIXTURE_RELEASE
    @test isfile(joinpath(b.out, "rows.csv")) && isfile(joinpath(b.out, "rows_anchor.csv"))
    # THE BACKGROUND POOL, AT THE DRIVER: the seed read out of the generator the run found in
    # its own tree, cited by line, and the pool's two tables found through the registry and
    # described by digest. The expected seed and line are read from the same file with the
    # same reader, so no number of the generator is typed here.
    seed = HE.pool_constant(generator, "SEED")
    @test get(b.summary, "background_pool_provenance_state", "") == "generator_source"
    @test get(b.summary, "background_pool_seed", "") == string(seed.value)
    @test get(b.summary, "background_pool_seed_origin", "") ==
          "scripts/prokaryote_background.jl:$(seed.line) (const SEED)"
    @test get(b.summary, "background_pool_generator_sha256", "") ==
          bytes2hex(sha256(read(generator)))
    @test get(b.summary, "background_pool_file", "") == "prokaryote_taxa.csv"
    @test get(b.summary, "background_pool_sha256", "") == bytes2hex(sha256(read(pool_taxa)))
    @test get(b.summary, "background_pool_taxa", "") == "3"
    @test get(b.summary, "background_frame_source_file", "") == "prokaryote_background.csv"
    @test get(b.summary, "background_frame_source_taxa", "") == "2"
    # three rows are not the generator's subsample size, and the run says so rather than
    # letting the seed vouch for a subsample that was not taken
    @test get(b.summary, "background_pool_subsample_applied", "") == "false"
    @test a.code == 0
    @test get(a.summary, "background_licence_code", "") == "strong"
    @test get(a.summary, "background_scan_evidence_ok", "") == "true"
    @test get(a.summary, "background_scan_evidence_reasons", "x") == ""
    @test get(a.summary, "background_removal_incomplete_reasons", "x") == ""
    @test get(a.summary, "background_carriers_left_hmm_hit", "") == "0"
    @test get(a.summary, "background_scan_frame_matches_analysis_frame", "") == "true"
    @test get(a.summary, "bg_removed_by_scan", "") == "2"
    # The disclosed arm's rows reach the plain row counters and never the registered ones.
    @test get(a.summary, "caonly_arm_taxa", "") == "1"
    @test parse(Int, get(a.summary, "rows_disclosed_contrasts", "0")) > 0
    @test parse(Int, get(a.summary, "rows_total", "0")) ==
          parse(Int, get(a.summary, "rows_registered_contrasts", "0")) +
          parse(Int, get(a.summary, "rows_disclosed_contrasts", "0"))
    @test parse(Int, get(a.summary, "rows_continuity_corrected", "0")) > 0
    @test get(a.summary, "rows_continuity_corrected_registered", "") == "0"
    @test parse(Int, get(a.summary, "rows_with_a_zero_arm_total", "0")) > 0
    @test get(a.summary, "rows_with_a_zero_arm_total_registered", "") == "0"
    @test parse(Int, get(a.summary, "rows_with_a_zero_stratum_cell_registered", "-1")) >= 0
    # The same two runs over the same table computed the same rows: the recall record changes
    # the licence and nothing else.
    @test read(joinpath(a.out, "rows.csv")) == read(joinpath(b.out, "rows.csv"))
    b.code == 0 && a.code == 0 || println(b.log, "\n", a.log)

    # THE THIRD RUN IS THE LEVEL THAT MOVES TAXA IN, over the same table and the same records
    # as the second, so every difference below is the level's doing and nothing else's.
    @test c.code == 0
    @test get(c.summary, "bg_exclusion_mode", "") == HE.MOVEIN_MODE
    @test get(c.summary, "bg_exclusion_mode_role", "") == "sensitivity"
    # The pool, as data: which rule ran, where it came from, and what it counted. None of
    # these numbers is in the script; each is a count of the file the run was handed.
    @test get(c.summary, "movein_pool_rule", "") == HE.REGISTERED_MOVEIN_POOL
    @test get(c.summary, "movein_pool_rule_role", "") == "registered"
    @test get(c.summary, "movein_pool_rule_source", "") == "explicit"
    @test get(c.summary, "movein_pool_source", "") == "pool.csv"
    @test get(c.summary, "movein_pool_sha256", "") == bytes2hex(sha256(read(poolfile)))
    @test get(c.summary, "movein_pool_rows_read", "") == "3"
    @test get(c.summary, "movein_pool_class_rows", "") == "2"      # the calcium row is not one
    @test get(c.summary, "movein_pool_taxa", "") == "2"
    @test get(c.summary, "movein_pool_taxa_in_frame", "") == "2"
    @test get(c.summary, "movein_pool_taxa_outside_frame", "") == "0"
    @test get(c.summary, "movein_moved_taxa", "") == "2"
    @test get(c.summary, "movein_moved_taxa_annotated", "") == "2"
    # The two rules disagree about which frame taxa are lanthanide-type, and the run says by
    # how much rather than leaving the reader to assume they are one set.
    @test get(c.summary, "movein_pool_taxa_under_the_scan_rule", "") == "2"
    @test get(c.summary, "movein_pool_taxa_in_both_rules", "") == "1"
    # n1 grows by exactly the annotated pool, and the arms say so from the table alone.
    @test get(c.summary, "n1_taxa", "") == "8" && get(a.summary, "n1_taxa", "") == "6"
    @test parse(Int, get(c.summary, "n1_taxa", "0")) ==
          parse(Int, get(c.summary, "movein_n1_before_the_move", "0")) +
          parse(Int, get(c.summary, "movein_moved_taxa_annotated", "0"))
    # n2 loses only the moved taxon the exclusion had left in; the other moved taxon was
    # already out of the background, and the run separates the two.
    @test get(a.summary, "n2_taxa", "") == "4" && get(c.summary, "n2_taxa", "") == "3"
    @test get(c.summary, "movein_background_taxa_the_move_removed", "") == "1"
    @test get(c.summary, "movein_moved_taxa_also_excluded_by_the_scan", "") == "1"
    @test get(c.summary, "movein_frame_size_before_the_move", "") == "6"
    @test get(c.summary, "frame_size", "") == "4"
    # THE PER-SEQUENCE UNIT, AT THE DRIVER. The moved taxa carry no census sequence, so the
    # sequence series is the one the unmoved run computed, and the run emits that rather than
    # leaving it to be inferred from an n1 that grew and a sequence count that did not.
    @test get(c.summary, "movein_moved_census_sequences", "") == "0"
    @test get(c.summary, "movein_sequence_units_from_moved_taxa", "") == "0"
    @test get(c.summary, "n_sequences_fg", "") == get(a.summary, "n_sequences_fg", "")
    # The moved taxa carry a domain and no kingdom, which is why the kingdom-keyed arm counts
    # above cannot see them; the run reports their domains instead of appearing to lose them.
    @test get(c.summary, "movein_moved_taxa_with_a_kingdom_label", "") == "0"
    @test get(c.summary, "movein_moved_taxa_by_domain", "") == "Bacteria:2"
    # THE OTHER MODES ARE UNTOUCHED BY THE NEW PATH: a run that moves nobody emits none of
    # these keys, tallies no move-in background size, and writes the rows it always wrote.
    for k in ("movein_pool_rule", "movein_moved_taxa", "movein_n1_before_the_move",
              "sensitivity_n2_$(HE.MOVEIN_MODE)")
        @test !haskey(a.summary, k) && !haskey(b.summary, k)
    end
    @test haskey(c.summary, "sensitivity_n2_hmmsearch")
    @test !haskey(c.summary, "sensitivity_n2_$(HE.MOVEIN_MODE)")
    @test read(joinpath(c.out, "rows.csv")) != read(joinpath(a.out, "rows.csv"))
    # The scan's coverage is still judged against the frame the scan was asked about, so the
    # level does not report its own definition as a hole in the table.
    @test get(c.summary, "background_scan_frame_matches_analysis_frame", "") == "true"
    @test get(c.summary, "background_carriers_left_hmm_hit_definitional", "") == "true"
    @test get(c.summary, "background_carriers_left_hmm_hit", "") == "0"
    c.code == 0 || println(c.log)
end

# ---------------------------------------------------------------- the command line, from outside
# The script copied into an EMPTY tree: its REPO holds nothing, so any stage past the command
# line would fail and the first write would create a directory. The copy is started as its
# own program, because the guard under test only runs then.
@testset "--help and a mistyped option exit before anything is loaded or written" begin
    tree = mktempdir()
    mkpath(joinpath(tree, "scripts"))
    mkpath(joinpath(tree, "src"))
    script = joinpath(tree, "scripts", "habitat_enrichment.jl")
    cp(joinpath(@__DIR__, "..", "scripts", "habitat_enrichment.jl"), script)
    # The registry the option defaults are named in. It is loaded AFTER --help and after the
    # option reader, so the two guards below still answer with nothing but the script; it is
    # here because the third case gets as far as the typed options, which resolve through it.
    cp(joinpath(@__DIR__, "..", "src", "ModuleDataInputs.jl"),
       joinpath(tree, "src", "ModuleDataInputs.jl"))
    function run_copy(args)
        out, err = IOBuffer(), IOBuffer()
        p = run(pipeline(ignorestatus(`$(Base.julia_cmd()) --startup-file=no $script $args`);
                         stdout = out, stderr = err))
        return (code = p.exitcode, out = String(take!(out)), err = String(take!(err)))
    end
    h = run_copy(["--help"])
    @test h.code == 0
    @test occursin("--bg-exclusion", h.out) && occursin("--anchor-reference", h.out)
    @test !occursin("Activating", h.out * h.err)            # no package was loaded
    u = run_copy(["--bg-exlusion", "lntype"])
    @test u.code == 1 && occursin("unknown option '--bg-exlusion'", u.err)
    v = run_copy(["--boot", "-1"])
    @test v.code == 1 && occursin("--boot", v.err)
    @test !occursin("Activating", u.err * v.err)
    @test sort(readdir(tree)) == ["scripts", "src"]         # nothing was written
    @test readdir(joinpath(tree, "scripts")) == ["habitat_enrichment.jl"]
    @test readdir(joinpath(tree, "src")) == ["ModuleDataInputs.jl"]
end

end
