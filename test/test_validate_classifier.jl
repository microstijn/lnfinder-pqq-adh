#=
    Offline end-to-end test for the control-set validation harness.

    No HMMER and no network are required: we hand-craft tiny AFA fixtures whose
    match-state columns are known, then drive the REAL production functions
    (calibrate_active_site, parse_hmmalign_results) plus the harness logic
    (diagnose_calibration, evaluate_controls). This proves the scoring + calibration
    diagnostics are correct before pointing them at real hmmalign output.

    Run:  julia --project=. test/test_validate_classifier.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
using DataFrames
using CSV

include(joinpath(@__DIR__, "..", "src", "validate_classifier.jl"))
using .ValidateClassifier
using .ValidateClassifier: calibrate_active_site, parse_hmmalign_results
using .ValidateClassifier: _inactive_uniparc_id, _inactive_reason, _rewrite_header,
                           _first_fasta_record, _source_map
using .ValidateClassifier: _hmm_length, _afa_match_length, _DEFAULT_HMM
using .ValidateClassifier: _search_tiers, _rank_candidates, _is_fragment
using .ValidateClassifier: _clopper_pearson, _LABELLED_COLUMNS, REFERENCE_ROLE, CONTESTED_TIER

# Helper: write records (id => raw aligned seq) to an AFA file.
function write_afa(path, recs)
    open(path, "w") do io
        for (id, seq) in recs
            write(io, ">", id, "\n", seq, "\n")
        end
    end
    return path
end

@testset "validate_classifier" begin

    tmp = mktempdir()

    # ---- Fixture A: unique catalytic column at model index 4 -------------------
    # Refs carry a lowercase insert ('q') after model pos 2 to exercise insert-stripping.
    # Model coords:  M D A D Y D E K L V  (positions 1..10)
    #   Ln ref column 4..6 = D Y D  -> D-x-D
    #   Ca ref column 4..6 = D Y A  -> D-x-[ATS]
    afaA = write_afa(joinpath(tmp, "A.afa"), [
        "LNREF" => "MDqADYDEKLV",   # Ln reference  (record 1)
        "CAREF" => "MDqADYAEKLV",   # Ca reference  (record 2)
        "LNT1"  => "GGGDFD",        # col4..6 D F D -> Lanthanide
        "LNT2"  => "AAADWD",        # col4..6 D W D -> Lanthanide
        "CAT1"  => "GGGDFT",        # col4..6 D F T -> Calcium
        "CAT2"  => "AAADWS",        # col4..6 D W S -> Calcium
    ])

    @testset "calibration picks the unique catalytic column" begin
        diag = diagnose_calibration(afaA)
        @test diag.collisions == [4]
        @test diag.chosen == 4
        @test calibrate_active_site(afaA) == 4
    end

    @testset "perfect control set scores 1.0" begin
        idx = calibrate_active_site(afaA)
        preds = parse_hmmalign_results(afaA, idx)
        # sanity: the calls we expect
        d = Dict(preds.accession .=> preds.classification)
        @test d["LNT1"] == :Lanthanide
        @test d["CAT2"] == :Calcium

        labels = DataFrame(
            accession = ["LNT1", "LNT2", "CAT1", "CAT2"],
            cofactor  = ["Lanthanide", "Lanthanide", "Calcium", "Calcium"],
        )
        ev = evaluate_controls(preds, labels)
        @test ev.sensitivity == 1.0
        @test ev.specificity == 1.0
        @test ev.accuracy == 1.0
    end

    @testset "a missed (Unknown) call lowers sensitivity and is reported" begin
        # Add a known-Ln protein whose column-6 residue is E -> classifier says :Unknown.
        afaA2 = write_afa(joinpath(tmp, "A2.afa"), [
            "LNREF" => "MDqADYDEKLV",
            "CAREF" => "MDqADYAEKLV",
            "LNT1"  => "GGGDFD",     # Lanthanide (correct)
            "UNK1"  => "GGGDFE",     # known Ln but col6=E -> :Unknown (a miss)
        ])
        idx = calibrate_active_site(afaA2)
        preds = parse_hmmalign_results(afaA2, idx)
        labels = DataFrame(accession = ["LNT1", "UNK1"],
                           cofactor  = ["Lanthanide", "Lanthanide"])
        ev = evaluate_controls(preds, labels)
        @test ev.tp == 1
        @test ev.fn == 1               # UNK1 counted against Ln recall
        @test ev.sensitivity == 0.5
    end

    # ---- Fixture B: AMBIGUOUS calibration (decoy + true collision) -------------
    # Ln ref:  D A D D Y D   -> D-x-D at col 1 (decoy) AND col 4 (true)
    # Ca ref:  D A S D Y A   -> D-x-[S] at col 1     AND D-x-[A] at col 4
    # Production calibrate_active_site returns the FIRST (col 1 = the decoy).
    @testset "ambiguous calibration is detected (the first-match flaw)" begin
        afaB = write_afa(joinpath(tmp, "B.afa"), [
            "LNREF" => "DADDYD",
            "CAREF" => "DASDYA",
        ])
        diag = diagnose_calibration(afaB)
        @test diag.collisions == [1, 4]      # two legitimate columns exist
        @test diag.chosen == 1               # production silently takes the decoy
        @test calibrate_active_site(afaB) == 1
    end

    # ---- Fetch hardening: a control row must never disappear quietly -------------
    #= The failure these suites exist to prevent: a UniProtKB entry that has been
       withdrawn still answers its .fasta URL with HTTP 200 and a ZERO-BYTE body. A
       harness that warns and skips that accession drops the labelled protein out
       of the inner join and out of every denominator, with nothing failing. Three
       things have to hold for that to be impossible: the UniParc id is readable from
       the entry JSON, the archived record is re-keyed to the control-set accession so
       it can still join to its label, and a label with no prediction stops the run. =#

    @testset "an inactive UniProtKB entry exposes its UniParc id" begin
        inactive = """{"entryType":"Inactive","primaryAccession":"A0ACD6B9F2","inactiveReason":{"inactiveReasonType":"DELETED"},"extraAttributes":{"uniParcId":"UPI0018C0EF6D"}}"""
        @test _inactive_uniparc_id(inactive) == "UPI0018C0EF6D"
        @test _inactive_reason(inactive) == "DELETED"

        active = """{"entryType":"UniProtKBReviewed","primaryAccession":"P16027","sequence":{"value":"MDA"}}"""
        @test _inactive_uniparc_id(active) === nothing
        @test _inactive_uniparc_id("") === nothing              # the zero-byte body itself
        @test _inactive_uniparc_id("not json at all") === nothing
    end

    @testset "an archived record is re-keyed to the control accession" begin
        hdr = _rewrite_header("A0ACD6B9F2", "UPI0018C0EF6D", "uniparc (UniProtKB deleted: DELETED)")
        @test hdr == ">up|A0ACD6B9F2|UPI0018C0EF6D uniparc (UniProtKB deleted: DELETED)"

        afaC = write_afa(joinpath(tmp, "C.afa"), [
            "LNREF" => "MDqADYDEKLV",
            "CAREF" => "MDqADYAEKLV",
            lstrip(hdr, '>') => "GGGDFD",       # the rescued record, as hmmalign sees it
        ])
        idx = calibrate_active_site(afaC)
        preds = parse_hmmalign_results(afaC, idx)
        # Keyed by the accession, NOT by the UPI — otherwise it could never join.
        @test "A0ACD6B9F2" in preds.accession
        @test !any(startswith.(preds.accession, "UPI"))
        @test Dict(preds.accession .=> preds.classification)["A0ACD6B9F2"] == :Lanthanide

        labels = DataFrame(accession = ["A0ACD6B9F2"], cofactor = ["Lanthanide"])
        ev = evaluate_controls(preds, labels)
        @test nrow(ev.table) == 1               # the deleted entry is still scored
        @test ev.missing == String[]
    end

    @testset "a label with no prediction cannot shrink n in silence" begin
        idx = calibrate_active_site(afaA)
        preds = parse_hmmalign_results(afaA, idx)
        labels = DataFrame(
            accession = ["LNT1", "CAT1", "MISSING1"],
            cofactor  = ["Lanthanide", "Calcium", "Calcium"],
        )
        @test_throws ArgumentError evaluate_controls(preds, labels)          # strict by default
        ev = evaluate_controls(preds, labels; strict = false)
        @test ev.missing == ["MISSING1"]
        @test nrow(ev.table) == 2
        @test evaluate_controls(preds, labels[1:2, :]).missing == String[]   # full coverage, no throw
    end

    @testset "fetched bodies are parsed and routed the way the manifest assumes" begin
        header, seq = _first_fasta_record(">UPI0018C0EF6D status=active\nMDAD\nYDEK\n")
        @test header == "UPI0018C0EF6D status=active"
        @test seq == "MDADYDEK"                      # unwrapped: the manifest length is exact
        @test _first_fasta_record("") == ("", "")     # a body with no record at all

        # A `source` column (or kwarg) routes single rows away from UniProtKB.
        @test _source_map(["P16027", "WP_1"], ["uniprot", "NCBI"])["WP_1"] == "ncbi"
        @test get(_source_map(["P16027"], nothing), "P16027", "uniprot") == "uniprot"
    end

    # ---- Accession resolver: which entry is proposed for (gene, organism) --------
    #= The failure these suites exist to prevent: the resolver proposes a short, unreviewed
       FRAGMENT where a curated full-length entry for the same protein and organism exists.
       Two independent defects are enough to produce that. A search that stops at the first
       tier returning anything, whose first tier queries the GENE symbol, never lets a
       Swiss-Prot entry with an empty `gene_names` field into the candidate list at all, so
       it cannot be ranked however good it is. And a ranking that knows only `reviewed` and
       the organism name cannot tell a short fragment from a complete sequence: to it they
       are the same kind of answer. The damage is silent and downstream: a fragment need not
       span the catalytic column, so the label of a control row can end up attached to a
       sequence that cannot be read at the position the whole classifier turns on. Two
       things must hold. A reviewed protein-name tier must be tried FIRST, so a gene-less
       curated entry is among the candidates; and the ranking must penalise fragments and
       short entries hard enough that the curated entry wins even when both match the
       organism. =#
    @testset "a gene-less curated entry can still enter the candidate list" begin
        tiers = _search_tiers("mxaF", "Methylophilus methylotrophus W3A1",
                              "methanol dehydrogenase")
        @test first.(tiers) == ["name_reviewed", "gene", "name_reviewed+genus",
                                "gene+genus", "name", "name+genus"]

        # The first query asks Swiss-Prot by protein name: the one question that does not
        # go through `gene_names`, and so the only one a gene-less entry can answer.
        q1 = last(tiers[1])
        @test occursin("reviewed:true", q1)
        @test occursin("protein_name:\"methanol dehydrogenase\"", q1)
        @test occursin("organism_name:\"Methylophilus methylotrophus W3A1\"", q1)
        @test !occursin("gene:", q1)

        # UniProt keeps the strain out of the organism name, so a control row that names
        # one matches only the genus/species form. The reviewed tier needs both, or it
        # never fires for the rows it was added for — and it must still come before the
        # gene question at the same organism specificity.
        q3 = last(tiers[3])
        @test occursin("reviewed:true", q3)
        @test occursin("organism_name:\"Methylophilus methylotrophus\"", q3)
        @test !occursin("W3A1", q3)
        @test occursin("W3A1", last(tiers[2]))             # the gene tier keeps the strain

        # With no protein-name keyword the gene tiers are all there is.
        @test first.(_search_tiers("mxaF", "Methylophilus methylotrophus W3A1", nothing)) ==
              ["gene", "gene+genus"]
    end

    @testset "a reviewed full-length entry outranks a short fragment" begin
        # The three entries one organism really offers for this gene: an unreviewed
        # fragment of the catalytic subunit, the curated full-length catalytic subunit with
        # NO gene name, and the curated small subunit. All three name the alpha or the beta
        # chain of the same enzyme in the same organism, so the organism bonus cannot
        # separate them and only the subunit bias, the fragment flag and the length can.
        hits = DataFrame(
            accession    = ["B7SBF8", "P38539", "P38540"],
            protein_name = ["Methanol dehydrogenase alpha subunit",
                            "Methanol dehydrogenase [cytochrome c] subunit 1 (MDH large subunit alpha)",
                            "Methanol dehydrogenase [cytochrome c] subunit 2 (MDH small subunit beta)"],
            gene_names   = ["mxaF", "", "moxI"],
            organism     = fill("Methylophilus methylotrophus (Bacterium W3A1)", 3),
            reviewed     = ["unreviewed", "reviewed", "reviewed"],
            length       = [169, 573, 91],
            fragment     = ["fragment", "", ""],
        )
        ranked = _rank_candidates(hits, "Methylophilus methylotrophus W3A1")
        @test ranked.accession[1] == "P38539"
        @test ranked.accession == ["P38539", "P38540", "B7SBF8"]
        @test issorted(ranked.score; rev = true)
        @test hits.accession == ["B7SBF8", "P38539", "P38540"]   # pure: input untouched
        @test !("score" in names(hits))

        # Each penalty stands on its own: neither is carrying the result alone.
        one(; kw...) = only(_rank_candidates(DataFrame(; accession = ["Q"],
                                                       protein_name = ["Methanol dehydrogenase"],
                                                       organism = ["Methylophilus methylotrophus"],
                                                       reviewed = ["unreviewed"], kw...),
                                             "Methylophilus methylotrophus").score)
        intact = one(length = [600], fragment = [""])
        @test one(length = [600], fragment = ["fragment"]) == intact - 5
        @test one(length = [169], fragment = [""])         == intact - 2
        @test one(length = [400], fragment = [""])         == intact   # the bar is exclusive

        # An unknown length is not a short one, and a table without the two columns at all
        # still ranks — a minimal fixture must not throw.
        @test one(length = [missing], fragment = [""]) == intact
        @test only(_rank_candidates(DataFrame(accession = ["Q"], reviewed = ["reviewed"]),
                                    "Any organism").score) == 2

        @test _is_fragment("fragment")
        @test !_is_fragment("")
        @test !_is_fragment(missing)
    end

    # ---- Coordinate guard: an index is only meaningful in model coordinates ------
    #= The failure these suites exist to prevent: the active-site column is an INDEX
       into the match states of the profile. Hand the scorer an alignment built against
       a different model, or one whose columns were re-derived after the fact, and every
       index still resolves — at a different residue in every sequence. No error, no
       change in n, and every call silently wrong. Two things must therefore hold: the
       model's own match-state count is readable from its header, and the alignment's
       match-column count is asserted against it before any column is read. =#

    @testset "the model's match-state count is read from its header" begin
        hmm3 = joinpath(tmp, "fixture.hmm")
        write(hmm3, "HMMER3/f\nNAME  x\nLENG  589\n")
        @test _hmm_length(hmm3) == 589

        # A header without LENG must raise: a silent 0 would disable the guard.
        nolen = joinpath(tmp, "noleng.hmm")
        write(nolen, "HMMER3/f\nNAME  x\nALPH  amino\n")
        @test_throws ErrorException _hmm_length(nolen)

        # Anything past the start of the emission table is model body, not header.
        late = joinpath(tmp, "late.hmm")
        write(late, "HMMER3/f\nNAME  x\nHMM          A        C\nLENG  10\n")
        @test_throws ErrorException _hmm_length(late)

        # The profile the harness actually aligns against parses the same way.
        if isfile(_DEFAULT_HMM)
            @test _hmm_length(_DEFAULT_HMM) == 589
        end
    end

    @testset "an alignment is scored only in the model's coordinates" begin
        # Fixture A's references span 10 match states ('q' is an insert, so stripped).
        @test _afa_match_length(afaA) == 10

        ctsv = joinpath(tmp, "controls.tsv")
        open(ctsv, "w") do io
            println(io, "# offline fixture — no network, no HMMER")
            println(io, "accession\tcofactor")
            for (a, c) in (("LNT1", "Lanthanide"), ("LNT2", "Lanthanide"),
                           ("CAT1", "Calcium"), ("CAT2", "Calcium"))
                println(io, a, "\t", c)
            end
        end
        score(dir; kw...) = validate_from_afa(; control_afa = afaA, control_tsv = ctsv,
                                               ln_ref_acc = "LNREF", ca_ref_acc = "CAREF",
                                               work_dir = joinpath(tmp, dir), kw...)

        res = score("guard_ok"; expected_len = 10)
        @test res.calibrated_index == 4
        @test res.evaluation.accuracy == 1.0

        # One column short of the model: the index would land on the wrong residue.
        @test_throws ErrorException score("guard_short"; expected_len = 9)
        @test_throws ErrorException score("guard_long"; expected_len = 11)

        # The count may also be taken straight from a profile header.
        hmm10 = joinpath(tmp, "len10.hmm")
        write(hmm10, "HMMER3/f\nNAME  fixture\nLENG  10\n")
        @test score("guard_hmm"; hmm = hmm10).evaluation.accuracy == 1.0

        # Given neither, the harness falls back to its own profile and warns — an
        # alignment in foreign coordinates is refused rather than scored.
        if isfile(_DEFAULT_HMM)
            @test_throws ErrorException score("guard_default")
        end
    end

    # ---- Tiered, stratified scoring ---------------------------------------------
    #= The failure these suites exist to prevent: the control set is not one population.
       Two of its rows ARE the sequences the active-site column was read off, one row's
       published cofactor is disputed, and the rest differ in how directly the cofactor was
       established. Pooled into a single accuracy, the calibration rows grade the rule on
       the sequences that set the rule — circular — and a disputed label moves a digit it
       cannot support. Three things have to hold. Held out means exactly role != reference,
       so a reference row cannot enter a held-out stratum through any tier it happens to
       carry. The contested row is listed with the call it received and never reaches a
       numerator or a denominator. And every rate travels with an exact (Clopper–Pearson)
       interval, because an all-correct stratum of a handful of rows has a lower bound far
       below 1, which the zero-width normal approximation at k == n hides completely. =#

    # accession, cofactor, role, tier, domain, organism, aligned sequence.
    # LNREF is the trap: role reference AND tier A, so it is eligible for heldout_A on
    # every test except the role test. A3 and C2 carry the same sequence, so the
    # redundancy count has something to find.
    tiered = [
        ("LNREF", "Lanthanide", "reference", "A",   "Bacteria",  "Ref organism L", "MDqADYDEKLV"),
        ("CAREF", "Calcium",    "reference", "ref", "Bacteria",  "Ref organism C", "MDqADYAEKLV"),
        ("A1",    "Lanthanide", "control",   "A",   "Bacteria",  "Organism A1",    "GGGDFD"),
        ("A2",    "Lanthanide", "control",   "A",   "Bacteria",  "Organism A2",    "AAADWD"),
        ("A3",    "Lanthanide", "control",   "A",   "Bacteria",  "Organism A3",    "GGGDYD"),
        ("A4",    "Lanthanide", "control",   "A",   "Bacteria",  "Organism A4",    "CCCDMD"),
        ("A5",    "Calcium",    "control",   "A",   "Bacteria",  "Organism A5",    "GGGDFA"),
        ("A6",    "Calcium",    "control",   "A",   "Bacteria",  "Organism A6",    "AAADWS"),
        ("A7",    "Calcium",    "control",   "A",   "Bacteria",  "Organism A7",    "GGGDYT"),
        ("A8",    "Calcium",    "control",   "A",   "Bacteria",  "Organism A8",    "CCCDMA"),
        ("B1",    "Lanthanide", "control",   "B",   "Bacteria",  "Organism B1",    "TTTDFD"),
        ("B2",    "Lanthanide", "control",   "B",   "Bacteria",  "Organism B2",    "VVVDWD"),
        ("B3",    "Calcium",    "control",   "B",   "Bacteria",  "Organism B3",    "TTTDFT"),
        ("B4",    "Calcium",    "control",   "B",   "Bacteria",  "Organism B4",    "VVVDWS"),
        ("C1",    "Lanthanide", "control",   "C",   "Bacteria",  "Organism C1",    "EEEDFD"),
        ("C2",    "Lanthanide", "control",   "C",   "Bacteria",  "Organism C2",    "GGGDYD"),
        ("C3",    "Calcium",    "control",   "C",   "Bacteria",  "Organism C3",    "EEEDFA"),
        ("C4",    "Calcium",    "control",   "C",   "Eukaryota", "Organism C4",    "KKKDWA"),
        # Contested: labelled Ln, called Ca. Any stratum that counts it drops below 1.0.
        ("X1",    "Lanthanide", "control",   "X",   "Bacteria",  "Organism X1",    "GGGDFS"),
    ]

    afaT = write_afa(joinpath(tmp, "T.afa"), [r[1] => r[7] for r in tiered])
    tsvT = joinpath(tmp, "controls_tiered.tsv")
    open(tsvT, "w") do io
        println(io, "# offline fixture — no network, no HMMER")
        println(io, "accession\tcofactor\trole\ttier\tdomain\torganism")
        for r in tiered
            println(io, join(r[1:6], '\t'))
        end
    end

    workT = joinpath(tmp, "tiered")
    resT = validate_from_afa(; control_afa = afaT, control_tsv = tsvT,
                               ln_ref_acc = "LNREF", ca_ref_acc = "CAREF",
                               work_dir = workT, expected_len = 10)
    S = resT.summary
    function stratum(df, name)
        hits = df[df.stratum .== name, :]
        @test nrow(hits) == 1
        return first(eachrow(hits))
    end

    @testset "held out means role != reference, whatever tier the row carries" begin
        @test is_heldout("control")
        @test is_heldout("")                       # no role recorded: nothing was demoted
        @test !is_heldout(REFERENCE_ROLE)
        @test !is_heldout(" Reference ")           # case and padding cannot smuggle one in

        # LNREF is a reference row carrying tier A. Tier A holds the eight control rows.
        @test stratum(S, "heldout_A").n == 8
        @test stratum(S, "heldout_A_Ln").n == 4
        @test stratum(S, "heldout_A_Ca").n == 4
        @test stratum(S, "reference").n == 2
        @test !("LNREF" in resT.evaluation_heldout.table.accession)
        @test !("CAREF" in resT.evaluation_heldout.table.accession)
        @test nrow(resT.evaluation_heldout.table) == 16

        # ...and it cannot sneak in even when it is the row that would change the answer.
        poisoned = DataFrame(
            accession   = ["R1", "H1", "H2"],
            cofactor    = ["Lanthanide", "Lanthanide", "Calcium"],
            role        = ["Reference", "control", "control"],
            tier        = ["A", "A", "A"],
            active_site = ["DYA", "DYD", "DYA"],
            pred        = [:Calcium, :Lanthanide, :Calcium],
            known       = [:Lanthanide, :Lanthanide, :Calcium],
        )
        poisoned.correct = poisoned.pred .== poisoned.known
        P = stratified_summary(poisoned; quiet = true)
        @test stratum(P, "heldout_A").n == 2
        @test stratum(P, "heldout_A").accuracy == 1.0     # the wrong reference row is out
        @test stratum(P, "heldout_all").n == 2
        @test stratum(P, "reference").n == 1
        @test stratum(P, "reference").accuracy == 0.0     # it is reported, just apart
    end

    @testset "a contested row is listed and never counted" begin
        x = stratum(S, "tier_$CONTESTED_TIER")
        @test x.n == 1
        @test x.correct === missing
        @test x.accuracy === missing                      # no numerator, no denominator
        @test occursin("Calcium", x.note)                 # listed with the call it got
        @test "X1" in resT.labelled.accession             # present in the labelled table

        # X1's call contradicts its label, so any accuracy it entered would be < 1.
        @test resT.evaluation_all.accuracy == 1.0
        @test resT.evaluation_heldout.accuracy == 1.0
        @test !("X1" in resT.evaluation_all.table.accession)
        @test !("X1" in resT.evaluation_heldout.table.accession)
        @test all(r -> r.accuracy === missing || r.accuracy == 1.0, eachrow(S))

        # The references are scored in `all` but not in `heldout`; both exclude X1.
        @test nrow(resT.evaluation_all.table) == 18
        @test stratum(S, "heldout_all").n == 16
    end

    @testset "every stratum carries an exact binomial interval" begin
        # Anchors: Clopper–Pearson lower bounds at 95% for k == n.
        @test isapprox(_clopper_pearson(8, 8)[1],   0.631; atol = 0.001)
        @test isapprox(_clopper_pearson(16, 16)[1], 0.794; atol = 0.001)
        @test _clopper_pearson(8, 8)[2] == 1.0
        @test _clopper_pearson(0, 0) === (missing, missing)

        @test isapprox(stratum(S, "heldout_A").ci_lo,   0.631; atol = 0.001)
        @test isapprox(stratum(S, "heldout_all").ci_lo, 0.794; atol = 0.001)
        @test stratum(S, "heldout_all").accuracy == 1.0
        @test occursin("mixed evidence", stratum(S, "heldout_all").note)
        @test all(r -> r.ci_lo === missing || r.ci_lo <= r.accuracy <= r.ci_hi, eachrow(S))
    end

    @testset "the set's composition is printed, zeros included" begin
        dom = domain_counts(resT.labelled)
        @test dom[dom.domain .== "Bacteria", :n] == [18]
        @test dom[dom.domain .== "Archaea", :n] == [0]        # stated, not left out
        @test dom[dom.domain .== "Eukaryota", :n] == [1]

        sw = switch_residue_counts(resT.labelled)
        @test Dict(sw.residue .=> sw.n) == Dict("A" => 5, "S" => 2, "T" => 2)

        printed = joinpath(tmp, "summary_stdout.txt")
        open(printed, "w") do io
            redirect_stdout(io) do
                stratified_summary(resT.labelled)
            end
        end
        out = read(printed, String)
        @test occursin(r"Archaea\s+0", out)
        @test occursin("tier_$CONTESTED_TIER", out)
    end

    @testset "near-identical rows are counted, so n_eff < n is visible" begin
        @test resT.redundancy.n_eff == 15                  # A3 and C2 are one draw, not two
        @test ["A3", "C2"] in resT.redundancy.clusters
        @test sum(length, resT.redundancy.clusters) == 16
        # Raise the bar and nothing changes; drop it and the whole set collapses.
        @test identity_clusters(afaT, ["A1", "A2"]).n_eff == 2
        @test identity_clusters(afaT, ["A1", "A2"]; threshold = 0.0).n_eff == 1
    end

    @testset "labelled rows and strata are written for a reader to check" begin
        labelled_csv = joinpath(workT, "control_predictions_labelled.csv")
        summary_csv  = joinpath(workT, "control_summary.csv")
        @test isfile(labelled_csv) && isfile(summary_csv)

        L = CSV.read(labelled_csv, DataFrame)
        @test names(L) == collect(_LABELLED_COLUMNS)
        @test nrow(L) == 19                                # every labelled row, X1 included
        @test L[L.accession .== "X1", :correct] == [false]

        W = CSV.read(summary_csv, DataFrame)
        @test nrow(W) == nrow(S)
        @test "ci_lo" in names(W) && "ci_hi" in names(W)
        @test isapprox(only(W[W.stratum .== "heldout_all", :ci_lo]), 0.794; atol = 0.001)
    end

    @testset "include_tiers narrows the scored set without touching calibration" begin
        resB = validate_from_afa(; control_afa = afaT, control_tsv = tsvT,
                                   ln_ref_acc = "LNREF", ca_ref_acc = "CAREF",
                                   work_dir = joinpath(tmp, "tiered_B"),
                                   include_tiers = ["B"], expected_len = 10)
        @test resB.calibrated_index == 4                   # the refs still set the column
        @test stratum(resB.summary, "heldout_all").n == 4
        @test stratum(resB.summary, "reference").n == 0
        @test stratum(resB.summary, "heldout_B").n == 4
        @test nrow(resB.summary[resB.summary.stratum .== "heldout_A", :]) == 0
        @test resB.evaluation_heldout.accuracy == 1.0

        # exclude_roles is the older switch; it still drops the references, and the two
        # evaluations then describe the same rows.
        resE = validate_from_afa(; control_afa = afaT, control_tsv = tsvT,
                                   ln_ref_acc = "LNREF", ca_ref_acc = "CAREF",
                                   work_dir = joinpath(tmp, "tiered_E"),
                                   exclude_roles = ["reference"], expected_len = 10)
        @test nrow(resE.evaluation_all.table) == 16
        @test resE.evaluation_heldout === resE.evaluation_all
        @test stratum(resE.summary, "reference").n == 0
    end
end
