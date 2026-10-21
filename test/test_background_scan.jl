#=
    Offline self-test for the background PQQ-ADH scan. No network and no HMMER: the
    script's pure half is included behind its own `PROGRAM_FILE` guard, so including this
    file defines the helpers without downloading anything or shelling out to WSL. One
    testset starts a copy of the script as a subprocess with a command line only, because
    what happens before the guard, when the file is the program, can only be measured
    from outside.

    Fixtures use arbitrary sizes and arbitrary values. None of them is a quantity measured
    on the data this script scans: a suite that encodes a run's own numbers passes only
    until the next run, and it publishes what the run found. The one exception is the
    --tblout table below, and it is not one in substance: its three targets are public
    Swiss-Prot entries and its scores are what the profile HMM under dat/models/ gives for
    them, so every number in it is recomputable from this repository plus three public
    sequences. It is copied rather than invented for the reason given under 1.

    The failures this suite exists to prevent:

    1. FIXED-WIDTH PARSING OF --tblout. hmmsearch's table is padded to look like columns,
       but the padding grows with the longest target name in the run. Slicing it by
       character offset works on one run and silently shifts a field on the next, which
       turns a bit score into an E-value. The fixture below is therefore copied verbatim
       from a real run of the census HMM, long target name and all — a hand-built table
       would have been padded evenly and would have agreed with a broken parser.

    2. THE TARGET NAME IS NOT AN ACCESSION. UniProt FASTA names a record `db|ACC|ENTRY`.
       Joining the raw target name against a table keyed by plain accession matches
       nothing at all, and an empty join here reads as "no background taxon carries the
       enzyme" — exactly the false negative this whole step exists to remove.

    3. AN EXCLUSION THAT DROPS THE WRONG SIDE. The scan's output is a per-taxon table over
       the WHOLE frame, zero rows included. An exclusion that keeps `n_hits > 0` instead
       of dropping it, or that inner-joins and thereby discards every taxon the scan never
       saw, empties or inverts the comparator arm without erroring.

    4. A THRESHOLD THAT DEFINES ITS OWN ANSWER. The recovery denominator must be every
       census accession present in the family, not just the ones that scored: counting
       only the scored ones makes an unrecovered accession disappear instead of failing.

    5. A MEMBERSHIP COLUMN DERIVED FROM HITS. Read off the taxon's hit rows,
       `in_ipr017512` is false for every zero-hit taxon whatever UniProt says about the
       organism, and an acceptance check built on it cannot fail — while the one carrier
       class that needs no HMM to find, a taxon InterPro already files in the
       methanol/ethanol dehydrogenase family, is exactly the one it is meant to catch.
       The roll-up is therefore tested with a taxon that has a family entry and NO hits.

    6. A RECALL THAT PASSES ON NOTHING. If the full-proteome arm flags no taxon at all,
       there is nothing to recover and the rate is undefined. Read as 1.0 it certifies the
       prefilter on an empty measurement, which is what a failed proteome download looks
       like. For the same reason a perfect rate is not a bound: the exact one-sided lower
       bound is computed from the same two numbers and is tested against its own defining
       equation and against a second, independent route to it, so a small denominator
       cannot be read as a tight result. The bound is recorded in recall.txt at full
       precision, because a consumer decides the stronger background statement on it; the
       scan's own gate stays on the rate, and the suite pins both halves: a perfect rate
       whose bound is below the bar passes the gate, and the bound row never gains a bar.

    7. AN ACCEPTANCE DENOMINATOR DRAWN FROM THE SEARCH SPACE. Counting "family taxa
       flagged" over the prefilter's own candidate rows can only ever count taxa the
       prefilter reached, so it reads 100 % precisely when the search space has a hole.
       The check has to run over the whole-frame roll-up, and the suite feeds it an
       unflagged family taxon to prove the MISS appears.

    8. THE SAME RECORD TWICE. The search space is a union of two overlapping downloads.
       A subset FASTA that emits a shared record once per page set doubles that sequence's
       taxon hit count.

    9. CACHES AND DOWNLOADS TRUSTED BY NAME. A hmmsearch table says in its own trailer
       what -T produced it, and a download can be compared with the total the server
       announced. Read neither, and a table from another cutoff and a page lost to a dropped
       connection both pass as current. The concatenated recall proteomes are the same hazard
       under a name: trusted by name, a changed sample scores the recall on the previous
       sample's sequences while the record still gives the requested size.

   10. A FAILING RUN THAT LEAVES A CONSUMABLE TABLE. Written in run order, the per-taxon
       table and its three records leave a run that MISSes an acceptance check having
       already overwritten the cutoff and the recall record, with the previous run's table
       still in place and nothing on disk saying the four do not belong together. The table
       and the records are published together, each record carrying the sha256 of the
       table it describes, or not at all.

   11. AN ACCEPTANCE VERDICT THAT DISAGREES WITH ITS OWN NUMBERS, and a bar remembered
       from an earlier run. Each gated row carries value, denominator and bar, and its
       verdict is derived from them, so a consumer can re-derive it; a row with nothing
       measured (an empty denominator, a NaN value) is a MISS rather than a pass.

   12. A RECORD THAT CANNOT BE READ BACK. A record is only a guard if a consumer can parse
       it; a reader that fills an absent `pass` with a default turns "not recorded" into
       "checked". Every record is therefore round-tripped, and missing or malformed fields
       come back as `nothing`.

   13. A RECALL SAMPLE REDRAWN INSTEAD OF EXTENDED. The sample is the first n draws of one
       seeded permutation, so enlarging it must keep every taxon the smaller sample had.
       A fresh subset per size would silently replace the taxa the earlier recall measured.
       And a failing run that displaced the earlier table must not leave the next run with
       nothing to compare its sample against.

   14. AN OPTION THAT IS IGNORED. An unknown option silently skipped lets a misspelt flag
       run the default scan while every record describes a run nobody asked for, and a
       `--help` read after the work starts begins the downloads.

    Run:  julia --project=. test/test_background_scan.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Test, DataFrames, CSV, SHA

include(joinpath(@__DIR__, "..", "scripts", "background_pqq_scan.jl"))

# `hmmsearch --tblout` (HMMER 3.4): two banner lines, one ruler line, then one row per
# target, then a trailing comment block. The first target name is four characters longer
# than the others and pushes every later field right — that is the point of the fixture.
const TBLOUT_FIXTURE = """
#                                                                      --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name            accession  query name              accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#    ------------------- ----------    -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
sp|A0A3F2YLY8|XOXF_METBY -          pqq_xoxF_mxaF_alignment -           1.1e-296  974.7  24.1  1.4e-296  974.4  24.1   1.0   1   0   0   1   1   1   1 Lanthanide-dependent methanol dehydrogenase XoxF OS=Methylotuvimicrobium buryatense OX=95641 GN=xoxF PE=1 SV=1
sp|A8R3S4|QEDH_PSEPU     -          pqq_xoxF_mxaF_alignment -           7.3e-143  466.8  10.4  9.8e-143  466.3  10.4   1.0   1   0   0   1   1   1   1 Quinoprotein ethanol dehydrogenase OS=Pseudomonas putida OX=303 GN=qedA PE=1 SV=1
sp|P54741|AFSK_STRCO     -          pqq_xoxF_mxaF_alignment -            1.5e-08   23.2   4.9   0.00018    9.8   0.4   4.3   1   1   2   4   4   4   2 Serine/threonine-protein kinase AfsK OS=Streptomyces coelicolor (strain ATCC BAA-471 / A3(2) / M145) OX=100226 GN=afsK PE=1 SV=2
#
# Program:         hmmsearch
# Version:         3.4 (Aug 2023)
# [ok]
"""

@testset "background_pqq_scan" begin

    @testset "tblout parsing" begin
        path = joinpath(mktempdir(), "fixture.tblout")
        write(path, TBLOUT_FIXTURE)
        df = parse_tblout(path)

        # Three data rows; the banner, ruler and trailer are comments.
        @test nrow(df) == 3
        @test names(df) == ["target", "accession", "full_evalue", "full_bits",
                            "best_dom_evalue", "best_dom_bits"]

        # Whitespace splitting, not character offsets: the first row's target name is
        # four characters longer than the others and its fields must still line up.
        @test df.accession == ["A0A3F2YLY8", "A8R3S4", "P54741"]
        @test df.full_bits ≈ [974.7, 466.8, 23.2]
        @test df.full_evalue ≈ [1.1e-296, 7.3e-143, 1.5e-08]
        @test df.best_dom_bits ≈ [974.4, 466.3, 9.8]
        @test df.best_dom_evalue ≈ [1.4e-296, 9.8e-143, 0.00018]

        # The last row is a serine/threonine kinase: a real low-scoring member of the
        # same beta-propeller superfamily. It is in the fixture because a table always
        # carries rows like it, which is the whole reason a cutoff exists.
        @test argmin(df.full_bits) == 3

        # A run that produced no hit at all is an empty table, never an error.
        empty_path = joinpath(mktempdir(), "empty.tblout")
        write(empty_path, "# target name\n#\n")
        @test nrow(parse_tblout(empty_path)) == 0
        @test nrow(parse_tblout(joinpath(mktempdir(), "absent.tblout"))) == 0
    end

    @testset "accession and TaxID extraction" begin
        @test uniprot_accession("sp|P16027|DHM1_PARDE") == "P16027"
        @test uniprot_accession("tr|I0JWN7|I0JWN7_9BACT Methanol dehydrogenase") == "I0JWN7"
        @test uniprot_accession("A0A2T6NCX9") == "A0A2T6NCX9"        # locally written header
        @test uniprot_accession("  Q5SLL9  extra text ") == "Q5SLL9"
        @test uniprot_accession("") == ""

        @test header_taxid("sp|P16027|DHM1_PARDE Methanol dehydrogenase OS=Paracoccus denitrificans OX=266 GN=mxaF PE=1 SV=2") == 266
        @test header_taxid("tr|X|X_9ZZZZ no organism field") === nothing
        # OX must not be read out of the middle of another token (OXA=1234).
        @test header_taxid("tr|X|X_9ZZZZ thing OXA=1234 OX=4321 PE=4") == 4321
    end

    @testset "frame: pool minus census, exact TaxID" begin
        pool = [266, 591154, 511746, 83333]
        census = [83333, 266, 99999]
        @test frame_taxids(pool, census) == [511746, 591154]
        # Exact matching only: a strain TaxID does not remove its species and vice versa.
        @test frame_taxids([1, 2], [10, 20]) == [1, 2]
        @test isempty(frame_taxids(Int[], [1, 2]))
        # Duplicates in the pool collapse; the result is sorted and unique.
        @test frame_taxids([5, 5, 3], Int[]) == [3, 5]
    end

    @testset "carrier exclusion join" begin
        frame = DataFrame(taxid = [100, 200, 300, 400], phylum = ["A", "B", "A", "C"])
        scan  = DataFrame(taxid = [100, 200, 300], n_hits = [0, 3, 0])

        kept, removed = exclude_carriers(frame, scan)
        @test kept.taxid == [100, 300, 400]      # 400 is absent from the scan and stays
        @test removed.taxid == [200]
        @test kept.phylum == ["A", "A", "C"]     # other columns travel with the rows
        @test nrow(kept) + nrow(removed) == nrow(frame)

        # Nothing flagged: the frame passes through untouched.
        k2, r2 = exclude_carriers(frame, DataFrame(taxid = [100, 200], n_hits = [0, 0]))
        @test nrow(k2) == 4 && nrow(r2) == 0
        # Everything flagged: the arm empties rather than silently keeping carriers.
        k3, r3 = exclude_carriers(frame, DataFrame(taxid = [100, 200, 300, 400], n_hits = [1, 1, 1, 1]))
        @test nrow(k3) == 0 && nrow(r3) == 4
    end

    @testset "empirical threshold and recovery" begin
        scores = collect(100.0:1.0:199.0)          # 100 calibration scores
        @test calibration_threshold(scores; rule = :min) == 100.0
        @test calibration_threshold(scores; rule = :p1) ≈ 100.99

        # The minimum keeps every scored accession by construction …
        @test recovery(scores, calibration_threshold(scores; rule = :min), length(scores)) == 1.0
        # … the 1st percentile gives up at most the bottom 1 %.
        @test recovery(scores, calibration_threshold(scores; rule = :p1), length(scores)) >= 0.99

        # Accessions that scored nothing at all stay in the denominator.
        @test recovery(scores, 100.0, 200) == 0.5
        @test recovery(Float64[], 10.0, 0) == 0.0
        @test_throws ErrorException calibration_threshold(Float64[])
        @test_throws ErrorException calibration_threshold(scores; rule = :median)
    end

    # The 100-score fixture above satisfies the >= 99 % criterion for a reason that has
    # nothing to do with the rule: at n = 100 the 1st percentile can drop at most one
    # score. Wherever 0.01·(n − 1) is not an integer the same rule drops more — at
    # n = 250 it drops 3, for a recovery of 0.988 — so the criterion is met only when
    # scores TIE at the cutoff. Both directions are asserted here, otherwise the suite
    # certifies a property the gate lacks.
    @testset "the recovery criterion is not satisfied by construction" begin
        distinct = collect(1.0:1.0:250.0)
        thr = calibration_threshold(distinct; rule = :p1)
        @test recovery(distinct, thr, length(distinct)) ≈ 0.988
        @test recovery(distinct, thr, length(distinct)) < 0.99

        # Ten scores tied at the cutoff: the 1st percentile lands inside the tie and
        # recovers everything, which is the shape a gate like this actually passes on.
        tied = vcat(fill(40.0, 10), collect(100.0:1.0:289.0))
        @test calibration_threshold(tied; rule = :p1) == 40.0
        @test recovery(tied, calibration_threshold(tied; rule = :p1), length(tied)) == 1.0
    end

    # One row per UniProtKB entry, as the union download delivers them. Taxon 10 is the
    # case a hit-derived membership column cannot express: a family carrier with no hit at
    # all. Taxon 30's entry carries no superfamily annotation, so a superfamily-only
    # prefilter never downloaded it — which is how such a taxon comes to have no hit.
    entries() = DataFrame(
        accession = ["A1", "A2", "A3", "A4", "A5", "A6", "A7"],
        taxid     = [ 10,   10,   20,   30,   40,   50,   60 ],
        len       = [ 600,  499,  700,  562,  499,  500,  800],
        interpro  = ["IPR011047;IPR017512;",   # full-length family member, taxon 10
                     "IPR011047;IPR017512;",   # 499 aa, same taxon: below the floor
                     "IPR011047;",             # superfamily only: not a family member
                     "IPR018391;IPR017512;",   # family member the prefilter cannot see
                     "IPR011047;IPR017512;",   # taxon 40 has nothing but a 499 aa entry
                     "IPR011047;IPR017512;",   # taxon 50 sits exactly on the floor
                     "IPR011047;IPR017512;"])  # taxon 60 is outside the frame
    frame_fx = [10, 20, 30, 40, 50]

    @testset "family membership is a property of the taxon, not of its hits" begin
        r = ipr017512_by_taxon(entries(), frame_fx)
        @test r.taxid == frame_fx                      # one row per frame taxon, in order
        @test names(r) == ["taxid", "n_ipr017512_entries", "n_ipr017512_ge500",
                            "in_ipr017512", "in_ipr017512_ge500",
                            "ipr017512_outside_prefilter"]

        # Taxon 10 has two family entries and NO hits anywhere in this fixture: the
        # roll-up never sees a hit table, so it says so regardless.
        t10 = r[r.taxid .== 10, :]
        @test t10.n_ipr017512_entries == [2]
        @test t10.in_ipr017512 == [true]

        # The >= 500 aa floor: the 499 aa entry does not count, the 600 aa one does.
        @test t10.n_ipr017512_ge500 == [1]
        @test r.in_ipr017512_ge500[r.taxid .== 40] == [false]   # only a 499 aa entry
        @test r.in_ipr017512[r.taxid .== 40] == [true]
        @test r.in_ipr017512_ge500[r.taxid .== 50] == [true]    # exactly 500 aa
    end

    @testset "family entries the superfamily prefilter cannot see" begin
        r = ipr017512_by_taxon(entries(), frame_fx)
        @test r.ipr017512_outside_prefilter[r.taxid .== 30] == [true]
        @test r.ipr017512_outside_prefilter[r.taxid .== 10] == [false]
        @test count(r.ipr017512_outside_prefilter) == 1
        # A taxon outside the frame contributes nothing, however it is annotated.
        @test !(60 in r.taxid)
        # A taxon with no family entry at all gets a zero row rather than no row, so the
        # per-taxon table covers the frame exactly.
        t20 = r[r.taxid .== 20, :]
        @test t20.n_ipr017512_entries == [0] && t20.in_ipr017512 == [false]
    end

    @testset "acceptance: an unflagged family taxon is a MISS" begin
        r = ipr017512_by_taxon(entries(), frame_fx)
        # Taxa 10, 30 and 50 carry a full-length family member. Flag two of them.
        @test ipr017512_unflagged(r, [10, 50]) == [30]
        @test verdict_of(isempty(ipr017512_unflagged(r, [10, 50]))) == "MISS"
        @test verdict_of(isempty(ipr017512_unflagged(r, [10, 30, 50]))) == "PASS"
        # A flagged taxon that carries nothing of the family neither helps nor hurts.
        @test ipr017512_unflagged(r, [10, 20, 50]) == [30]
        @test verdict_of(nothing) == "DISCLOSED"

        # The gated row the check feeds carries its own numbers: flagged of family taxa,
        # with the bar computed by the run rather than remembered from an earlier one.
        n_family = count(r.in_ipr017512_ge500)
        missed = ipr017512_unflagged(r, [10, 50])
        @test floor_verdict(n_family - length(missed), n_family, n_family) === false
        @test floor_verdict(n_family, n_family, n_family) === true
    end

    @testset "the specificity disclosure counts flagged taxa only" begin
        s = specificity_disclosure([0, 3, 1, 2, 0], [0, 0, 1, 2, 0])
        @test (s.n_flagged, s.n_without_ln_or_ca) == (3, 1)
        @test s.share ≈ 1 / 3
        @test isnan(specificity_disclosure([0, 0], [0, 0]).share)
    end

    @testset "a cached table names the -T it was produced at" begin
        dir = mktempdir()
        p = joinpath(dir, "candidates.tblout")
        write(p, TBLOUT_FIXTURE *
              "# Option settings: hmmsearch -o /dev/null --tblout x.tblout -T 12.5000 --incT 12.5000 --cpu 8 model.hmm seqs.fasta \n" *
              "# [ok]\n")
        @test tblout_threshold(p) ≈ 12.5
        # No option line (an older or hand-made table): unknown, never "the right one".
        q = joinpath(dir, "bare.tblout")
        write(q, TBLOUT_FIXTURE)
        @test tblout_threshold(q) === nothing
        @test tblout_threshold(joinpath(dir, "absent.tblout")) === nothing
        # The threshold must not be read out of another option's value.
        s = joinpath(dir, "other.tblout")
        write(s, "# Option settings: hmmsearch --cpu 8 --domT 55.0 model.hmm seqs.fasta \n")
        @test tblout_threshold(s) === nothing
    end

    # Two page sets that overlap, exactly as the two InterPro downloads do: P1 carries a
    # shared entry and a superfamily-only one, P2 carries the shared entry again and one
    # the superfamily prefilter never sees.
    function page_files(dir)
        mkpath(dir)
        p1 = joinpath(dir, "page_0001.fasta")
        p2 = joinpath(dir, "page_0002.fasta")
        write(p1, ">tr|SHARED|SHARED_9ZZZZ x OX=10\nMAAA\nMBBB\n>tr|SFONLY|SFONLY_9ZZZZ x OX=20\nMCCC\n")
        write(p2, ">tr|SHARED|SHARED_9ZZZZ x OX=10\nMAAA\nMBBB\n>tr|FAMONLY|FAMONLY_9ZZZZ x OX=30\nMDDD\n")
        return [p1, p2]
    end

    @testset "the union writes each record once" begin
        files = page_files(mktempdir())
        out = joinpath(mktempdir(), "candidates.fasta")
        n = filter_fasta(files, Set(["SHARED", "FAMONLY"]), out)
        # SHARED appears in both page sets; twice in the FASTA is twice in the tblout,
        # and its taxon's hit count would double.
        @test n == 2
        @test count(l -> startswith(l, '>'), readlines(out)) == 2
        @test fasta_accessions(out) == Set(["SHARED", "FAMONLY"])
        # Sequence lines travel with their header, and nothing unrequested is emitted.
        @test count(==("MCCC"), readlines(out)) == 0
        @test count(==("MDDD"), readlines(out)) == 1
    end

    @testset "a subset FASTA is reused only for the set it holds" begin
        files = page_files(mktempdir())
        out = joinpath(mktempdir(), "candidates.fasta")
        a = sync_subset_fasta(out, files, Set(["SHARED"]); label = "test")
        @test (a.n, a.rebuilt) == (1, true)
        b = sync_subset_fasta(out, files, Set(["SHARED"]); label = "test")
        @test (b.n, b.rebuilt) == (1, false)          # unchanged request, no rebuild
        c = sync_subset_fasta(out, files, Set(["SHARED", "SFONLY"]); label = "test")
        @test (c.n, c.rebuilt) == (2, true)           # the set changed: rebuilt
        @test fasta_accessions(out) == Set(["SHARED", "SFONLY"])
        # An accession the pages do not contain must stop the run, not shrink the search.
        @test_throws ErrorException sync_subset_fasta(out, files, Set(["GHOST"]); label = "test")
    end

    @testset "downloads are checked against the total the server announced" begin
        dir = mktempdir()
        t1 = joinpath(dir, "page_0001.tsv")
        write(t1, "Entry\tOrganism (ID)\tLength\tInterPro\n" *
                  "SHARED\t10\t600\tIPR011047;IPR017512;\n" *
                  "SFONLY\t20\t700\tIPR011047;\n")
        t2 = joinpath(dir, "page_0002.tsv")
        write(t2, "Entry\tOrganism (ID)\tLength\tInterPro\n" *
                  "FAMONLY\t30\t562\tIPR018391;IPR017512;\n")
        rows = parse_entry_tsv([t1, t2])
        @test nrow(rows) == 3                        # the header line is not a row
        @test rows.accession == ["SHARED", "SFONLY", "FAMONLY"]
        @test rows.len == [600, 700, 562]
        @test rows.taxid == [10, 20, 30]
    end

    @testset "a page short of the announced total stops the run" begin
        tdir = mktempdir()
        t1 = joinpath(tdir, "page_0001.tsv")
        write(t1, "Entry\tOrganism (ID)\tLength\tInterPro\n" *
                  "SHARED\t10\t600\tIPR011047;IPR017512;\n" *
                  "SFONLY\t20\t700\tIPR011047;\n" *
                  "FAMONLY\t30\t562\tIPR018391;IPR017512;\n")
        rows = parse_entry_tsv([t1])
        fdir = mktempdir()
        fa = joinpath(fdir, "page_0001.fasta")
        write(fa, ">tr|SHARED|S x\nMA\n>tr|SFONLY|S x\nMB\n>tr|FAMONLY|F x\nMC\n")
        @test verify_download("test", rows, [fa], 3) == 3
        # One page lost to a dropped connection: fewer rows than the server announced.
        @test_throws ErrorException verify_download("test", rows, [fa], 4)
        # TSV and FASTA out of step with each other.
        write(fa, ">tr|SHARED|S x\nMA\n>tr|SFONLY|S x\nMB\n>tr|OTHER|O x\nMC\n")
        @test_throws ErrorException verify_download("test", rows, [fa], 3)
    end

    @testset "a zero row can mean never searched" begin
        # 300 has no candidate entry, so its zero is "UniProtKB annotates no member of
        # either InterPro entry for this TaxID", not "the HMM found nothing". Both kinds
        # of zero are kept in the background: this step removes what it found.
        frame = DataFrame(taxid = [100, 200, 300])
        scan = DataFrame(taxid = [100, 200, 300], n_hits = [0, 2, 0],
                         n_candidates = [7, 3, 0])
        kept, removed = exclude_carriers(frame, scan)
        @test kept.taxid == [100, 300]
        @test removed.taxid == [200]
    end

    @testset "a gated row's verdict follows its own numbers" begin
        @test floor_verdict(4, 4, 4) === true
        @test floor_verdict(3, 4, 4) === false
        @test floor_verdict(0.99, 250, 0.98) === true
        @test floor_verdict(0.97, 250, 0.98) === false
        @test floor_verdict(2, missing, 1) === true      # a floor with no population
        # Nothing measured is never a pass: an empty denominator, an undefined rate, a
        # control taxon with no row at all.
        @test floor_verdict(0, 0, 0) === false
        @test floor_verdict(NaN, 0, 0.98) === false
        @test floor_verdict(missing, missing, 1) === false
        # No bar: a disclosure, whatever the numbers.
        @test floor_verdict(2, 4, missing) === nothing
        @test floor_verdict(0, 0, missing) === nothing

        tbl = acceptance_table()
        accept!(tbl, "gated, holds"; value = 7, denominator = 9, bar = 7,
                note = "bar computed; 2 of 9 short of all; 7 %")
        accept!(tbl, "gated, measured nothing"; value = NaN, denominator = 0, bar = 0.98)
        accept!(tbl, "disclosed"; value = 3, denominator = 9, note = "no bar")
        @test names(tbl) == ["check", "value", "denominator", "bar", "verdict", "note",
                             "run_id", "scan_sha256"]
        @test tbl.verdict == ["PASS", "MISS", "DISCLOSED"]
        @test (tbl.value, tbl.denominator, tbl.bar) ==
              (["7", "NaN", "3"], ["9", "0", "9"], ["7", "0.98", ""])
        @test acceptance_blocks(tbl)
        @test !acceptance_blocks(tbl[tbl.verdict .!= "MISS", :])
        # A share is written so it reads back as the same Float64: a value rounded for
        # display can print as its bar while failing it, and the re-derived verdict flips.
        accept!(tbl, "a share just under its bar"; value = 0.97999999, denominator = 3,
                bar = 0.98)
        @test tbl.verdict[end] == "MISS" && parse(Float64, tbl.value[end]) == 0.97999999
        tbl = tbl[1:3, :]
        # A cell with a tab or a newline would split the row it is written into, and one
        # with a double quote would be written quoted, which a tab-splitting reader keeps.
        @test_throws ErrorException accept!(tbl, "tabbed"; note = "a\tb")
        @test_throws ErrorException accept!(tbl, "broken"; note = "a\nb")
        @test_throws ErrorException accept!(tbl, "quoted \"check\""; note = "x")
        @test_throws ErrorException accept!(tbl, "quoted note"; note = "the \"bar\"")
        @test nrow(tbl) == 3

        # Written and read back: numbers as numbers, blanks as missing, text verbatim.
        p = joinpath(mktempdir(), "acceptance.tsv")
        write_acceptance(p, tbl)
        back = read_acceptance(p)
        @test back.check == tbl.check && back.verdict == tbl.verdict && back.note == tbl.note
        @test back.value[1] == 7.0 && isnan(back.value[2]) && back.value[3] == 3.0
        @test back.denominator == [9.0, 0.0, 9.0]
        @test back.bar[1] == 7.0 && back.bar[2] == 0.98 && ismissing(back.bar[3])
        # Every gated verdict can be re-derived from the row as read back.
        for r in eachrow(back)
            r.verdict == "DISCLOSED" && continue
            @test verdict_of(floor_verdict(r.value, r.denominator, r.bar)) == r.verdict
        end
        # No file is no verdicts, never "nothing failed"; another shape is an error.
        @test nrow(read_acceptance(joinpath(mktempdir(), "absent.tsv"))) == 0
        q = joinpath(mktempdir(), "old.tsv")
        write(q, "check\tvalue\tbar\tverdict\nc\t1\t1\tPASS\n")
        @test_throws ErrorException read_acceptance(q)
    end

    @testset "a recall that measured nothing is a MISS, and the record re-derives its pass" begin
        s = recall_summary(n_sampled = 40, full_flagged = [11, 12, 13], pre_flagged = [12, 13, 14],
                           n_never_searched = 6, bar = 0.98, seed = 12345)
        @test (s.n_full_flagged, s.n_recovered, s.missed) == (3, 2, [11])
        @test s.rate ≈ 2 / 3
        @test !s.pass
        # Nothing flagged by the full scan: the rate is undefined and the verdict is MISS.
        # Read as "everything recovered" this is how an empty proteome download passes.
        e = recall_summary(n_sampled = 40, full_flagged = Int[], pre_flagged = [1, 2],
                           n_never_searched = 6, bar = 0.98, seed = 12345)
        @test isnan(e.rate) && !e.pass && verdict_of(e.pass) == "MISS"
        @test read_keyvals(write_keyvals(joinpath(mktempdir(), "r.txt"),
            recall_keyvals(e; threshold_bits = 12.5, uniprot_release = "REL_A")))["rate"] == "NaN"
        # Every taxon recovered on a non-empty sample passes; a recall at the bar passes too.
        p = recall_summary(n_sampled = 40, full_flagged = [1, 2], pre_flagged = [1, 2, 3],
                           n_never_searched = 0, bar = 0.98, seed = 12345)
        @test p.pass && p.rate == 1.0 && isempty(p.missed)
        at_bar = recall_summary(n_sampled = 40, full_flagged = 1:4, pre_flagged = 1:3,
                                n_never_searched = 0, bar = 0.75, seed = 12345)
        @test at_bar.pass
        # The record carries the counts its pass is decided on, so a reader can re-derive
        # the verdict as n_full_flagged > 0 and n_recovered / n_full_flagged >= bar, and a
        # written pass that disagrees with its own numbers can be refused.
        for x in (s, e, p, at_bar)
            path = write_keyvals(joinpath(mktempdir(), "recall.txt"),
                                 recall_keyvals(x; threshold_bits = 12.5, uniprot_release = "REL_A"))
            r = read_recall(path)
            @test r.pass === x.pass
            @test r.pass === (r.n_full_flagged > 0 && r.n_recovered / r.n_full_flagged >= r.bar)
            # The bound comes back as the Float64 that was written, and a reader can
            # recompute it from the two counts beside it: NaN exactly when nothing was flagged.
            @test isequal(r.lower_bound_95, x.lower_bound_95)
            @test isequal(r.lower_bound_95, clopper_pearson_lower(r.n_recovered, r.n_full_flagged))
            @test isnan(r.lower_bound_95) == (r.n_full_flagged == 0)
        end
        # The gate is the rate, not the bound: every taxon recovered on a small sample passes
        # the scan's own check while its bound sits far below the same bar.
        @test p.pass && p.lower_bound_95 < p.bar
        @test p.lower_bound_95 == 0.05^(1 / 2)
        @test isnan(e.lower_bound_95)
        @test s.lower_bound_95 == clopper_pearson_lower(2, 3)
    end

    @testset "a perfect rate is not a bound" begin
        @test isnan(clopper_pearson_lower(0, 0))          # an empty sample bounds nothing
        @test clopper_pearson_lower(0, 8) == 0.0
        for n in (7, 33, 120)
            @test clopper_pearson_lower(n, n) ≈ 0.05^(1 / n) atol = 1e-12
        end
        # The defining equation, checked in exact arithmetic: at the bound, k or more
        # successes out of n have probability alpha.
        tail(k, n, p) = sum(binomial(big(n), big(i)) * big(p)^i * (1 - big(p))^(n - i) for i in k:n)
        for (k, n) in ((3, 5), (21, 25), (38, 41), (117, 120))
            @test Float64(tail(k, n, clopper_pearson_lower(k, n))) ≈ 0.05 atol = 1e-9
        end

        # A second route to the same number, sharing nothing with the script's: the binomial
        # upper tail equals the regularised incomplete beta function I_p(k, n - k + 1), whose
        # integrand is expanded here into a polynomial with exact integer coefficients and
        # solved by bisection in 256-bit arithmetic. The script sums binomial terms in log
        # space in Float64; if either route were wrong, they would not meet to 1e-12.
        function beta_route(k, n; alpha = 0.05)
            setprecision(BigFloat, 256) do
                c = big(n) * binomial(big(n - 1), big(k - 1))       # n! / ((k-1)! (n-k)!)
                terms = [(binomial(big(n - k), big(j)) * (-1)^j, k + j) for j in 0:(n - k)]
                cdf(p) = c * sum(b * p^e / e for (b, e) in terms)
                lo, hi = big(0.0), big(1.0)
                for _ in 1:200
                    mid = (lo + hi) / 2
                    cdf(mid) < alpha ? (lo = mid) : (hi = mid)
                end
                Float64((lo + hi) / 2)
            end
        end
        for (k, n) in ((3, 5), (21, 25), (38, 41), (57, 60), (117, 120))
            @test clopper_pearson_lower(k, n) ≈ beta_route(k, n) atol = 1e-12
        end
        # At k == n the closed form, the script and the second route are one number.
        for n in (7, 33)
            @test clopper_pearson_lower(n, n) ≈ beta_route(n, n) atol = 1e-12
        end
        # The recall summary records exactly this routine's 95 % bound on its own counts,
        # whatever the routine's default level is set to.
        for (full, pre) in ((1:25, 1:21), (1:12, 1:12), (Int[], [1]), (1:41, 4:41))
            x = recall_summary(n_sampled = 90, full_flagged = full, pre_flagged = pre,
                               n_never_searched = 0, bar = 0.98, seed = 12345)
            @test isequal(x.lower_bound_95,
                          clopper_pearson_lower(x.n_recovered, x.n_full_flagged; alpha = 0.05))
        end
        x = recall_summary(n_sampled = 90, full_flagged = 1:25, pre_flagged = 1:21,
                           n_never_searched = 0, bar = 0.98, seed = 12345)
        @test x.lower_bound_95 ≈ beta_route(21, 25) atol = 1e-12
        # Written to recall.txt and read back, the bound is the same Float64, not the rate's
        # six decimals: a consumer that recomputes it must find no rounding to forgive.
        back = read_recall(write_keyvals(joinpath(mktempdir(), "recall.txt"),
            recall_keyvals(x; threshold_bits = 12.5, uniprot_release = "REL_A")))
        @test back.lower_bound_95 === x.lower_bound_95
        @test round(back.lower_bound_95; digits = 6) != back.lower_bound_95
        # More evidence tightens it; a bound never exceeds the observed rate.
        @test clopper_pearson_lower(12, 12) < clopper_pearson_lower(30, 30)
        @test clopper_pearson_lower(9, 10) < 0.9
        @test clopper_pearson_lower(9, 10; alpha = 0.2) > clopper_pearson_lower(9, 10)
        @test_throws ErrorException clopper_pearson_lower(5, 4)
        @test_throws ErrorException clopper_pearson_lower(1, 4; alpha = 1.5)
    end

    @testset "the recall rows: the gate on the rate, the bound disclosed" begin
        gate = "prefilter recall"
        bound = "prefilter recall, exact one-sided 95 % lower bound"
        rows(s) = recall_acceptance!(acceptance_table(), s)
        rsum(full, pre) = recall_summary(n_sampled = 60, full_flagged = full,
                                         pre_flagged = pre, n_never_searched = 3,
                                         bar = 0.98, seed = 12345)

        # Every flagged taxon recovered on a small sample: the rate clears the bar, the
        # bound does not. The scan's gate passes and nothing blocks publication, because the
        # prefilter did find every carrier; the bound is on the record for its consumer.
        perfect = rsum(1:9, 1:12)
        @test perfect.lower_bound_95 < perfect.bar
        t = rows(perfect)
        @test t.check == [gate, bound]
        @test t.verdict == ["PASS", "DISCLOSED"]
        @test !acceptance_blocks(t)
        @test t.bar == ["0.98", ""]                         # the bound row has no bar at all
        @test t.value == [num_cell(1.0), num_cell(perfect.lower_bound_95)]
        @test t.denominator == ["9", "9"]
        @test occursin("disclosed", t.note[2])

        # A carrier the prefilter missed fails the gate; the bound row still never blocks.
        missed = rsum(1:9, 2:9)
        m = rows(missed)
        @test m.verdict == ["MISS", "DISCLOSED"]
        # Nothing flagged: the rate is undefined and the gate fails, the bound is NaN and
        # disclosed, so the one MISS comes from the gate and not from the bound.
        empty_scan = rows(rsum(Int[], 1:3))
        @test empty_scan.verdict == ["MISS", "DISCLOSED"]
        @test empty_scan.value == ["NaN", "NaN"] && empty_scan.denominator == ["0", "0"]

        # Read back, each check is present exactly once, the gate re-derives from its row,
        # and the bound cell is the recorded Float64 a consumer can compare with recall.txt.
        p = joinpath(mktempdir(), "acceptance.tsv")
        write_acceptance(p, vcat(t, m))
        back = read_acceptance(p)
        @test count(==(gate), back.check) == 2 && count(==(bound), back.check) == 2
        g = back[1, :]
        @test verdict_of(floor_verdict(g.value, g.denominator, g.bar)) == g.verdict == "PASS"
        @test back.value[2] === perfect.lower_bound_95
        @test all(ismissing, back.bar[back.check .== bound])
        once = read_acceptance(write_acceptance(joinpath(mktempdir(), "a.tsv"), t))
        @test count(==(gate), once.check) == 1 && count(==(bound), once.check) == 1
    end

    @testset "the recall sample is a prefix of one seeded permutation" begin
        frame = collect(3001:3250)                     # 250 arbitrary TaxIDs
        small = recall_sample(frame, 30, 4242)
        large = recall_sample(frame, 90, 4242)
        @test length(small.taxa) == 30 && length(large.taxa) == 90
        @test small.taxa == sort(small.draw_order)
        # Enlarging the sample adds draws after the earlier ones and swaps none out.
        @test large.draw_order[1:30] == small.draw_order
        @test issubset(Set(small.taxa), Set(large.taxa))
        # The seed, not luck, fixes the sample; beyond the frame it is the whole frame.
        @test recall_sample(frame, 30, 4243).taxa != small.taxa
        @test recall_sample(frame, 400, 4242).taxa == frame
        @test_throws ErrorException recall_sample(frame, 0, 4242)

        e = sample_extends(small.taxa, large.draw_order)
        @test (e.n_earlier, e.n_matched, e.ok) == (30, 30, true)
        # Same size, one taxon different: a frame or seed that moved between runs.
        swapped = copy(small.taxa)
        swapped[1] = first(setdiff(frame, large.taxa))
        @test sample_extends(swapped, large.draw_order).ok === false
        # The earlier 30 are inside the 90 but are not its FIRST 30 draws.
        @test sample_extends(large.draw_order[31:60], large.draw_order).ok === false
        # An earlier sample larger than this run's cannot be extended by it.
        @test sample_extends(large.taxa, small.draw_order).ok === false
        @test sample_extends(Int[], large.draw_order).ok === nothing
        # The acceptance row re-derives the same verdict from the two numbers it carries.
        for earlier in (small.taxa, swapped, large.draw_order[31:60], large.taxa)
            x = sample_extends(earlier, large.draw_order)
            y = sample_extends(earlier, small.draw_order)
            @test floor_verdict(x.n_matched, x.n_earlier, x.n_earlier) === x.ok
            @test floor_verdict(y.n_matched, y.n_earlier, y.n_earlier) === y.ok
        end
    end

    @testset "the earlier sample is read from a scan table, or the check stops" begin
        dir = mktempdir()
        p = joinpath(dir, "scan.csv")
        CSV.write(p, DataFrame(taxid = [7, 3, 5, 9], n_hits = [0, 1, 0, 2],
                               in_prefilter_recall_sample = [true, false, true, true]))
        @test csv_recall_sample(p; required = true) == [5, 7, 9]
        # No earlier table: nothing to extend, unless one was named on the command line.
        absent = joinpath(dir, "absent.csv")
        @test csv_recall_sample(absent; required = false) == Int[]
        @test_throws ErrorException csv_recall_sample(absent; required = true)
        # A table without the column, or with no sampled taxon, is not a scan table.
        q = joinpath(dir, "other.csv")
        CSV.write(q, DataFrame(taxid = [7, 3], n_hits = [0, 1]))
        @test_throws ErrorException csv_recall_sample(q; required = false)
        r = joinpath(dir, "none.csv")
        CSV.write(r, DataFrame(taxid = [7, 3], in_prefilter_recall_sample = [false, false]))
        @test_throws ErrorException csv_recall_sample(r; required = false)

        # A failing run displaces the published table to .stale. The next run still checks
        # its sample against that table instead of finding nothing and comparing nothing.
        pub = joinpath(mktempdir(), "scan.csv")
        @test default_recall_reference(pub) == pub                 # neither: no earlier sample
        @test csv_recall_sample(default_recall_reference(pub); required = false) == Int[]
        cp(p, pub * ".stale")
        @test default_recall_reference(pub) == pub * ".stale"
        @test csv_recall_sample(default_recall_reference(pub); required = false) == [5, 7, 9]
        cp(q, pub)                                                  # a published table wins
        @test default_recall_reference(pub) == pub
    end

    # Three tiny proteome downloads. Taxon 22's second entry names another TaxID, the case
    # the recall has to count rather than assume away.
    function proteome_dir()
        prot = joinpath(mktempdir(), "proteomes")
        mkpath(prot)
        write(joinpath(prot, "tax_11.fasta"), ">tr|P11A|P11A_9ZZZZ x OX=11\nMAAA\n")
        write(joinpath(prot, "tax_22.fasta"),
              ">tr|P22A|P22A_9ZZZZ x OX=22\nMBBB\n>tr|P22B|P22B_9ZZZZ x OX=2222\nMBBC\n")
        write(joinpath(prot, "tax_33.fasta"), ">tr|P33A|P33A_9ZZZZ x OX=33\nMCCC\n")
        return prot
    end

    @testset "the recall concatenation is reused only for the sample it was built from" begin
        prot = proteome_dir()
        out = joinpath(dirname(prot), "recall_proteomes.fasta")
        a = sync_recall_fasta(out, prot, [11, 22])
        @test (a.n, a.rebuilt) == (2, true)
        @test fasta_accessions(out) == Set(["P11A", "P22A", "P22B"])
        @test !sync_recall_fasta(out, prot, [11, 22]).rebuilt        # same sample, reused
        # The sample grew: the old concatenation would score the recall on too few taxa.
        @test sync_recall_fasta(out, prot, [11, 22, 33]).rebuilt
        @test fasta_accessions(out) == Set(["P11A", "P22A", "P22B", "P33A"])
        # A proteome replaced on disk changes the source, so the concatenation is stale.
        write(joinpath(prot, "tax_33.fasta"), ">tr|P33A|P33A_9ZZZZ x OX=33\nMCCCDDDD\n")
        @test sync_recall_fasta(out, prot, [11, 22, 33]).rebuilt
        # A concatenation with no manifest beside it at all.
        rm(out * ".manifest")
        @test sync_recall_fasta(out, prot, [11, 22, 33]).rebuilt
        @test_throws ErrorException sync_recall_fasta(out, prot, [11, 44])
        # A download that lacks its final newline must not glue a sequence to a header.
        write(joinpath(prot, "tax_11.fasta"), ">tr|P11A|P11A_9ZZZZ x OX=11\nMAAA")
        sync_recall_fasta(out, prot, [11, 33])
        @test count(l -> startswith(l, '>'), readlines(out)) == 2
        @test "MAAA" in readlines(out)
    end

    @testset "a full-proteome hit counts for the TaxID its entry names" begin
        prot = proteome_dir()
        m = proteome_taxon_map(prot, [11, 22, 33])
        @test (m.n_entries, m.n_other_taxid) == (4, 1)
        @test m.acc2tax["P22B"] == 2222 && m.acc2tax["P22A"] == 22
        # P22B is filed under taxon 22 but names 2222: the prefilter would attribute it to
        # 2222 too, so it flags neither sampled taxon and is counted as set aside.
        f = recall_flagged(["P11A", "P22B"], m.acc2tax, Set([11, 22, 33]))
        @test f.flagged == Set([11]) && f.n_outside == 1
        g = recall_flagged(["P11A", "P22A", "P22B", "P33A"], m.acc2tax, Set([11, 22, 33]))
        @test g.flagged == Set([11, 22, 33]) && g.n_outside == 1
        # An accession no sampled proteome holds is set aside, never guessed into a taxon.
        @test recall_flagged(["GHOST"], m.acc2tax, Set([11])).n_outside == 1
        @test m.n_no_taxid == 0

        # A header without an OX= token has no TaxID of its own: it falls back to the
        # taxon downloaded, and that fallback is counted instead of passing unseen.
        write(joinpath(prot, "tax_33.fasta"),
              ">tr|P33A|P33A_9ZZZZ x OX=33\nMCCC\n>tr|P33N|P33N_9ZZZZ no organism token\nMCCD\n")
        n = proteome_taxon_map(prot, [11, 22, 33])
        @test (n.n_entries, n.n_other_taxid, n.n_no_taxid) == (5, 1, 1)
        @test n.acc2tax["P33N"] == 33
        @test proteome_taxon_map(prot, [11, 22]).n_no_taxid == 0
    end

    @testset "each proteome download is checked against the search-space release" begin
        prot = proteome_dir()
        write_keyvals(release_record_file(prot, 11), ["X-UniProt-Release" => "REL_B",
                                                      "proteome" => "UP_FIXTURE_1"])
        write_keyvals(release_record_file(prot, 22), ["X-UniProt-Release" => "REL_B",
                                                      "proteome" => ""])
        a = proteome_release_audit(prot, [11, 22, 33], "REL_B")
        @test (a.n_matching, a.n_other, a.n_unrecorded, a.n_recorded) == (2, 0, 1, 2)
        @test a.releases == ["REL_B"]
        # Taxon 22 had no reference proteome, so its download is every entry of the TaxID;
        # taxon 33 has no record, so it cannot be counted either way.
        @test a.n_no_reference_proteome == 1
        # The gated row: recorded downloads from the search-space release, of all recorded.
        @test floor_verdict(a.n_matching, a.n_recorded, a.n_recorded) === true
        write_keyvals(release_record_file(prot, 33), ["X-UniProt-Release" => "REL_C"])
        b = proteome_release_audit(prot, [11, 22, 33], "REL_B")
        @test (b.n_matching, b.n_other, b.n_unrecorded) == (2, 1, 0)
        @test floor_verdict(b.n_matching, b.n_recorded, b.n_recorded) === false
        # A record that does not name a proteome at all cannot vouch for one.
        @test b.n_no_reference_proteome == 2
        # A release record left blank is unrecorded, and an unknown search-space release
        # matches nothing.
        write_keyvals(release_record_file(prot, 33), ["X-UniProt-Release" => ""])
        @test proteome_release_audit(prot, [33], "REL_B").n_unrecorded == 1
        @test proteome_release_audit(prot, [11], "").n_matching == 0
    end

    @testset "one run, one identity: four files promoted together or not at all" begin
        dir = mktempdir()
        csvp = joinpath(dir, "scan.csv")
        thr = joinpath(dir, "threshold.txt")
        rec = joinpath(dir, "recall.txt")
        acc = joinpath(dir, "acceptance.tsv")
        four = (thr, rec, acc, csvp)
        # An earlier run's four files, which a failing run must not leave consumable.
        earlier_csv = "taxid,n_hits\n7,1\n"
        write(csvp, earlier_csv)
        write(thr, "threshold_bits\t12.5000\n")
        write(rec, "pass\ttrue\n")
        write(acc, "check\tvalue\tdenominator\tbar\tverdict\tnote\trun_id\tscan_sha256\n")

        s = recall_summary(n_sampled = 40, full_flagged = [11, 12, 13],
                           pre_flagged = [12, 13, 14], n_never_searched = 6,
                           bar = 0.98, seed = 12345)
        rpairs = recall_keyvals(s; threshold_bits = 12.5, uniprot_release = "REL_A")
        tpairs = ["rule" => "fixture", "threshold_bits" => "12.5000",
                  "active_site_column" => "17"]

        # ---- a MISS: nothing promoted, the earlier table displaced, pending files kept
        failing = acceptance_table()
        accept!(failing, "a floor that holds"; value = 4, denominator = 4, bar = 4)
        accept!(failing, "a recall below its bar"; value = s.rate,
                denominator = s.n_full_flagged, bar = s.bar)
        accept!(failing, "a disclosure"; value = 1, denominator = 4)
        miss_csv = Vector{UInt8}(codeunits("taxid,n_hits\n9,0\n"))
        m = publish_scan(; out_csv = csvp, csv = miss_csv, threshold_path = thr,
                           threshold_pairs = tpairs, recall_path = rec, recall_pairs = rpairs,
                           acceptance_path = acc, acceptance = failing, run_id = "RUN-A")
        @test !m.promoted
        @test !any(isfile, four)                              # no final name holds anything
        @test all(p -> isfile(p * ".pending"), four)          # only this run's pending files
        @test read(csvp * ".stale", String) == earlier_csv    # the earlier table, unconsumable
        @test all(p -> isfile(p * ".stale"), four)
        @test sort(basename.(m.staled)) == sort(basename.(collect(four)) .* ".stale")
        # The pending records name the pending table, not the one that was displaced.
        miss_sha = bytes2hex(sha256(read(csvp * ".pending")))
        @test m.scan_sha256 == miss_sha
        @test read_recall(rec * ".pending").scan_sha256 == miss_sha
        @test read_threshold(thr * ".pending").scan_sha256 == miss_sha
        @test unique(read_acceptance(acc * ".pending").scan_sha256) == [miss_sha]
        @test read_acceptance(acc * ".pending").verdict == ["PASS", "MISS", "DISCLOSED"]

        # ---- a PASS: all four promoted, every record naming the promoted table's hash
        passing = acceptance_table()
        accept!(passing, "a floor that holds"; value = 4, denominator = 4, bar = 4)
        accept!(passing, "a disclosure"; value = 1, denominator = 4)
        pass_csv = Vector{UInt8}(codeunits("taxid,n_hits\n9,0\n11,2\n"))
        p = publish_scan(; out_csv = csvp, csv = pass_csv, threshold_path = thr,
                           threshold_pairs = tpairs, recall_path = rec, recall_pairs = rpairs,
                           acceptance_path = acc, acceptance = passing, run_id = "RUN-B")
        @test p.promoted && isempty(p.staled)
        @test all(isfile, four)
        @test !any(q -> isfile(q * ".pending"), four)
        @test read(csvp) == pass_csv
        sha = bytes2hex(sha256(read(csvp)))                   # what a consumer computes
        @test p.scan_sha256 == sha
        @test read_threshold(thr).scan_sha256 == sha
        @test read_recall(rec).scan_sha256 == sha
        @test unique(read_acceptance(acc).scan_sha256) == [sha]
        @test read_threshold(thr).run_id == read_recall(rec).run_id == "RUN-B"
        @test unique(read_acceptance(acc).run_id) == ["RUN-B"]
        @test read_recall(rec).derived_at == "RUN-B"
        # The earlier failure's displaced table stays displaced; nothing reads it.
        @test read(csvp * ".stale", String) == earlier_csv

        # ---- the recall record: every key, in order, and back again
        @test [first(split(l, '\t')) for l in eachline(rec)] ==
              ["n_sampled", "n_full_flagged", "n_recovered", "rate", "bar", "pass",
               "n_never_searched_in_sample", "missed_taxids", "seed", "threshold_bits",
               "uniprot_release", "lower_bound_95", "run_id", "scan_sha256", "derived_at"]
        r = read_recall(rec)
        @test (r.n_sampled, r.n_full_flagged, r.n_recovered) == (40, 3, 2)
        @test r.rate ≈ 2 / 3 atol = 1e-6
        @test (r.bar, r.pass, r.n_never_searched_in_sample) == (0.98, false, 6)
        @test (r.missed_taxids, r.seed, r.threshold_bits, r.uniprot_release) ==
              ([11], 12345, 12.5, "REL_A")
        # The promoted record carries the bound of its own two counts, bit for bit.
        @test r.lower_bound_95 === clopper_pearson_lower(2, 3) === s.lower_bound_95
        @test read_recall(rec * ".stale").lower_bound_95 === nothing   # an older record's shape
        t = read_threshold(thr)
        @test (t.rule, t.threshold_bits, t.active_site_column) == ("fixture", 12.5, 17)
        @test [first(split(l, '\t')) for l in eachline(thr)][end-2:end] ==
              ["run_id", "scan_sha256", "derived_at"]

        # A run that failed after an earlier success displaces the promoted four again.
        m2 = publish_scan(; out_csv = csvp, csv = miss_csv, threshold_path = thr,
                            threshold_pairs = tpairs, recall_path = rec, recall_pairs = rpairs,
                            acceptance_path = acc, acceptance = failing, run_id = "RUN-C")
        @test !m2.promoted && !any(isfile, four)
        @test read(csvp * ".stale") == pass_csv
        @test_throws ErrorException publish_scan(; out_csv = csvp, csv = pass_csv,
            threshold_path = thr, threshold_pairs = tpairs, recall_path = rec,
            recall_pairs = rpairs, acceptance_path = acc, acceptance = passing, run_id = "")
    end

    @testset "the CSV bytes that are hashed are the bytes the file holds" begin
        df = DataFrame(taxid = [5, 6], best_bits = [missing, 40.0], hit_classes = ["", "Gap;Lanthanide"],
                       in_prefilter_recall_sample = [true, false])
        p = joinpath(mktempdir(), "t.csv")
        CSV.write(p, df)
        @test csv_bytes(df) == read(p)
    end

    @testset "an empty, absent or malformed record field reads as nothing" begin
        dir = mktempdir()
        p = joinpath(dir, "recall.txt")
        write_keyvals(p, ["n_sampled" => "40", "n_full_flagged" => "", "rate" => "NaN",
                          "pass" => "yes", "missed_taxids" => "11;x", "scan_sha256" => "abc123",
                          "lower_bound_95" => "0.9x"])
        r = read_recall(p)
        @test r.n_sampled == 40
        @test r.n_full_flagged === nothing             # empty
        @test r.n_recovered === nothing                # absent
        @test isnan(r.rate)                            # an undefined rate stays undefined
        @test r.lower_bound_95 === nothing             # an unreadable bound is no bound
        @test isnan(read_recall(write_keyvals(joinpath(dir, "nan.txt"),
                                              ["lower_bound_95" => "NaN"])).lower_bound_95)
        @test r.pass === nothing                       # only true / false are verdicts
        @test r.missed_taxids === nothing              # a partly unreadable list is no list
        @test r.scan_sha256 === nothing                # not a sha256, cannot match a table
        a = read_recall(joinpath(dir, "absent.txt"))
        @test all(isnothing, (a.pass, a.n_full_flagged, a.run_id, a.scan_sha256,
                              a.lower_bound_95))
        @test read_recall(write_keyvals(joinpath(dir, "r2.txt"),
                                        ["missed_taxids" => ""])).missed_taxids == Int[]
        @test read_threshold(joinpath(dir, "absent.txt")).threshold_bits === nothing
    end

    @testset "an option the script does not define stops the run" begin
        spec = Dict("--alpha" => :value, "--beta" => :flag)
        @test parse_options(["--alpha", "7", "--beta"], spec) ==
              Dict("--alpha" => "7", "--beta" => "")
        @test parse_options(String[], spec) == Dict{String,String}()
        @test_throws ErrorException parse_options(["--alpa", "7"], spec)       # misspelt
        @test_throws ErrorException parse_options(["--alpha"], spec)           # no value
        @test_throws ErrorException parse_options(["--alpha", "--beta"], spec) # value forgotten
        @test_throws ErrorException parse_options(["--beta", "--beta"], spec)  # twice
        @test_throws ErrorException parse_options(["7"], spec)                 # a stray word
        # Every option the script defines is one its usage text names.
        @test all(o -> occursin(o, USAGE), keys(OPTION_SPEC))
        @test parse_options(["--recall-n", "40", "--threshold", "12.5"], OPTION_SPEC) ==
              Dict("--recall-n" => "40", "--threshold" => "12.5")

        cli = Dict("--threshold" => "12.5", "--recall-n" => "40")
        @test option_number(cli, "--threshold", Float64, nothing; what = "a bit score") == 12.5
        @test option_number(cli, "--recall-n", Int, 150; min_value = 1) == 40
        @test option_number(Dict{String,String}(), "--recall-n", Int, 150) == 150
        @test option_number(Dict{String,String}(), "--threshold", Float64, nothing) === nothing
        # A bad value is an error, never the default under the requested name.
        @test_throws ErrorException option_number(Dict("--recall-n" => "0"), "--recall-n",
                                                  Int, 150; min_value = 1)
        @test_throws ErrorException option_number(Dict("--recall-n" => "4.5"), "--recall-n", Int, 150)
        @test_throws ErrorException option_number(Dict("--threshold" => "abc"), "--threshold",
                                                  Float64, nothing)
        @test_throws ErrorException option_number(Dict("--threshold" => "Inf"), "--threshold",
                                                  Float64, nothing)
    end

    # The script copied into an empty tree: its REPO is that tree, which holds no data, so
    # any stage past the command line would fail and the first write would create `dat`.
    # The copy is started as its own program, because the guard under test only runs then.
    @testset "--help and a mistyped option exit before anything is loaded or written" begin
        tree = mktempdir()
        mkpath(joinpath(tree, "scripts"))
        script = joinpath(tree, "scripts", "background_pqq_scan.jl")
        cp(joinpath(@__DIR__, "..", "scripts", "background_pqq_scan.jl"), script)
        function run_copy(args)
            out, err = IOBuffer(), IOBuffer()
            p = run(pipeline(ignorestatus(`$(Base.julia_cmd()) --startup-file=no $script $args`);
                             stdout = out, stderr = err))
            return (code = p.exitcode, out = String(take!(out)), err = String(take!(err)))
        end
        h = run_copy(["--help"])
        @test h.code == 0
        @test occursin("--recall-n", h.out) && occursin("--recall-reference", h.out)
        @test !occursin("Activating", h.out * h.err)           # no package was loaded
        u = run_copy(["--recall-nn", "300"])
        @test u.code == 1
        @test occursin("unknown option '--recall-nn'", u.err)
        v = run_copy(["--recall-n", "0"])
        @test v.code == 1 && occursin("--recall-n", v.err)
        @test !occursin("Activating", u.err * v.err)
        @test readdir(tree) == ["scripts"]                     # nothing was written
        @test readdir(joinpath(tree, "scripts")) == ["background_pqq_scan.jl"]
    end
end
