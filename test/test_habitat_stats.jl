#=
    Offline self-test for the habitat enrichment estimators. No network, no data files:
    every fixture is small enough to be checked by hand on paper, which is the point — an
    estimator that is only ever run on 10^3 rows can be wrong by a factor and look
    plausible.

    Four failures this suite exists to prevent:

    1. AN ASYMMETRIC UNIT. Aggregating the foreground to unique genus-by-category PAIRS
       while the background stays one row per taxon lets a genus that occupies four
       categories contribute four foreground units and one background unit: n1 and n2 are
       then not in the same currency, and the contrast cannot be read as a conservative
       bound. `collapse_to_genus` puts both arms on genera, and the test below pins n1 to
       the number of genera — which is strictly smaller than the number of pairs in the
       fixture, so an asymmetric unit cannot pass it.

    2. SILENT DRIFT BETWEEN THE TWO WEIGHT FORMS. One numeric path serves both weight
       representations, so a Bool hit vector and the equivalent 0/1 Float vector must give
       bitwise identical output through every path — if generalising to fractional weights
       perturbed the arithmetic, results would shift while nothing in the shape of the
       code showed it.

    3. AN UNCONDITIONAL CONTINUITY CORRECTION. Adding 0.5 to every cell (rather than only
       to empty ones) shrinks every estimate towards zero and would quietly damp real
       contrasts. The correction has to fire on an empty cell and nowhere else.

    4. A STANDARDISED INTERVAL THAT RESTS ON AN EMPTY STRATUM CELL, SILENTLY. When one
       stratum carries most of the foreground weight, the background hit count inside it is
       what decides the interval — and an empty cell there adds nothing to p0 and nothing to
       v0, so the estimate and its interval look exactly as they would if that stratum had
       been absent. The row totals cannot show it: `a` and `b` can be comfortably above
       any gating threshold while a contributing stratum cell is zero. `n_zero_stratum_cells` and
       `min_bg_stratum_cell` are that disclosure, and the tests below pin what counts as a
       cell (contributing strata only, after the collapse) and that a fractional cell is
       reported as it is rather than rounded into or out of emptiness.

    Run:  julia --project=. test/test_habitat_stats.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
include(joinpath(@__DIR__, "..", "src", "ModuleHabitatStats.jl"))
using .HabitatStats

# ---------------------------------------------------------------------------
# The worked example. Foreground: 8 taxa of phylum P (4 hits) and 2 of phylum B
# (1 hit). Background: 20 of P (5 hits) and 20 of B (10 hits).
#
#   p1      = 5/10                              = 0.5
#   p0_std  = 0.8*(5/20) + 0.2*(10/20)          = 0.30
#   log2 RR = log2(0.5/0.30)                    = 0.7369655941662...
#   var     = (1/5 - 1/10) + v0/0.30^2,
#             v0 = 0.8^2*5*15/20^3 + 0.2^2*10*10/20^3 = 0.0065
#                                               = 0.1722222222222...
#
# The thresholds are lowered to 2 so that P and B stay separate strata; at the module
# defaults this fixture would collapse into a single "Other phyla" stratum, which is
# checked separately below.
# ---------------------------------------------------------------------------
const FG_PHY = vcat(fill("P", 8), fill("B", 2))
const FG_HIT = vcat(trues(4), falses(4), [true], [false])
const BG_PHY = vcat(fill("P", 20), fill("B", 20))
const BG_HIT = vcat(trues(5), falses(15), trues(10), falses(10))

@testset "standardised — hand-computed two-stratum example" begin
    st = standardised(FG_PHY, FG_HIT, BG_PHY, BG_HIT; min_fg_k = 2, min_bg_k = 2)

    @test st.n_strata == 2
    @test st.p1 == 0.5
    @test st.p0 ≈ 0.8 * (5 / 20) + 0.2 * (10 / 20) atol = 1e-12
    @test st.p0 ≈ 0.30 atol = 1e-12
    @test st.log2_rr ≈ 0.736965594166206 atol = 1e-12
    @test round(st.log2_rr, digits = 3) == 0.737

    # the variance formula, written out by hand rather than re-using the module's
    v0_hand   = 0.8^2 * 5 * (20 - 5) / 20^3 + 0.2^2 * 10 * (20 - 10) / 20^3
    var_hand  = (1 / 5 - 1 / 10) + v0_hand / 0.30^2
    @test v0_hand ≈ 0.0065 atol = 1e-12
    @test st.var ≈ var_hand atol = 1e-10
    @test st.var ≈ 0.17222222222222222 atol = 1e-10

    # the interval is that variance carried to log2, 1.96 SE either side — so the
    # reported bounds must reconstruct the reported variance
    se = sqrt(st.var) / log(2)
    @test st.ci_lo ≈ st.log2_rr - 1.96se atol = 1e-12
    @test st.ci_hi ≈ st.log2_rr + 1.96se atol = 1e-12
    @test ((st.ci_hi - st.log2_rr) / 1.96 * log(2))^2 ≈ st.var atol = 1e-10

    # standardisation moves the estimate: the crude ratio uses the pooled background
    # rate 15/40, the standardised one the reweighted 0.30, and they must differ here
    @test crude(5, 10, 15, 40).log2_rr ≈ log2(0.5 / (15 / 40)) atol = 1e-12
    @test !isapprox(st.log2_rr, crude(5, 10, 15, 40).log2_rr; atol = 1e-6)
end

@testset "Bool and the equivalent 0/1 Float weights agree exactly" begin
    fgf = Float64.(FG_HIT); bgf = Float64.(BG_HIT)

    sb = standardised(FG_PHY, FG_HIT, BG_PHY, BG_HIT; min_fg_k = 2, min_bg_k = 2)
    sf = standardised(FG_PHY, fgf,    BG_PHY, bgf;    min_fg_k = 2, min_bg_k = 2)
    # `===` on Float64 is bitwise: no tolerance is allowed to hide drift here
    for k in keys(sb)
        @test getfield(sb, k) === getfield(sf, k)
    end

    cb = crude(sum(FG_HIT), 10, sum(BG_HIT), 40)
    cf = crude(sum(fgf),    10, sum(bgf),    40)
    for k in keys(cb)
        @test getfield(cb, k) === getfield(cf, k)
    end
    @test cb.rounded_for_fisher == false      # 0/1 Floats are integral: nothing to round

    # …and through the row builder, where the q-values and the gating also have to match
    rb = contrast_rows(FG_PHY, Dict("H" => FG_HIT), BG_PHY, Dict("H" => BG_HIT);
                       contrast = "fixture", unit = "taxa", min_fg_k = 2, min_bg_k = 2)
    rf = contrast_rows(FG_PHY, Dict("H" => fgf),    BG_PHY, Dict("H" => bgf);
                       contrast = "fixture", unit = "taxa", min_fg_k = 2, min_bg_k = 2)
    @test length(rb) == length(rf) == 1
    @test isequal(rb[1], rf[1])               # isequal, so a NaN q-value compares equal
    @test rb[1].a == 5 && rb[1].n1 == 10 && rb[1].b == 15 && rb[1].n2 == 40
end

@testset "thin phyla collapse into 'Other phyla'" begin
    # At the module defaults P (8) is under MIN_FG_K and B (2) is under both thresholds,
    # so every unit lands in one stratum …
    st1 = standardised(FG_PHY, FG_HIT, BG_PHY, BG_HIT)
    @test st1.n_strata == 1
    # … and a single stratum is arithmetically the same as not standardising at all,
    # which is the check that the collapse relabels BOTH arms and not just one.
    @test st1.p0 ≈ 15 / 40 atol = 1e-12
    @test st1.log2_rr ≈ crude(5, 10, 15, 40).log2_rr atol = 1e-12

    # relabelling by hand must reproduce it exactly
    st2 = standardised(fill("Other phyla", 10), FG_HIT, fill("Other phyla", 40), BG_HIT)
    @test st2.log2_rr === st1.log2_rr && st2.var === st1.var

    # A phylum thick enough in the foreground but thin in the BACKGROUND collapses too:
    # the alternative is a full stratum weight resting on a handful of background units,
    # here a 3/3 rate that would be carried into p0 as if it were an estimate.
    fg   = vcat(fill("P", 20), fill("Rare1", 10), fill("Rare2", 10))
    bg   = vcat(fill("P", 40), fill("Rare1", 3),  fill("Rare2", 5))
    fgh  = vcat(trues(10), falses(10), trues(5), falses(5), trues(2), falses(8))
    bgh  = vcat(trues(10), falses(30), trues(3),            falses(5))
    st3  = standardised(fg, fgh, bg, bgh)
    @test st3.n_strata == 2                        # P, and the two thin ones merged
    @test st3.p0 ≈ 0.5 * (10 / 40) + 0.5 * (3 / 8) atol = 1e-12
    lumped = vcat(fill("P", 20), fill("Other phyla", 20))
    st3h = standardised(lumped, fgh, vcat(fill("P", 40), fill("Other phyla", 8)), bgh)
    @test st3h.log2_rr === st3.log2_rr && st3h.var === st3.var
    # and the collapse is load-bearing: keeping the thin strata gives a different answer
    st3k = standardised(fg, fgh, bg, bgh; min_fg_k = 1, min_bg_k = 1)
    @test st3k.n_strata == 3
    @test st3k.p0 ≈ 0.5 * (10 / 40) + 0.25 * 1.0 + 0.25 * 0.0 atol = 1e-12
    @test !isapprox(st3k.p0, st3.p0; atol = 1e-6)

    # A stratum with foreground units but NO background units is dropped rather than
    # guessed at. Pinned deliberately: the remaining weights are not renormalised, so the
    # standardised background share is computed over less than the whole foreground.
    st4 = standardised(vcat(fill("P", 10), fill("X", 10)), vcat(trues(5), falses(5), trues(5), falses(5)),
                       fill("P", 40), vcat(trues(10), falses(30)); min_fg_k = 1, min_bg_k = 0)
    @test st4.n_strata == 1
    @test st4.p0 ≈ 0.5 * (10 / 40) atol = 1e-12
end

# ---------------------------------------------------------------------------
# The zero-cell example. The foreground is the worked example's; the background keeps its
# 20 P and 20 B units but has NO hits in B. Stratum B still contributes (it has background
# units), and its cell is exactly zero.
#
#   p0_std  = 0.8*(5/20) + 0.2*(0/20)              = 0.20
#   log2 RR = log2(0.5/0.20) = log2(2.5)
#   var     = (1/5 - 1/10) + v0/0.20^2,
#             v0 = 0.8^2*5*15/20^3 + 0.2^2*0*20/20^3 = 0.006 + 0
#                                                  = 0.25
#
# The second term of v0 is the point: the empty cell adds nothing to the variance.
# ---------------------------------------------------------------------------
const BG_HIT_Z = vcat(trues(5), falses(15), falses(20))

@testset "standardised — an empty background stratum cell is disclosed" begin
    st = standardised(FG_PHY, FG_HIT, BG_PHY, BG_HIT_Z; min_fg_k = 2, min_bg_k = 2)

    @test st.n_strata == 2
    @test st.n_zero_stratum_cells == 1
    @test st.min_bg_stratum_cell === 0.0             # an exact zero, carried as a Float64
    @test st.p0 ≈ 0.20 atol = 1e-12
    @test st.log2_rr ≈ log2(2.5) atol = 1e-12
    @test st.var ≈ (1 / 5 - 1 / 10) + (0.8^2 * 5 * 15 / 20^3) / 0.20^2 atol = 1e-12
    @test st.var ≈ 0.25 atol = 1e-12
    @test isfinite(st.ci_lo) && isfinite(st.ci_hi)   # an interval IS reported

    # WHY IT MUST BE DISCLOSED. Take stratum B's background away altogether, so B is
    # dropped rather than contributing an empty cell: the estimate and the variance are
    # bitwise the same. Nothing in the estimate distinguishes "no hits among 20 units" from
    # "no units at all", so only the disclosure can.
    gone = standardised(FG_PHY, FG_HIT, fill("P", 20), BG_HIT_Z[1:20]; min_fg_k = 2, min_bg_k = 0)
    @test gone.log2_rr === st.log2_rr && gone.var === st.var
    @test gone.n_strata == 1
    # a dropped stratum is not a cell, so it is not an empty one either
    @test gone.n_zero_stratum_cells == 0
    @test gone.min_bg_stratum_cell === 5.0

    # The row totals do not show it. `b` is non-zero, the row is tested, and no continuity
    # correction fires, while a contributing stratum cell is empty.
    rows = contrast_rows(FG_PHY, Dict("H" => FG_HIT), BG_PHY, Dict("H" => BG_HIT_Z);
                         contrast = "fixture", unit = "taxa", min_fg_k = 2, min_bg_k = 2)
    r = only(rows)
    @test r.b == 5 && r.tested && !r.continuity_corrected
    @test r.n_zero_stratum_cells == 1
    @test r.min_bg_stratum_cell === 0.0

    # Only strata the FOREGROUND occupies contribute. Twenty extra background units in a
    # phylum the foreground never reaches, all without hits, are not cells of this estimate.
    bg_extra = vcat(BG_PHY, fill("Z", 20)); bgh_extra = vcat(BG_HIT_Z, falses(20))
    ex = standardised(FG_PHY, FG_HIT, bg_extra, bgh_extra; min_fg_k = 2, min_bg_k = 2)
    @test ex.n_strata == 2
    @test ex.n_zero_stratum_cells == 1
    @test ex.min_bg_stratum_cell === 0.0

    # When every contributing cell is empty the estimate is undefined, and the disclosure
    # still says why rather than being blanked along with it.
    allzero = standardised(FG_PHY, FG_HIT, BG_PHY, falses(40); min_fg_k = 2, min_bg_k = 2)
    @test isnan(allzero.log2_rr) && isnan(allzero.ci_lo)
    @test allzero.n_strata == 2
    @test allzero.n_zero_stratum_cells == 2
    @test allzero.min_bg_stratum_cell === 0.0

    # With no foreground hits no stratum is visited: nothing contributed, nothing to report.
    nofg = standardised(FG_PHY, falses(10), BG_PHY, BG_HIT_Z; min_fg_k = 2, min_bg_k = 2)
    @test nofg.n_strata == 0
    @test nofg.n_zero_stratum_cells == 0
    @test isnan(nofg.min_bg_stratum_cell)

    # the 0/1 Float equivalent agrees bitwise on the new fields as well, on the fixture
    # where they are non-trivial, and through the row builder
    sf = standardised(FG_PHY, Float64.(FG_HIT), BG_PHY, Float64.(BG_HIT_Z);
                      min_fg_k = 2, min_bg_k = 2)
    @test sf.n_zero_stratum_cells === st.n_zero_stratum_cells
    @test sf.min_bg_stratum_cell === st.min_bg_stratum_cell
    rf = contrast_rows(FG_PHY, Dict("H" => Float64.(FG_HIT)), BG_PHY,
                       Dict("H" => Float64.(BG_HIT_Z));
                       contrast = "fixture", unit = "taxa", min_fg_k = 2, min_bg_k = 2)
    @test only(rf).n_zero_stratum_cells === r.n_zero_stratum_cells
    @test only(rf).min_bg_stratum_cell === r.min_bg_stratum_cell
    @test isequal(only(rf), r)
end

@testset "standardised — a thin fractional cell is reported as it is, not rounded" begin
    # Genus-collapsed weights are shares, so a stratum cell can hold a third of a unit.
    # Stratum P's background has one unit carrying 1/3 and nineteen carrying nothing;
    # stratum Q's background is fully hit. The row total b is over twenty, the P cell is
    # a third of a unit.
    fg  = vcat(fill("P", 10), fill("Q", 10))
    fgh = vcat(trues(5), falses(5), trues(5), falses(5))
    bg  = vcat(fill("P", 20), fill("Q", 20))
    bgh = vcat([1 / 3], zeros(19), ones(20))
    st  = standardised(fg, fgh, bg, bgh; min_fg_k = 2, min_bg_k = 2)

    @test st.n_strata == 2
    @test st.p0 ≈ 0.5 * ((1 / 3) / 20) + 0.5 * (20 / 20) atol = 1e-12
    @test st.min_bg_stratum_cell ≈ 1 / 3 atol = 1e-12
    @test !isinteger(st.min_bg_stratum_cell)
    # rounding either way would misstate the cell: to 0 it becomes an empty cell that is
    # not empty, and to 1 a whole unit that does not exist
    @test st.min_bg_stratum_cell != round(1 / 3)
    @test st.min_bg_stratum_cell < 1
    @test st.n_zero_stratum_cells == 0
    @test sum(bgh) > 20                               # the total gives no hint of it

    # a thin but whole cell comes back as the same value, just as a Float64
    bgh2 = vcat(trues(2), falses(18), trues(20))
    st2 = standardised(fg, fgh, bg, bgh2; min_fg_k = 2, min_bg_k = 2)
    @test st2.min_bg_stratum_cell === 2.0
    @test st2.n_zero_stratum_cells == 0
end

@testset "the 'Other phyla' collapse is one cell in the disclosure" begin
    # The collapse fixture from above: Rare1 and Rare2 are thick in the foreground but thin
    # in the background, so at the module defaults they merge into one "Other phyla"
    # stratum. The disclosure has to count that merged stratum as ONE cell.
    fg  = vcat(fill("P", 20), fill("Rare1", 10), fill("Rare2", 10))
    bg  = vcat(fill("P", 40), fill("Rare1", 3),  fill("Rare2", 5))
    fgh = vcat(trues(10), falses(10), trues(5), falses(5), trues(2), falses(8))

    # Rare1 fully hit (3 of 3), Rare2 empty (0 of 5). Merged, the cell holds 3 hits and is
    # not empty; kept apart, Rare2 is an empty cell of its own.
    bgh = vcat(trues(10), falses(30), trues(3), falses(5))
    merged = standardised(fg, fgh, bg, bgh)
    @test merged.n_strata == 2
    @test merged.n_zero_stratum_cells == 0
    @test merged.min_bg_stratum_cell === 3.0
    apart = standardised(fg, fgh, bg, bgh; min_fg_k = 1, min_bg_k = 1)
    @test apart.n_strata == 3
    @test apart.n_zero_stratum_cells == 1
    @test apart.min_bg_stratum_cell === 0.0

    # Both thin strata empty. Merged they are one empty cell, not two.
    bgh0 = vcat(trues(10), falses(30), falses(3), falses(5))
    merged0 = standardised(fg, fgh, bg, bgh0)
    @test merged0.n_strata == 2
    @test merged0.n_zero_stratum_cells == 1
    @test merged0.min_bg_stratum_cell === 0.0
    apart0 = standardised(fg, fgh, bg, bgh0; min_fg_k = 1, min_bg_k = 1)
    @test apart0.n_zero_stratum_cells == 2

    # relabelling by hand gives the same disclosure, so the counts follow the collapse in
    # BOTH arms rather than being taken over the raw labels
    lumped_fg = vcat(fill("P", 20), fill("Other phyla", 20))
    lumped_bg = vcat(fill("P", 40), fill("Other phyla", 8))
    for h in (bgh, bgh0)
        auto = standardised(fg, fgh, bg, h)
        hand = standardised(lumped_fg, fgh, lumped_bg, h)
        @test hand.n_zero_stratum_cells === auto.n_zero_stratum_cells
        @test hand.min_bg_stratum_cell === auto.min_bg_stratum_cell
    end
end

@testset "crude — Haldane correction only on an empty cell" begin
    both = crude(5, 10, 15, 40)
    @test both.continuity_corrected == false
    @test both.log2_rr ≈ log2((5 / 10) / (15 / 40)) atol = 1e-12
    # the uncorrected Katz standard error, by hand
    @test (both.ci_hi - both.log2_rr) / 1.96 ≈
          sqrt(1 / 5 + 1 / 15 - 1 / 10 - 1 / 40) / log(2) atol = 1e-12

    for (a, b) in ((0, 15), (5, 0))
        c = crude(a, 10, b, 40)
        @test c.continuity_corrected == true
        @test isfinite(c.log2_rr) && isfinite(c.ci_lo) && isfinite(c.ci_hi)
    end
    zero_fg = crude(0, 10, 15, 40)
    @test zero_fg.log2_rr ≈ log2((0.5 / 11) / (15.5 / 41)) atol = 1e-12
    # a small non-zero cell must NOT be corrected — that is what would damp real contrasts
    @test crude(1, 10, 15, 40).continuity_corrected == false
    @test crude(1, 10, 15, 40).log2_rr ≈ log2((1 / 10) / (15 / 40)) atol = 1e-12

    # p is Fisher's exact on the 2x2 table, and the margins stay consistent
    @test 0.0 <= both.fisher_p <= 1.0
    @test crude(5, 10, 5, 10).log2_rr == 0.0
    @test crude(5, 10, 5, 10).fisher_p == 1.0
end

@testset "crude — fractional hits are rounded for Fisher only" begin
    frac = crude(2.5, 10, 3.0, 40)
    @test frac.rounded_for_fisher == true
    # the point estimate and the interval use the UNROUNDED sums
    @test frac.log2_rr ≈ log2((2.5 / 10) / (3.0 / 40)) atol = 1e-12
    @test (frac.ci_hi - frac.log2_rr) / 1.96 ≈
          sqrt(1 / 2.5 + 1 / 3.0 - 1 / 10 - 1 / 40) / log(2) atol = 1e-12
    # and the p-value is the one for the rounded table
    @test frac.fisher_p == crude(2, 10, 3, 40).fisher_p
    @test crude(3.0, 10, 3.0, 40).rounded_for_fisher == false   # integral Floats
    @test crude(3.0, 10, 3.0, 40).fisher_p == crude(3, 10, 3, 40).fisher_p
end

@testset "bh — q-values on a known vector" begin
    # n = 5; sorted p*n/k = [0.005, 0.02, 0.065, 0.05125, 0.9], and the running minimum
    # from the top pulls the third value down to the fourth.
    ps = [0.001, 0.008, 0.039, 0.041, 0.9]
    q  = bh(ps)
    @test q ≈ [0.005, 0.02, 0.05125, 0.05125, 0.9] atol = 1e-12
    @test issorted(q)                                    # monotone for sorted input
    # the order of the input must not matter: permute, and the q-values follow
    perm = [4, 1, 5, 3, 2]
    @test bh(ps[perm]) ≈ q[perm] atol = 1e-12
    @test bh([0.04]) ≈ [0.04] atol = 1e-12
    @test all(bh(fill(0.5, 4)) .≈ 0.5)
end

# ---------------------------------------------------------------------------
# Genus fixture. Six foreground taxa in three genera plus one row with no genus;
# the genera occupy FIVE distinct genus-by-category pairs, so the pair count and the
# genus count cannot be confused with each other.
# ---------------------------------------------------------------------------
const G_GEN = ["Methylobacterium", "Methylobacterium", "Methylobacterium",
               "Bradyrhizobium", "Bradyrhizobium", "Hyphomicrobium", missing]
const G_PHY = ["Pseudomonadota", "Pseudomonadota", "Pseudomonadota",
               "Pseudomonadota", "Bacillota", "Pseudomonadota", "Pseudomonadota"]
const G_CAT = ["Terrestrial", "Terrestrial", "Aquatic",
               "Terrestrial", "Host-plant", "Aquatic", "Terrestrial"]

@testset "collapse_to_genus — one unit per genus, on both arms" begin
    col = collapse_to_genus(G_GEN, G_PHY, flag_weights(G_CAT))

    @test col.genus == ["Bradyrhizobium", "Hyphomicrobium", "Methylobacterium"]
    @test col.n_taxa == [2, 1, 3]
    @test col.n_dropped == 1                       # the row with no genus
    # modal phylum: a 2-1 majority, and a 1-1 tie resolved alphabetically
    @test col.phylum == ["Bacillota", "Pseudomonadota", "Pseudomonadota"]
    @test collapse_to_genus(["g", "g"], ["Bacillota", "Pseudomonadota"],
                            Dict("H" => [true, false])).phylum == ["Bacillota"]

    # a genus weight is the SHARE of its taxa in the category …
    @test col.weights["Terrestrial"] ≈ [1 / 2, 0.0, 2 / 3] atol = 1e-12
    @test col.weights["Aquatic"]     ≈ [0.0, 1.0, 1 / 3] atol = 1e-12
    @test col.weights["Host-plant"]  ≈ [1 / 2, 0.0, 0.0] atol = 1e-12
    # … so every genus is worth exactly one unit in total, however many taxa it has
    for g in eachindex(col.genus)
        @test sum(col.weights[c][g] for c in keys(col.weights)) ≈ 1.0 atol = 1e-12
    end
end

@testset "the genus unit is symmetric: n1 and n2 both count genera" begin
    bg_gen = ["Escherichia", "Escherichia", "Bacillus", "Bacillus", "Bacillus", "Vibrio"]
    bg_phy = ["Pseudomonadota", "Pseudomonadota", "Bacillota", "Bacillota", "Bacillota",
              "Pseudomonadota"]
    bg_cat = ["Aquatic", "Terrestrial", "Terrestrial", "Terrestrial", "Aquatic", "Aquatic"]

    fg = collapse_to_genus(G_GEN, G_PHY, flag_weights(G_CAT))
    bg = collapse_to_genus(bg_gen, bg_phy, flag_weights(bg_cat))
    rows = contrast_rows(fg.phylum, fg.weights, bg.phylum, bg.weights;
                         contrast = "fixture", unit = "genera", min_fg_k = 1, min_bg_k = 1)

    # THE REGRESSION. Unique genus-by-category pairs number 5 in the foreground and 4 in
    # the background; genera number 3 and 3. Both arms must report the genus count.
    n_fg_pairs = length(unique(collect(zip(G_GEN, G_CAT)))) - 1   # minus the no-genus row
    @test n_fg_pairs == 5
    @test all(r -> r.n1 == 3, rows)
    @test all(r -> r.n2 == 3, rows)
    @test all(r -> r.n1 != n_fg_pairs, rows)
    @test all(r -> r.unit == "genera", rows)

    # the hit totals are fractional, and every unit is still worth one
    terr = only(filter(r -> r.habitat == "Terrestrial", rows))
    @test terr.a ≈ 1 / 2 + 2 / 3 atol = 1e-12
    @test terr.b ≈ 1 / 2 + 2 / 3 atol = 1e-12     # Escherichia 1/2, Bacillus 2/3
    @test terr.log2_rr ≈ 0.0 atol = 1e-12         # identical shares in identical arms
    @test terr.rounded_for_fisher == true
    @test sum(r.a for r in rows) ≈ 3.0 atol = 1e-12
    @test sum(r.b for r in rows) ≈ 3.0 atol = 1e-12
end

@testset "contrast_rows — categories, gating and the family contrast" begin
    fg_cat = vcat(fill("Aquatic", 6), fill("Terrestrial", 4))
    bg_cat = vcat(fill("Aquatic", 10), fill("Terrestrial", 28), fill("Geological", 2))
    rows = contrast_rows(FG_PHY, flag_weights(fg_cat), BG_PHY, flag_weights(bg_cat);
                         contrast = "fixture", unit = "taxa", min_fg_k = 2, min_bg_k = 2)

    @test [r.habitat for r in rows] == ["Aquatic", "Geological", "Terrestrial"]
    @test [r.contrast for r in rows] == fill("fixture", 3)
    geo = only(filter(r -> r.habitat == "Geological", rows))
    @test geo.a == 0 && geo.b == 2                 # absent from one arm, still reported
    @test geo.tested == false                      # a+b = 2 < MIN_CELL
    @test isnan(geo.q_bh)                          # untested rows carry no q-value
    @test all(r -> r.tested == (r.a + r.b >= 5), rows)
    @test all(r -> r.pct_fg ≈ 100r.a / r.n1, rows)

    # min_cell is a kwarg, and it compares the SUM, so fractional hits gate too
    lifted = contrast_rows(FG_PHY, flag_weights(fg_cat), BG_PHY, flag_weights(bg_cat);
                           contrast = "fixture", unit = "taxa", min_cell = 100,
                           min_fg_k = 2, min_bg_k = 2)
    @test all(r -> r.tested == false, lifted)
    @test all(r -> isnan(r.q_bh), lifted)

    # non-exclusive flags: a unit carrying "A+B" belongs to both categories, which is what
    # gives the chemistry rows a family contrast to be ticked against
    combo = flag_weights(["Hypersaline+Thermal", "Thermal", "", missing]; splitcombo = true)
    @test sort(collect(keys(combo))) == ["Hypersaline", "Thermal"]
    @test combo["Hypersaline"] == [true, false, false, false]
    @test combo["Thermal"]     == [true, true, false, false]
    @test !haskey(flag_weights(["Hypersaline+Thermal"]), "Thermal")   # off by default

    # the same weights run against a second foreground arm give the family rows: the
    # contrast label is what separates them in the concatenated table
    fam = contrast_rows(FG_PHY, flag_weights(fg_cat), BG_PHY, flag_weights(bg_cat);
                        contrast = "family", unit = "taxa (chemistry)",
                        min_fg_k = 2, min_bg_k = 2)
    @test all(r -> r.contrast == "family" && r.unit == "taxa (chemistry)", fam)
    @test [r.habitat for r in fam] == [r.habitat for r in rows]

    # a weight vector of the wrong length is an error, not a silent recycle
    @test_throws ErrorException contrast_rows(FG_PHY, Dict("A" => trues(3)),
                                              BG_PHY, Dict("A" => trues(40));
                                              contrast = "x", unit = "taxa")
end

@testset "contrast_rows — disclosure columns appended, earlier columns unmoved" begin
    # The table these rows become has consumers that were written against its earlier
    # layout, so the two disclosure columns go after every existing one and nothing before
    # them may move or be renamed.
    earlier = (:contrast, :unit, :habitat, :a, :n1, :b, :n2, :pct_fg, :pct_bg,
               :log2_rr, :ci_lo, :ci_hi, :log2_rr_std, :ci_lo_std, :ci_hi_std, :n_strata,
               :fisher_p, :q_bh, :tested, :continuity_corrected, :rounded_for_fisher)
    fg_cat = vcat(fill("Aquatic", 6), fill("Terrestrial", 4))
    bg_cat = vcat(fill("Aquatic", 10), fill("Terrestrial", 28), fill("Geological", 2))
    rows = contrast_rows(FG_PHY, flag_weights(fg_cat), BG_PHY, flag_weights(bg_cat);
                         contrast = "fixture", unit = "taxa", min_fg_k = 2, min_bg_k = 2)
    for r in rows
        @test keys(r) == (earlier..., :n_zero_stratum_cells, :min_bg_stratum_cell)
    end

    # and every row carries exactly what standardised() reports for its category. Here
    # the background's Aquatic units all sit in stratum P, so Aquatic has an empty B cell
    # beside a background total of 10; Geological has no foreground hits and so no
    # contributing strata at all.
    fw = flag_weights(fg_cat); bw = flag_weights(bg_cat)
    for r in rows
        st = standardised(FG_PHY, get(fw, r.habitat, falses(10)), BG_PHY,
                          get(bw, r.habitat, falses(40)); min_fg_k = 2, min_bg_k = 2)
        @test r.n_zero_stratum_cells === st.n_zero_stratum_cells
        @test isequal(r.min_bg_stratum_cell, st.min_bg_stratum_cell)
    end
    aq = only(filter(r -> r.habitat == "Aquatic", rows))
    @test aq.b == 10 && aq.n_zero_stratum_cells == 1 && aq.min_bg_stratum_cell === 0.0
    geo = only(filter(r -> r.habitat == "Geological", rows))
    @test geo.n_strata == 0 && geo.n_zero_stratum_cells == 0 && isnan(geo.min_bg_stratum_cell)
    ter = only(filter(r -> r.habitat == "Terrestrial", rows))
    @test ter.n_zero_stratum_cells == 0 && ter.min_bg_stratum_cell === 10.0
end

@testset "well_supported_mask — the thinly annotated tail" begin
    n_distinct = [1, 2, 3, 7, missing, 0]
    @test well_supported_mask(n_distinct) == BitVector([0, 0, 1, 1, 0, 0])
    @test well_supported_mask(n_distinct; min_distinct = 2) == BitVector([0, 1, 1, 1, 0, 0])
    @test well_supported_mask(n_distinct; min_distinct = 1) == BitVector([1, 1, 1, 1, 0, 0])
    @test count(well_supported_mask(Int[])) == 0
    # it is a mask, so it composes with the arm vectors by indexing
    keep = well_supported_mask([3, 1, 4])
    @test ["a", "b", "c"][keep] == ["a", "c"]
end
