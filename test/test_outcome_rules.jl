#=
    Offline self-test for the recorded outcome rules. No network, no data files and no
    estimator: every input is a bare bound pair, a ratio, a list of category labels or a
    record's fields, so each case can be checked by reading it.

    EVERY FIXTURE IS ABSTRACT. Category labels are role names such as "kept-u1", digests are
    repeated hex patterns, and counts are small arbitrary integers. None is taken from a run
    of any analysis, so each case pins a logical property of a rule and says nothing about
    what any run found.

    The failures this suite exists to prevent:

    1. AN UNCOMPUTED ESTIMATE REPORTED AS A NULL RESULT. A Bool interval test reads a NaN
       interval as "includes zero", and a consumer prints that as a finding. Undecidable has
       to survive as its own value, and every rule that consumes an interval has to carry it
       through instead of collapsing it into one side — including the code mapper a consumer
       uses for rules keyed on a single interval.

    2. "NEW" DECIDED AGAINST A CONSTANT. A comparison made against a hard-coded tuple of
       names rather than against the anchor lists a category that already excluded zero in
       the anchor run as new. A category that excludes zero in both runs and
       that no other rule governs must come out pre-existing; an implementation that ignores
       the anchor files it under `legend_only` and fails here, by name.

    3. A LOSS THAT WAS NEVER MEASURED. A category that excluded zero in the anchor and has no
       computed primary interval must not be reported as ceasing to exclude zero. It has to
       come back as undecidable in the primary run, apart from the categories whose primary
       interval was computed and includes zero.

    4. A STRONG ROBUSTNESS TIER WITHOUT ITS GRID. A depletion whose primary interval
       excludes zero is only a candidate for the two strong tiers; they are defined by the
       sensitivity grid. With no grid the answer must be undecidable, not whichever tier the
       primary row resembles.

    5. THE STRONGEST GRADE FROM AN UNBOUNDED RATIO. A ratio above one, a negative ratio, and
       the huge ratio a near-zero denominator produces all clear a bare threshold and so all
       return the strongest sharing grade. Each is pinned to :none here, labelled as a
       refusal rather than as the grading form that happened to see it.

    6. A LICENCE THAT IS REALLY A FLAG. A predicate that reduces to "the exclusion flag was
       set" grants the licence on no measurement at all, because both counts it reads are
       zero by construction. Below, each input is spoiled on its own — to a failing value
       and to an unknown one — and each spoilage alone must withdraw the licence.

    7. A LICENCE READ OFF A RECORD NOBODY CHECKED. A scan table that covers the frame but
       reports no carrier, a recall record that says pass while measuring nothing, a record
       from another scan, and a frame the scan does not cover each license the strong code
       where these guards are absent. They are rebuilt here from the pure guards, end to end,
       so that each stays a failing case of this suite.

    8. A RECORD'S OTHER NUMBERS TAKEN AS WRITTEN. With every digest intact, a rate that
       disagrees with its own counts, a bar lowered in the record whose counts it judges,
       more flagged taxa than were sampled, an acceptance table with its gated rows cut out,
       and a flagged count the table read does not have would each license the strong code too.
       Each is spoiled on its own below and must be refused for its own reason.

    9. A RECALL BAR MET BY A POINT ESTIMATE. A sample that recovered every taxon it flagged
       has a point estimate of one at any size. The verdict is taken on the exact lower
       confidence bound, so a clean sample too small for its bound to reach the bar must not
       pass, and the same counts under a bar the bound does reach must.

    Sets are compared as sorted vectors throughout, never in hash order.

    Run:  julia --project=. test/test_outcome_rules.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
include(joinpath(@__DIR__, "..", "src", "ModuleOutcomeRules.jl"))
using .OutcomeRules

@testset "outcome rules" begin

# The keys the reason strings are prefixed with, sorted, and the check that there is exactly
# one reason and it concerns `key`. Local to this suite, so the aggregate runner, which
# includes every suite into one process, cannot pick up a different definition.
reason_keys(res) = sort([String(first(split(r, ':'))) for r in res.reasons])
only_reason(res, key) = length(res.reasons) == 1 && startswith(res.reasons[1], key * ":")

@testset "interval_excludes_zero — three states, never two" begin
    @test interval_excludes_zero(-1.2, -0.1) === :excludes_zero     # depleted
    @test interval_excludes_zero(0.1, 2.0)   === :excludes_zero     # enriched
    @test interval_excludes_zero(-0.5, 0.5)  === :includes_zero
    @test interval_excludes_zero(1, 2)       === :excludes_zero     # integer bounds are bounds
    @test interval_excludes_zero(0.3, Inf)   === :excludes_zero
    @test interval_excludes_zero(-Inf, Inf)  === :includes_zero

    # a bound that reaches zero has not excluded it
    @test interval_excludes_zero(-0.4, 0.0) === :includes_zero
    @test interval_excludes_zero(0.0, 0.4)  === :includes_zero
    @test interval_excludes_zero(0.0, 0.0)  === :includes_zero

    # THE REGRESSION. An estimate that was never computed is neither side. One NaN bound
    # is enough: (NaN, -0.1) would read as "excludes zero" to a test that only looked at hi.
    for (lo, hi) in ((NaN, 0.5), (-0.5, NaN), (NaN, NaN), (NaN, -0.1), (0.1, NaN),
                     (missing, 0.5), (0.1, nothing), (nothing, nothing), (missing, missing))
        @test interval_excludes_zero(lo, hi) === :undecidable
    end

    # a pair handed over the wrong way round is not an interval, whatever its signs
    @test interval_excludes_zero(0.5, -0.5)  === :undecidable
    @test interval_excludes_zero(-0.1, -1.2) === :undecidable
    # nor is a flag in a bound slot, or an unparsed text column
    @test interval_excludes_zero(true, true)   === :undecidable
    @test interval_excludes_zero("0.1", "0.2") === :undecidable
end

@testset "interval_code — two caller codes, and a third that is not the caller's" begin
    ic(lo, hi) = interval_code(lo, hi; excludes = :code_x, includes = :code_y)

    @test ic(-1.2, -0.1) === :code_x          # excludes zero on the depleted side
    @test ic(0.1, 2.0)   === :code_x          # and on the enriched side
    @test ic(-0.5, 0.5)  === :code_y
    @test ic(-0.4, 0.0)  === :code_y          # a bound reaching zero includes it
    @test ic(0, 0)       === :code_y

    # THE REGRESSION. An uncomputed or malformed interval maps to neither caller code.
    for (lo, hi) in ((NaN, 0.5), (-0.5, NaN), (NaN, NaN), (NaN, -0.1), (0.1, NaN),
                     (missing, 0.5), (0.1, nothing), (0.5, -0.5), (true, true), ("0.1", "0.2"))
        @test ic(lo, hi) === :undecidable
    end

    # The two rule shapes a consumer maps — one acts when the interval includes zero, the
    # other when it excludes zero — are the same mapper with the codes swapped, and they
    # agree with the three-valued test case by case.
    table = [(-1.0, -0.2) => (:excludes_zero, :off, :keep),
             (0.3, 0.9)   => (:excludes_zero, :off, :keep),
             (-0.3, 0.3)  => (:includes_zero, :on,  :drop),
             (NaN, NaN)   => (:undecidable,   :undecidable, :undecidable),
             (0.2, 0.1)   => (:undecidable,   :undecidable, :undecidable)]
    for ((lo, hi), (status, act_on_include, act_on_exclude)) in table
        @test interval_excludes_zero(lo, hi) === status
        @test interval_code(lo, hi; excludes = :off,  includes = :on)   === act_on_include
        @test interval_code(lo, hi; excludes = :keep, includes = :drop) === act_on_exclude
    end

    # A MAPPING WITH A HOLE IS AN ERROR, whatever the interval. Both codes are required ...
    @test_throws UndefKeywordError interval_code(0.1, 0.2; excludes = :code_x)
    @test_throws UndefKeywordError interval_code(0.1, 0.2; includes = :code_y)
    # ... they must differ, or the interval decides nothing ...
    @test_throws ArgumentError interval_code(0.1, 0.2; excludes = :code_x, includes = :code_x)
    # ... neither may be the reserved third value, or uncomputed reads as decided ...
    @test_throws ArgumentError interval_code(0.1, 0.2; excludes = :undecidable, includes = :code_y)
    @test_throws ArgumentError interval_code(-0.1, 0.2; excludes = :code_x, includes = :undecidable)
    # ... and both are Symbols, not flags or text
    @test_throws ArgumentError interval_code(0.1, 0.2; excludes = true, includes = false)
    @test_throws ArgumentError interval_code(0.1, 0.2; excludes = "code_x", includes = :code_y)
    # the codes are checked before the interval, so an undecidable interval cannot hide a bad call
    @test_throws ArgumentError interval_code(NaN, NaN; excludes = :code_x, includes = :code_x)
    @test_throws ArgumentError interval_code(missing, missing; excludes = :code_x, includes = :undecidable)
end

@testset "sharing_wording_code — point, band, interval and the comparator override" begin
    sw(r, lo = NaN, hi = NaN; c = :excludes_zero) =
        sharing_wording_code(r, lo, hi; comparator_excludes_zero = c)

    # point form, clear of both bands
    @test sw(0.90) == (code = :large,   rule_used = :point)
    @test sw(0.50) == (code = :partial, rule_used = :point)
    @test sw(0.10) == (code = :none,    rule_used = :point)
    @test sw(1.0)  == (code = :large,   rule_used = :point)     # the whole effect is a share

    # the upper band: 0.70–0.80 takes the middle grade, from either side of the boundary
    @test sw(0.78) == (code = :partial, rule_used = :band)      # the point form says :large
    @test sw(0.75) == (code = :partial, rule_used = :band)
    @test sw(0.72) == (code = :partial, rule_used = :band)
    @test sw(0.70) == (code = :partial, rule_used = :band)      # both edges are inside
    @test sw(0.80) == (code = :partial, rule_used = :band)
    @test sw(nextfloat(0.80)) == (code = :large,   rule_used = :point)
    @test sw(prevfloat(0.70)) == (code = :partial, rule_used = :point)

    # the lower band: 0.20–0.30 takes no grade at all
    @test sw(0.28) == (code = :none, rule_used = :band)         # the point form says :partial
    @test sw(0.25) == (code = :none, rule_used = :band)
    @test sw(0.22) == (code = :none, rule_used = :band)
    @test sw(0.20) == (code = :none, rule_used = :band)
    @test sw(0.30) == (code = :none, rule_used = :band)
    @test sw(nextfloat(0.30)) == (code = :partial, rule_used = :point)
    @test sw(prevfloat(0.20)) == (code = :none,    rule_used = :point)

    # the interval form supersedes the band: an interval wholly above the boundary resolves
    # it, even though the point sits inside the band
    @test sw(0.78, 0.76, 0.79) == (code = :large, rule_used = :interval)
    # an interval straddling a boundary has not resolved it, and takes the weaker grade
    @test sw(0.90, 0.60, 0.95) == (code = :partial, rule_used = :interval)
    @test sw(0.40, 0.20, 0.60) == (code = :none,    rule_used = :interval)
    @test sw(0.50, 0.30, 0.70) == (code = :partial, rule_used = :interval)
    @test sw(0.10, 0.05, 0.20) == (code = :none,    rule_used = :interval)
    @test sw(0.90, 0.80, 1.0)  == (code = :large,   rule_used = :interval)   # one is inside
    # an interval decides even where the point ratio is missing
    @test sw(NaN, 0.30, 0.70) == (code = :partial, rule_used = :interval)

    # an interval that exists on one side only is no interval: the band and point govern
    @test sw(0.78, 0.76, NaN)     == (code = :partial, rule_used = :band)
    @test sw(0.90, NaN, 0.95)     == (code = :large,   rule_used = :point)
    @test sw(0.50, missing, 0.90) == (code = :partial, rule_used = :point)
    @test sw(0.50, 0.40, nothing) == (code = :partial, rule_used = :point)
    # and neither is a pair handed over the wrong way round
    @test sw(0.50, 0.90, 0.30)    == (code = :partial, rule_used = :point)

    # OUTSIDE THE DOMAIN. A bare threshold returns the strongest grade for each of these;
    # here each is a :none that says it refused, never a label a grading form also uses.
    @test sw(1.4)  == (code = :none, rule_used = :out_of_domain)   # comparator more enriched
    @test sw(-0.6) == (code = :none, rule_used = :out_of_domain)   # arms move in opposite directions
    @test sw(0.0)  == (code = :none, rule_used = :out_of_domain)   # comparator does not move
    # a reference effect near zero: 0.45 over 0.05
    near_zero_den = 0.45 / 0.05
    @test sw(near_zero_den) == (code = :none, rule_used = :out_of_domain)
    @test sw(near_zero_den, 0.9, 400.0) == (code = :none, rule_used = :out_of_domain)
    @test sw(Inf)  == (code = :none, rule_used = :out_of_domain)   # a denominator of exactly zero
    @test sw(-Inf) == (code = :none, rule_used = :out_of_domain)
    @test sw(0.95, 0.85, 1.30) == (code = :none, rule_used = :out_of_domain)
    @test sw(2.0, 1.5, 3.0)    == (code = :none, rule_used = :out_of_domain)
    # an interval reaching zero or below is outside the domain at its lower end, as a point
    # ratio at or below zero is
    @test sw(0.30, -0.4, 0.90) == (code = :none, rule_used = :out_of_domain)
    @test sw(0.60, 0.0, 0.90)  == (code = :none, rule_used = :out_of_domain)
    # an interval whose lower bound sits in the middle grade but whose upper bound runs past
    # one: the grading form would have said :partial, and the refusal is told apart from it
    res = sw(0.90, 0.72, 1.05)
    @test res == (code = :none, rule_used = :out_of_domain)
    @test res.rule_used !== :interval
    # one value outside is enough: a point beside an interval inside the domain, a lone bound
    @test sw(1.4, 0.30, 0.90)  == (code = :none, rule_used = :out_of_domain)
    @test sw(0.90, NaN, 1.30)  == (code = :none, rule_used = :out_of_domain)
    @test sw(NaN, -0.2, missing) == (code = :none, rule_used = :out_of_domain)
    @test sw(0.50, 1.30, 0.30) == (code = :none, rule_used = :out_of_domain)   # swapped, and outside
    # inside the domain a small fraction is still graded, by the form that graded it
    @test sw(0.10).rule_used === :point && sw(0.10, 0.05, 0.20).rule_used === :interval

    # nothing to grade at all
    @test sw(NaN)     == (code = :none, rule_used = :undecidable)
    @test sw(missing) == (code = :none, rule_used = :undecidable)
    @test sw(0.0 / 0.0) == (code = :none, rule_used = :undecidable)

    # THE COMPARATOR OVERRIDE beats every form, including an interval that says :large
    @test sw(0.90, 0.80, 0.95; c = :includes_zero) == (code = :none, rule_used = :comparator_null)
    @test sw(0.90; c = false) == (code = :none, rule_used = :comparator_null)
    @test sw(0.90; c = true)  == (code = :large, rule_used = :point)
    @test sw(0.90; c = interval_excludes_zero(-0.2, 0.4)) ==
          (code = :none, rule_used = :comparator_null)
    # ... and the domain refusal: a null comparator is the more basic fact
    @test sw(1.4, 1.1, 1.9; c = :includes_zero) == (code = :none, rule_used = :comparator_null)
    # a comparator whose own estimate was never computed is NOT reported as including zero
    for c in (:undecidable, missing, nothing, interval_excludes_zero(NaN, 0.4), :excludes)
        res = sw(0.90, 0.80, 0.95; c = c)
        @test res == (code = :none, rule_used = :undecidable)
        @test res.rule_used !== :comparator_null
        @test sw(1.4; c = c) == (code = :none, rule_used = :undecidable)
    end
    # the override cannot be forgotten: the keyword has no default
    @test_throws UndefKeywordError sharing_wording_code(0.9, NaN, NaN)

    # every answer stays inside the contracted vocabulary
    for r in (-1.0, 0.0, 0.1, 0.22, 0.5, 0.78, 0.9, 1.0, 1.5, NaN),
        (lo, hi) in ((NaN, NaN), (0.3, 0.7), (0.8, 0.9), (0.9, 400.0), (-0.2, 0.6)),
        c in (:excludes_zero, :includes_zero, :undecidable)
        res = sw(r, lo, hi; c = c)
        @test res.code in (:none, :partial, :large)
        @test res.rule_used in (:interval, :band, :point, :comparator_null, :undecidable,
                                :out_of_domain)
        # a refusal never carries a grade
        res.rule_used === :out_of_domain && @test res.code === :none
    end
end

@testset "depletion_tier — the strong tiers need the grid" begin
    depleted = (-1.0, -0.1)
    clean = ["variant A" => (-1.3, -0.2), "variant B" => (-1.1, -0.05), "variant C" => (-0.9, -0.01)]

    @test depletion_tier(depleted..., clean) == (tier = :t1, breakers = String[])

    # one variant reaching zero breaks it, and is named
    one = ["variant A" => (-1.3, -0.2), "variant B" => (-0.9, 0.03)]
    @test depletion_tier(depleted..., one) == (tier = :t2, breakers = ["variant B"])
    @test depletion_tier(depleted..., ["variant A" => (-1.0, 0.0)]).tier === :t2   # hi == 0
    # several breakers come back sorted and once each, whatever order they were given in
    many = ["variant Z" => (-0.5, 0.4), "variant A" => (-0.4, 0.2), "variant M" => (-1.0, -0.3),
            "variant Z" => (-0.6, 0.1)]
    @test depletion_tier(depleted..., many) == (tier = :t2, breakers = ["variant A", "variant Z"])
    # a variant excluding zero on the ENRICHED side has not reproduced a depletion
    @test depletion_tier(depleted..., ["variant E" => (0.1, 0.9)]) ==
          (tier = :t2, breakers = ["variant E"])
    # a variant whose interval could not be computed has not reproduced it either
    @test depletion_tier(depleted..., vcat(clean, ["variant N" => (NaN, NaN)])) ==
          (tier = :t2, breakers = ["variant N"])

    # T3 is decided by the primary alone: no grid needed, and none consulted
    for grid in (nothing, missing, Pair[], Dict{String,Tuple{Float64,Float64}}(), one)
        @test depletion_tier(-0.8, 0.6, grid) == (tier = :t3, breakers = String[])
    end
    @test depletion_tier(-1.0, 0.0, nothing).tier === :t3          # a bound reaching zero
    @test depletion_tier(0.2, 0.9, nothing).tier  === :t3          # enriched is not depleted
    # ... even over a grid of depleted variants, which would otherwise pass for robust
    @test depletion_tier(0.2, 0.9, clean) == (tier = :t3, breakers = String[])

    # THE REGRESSION. The primary excludes zero, but there is no grid: undecidable, and in
    # particular neither of the tiers the grid alone can establish.
    for grid in (nothing, missing, Pair[], Dict{String,Tuple{Float64,Float64}}(), NamedTuple())
        res = depletion_tier(depleted..., grid)
        @test res == (tier = :undecidable, breakers = String[])
        @test res.tier !== :t1 && res.tier !== :t2
    end

    # A SHORTER GRID THAN THE ONE SET OUT. A variant that failed to run and was dropped has
    # not agreed with the primary, so the strong tier is withheld while it is missing ...
    planned = ["variant A", "variant B", "variant C", "variant D"]
    @test depletion_tier(depleted..., clean) == (tier = :t1, breakers = String[])
    res = depletion_tier(depleted..., clean; expected_variants = planned)
    @test res == (tier = :undecidable, breakers = String[])
    @test res.tier !== :t1
    # ... but a breaker that is present decides the middle tier without the missing one
    @test depletion_tier(depleted..., one; expected_variants = planned) ==
          (tier = :t2, breakers = ["variant B"])
    # the complete grid, named with stray padding, and extras beyond the plan, still reach T1
    @test depletion_tier(depleted..., vcat(clean, ["variant X" => (-0.8, -0.3)]);
                       expected_variants = [" variant A", "variant B ", "variant C"]).tier === :t1
    # and the plan changes nothing that the primary alone decides
    @test depletion_tier(-0.8, 0.6, clean; expected_variants = planned).tier === :t3
    @test depletion_tier(depleted..., nothing; expected_variants = planned).tier === :undecidable

    # an undecidable primary is undecidable even over a clean grid
    @test depletion_tier(NaN, -0.1, clean) == (tier = :undecidable, breakers = String[])
    @test depletion_tier(missing, missing, clean).tier === :undecidable

    # the grid may arrive as a Dict, a NamedTuple, or with NamedTuple bounds
    @test depletion_tier(depleted..., Dict("variant A" => (-1.0, -0.2), "variant B" => (-0.5, 0.1))) ==
          (tier = :t2, breakers = ["variant B"])
    @test depletion_tier(depleted..., (variant_a = (-1.0, -0.2), variant_b = [-0.8, -0.1])).tier === :t1
    @test depletion_tier(depleted..., ["variant A" => (lo = -1.0, hi = -0.2)]).tier === :t1

    # a malformed grid is an error, not a guess
    @test_throws ErrorException depletion_tier(depleted..., [(-1.0, -0.2)])
    @test_throws ErrorException depletion_tier(depleted..., ["variant A" => (-1.0, -0.5, -0.2)])
    @test_throws ErrorException depletion_tier(depleted..., ["variant A" => (low = -1.0, high = -0.2)])
end

@testset "category_sets — new against the anchor, not against a constant" begin
    # Role names, not categories. The prefix says where a label stands — in both runs
    # ("kept"), in the primary run only ("new"), in the anchor only ("lost"), in neither
    # ("idle") — and the suffix whether another rule governs it ("g") or not ("u"). Every
    # case the partition distinguishes appears at least once.
    primary  = ["kept-u1", "kept-u2", "kept-g", "new-u1", "new-u2", "new-g"]
    anchor   = ["kept-u1", "kept-u2", "kept-g", "lost-u", "lost-g"]
    governed = ["kept-g", "new-g", "lost-g", "idle-g"]
    cs(p, a, g) = category_sets(p, a, g; primary_undecidable = String[])

    s = cs(primary, anchor, governed)
    @test s.anchor_available
    @test s.newly       == ["new-g", "new-u1", "new-u2"]
    @test s.legend_only == ["new-u1", "new-u2"]
    @test s.preexisting == ["kept-u1", "kept-u2"]
    @test s.no_longer   == ["lost-g", "lost-u"]
    @test s.undecidable_in_primary == String[]
    # an ungoverned label that excluded zero in the anchor too is carried over, never new
    for c in ("kept-u1", "kept-u2")
        @test !(c in s.legend_only)
        @test !(c in s.newly)
    end
    # a governed label is new when it is new, and no ungoverned list claims it
    @test "new-g" in s.newly
    @test !("new-g" in s.legend_only) && !("kept-g" in s.preexisting)
    # a governed label in neither run appears nowhere
    @test !any(v -> "idle-g" in v, (s.newly, s.preexisting, s.no_longer, s.legend_only))

    # the anchor measured against itself: by construction nothing is new and nothing lost
    self = cs(anchor, anchor, governed)
    @test self.newly == String[] && self.legend_only == String[] && self.no_longer == String[]
    @test self.preexisting == ["kept-u1", "kept-u2", "lost-u"]

    # input order, container type, duplicates, padding and blank labels change nothing —
    # a split of a joined key gives SubStrings, and a split of an empty key gives ""
    s2 = category_sets(Set(reverse(primary)), Tuple(anchor),
                       split(" kept-g;new-g;lost-g;;idle-g;new-g", ';');
                       primary_undecidable = split("", ';'))
    @test s2 == s
    # a blank label kept as a category would surface as a "new" category named ""
    blank = cs(split("", ';'), anchor, governed)
    @test blank.newly == String[]
    @test blank.no_longer == sort(anchor)
    for v in (s.newly, s.preexisting, s.no_longer, s.legend_only, s.undecidable_in_primary)
        @test issorted(v) && allunique(v) && eltype(v) === String
    end
    # a lone label is one label, not a string of characters
    @test cs("lone-u", String[], governed).newly == ["lone-u"]

    # positional destructuring still works for a caller that wants the four sets only
    newly, preexisting, no_longer, legend_only = s
    @test (newly, preexisting, no_longer, legend_only) ==
          (["new-g", "new-u1", "new-u2"], ["kept-u1", "kept-u2"], ["lost-g", "lost-u"],
           ["new-u1", "new-u2"])

    # AN UNAVAILABLE ANCHOR IS UNDECIDABLE: every set empty, and the caller is told so.
    # It must never be read as an empty anchor, which would declare the whole primary new.
    for missing_anchor in (nothing, missing)
        u = category_sets(primary, missing_anchor, governed; primary_undecidable = ["lost-u"])
        @test u.anchor_available == false
        @test u.newly == String[] && u.preexisting == String[] &&
              u.no_longer == String[] && u.legend_only == String[] &&
              u.undecidable_in_primary == String[]
    end
    # ... whereas an EMPTY anchor is a real answer: the earlier run had nothing, so all is new
    e = cs(primary, String[], governed)
    @test e.anchor_available
    @test e.newly == sort(primary)
    @test e.legend_only == ["kept-u1", "kept-u2", "new-u1", "new-u2"]
    @test e.preexisting == String[]
end

@testset "category_sets — a loss is a measurement, an uncomputed interval is not" begin
    primary  = ["kept-u", "new-u"]
    anchor   = ["kept-u", "lost-u", "lost-g", "unmeasured-u", "unmeasured-g"]
    governed = ["lost-g", "unmeasured-g"]
    # undecidable in the primary run: two anchor labels, and one the anchor never had
    undecidable = ["unmeasured-u", "unmeasured-g", "unmeasured-elsewhere"]

    s = category_sets(primary, anchor, governed; primary_undecidable = undecidable)
    @test s.anchor_available
    @test s.no_longer == ["lost-g", "lost-u"]
    @test s.undecidable_in_primary == ["unmeasured-g", "unmeasured-u"]
    # THE REGRESSION: no label whose primary interval was never computed is reported lost
    @test isempty(intersect(s.no_longer, undecidable))
    # the two lists partition what the anchor had and the primary set lacks, exactly
    @test isempty(intersect(s.no_longer, s.undecidable_in_primary))
    @test sort(vcat(s.no_longer, s.undecidable_in_primary)) == sort(setdiff(anchor, primary))
    # an undecidable label the anchor never had was not lost, so neither list names it
    @test !("unmeasured-elsewhere" in s.undecidable_in_primary)
    @test !("unmeasured-elsewhere" in s.no_longer)
    # the undecidable list moves nothing else
    base = category_sets(primary, anchor, governed; primary_undecidable = String[])
    @test (s.newly, s.preexisting, s.legend_only) == (base.newly, base.preexisting, base.legend_only)
    # and with no undecidable interval the whole remainder is a measured loss
    @test base.no_longer == ["lost-g", "lost-u", "unmeasured-g", "unmeasured-u"]
    @test base.undecidable_in_primary == String[]

    # a label that excluded zero on one unit and is undecidable on another was measured
    # excluding zero, so it stays with the primary set and is neither lost nor unmeasured
    mixed = category_sets(primary, anchor, governed; primary_undecidable = ["kept-u", "new-u"])
    @test mixed.preexisting == ["kept-u"] && mixed.newly == ["new-u"]
    @test mixed.undecidable_in_primary == String[]
    @test mixed.no_longer == base.no_longer

    # a caller that cannot say which intervals were computed gets no loss at all
    for unknown in (nothing, missing)
        u = category_sets(primary, anchor, governed; primary_undecidable = unknown)
        @test u.anchor_available
        @test u.no_longer == String[]
        @test u.undecidable_in_primary == sort(setdiff(anchor, primary))
        @test (u.newly, u.preexisting, u.legend_only) == (base.newly, base.preexisting, base.legend_only)
    end

    # the keyword cannot be forgotten: it has no default
    @test_throws UndefKeywordError category_sets(primary, anchor, governed)

    # a lone padded label, SubStrings and duplicates are labels like any other
    @test category_sets(primary, anchor, governed;
                        primary_undecidable = " unmeasured-u ").undecidable_in_primary == ["unmeasured-u"]
    @test category_sets(primary, anchor, governed;
                        primary_undecidable = split("unmeasured-g;;unmeasured-g", ';')) ==
          merge(base, (no_longer = ["lost-g", "lost-u", "unmeasured-u"],
                       undecidable_in_primary = ["unmeasured-g"]))
    for v in (s.no_longer, s.undecidable_in_primary)
        @test issorted(v) && allunique(v) && eltype(v) === String
    end
end

@testset "recall_lower_bound — the exact one-sided bound, not the rate" begin
    # KNOWN VALUES. Recovering nothing is consistent with a recall of zero, and one of one is
    # bounded at the level itself.
    @test recall_lower_bound(0, 1) === 0.0
    @test recall_lower_bound(0, 37) === 0.0
    @test recall_lower_bound(1, 1) ≈ 0.05 atol = 1e-15
    @test recall_lower_bound(1, 1; alpha = 0.2) ≈ 0.2 atol = 1e-15
    @test recall_lower_bound(1, 1) isa Float64 && recall_lower_bound(22, 37) isa Float64
    # an empty sample bounds nothing, at any level
    @test isnan(recall_lower_bound(0, 0))
    @test isnan(recall_lower_bound(0, 0; alpha = 0.3))

    # THE CLOSED FORM AGAINST THE BISECTION, where both apply
    for n in (1, 2, 5, 12, 29, 64, 333, 4096), alpha in (0.01, 0.05, 0.2)
        @test recall_lower_bound(n, n; alpha = alpha) == alpha^(1 / n)
        @test OutcomeRules._bisect_lower_bound(n, n, alpha) ≈ alpha^(1 / n) atol = 1e-11
    end

    # A SECOND ROUTE for k < n. Each literal is scipy.stats.beta.ppf(alpha, k, n - k + 1), the
    # beta-quantile form of the same bound, computed outside this suite; each also agrees to
    # 1e-15 with a 50-digit bisection on the binomial tail. The pairs are arbitrary. The last
    # has n large enough that the log-factorials come from the series branch.
    @test recall_lower_bound(22, 37)             ≈ 0.4463441121359716  atol = 1e-11
    @test recall_lower_bound(3, 25; alpha = 0.1) ≈ 0.04491346394064166 atol = 1e-11
    @test recall_lower_bound(9990, 10000)        ≈ 0.9983043685463546  atol = 1e-11

    # MONOTONE THE WAY A BOUND MUST BE. More recovered out of the same sample raises it; the
    # same number recovered out of a larger sample lowers it; a stricter level lowers it; and
    # at the default level it never exceeds the point estimate. The k-to-n step crosses from
    # the bisection to the closed form, so the two branches are also checked against each other.
    for n in 1:40
        bounds = [recall_lower_bound(k, n) for k in 0:n]
        @test all(diff(bounds) .> 0)
        @test all(bounds[k + 1] <= k / n for k in 0:n)
        @test all(b -> 0 <= b < 1, bounds)
        @test all(recall_lower_bound(k, n + 1) < bounds[k + 1] for k in 1:n)
        @test all(recall_lower_bound(k, n; alpha = 0.01) < bounds[k + 1] for k in 1:n)
    end

    # THE SIZE A CLEAN SAMPLE NEEDS. With no miss the bound is alpha^(1/n), so a sample with no
    # miss reaches a bar only from the size where alpha^(1/n) does; one taxon fewer does not.
    for bar in (0.6, 0.9, 0.95)
        n_min = ceil(Int, log(0.05) / log(bar))
        @test recall_lower_bound(n_min, n_min) >= bar
        @test recall_lower_bound(n_min - 1, n_min - 1) < bar
    end

    # The log-factorials every tail term starts from are exact to double precision on both
    # sides of the switch from direct summation to the series. A dropped correction term moves
    # them by more than this just past the switch, while the bounds above would not notice.
    @test all(abs(OutcomeRules._logfact(m) - Float64(sum(log, big.(2:max(m, 1)); init = big(0.0)))) < 5e-14
              for m in 0:40)

    # a large sample is summed where the tail's mass is, not term by term over all of it
    # (it sits below the rate by less than three standard errors of that rate)
    big_rate = 999_000 / 1_000_000
    @test 0 < big_rate - recall_lower_bound(999_000, 1_000_000) <
          3 * sqrt(big_rate * (1 - big_rate) / 1_000_000)
    @test recall_lower_bound(big(22), big(37)) == recall_lower_bound(22, 37)

    # A CALL THAT DESCRIBES NO SAMPLE IS AN ERROR, not a number
    @test_throws ArgumentError recall_lower_bound(5, 4)             # more recovered than flagged
    @test_throws ArgumentError recall_lower_bound(-1, 4)
    @test_throws ArgumentError recall_lower_bound(0, -1)
    @test_throws ArgumentError recall_lower_bound(true, true)       # a flag is not a count
    @test_throws ArgumentError recall_lower_bound(1, 2^60)          # beyond exact double arithmetic
    for alpha in (0.0, 1.0, -0.1, 1.5, NaN)
        @test_throws ArgumentError recall_lower_bound(1, 2; alpha = alpha)
    end
    @test_throws MethodError recall_lower_bound(2.0, 3)             # counts are integers
end

@testset "scan_evidence_status — the records must vouch for the table read" begin
    sha_read  = repeat("0123456789abcdef", 4)     # the digest of the table that was read
    sha_other = repeat("fedcba9876543210", 4)     # the digest of some other table
    # The recall fields of a faithful record of k recovered out of n flagged: the rate and the
    # bound as its counts give them, the bar in both records that carry it, and a pass field
    # taken on the point estimate, the way the scan's own gate takes it.
    recall_of(k, n; bar = 0.8, pass = n > 0 && k / n >= bar) =
        (recall_n_full_flagged = n, recall_n_recovered = k,
         recall_rate_field = n == 0 ? NaN : k / n,
         recall_lower_bound_field = recall_lower_bound(k, n),
         recall_bar = bar, bar_in_other_records = ["acceptance record" => bar],
         recall_pass_field = pass)
    good = merge((scan_sha256 = sha_read,
                  record_sha256s = Dict("threshold" => sha_read, "recall" => sha_read,
                                        "acceptance" => [sha_read, sha_read, sha_read]),
                  acceptance_checks = ["gate-a", "gate-b", "note-c"],
                  acceptance_verdicts = ["PASS", "PASS", "DISCLOSED"],
                  required_checks = ["gate-a", "gate-b"],
                  flagged_row_value = 7.0, table_n_flagged = 7,
                  family_row_denominator = 5.0, table_n_family = 5,
                  recall_n_sampled = 90,
                  recall_threshold_bits = 12.5, recall_uniprot_release = "release-x",
                  threshold_bits_in_other_records = ["threshold record" => 12.5],
                  uniprot_release_in_other_records = ["threshold record" => "release-x"]),
                 recall_of(38, 40))
    st(kw = (;)) = scan_evidence_status(; merge(good, kw)...)
    one(key, value) = st(NamedTuple{(key,)}((value,)))

    ok = st()
    @test ok.ok === true && ok.recall_pass_derived === true && ok.reasons == String[]
    @test keys(ok) == (:ok, :recall_pass_derived, :reasons, :recall_rate, :recall_lower_bound)
    # the rate and bound handed back are the counts', for the caller to report
    @test ok.recall_rate == 38 / 40
    @test ok.recall_lower_bound == recall_lower_bound(38, 40)

    # the records may arrive keyed, as pairs, as a bare list or as one digest; digest case
    # and a CRLF left on the last column of a record are not part of the digest
    for recs in ((threshold = sha_read, recall = sha_read, acceptance = sha_read),
                 ["threshold" => sha_read, "recall" => sha_read],
                 [sha_read, sha_read, sha_read], (sha_read, sha_read), sha_read,
                 [uppercase(sha_read), " $(sha_read)\r"])
        @test st((record_sha256s = recs,)).ok
    end

    # 1. A RECORD FROM ANOTHER RUN. One foreign digest is enough, and the reason names it.
    stale = st((record_sha256s = Dict("threshold" => sha_read, "recall" => sha_other,
                                      "acceptance" => sha_read),))
    @test stale.ok == false
    @test only_reason(stale, "record_sha256s")
    @test occursin("recall", stale.reasons[1])
    @test stale.recall_pass_derived === true          # the recall counts are not what failed
    # ... including a single row of a record whose rows each carry a digest
    @test only_reason(st((record_sha256s = Dict("acceptance" => [sha_read, sha_other, sha_read]),)),
                      "record_sha256s")
    # records that agree with one another but not with the table read are foreign too
    @test st((record_sha256s = [sha_other, sha_other, sha_other],)).ok == false
    moved = st((scan_sha256 = sha_other,))
    @test moved.ok == false
    @test reason_keys(moved) == fill("record_sha256s", 3)

    # an unknown or malformed digest for the table read is its own reason
    for scan in (nothing, missing, "", first(sha_read, 12), sha_read * "0",
                 replace(sha_read, "a" => "g"), 12345, Symbol(sha_read))
        res = st((scan_sha256 = scan,))
        @test res.ok == false
        @test only_reason(res, "scan_sha256")
    end
    # ... and so is a record digest that is absent, empty, truncated or not text
    for recs in (nothing, missing, String[], Dict{String,String}(), NamedTuple(),
                 Dict("recall" => nothing), Dict("recall" => ""), Dict("recall" => first(sha_read, 40)),
                 Dict("acceptance" => String[]), Dict("acceptance" => Any[sha_read, "", sha_read]),
                 Dict("acceptance" => Any[sha_read, missing]), Any[sha_read, missing], [sha_read, 7])
        res = st((record_sha256s = recs,))
        @test res.ok == false
        @test only_reason(res, "record_sha256s")
    end

    # 2. A REFUSED TABLE. One MISS among passes withholds, and it is the only reason given.
    miss = st((acceptance_verdicts = ["PASS", "MISS", "DISCLOSED"],))
    @test miss.ok == false
    @test only_reason(miss, "acceptance_verdicts")
    @test miss.recall_pass_derived === true
    # no verdicts, a blank verdict, a lower-case one, one outside the vocabulary: none is a pass
    for v in (nothing, missing, String[], ["PASS", "PASS", ""], ["PASS", "pass", "DISCLOSED"],
              ["FAIL", "PASS", "DISCLOSED"], Any["PASS", "PASS", missing],
              Any["PASS", nothing, "DISCLOSED"])
        res = st((acceptance_verdicts = v,))
        @test res.ok == false
        @test only_reason(res, "acceptance_verdicts")
    end
    for v in ("OK", 1)
        res = st((acceptance_checks = "gate-a", acceptance_verdicts = v, required_checks = "gate-a"))
        @test res.ok == false && only_reason(res, "acceptance_verdicts")
    end
    # a MISS and an unknown verdict are two failures, and both are named
    @test length(st((acceptance_verdicts = ["MISS", "?", "DISCLOSED"],)).reasons) == 2
    # a lone row, padding and Symbols are checks and verdicts like any other
    @test st((acceptance_checks = "gate-a", acceptance_verdicts = "PASS", required_checks = "gate-a")).ok
    @test st((acceptance_checks = ["gate-a ", " gate-b", "note-c"],
              acceptance_verdicts = [" PASS", "PASS", "DISCLOSED\r"])).ok
    @test st((acceptance_verdicts = [:PASS, :PASS, :DISCLOSED],)).ok

    # A GATE LEFT OUT. A table with its gated rows removed has no MISS left to find, so each of
    # these would pass on verdicts alone; each is refused, on its own, for the gate it lacks.
    # Disclosures only, which pass on verdicts alone:
    disclosures = st((acceptance_checks = ["note-c"], acceptance_verdicts = ["DISCLOSED"]))
    @test disclosures.ok == false
    @test reason_keys(disclosures) == ["acceptance_checks", "acceptance_checks"]
    # one gated row cut out, and the reason names it
    absent = st((acceptance_checks = ["gate-a", "note-c"], acceptance_verdicts = ["PASS", "DISCLOSED"]))
    @test absent.ok == false && only_reason(absent, "acceptance_checks")
    @test occursin("gate-b", absent.reasons[1])
    # a gated row present twice: which of its verdicts stands is not for the reader to pick
    twice = st((acceptance_checks = ["gate-a", "gate-b", "gate-b", "note-c"],
                acceptance_verdicts = ["PASS", "PASS", "PASS", "DISCLOSED"]))
    @test twice.ok == false && only_reason(twice, "acceptance_checks")
    # a gated row demoted to a disclosure, whose bar was therefore never applied
    demoted = st((acceptance_verdicts = ["PASS", "DISCLOSED", "DISCLOSED"],))
    @test demoted.ok == false && only_reason(demoted, "acceptance_checks")
    # check names that do not line up with the verdicts tie no verdict to any check
    @test only_reason(st((acceptance_checks = ["gate-a", "gate-b"],)), "acceptance_checks")
    # names are matched whole and exactly, after padding; a repeated required name is one name
    @test st((required_checks = [" gate-a", "gate-b ", "gate-a"],)).ok
    @test only_reason(st((required_checks = ["gate", "gate-b"],)), "acceptance_checks")
    @test only_reason(st((required_checks = ["GATE-A", "gate-b"],)), "acceptance_checks")
    # no check names, and no required check at all, are reasons of their own
    for checks in (nothing, missing, String[])
        res = st((acceptance_checks = checks,))
        @test res.ok == false && only_reason(res, "acceptance_checks")
    end
    for req in (nothing, missing, String[], "", ["", " "])
        res = st((required_checks = req,))
        @test res.ok == false && only_reason(res, "required_checks")
    end

    # 3. THE ACCEPTANCE TABLE AGAINST THE TABLE READ. A record tied to the table by its digest
    # can still carry counts the table does not have; each side of each comparison fails alone.
    @test only_reason(st((flagged_row_value = 6.0,)), "flagged_row_value")
    @test only_reason(st((table_n_flagged = 8,)), "flagged_row_value")
    @test only_reason(st((family_row_denominator = 4.0,)), "family_row_denominator")
    @test only_reason(st((table_n_family = 6,)), "family_row_denominator")
    @test st((flagged_row_value = 7, family_row_denominator = 5)).ok      # an Int is a count too
    for (key, bad) in ((:flagged_row_value, nothing), (:flagged_row_value, NaN),
                       (:flagged_row_value, 6.5), (:table_n_flagged, missing),
                       (:table_n_flagged, "7"), (:family_row_denominator, -5.0),
                       (:table_n_family, true), (:table_n_family, Inf))
        res = one(key, bad)
        @test res.ok == false && only_reason(res, String(key))
    end

    # A WHOLE NUMBER TOO LARGE TO BE A COUNT is unknown, with a reason, and never an
    # InexactError: whole-valued floats beyond Int range, and beyond exact double arithmetic
    for (key, huge) in ((:recall_n_full_flagged, 1e20), (:recall_n_recovered, 2.0^63),
                        (:recall_n_sampled, 1e300), (:flagged_row_value, 1e19),
                        (:table_n_family, big(10)^30), (:table_n_flagged, 2.0^53 + 2))
        res = one(key, huge)
        @test res.ok == false
        @test only_reason(res, String(key))
    end

    # 4. THE RECALL RECORD AGAINST ITS OWN COUNTS. A rate that disagrees with the counts it
    # summarises, every digest and verdict intact: refused, although the counts pass.
    wrong_rate = st((recall_rate_field = 0.5,))
    @test wrong_rate.ok == false && only_reason(wrong_rate, "recall_rate_field")
    @test wrong_rate.recall_pass_derived === true
    @test wrong_rate.recall_rate == 38 / 40                # the counts' rate, not the record's
    # A rate is written to six decimals, so an honest one sits within half a unit of its last
    # decimal whichever way the half was rounded. Here the rounding alone puts it a hair past
    # the half unit once the decimal is read back, which is what the parsing slack is for.
    @test abs(0.007813 - 1 / 128) > 5e-7
    for written in (0.007812, 0.007813)
        @test st(merge(recall_of(1, 128; bar = 1e-4), (recall_rate_field = written,
                                                       recall_n_sampled = 200))).ok
    end
    @test only_reason(st(merge(recall_of(1, 128; bar = 1e-4), (recall_rate_field = 0.007811,
                                                               recall_n_sampled = 200))),
                      "recall_rate_field")
    # the bound the record writes agrees to 1e-9, and not beyond
    lb = recall_lower_bound(38, 40)
    @test st((recall_lower_bound_field = lb + 0.5e-9,)).ok
    @test only_reason(st((recall_lower_bound_field = lb + 2e-9,)), "recall_lower_bound_field")
    @test only_reason(st((recall_lower_bound_field = 0.95,)), "recall_lower_bound_field")  # the rate, in the bound's place
    # NaN beside counts that flagged taxa is a disagreement, for either field
    @test only_reason(st((recall_rate_field = NaN,)), "recall_rate_field")
    @test only_reason(st((recall_lower_bound_field = NaN,)), "recall_lower_bound_field")
    # over a sample that flagged nothing, NaN is right and a number is not
    empty_kw = recall_of(0, 0)
    @test reason_keys(st(empty_kw)) == ["recall_n_full_flagged"]
    @test reason_keys(st(merge(empty_kw, (recall_rate_field = 1.0,)))) ==
          ["recall_n_full_flagged", "recall_rate_field"]
    @test reason_keys(st(merge(empty_kw, (recall_lower_bound_field = 0.0,)))) ==
          ["recall_lower_bound_field", "recall_n_full_flagged"]
    for (key, bad) in ((:recall_rate_field, nothing), (:recall_rate_field, "0.95"),
                       (:recall_rate_field, true), (:recall_lower_bound_field, missing),
                       (:recall_lower_bound_field, "0.85"))
        res = one(key, bad)
        @test res.ok == false && only_reason(res, String(key))
    end
    # more taxa flagged than the sample held; all of it flagged is possible
    @test only_reason(st((recall_n_sampled = 39,)), "recall_n_sampled")
    @test st((recall_n_sampled = 40,)).ok
    for bad in (nothing, -1, 90.5, "90")
        @test only_reason(one(:recall_n_sampled, bad), "recall_n_sampled")
    end

    # THE BAR, THE CUTOFF AND THE RELEASE AGAINST THE OTHER RECORDS. A bar lowered in the recall
    # record, beside counts that fail the bar the acceptance table holds: under its own bar the
    # record passes, and it is refused because that bar is not the other record's.
    lowered = st(merge(recall_of(20, 40; bar = 0.3), (bar_in_other_records = ["acceptance record" => 0.8],)))
    @test recall_lower_bound(20, 40) >= 0.3 && 20 / 40 < 0.8
    @test lowered.ok == false && only_reason(lowered, "bar_in_other_records")
    @test only_reason(st((bar_in_other_records = ["acceptance record" => 0.8 + 1e-3],)),
                      "bar_in_other_records")
    @test st((bar_in_other_records = ["acceptance record" => 0.8 + 4e-7],)).ok     # written precision
    @test only_reason(st((threshold_bits_in_other_records = ["threshold record" => 12.6],)),
                      "threshold_bits_in_other_records")
    @test st((threshold_bits_in_other_records = ["threshold record" => 12.50004],)).ok
    @test only_reason(st((uniprot_release_in_other_records = ["threshold record" => "release-y"],)),
                      "uniprot_release_in_other_records")
    @test st((uniprot_release_in_other_records = (threshold = " release-x\r", acceptance = "release-x"),)).ok
    # every record named must agree, not just one of them
    @test only_reason(st((threshold_bits_in_other_records = ["threshold record" => 12.5,
                                                             "acceptance record" => 13.0],)),
                      "threshold_bits_in_other_records")
    # no other record, or one carrying nothing usable, has checked nothing
    for (key, bad) in ((:bar_in_other_records, nothing), (:bar_in_other_records, Pair[]),
                       (:bar_in_other_records, ["acceptance record" => NaN]),
                       (:threshold_bits_in_other_records, missing),
                       (:threshold_bits_in_other_records, Dict("threshold record" => "12.5")),
                       (:uniprot_release_in_other_records, ["threshold record" => ""]),
                       (:uniprot_release_in_other_records, NamedTuple()))
        res = one(key, bad)
        @test res.ok == false && only_reason(res, String(key))
    end
    # and the recall record's own copy has to be there to be compared
    for (key, bad) in ((:recall_threshold_bits, nothing), (:recall_threshold_bits, "12.5"),
                       (:recall_uniprot_release, " "), (:recall_uniprot_release, missing))
        res = one(key, bad)
        @test res.ok == false && only_reason(res, String(key))
    end

    # 5. THE RECALL VERDICT, ON THE LOWER BOUND.
    # A SAMPLE THAT MISSED NOTHING BUT IS TOO SMALL FOR ITS BOUND TO REACH THE BAR. The point
    # estimate is one and clears the bar; the bound does not. The recall has not passed, and
    # the record, honest about its point estimate, is refused for the bound alone.
    @test recall_lower_bound(12, 12) < 0.9 <= 12 / 12
    clean_small = st(recall_of(12, 12; bar = 0.9))
    @test clean_small.ok == false
    @test clean_small.recall_pass_derived === false
    @test only_reason(clean_small, "recall_n_recovered")
    @test clean_small.recall_rate == 1.0
    @test clean_small.recall_lower_bound == recall_lower_bound(12, 12)
    # a pass field set to agree with the bound instead contradicts the point estimate it is
    # taken on, so it is named as well
    @test reason_keys(st(recall_of(12, 12; bar = 0.9, pass = false))) ==
          ["recall_n_recovered", "recall_pass_field"]
    # the same counts under a bar their bound reaches pass
    @test recall_lower_bound(12, 12) >= 0.75
    clean_ok = st(recall_of(12, 12; bar = 0.75))
    @test clean_ok.ok && clean_ok.recall_pass_derived === true
    # ... and the higher bar is reached by a clean sample large enough for it
    @test st(recall_of(29, 29; bar = 0.9)).ok
    # the bar is inclusive on the bound: a bar exactly at it passes, the next float up does not
    at_bound = recall_lower_bound(27, 30)
    @test st(recall_of(27, 30; bar = at_bound)).ok
    above = st(recall_of(27, 30; bar = nextfloat(at_bound)))
    @test above.recall_pass_derived === false && only_reason(above, "recall_n_recovered")

    # A RECORD THAT MEASURED NOTHING BUT SAYS PASS. No flagged taxon, nothing to recover: the
    # derived verdict is false whatever the field reads, and the disagreement is named.
    hollow = st(recall_of(0, 0; pass = true))
    @test hollow.ok == false
    @test hollow.recall_pass_derived === false
    @test isnan(hollow.recall_rate) && isnan(hollow.recall_lower_bound)
    @test reason_keys(hollow) == ["recall_n_full_flagged", "recall_pass_field"]
    # ... and saying fail about it is honest, but it still measured nothing
    hollow_honest = st(recall_of(0, 0; pass = false))
    @test hollow_honest.ok == false && only_reason(hollow_honest, "recall_n_full_flagged")
    # A RECORD WHOSE COUNTS CONTRADICT ITS VERDICT
    contra = st(recall_of(14, 30; pass = true))
    @test contra.ok == false
    @test contra.recall_pass_derived === false
    @test reason_keys(contra) == ["recall_n_recovered", "recall_pass_field"]
    # an honest failing record is a failure too
    fails = st(recall_of(14, 30))
    @test fails.ok == false && fails.recall_pass_derived === false
    @test only_reason(fails, "recall_n_recovered")
    # and counts that pass beside a field that says fail disagree just the same
    modest = st((recall_pass_field = false,))
    @test modest.ok == false && modest.recall_pass_derived === true
    @test only_reason(modest, "recall_pass_field")
    # whole counts stored as floats are counts
    @test st((recall_n_full_flagged = 40.0, recall_n_recovered = 38.0)).ok
    # a record recovering more than it flagged is not a rate of anything
    over = st((recall_n_full_flagged = 10, recall_n_recovered = 12))
    @test over.ok == false && over.recall_pass_derived === nothing
    @test over.recall_rate === nothing && over.recall_lower_bound === nothing
    @test only_reason(over, "recall_n_recovered")

    # UNKNOWN IS NOT CLEAN, for each recall input on its own; the derived verdict is withheld
    # whenever an input it needs is unknown. A bar that is known but out of range is given in
    # both records, so that the range is the only thing wrong with it.
    for (kw, key, derivable) in
            (((recall_n_full_flagged = nothing,), "recall_n_full_flagged", false),
             ((recall_n_full_flagged = missing,), "recall_n_full_flagged", false),
             ((recall_n_full_flagged = "40",),    "recall_n_full_flagged", false),
             ((recall_n_full_flagged = 40.5,),    "recall_n_full_flagged", false),
             ((recall_n_full_flagged = -1,),      "recall_n_full_flagged", false),
             ((recall_n_full_flagged = true,),    "recall_n_full_flagged", false),
             ((recall_n_full_flagged = NaN,),     "recall_n_full_flagged", false),
             ((recall_n_recovered = nothing,),    "recall_n_recovered",    false),
             ((recall_n_recovered = Inf,),        "recall_n_recovered",    false),
             ((recall_bar = nothing,),            "recall_bar",            false),
             ((recall_bar = NaN,),                "recall_bar",            false),
             ((recall_bar = "0.9",),              "recall_bar",            false),
             ((recall_bar = 0.0,  bar_in_other_records = 0.0),  "recall_bar", false),  # passes anything
             ((recall_bar = -0.5, bar_in_other_records = -0.5), "recall_bar", false),
             ((recall_bar = 1.5,  bar_in_other_records = 1.5),  "recall_bar", false),  # passes nothing
             ((recall_pass_field = nothing,),     "recall_pass_field",     true),
             ((recall_pass_field = "true",),      "recall_pass_field",     true),
             ((recall_pass_field = 1,),           "recall_pass_field",     true))
        res = st(kw)
        @test res.ok == false
        @test only_reason(res, key)
        @test (res.recall_pass_derived !== nothing) == derivable
    end

    # leaving any one keyword out cannot earn `ok`
    for key in keys(good)
        kept = NamedTuple{Tuple(k for k in keys(good) if k !== key)}(good)
        @test scan_evidence_status(; kept...).ok == false
    end
    none_given = scan_evidence_status()
    @test none_given.ok == false
    @test none_given.recall_pass_derived === nothing
    @test reason_keys(none_given) == sort(["scan_sha256", "record_sha256s",
        "acceptance_verdicts", "acceptance_checks", "required_checks",
        "flagged_row_value", "table_n_flagged", "family_row_denominator", "table_n_family",
        "recall_n_sampled", "recall_n_full_flagged", "recall_n_recovered", "recall_bar",
        "recall_rate_field", "recall_lower_bound_field", "recall_threshold_bits",
        "recall_uniprot_release", "bar_in_other_records", "threshold_bits_in_other_records",
        "uniprot_release_in_other_records", "recall_pass_field"])

    # every failure is reported, not just the first one found
    all_bad = st((scan_sha256 = sha_other, acceptance_verdicts = ["MISS", "PASS", "DISCLOSED"],
                  required_checks = ["gate-z"], flagged_row_value = 6.0,
                  recall_n_full_flagged = 0, recall_n_recovered = 0))
    @test unique(reason_keys(all_bad)) == ["acceptance_checks", "acceptance_verdicts",
        "flagged_row_value", "recall_lower_bound_field", "recall_n_full_flagged",
        "recall_pass_field", "recall_rate_field", "record_sha256s"]
end

@testset "background_removal_complete — every input can withdraw it" begin
    clean = (left_hmm_hit = 0, left_ipr_ge500 = 0, recall_pass = true, frame_matches = true,
             scan_evidence_ok = true)
    ok = background_removal_complete(; clean...)
    @test ok.complete
    @test ok.reasons == String[]
    @test background_licence_code(ok) === :strong
    @test background_removal_complete(; merge(clean, (left_hmm_hit = 0.0, left_ipr_ge500 = 0.0))...).complete

    spoil(kw) = background_removal_complete(; merge(clean, kw)...)

    # a left-over family-signature carrier in the analysed arm
    ipr = spoil((left_ipr_ge500 = 3,))
    @test ipr.complete == false
    @test only_reason(ipr, "left_ipr_ge500")
    @test occursin("3", ipr.reasons[1])
    @test background_licence_code(ipr) === :fallback

    # a left-over HMM-hit carrier
    hmm = spoil((left_hmm_hit = 2,))
    @test hmm.complete == false
    @test only_reason(hmm, "left_hmm_hit")
    # a count too large to print as an integer is still a reason, not an InexactError raised
    # while the reason is being written
    huge = spoil((left_hmm_hit = 1e20,))
    @test huge.complete == false && only_reason(huge, "left_hmm_hit")

    # a recall check that did not pass
    rec = spoil((recall_pass = false,))
    @test rec.complete == false
    @test only_reason(rec, "recall_pass")
    @test background_licence_code(rec) === :fallback

    # a scan that does not cover the analysed frame
    frm = spoil((frame_matches = false,))
    @test frm.complete == false
    @test only_reason(frm, "frame_matches")
    @test background_licence_code(frm) === :fallback

    # a scan table its own persisted records do not vouch for
    evd = spoil((scan_evidence_ok = false,))
    @test evd.complete == false
    @test only_reason(evd, "scan_evidence_ok")
    @test background_licence_code(evd) === :fallback

    # UNKNOWN IS NOT CLEAN, for each input on its own
    for key in (:left_hmm_hit, :left_ipr_ge500, :recall_pass, :frame_matches, :scan_evidence_ok),
        unknown in (nothing, missing)
        res = spoil(NamedTuple{(key,)}((unknown,)))
        @test res.complete == false
        @test only_reason(res, String(key))
        @test background_licence_code(res) === :fallback
    end
    # a count that is not a count, and a verdict that is not a Bool — including the whole
    # status tuple handed over where only its `ok` belongs
    for kw in ((left_hmm_hit = NaN,), (left_ipr_ge500 = -1,), (left_hmm_hit = false,),
               (left_ipr_ge500 = "0",), (recall_pass = 1,), (recall_pass = "true",),
               (frame_matches = 1.0,), (scan_evidence_ok = "true",),
               (scan_evidence_ok = (ok = true, recall_pass_derived = true, reasons = String[]),))
        res = spoil(kw)
        @test res.complete == false
        @test only_reason(res, String(only(keys(kw))))
    end
    # leaving an input out cannot earn the licence
    @test background_removal_complete(left_hmm_hit = 0, left_ipr_ge500 = 0, recall_pass = true,
                                      frame_matches = true).complete == false
    @test background_removal_complete(left_hmm_hit = 0, recall_pass = true, frame_matches = true,
                                      scan_evidence_ok = true).complete == false
    none_given = background_removal_complete()
    @test none_given.complete == false
    @test length(none_given.reasons) == 5

    # every failure is reported, not just the first one found
    all_bad = background_removal_complete(left_hmm_hit = 7, left_ipr_ge500 = 11,
                                          recall_pass = false, frame_matches = false,
                                          scan_evidence_ok = false)
    @test reason_keys(all_bad) ==
          ["frame_matches", "left_hmm_hit", "left_ipr_ge500", "recall_pass", "scan_evidence_ok"]
end

@testset "end to end over the pure guards — each construction withholds the strong code" begin
    # The chain a consumer runs: the scan's records into scan_evidence_status, its verdicts
    # into background_removal_complete beside the counts, the result into the licence code.
    sha_read  = repeat("5a", 32)
    sha_other = repeat("c3", 32)
    recall_of(k, n; bar = 0.8) =
        (recall_n_full_flagged = n, recall_n_recovered = k,
         recall_rate_field = n == 0 ? NaN : k / n,
         recall_lower_bound_field = recall_lower_bound(k, n),
         recall_bar = bar, bar_in_other_records = (acceptance = bar,),
         recall_pass_field = n > 0 && k / n >= bar)
    evidence(kw = (;)) = scan_evidence_status(; merge(
        (scan_sha256 = sha_read,
         record_sha256s = (threshold = sha_read, recall = sha_read, acceptance = sha_read),
         acceptance_checks = ["gate-p", "gate-q", "gate-r", "note-s"],
         acceptance_verdicts = ["PASS", "PASS", "PASS", "DISCLOSED"],
         required_checks = ["gate-p", "gate-q", "gate-r"],
         flagged_row_value = 9.0, table_n_flagged = 9,
         family_row_denominator = 4.0, table_n_family = 4,
         recall_n_sampled = 66, recall_threshold_bits = 30.25, recall_uniprot_release = "rel-1",
         threshold_bits_in_other_records = (threshold = 30.25,),
         uniprot_release_in_other_records = (threshold = "rel-1",)),
        recall_of(24, 25), kw)...)
    removal(ev; kw...) = background_removal_complete(; merge(
        (left_hmm_hit = 0, left_ipr_ge500 = 0, recall_pass = ev.recall_pass_derived,
         frame_matches = true, scan_evidence_ok = ev.ok), (; kw...))...)
    licence(ev; kw...) = background_licence_code(removal(ev; kw...))
    # the construction a consumer that forwarded the record's own pass field would run
    forwarded(ev, field) = background_licence_code(background_removal_complete(left_hmm_hit = 0,
        left_ipr_ge500 = 0, recall_pass = field, frame_matches = true, scan_evidence_ok = ev.ok))

    # the control: honest records for the table read, nothing left, the frame covered
    @test licence(evidence()) === :strong

    # A TABLE THAT COVERS THE FRAME BUT REPORTS NOTHING. Every count taken from it is clean —
    # no carrier is left because none is reported, and its taxon set is the frame — so the
    # counts above stay at zero and frame_matches stays true. The scan's acceptance check is
    # what refuses such a table; its verdicts, tied to the table by digest, withhold the code.
    refused = evidence((acceptance_verdicts = ["MISS", "MISS", "PASS", "DISCLOSED"],))
    @test refused.recall_pass_derived === true
    @test only_reason(removal(refused), "scan_evidence_ok")
    @test licence(refused) === :fallback
    # ... and the same table read with no acceptance record at all
    @test licence(evidence((acceptance_verdicts = nothing,))) === :fallback
    # A TABLE LEFT FROM AN EARLIER RUN, beside the records of a later run that refused its
    # own table: the records describe a different digest than the one read
    stale = evidence((record_sha256s = (threshold = sha_other, recall = sha_other,
                                       acceptance = sha_other),))
    @test licence(stale) === :fallback

    # A RECALL RECORD THAT SAYS PASS WHILE MEASURING NOTHING
    hollow = evidence(merge(recall_of(0, 0), (recall_pass_field = true,)))
    @test hollow.recall_pass_derived === false
    @test reason_keys(removal(hollow)) == ["recall_pass", "scan_evidence_ok"]
    @test licence(hollow) === :fallback
    # a consumer that forwards the record's own field instead of the derived verdict is
    # stopped all the same, by the evidence leg
    @test forwarded(hollow, true) === :fallback
    # A RECALL RECORD WHOSE COUNTS CONTRADICT ITS PASS FIELD
    @test licence(evidence(merge(recall_of(14, 30), (recall_pass_field = true,)))) === :fallback

    # A RECALL RECORD WHOSE RATE DISAGREES WITH ITS OWN COUNTS. Everything else is honest and
    # the counts pass, so the evidence leg is the only one that withdraws the code.
    wrong_rate = evidence((recall_rate_field = 0.5,))
    @test wrong_rate.recall_pass_derived === true
    @test only_reason(removal(wrong_rate), "scan_evidence_ok")
    @test licence(wrong_rate) === :fallback
    # A BAR LOWERED IN THE RECALL RECORD beside counts the acceptance table's bar fails
    lowered = evidence(merge(recall_of(12, 25; bar = 0.2), (bar_in_other_records = (acceptance = 0.8,),)))
    @test lowered.recall_pass_derived === true
    @test licence(lowered) === :fallback
    # AN ACCEPTANCE TABLE WITH A GATED ROW CUT OUT, and one cut down to its disclosures
    @test licence(evidence((acceptance_checks = ["gate-p", "gate-q", "note-s"],
                            acceptance_verdicts = ["PASS", "PASS", "DISCLOSED"]))) === :fallback
    @test licence(evidence((acceptance_checks = ["note-s"], acceptance_verdicts = ["DISCLOSED"]))) ===
          :fallback
    # AN ACCEPTANCE TABLE WHOSE FLAGGED COUNT IS NOT THE PARSED TABLE'S, digests intact
    @test licence(evidence((table_n_flagged = 0,))) === :fallback

    # A CLEAN SAMPLE TOO SMALL FOR ITS BOUND. Every flagged taxon recovered, a point estimate of
    # one above the bar, every record honest: the bound is below the bar, so the recall has
    # not passed and the code is the fallback ...
    @test recall_lower_bound(25, 25) < 0.9
    small = evidence(recall_of(25, 25; bar = 0.9))
    @test small.recall_pass_derived === false
    @test reason_keys(removal(small)) == ["recall_pass", "scan_evidence_ok"]
    @test licence(small) === :fallback
    # ... also for a consumer that forwarded the record's pass field, which is true here
    @test small.recall_pass_derived === false && forwarded(small, true) === :fallback
    # ... while the same counts under a bar their bound reaches license the strong code
    @test recall_lower_bound(25, 25) >= 0.85
    @test licence(evidence(recall_of(25, 25; bar = 0.85))) === :strong

    # A FRAME MISMATCH. A carrier whose row is missing from the table reads as no carrier, so
    # the leftover counts stay at zero and the evidence is honest about the table it has;
    # only the frame comparison withdraws the code, and it is the only reason given.
    dropped_row = removal(evidence(); frame_matches = false)
    @test only_reason(dropped_row, "frame_matches")
    @test background_licence_code(dropped_row) === :fallback

    # A FAMILY CARRIER WITH NO HIT left in the analysed arm: one is enough
    @test only_reason(removal(evidence(); left_ipr_ge500 = 1), "left_ipr_ge500")
    @test licence(evidence(); left_ipr_ge500 = 1) === :fallback

    # AN UNUSABLE ANCHOR. Whatever the primary run found, and whatever it could not compute,
    # no partition is emitted and the caller is told there was no anchor.
    for unusable in (nothing, missing)
        p = category_sets(["cat-p", "cat-q"], unusable, ["cat-q"]; primary_undecidable = ["cat-r"])
        @test p.anchor_available == false
        @test all(isempty, (p.newly, p.preexisting, p.no_longer, p.legend_only,
                            p.undecidable_in_primary))
    end

    # AN UNCOMPUTED PRIMARY INTERVAL, from bound pairs to the partition: the category with no
    # interval is undecidable in the primary run, not lost, and maps to no caller code.
    rows = ["cat-a" => (0.2, 0.9), "cat-b" => (-0.3, 0.4), "cat-c" => (NaN, NaN)]
    sig = [c for (c, (lo, hi)) in rows if interval_excludes_zero(lo, hi) === :excludes_zero]
    und = [c for (c, (lo, hi)) in rows if interval_excludes_zero(lo, hi) === :undecidable]
    part = category_sets(sig, first.(rows), String[]; primary_undecidable = und)
    @test part.preexisting == ["cat-a"]
    @test part.no_longer == ["cat-b"]
    @test part.undecidable_in_primary == ["cat-c"]
    @test [interval_code(lo, hi; excludes = :keep, includes = :drop) for (_, (lo, hi)) in rows] ==
          [:keep, :drop, :undecidable]
end

@testset "background_licence_code — only exactly true is strong" begin
    @test background_licence_code(true)  === :strong
    @test background_licence_code(false) === :fallback
    for x in (nothing, missing, 1, "true", :strong, (complete = 1,), (complete = false,),
              (verdict = true,), (ok = true, recall_pass_derived = true, reasons = String[]))
        @test background_licence_code(x) === :fallback
    end
    @test background_licence_code((complete = true, reasons = String[])) === :strong
end

@testset "the module is pure and exports exactly its contract" begin
    @test Set(names(OutcomeRules)) == Set([:OutcomeRules, :interval_excludes_zero,
        :interval_code, :sharing_wording_code, :depletion_tier, :category_sets,
        :background_removal_complete, :background_licence_code, :scan_evidence_status,
        :recall_lower_bound])
    # No I/O and no dependencies: a predicate that reads a file is not a predicate, and the
    # two consumers must be able to call it on values they computed themselves. A digest is
    # compared here and computed by the caller, which has the file.
    src = read(joinpath(@__DIR__, "..", "src", "ModuleOutcomeRules.jl"), String)
    @test !occursin(r"^\s*(using|import)\s"m, src)
    @test !occursin(r"\b(open|read|write|readline|readlines|print|println|mkpath|rm|download|sha256)\(", src)
    @test !occursin(r"@(printf|info|warn|error|show)\b", src)
end

end
