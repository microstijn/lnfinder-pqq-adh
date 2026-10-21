module OutcomeRules

export interval_excludes_zero, interval_code, sharing_wording_code, depletion_tier,
       category_sets, background_removal_complete, background_licence_code,
       scan_evidence_status, recall_lower_bound

#=
    OutcomeRules — the recorded outcome rules, as pure predicates over an estimate.

    WHY THIS MODULE EXISTS. The rules that decide which claim a habitat-enrichment result
    licenses are pure predicates here rather than inline expressions in
    scripts/habitat_enrichment.jl, where a one-liner sitting next to the code that prints
    its answer is reachable by no test. Each entry below is a failure this module forecloses;
    the first four are rules that cannot fail, or cannot be stated at all, inline.

      * THE TWO-VALUED INTERVAL TEST. Answering "does this interval exclude zero?" with a
        Bool folds a NaN interval — an estimate that was never computed — into
        "includes zero", and the output then reports a decision about a number that does not
        exist. Three values are the remedy: excluding zero, including zero and undecidable
        are three different states, and only the first two are evidence.

      * THE NOVELTY TEST THAT COMPARES AGAINST A CONSTANT. Whether a category is NEW can
        only be asked against the earlier result it would be new relative to. Comparing a
        run against a hard-coded list of category names reports a category that had already
        excluded zero as new — including on the earlier run itself, where by construction
        nothing can be new. `category_sets` takes that earlier (anchor) set as an argument
        and has no constant to fall back on.

      * A ROBUSTNESS TIER CLAIMED WITHOUT THE EVIDENCE FOR IT. The two strong tiers are
        defined by what a one-at-a-time sensitivity grid does, so they cannot be read off
        the primary row. `depletion_tier` returns :undecidable when no grid is supplied, and
        the only tier it will decide without one is the tier that needs no grid: the one
        where the primary interval already includes zero.

      * A GRADED WORD TAKEN FROM AN UNBOUNDED RATIO. The sharing rule grades a comparator
        arm's effect as a fraction of the reference arm's. A fraction above one, a negative
        fraction, or one produced by a denominator near zero is not a fraction of anything,
        and a bare comparison returns the strongest grade for all three.

    Four more places where an absent measurement can pass for a present one — in the module
    itself and in what its consumer would otherwise decide on its own.

      * AN ABSENCE REPORTED AS A LOSS. A category that excluded zero in the anchor and is
        missing from the primary set is not thereby one that no longer excludes zero: its
        primary interval may have been computed and included zero, or may never have been
        computed at all, and only the first is a measurement. `category_sets` takes the
        categories whose primary interval is undecidable and returns them as a list of their
        own.

      * A MAPPING THAT LIVES IN THE CONSUMER. Rules that turn one interval into one of two
        caller codes are two-line ternaries when they sit in the script, which no test
        reaches and which the next consumer has to copy. `interval_code` is the one mapping,
        and its undecidable case is not the caller's to rename.

      * A REFUSAL THAT LOOKS LIKE A GRADE. Labelling the sharing rule's refusal of an
        out-of-domain fraction with the form that happened to see it leaves a table unable to
        tell "graded, and too small to count" from "not a fraction at all". The refusal
        carries its own label.

      * A LICENCE READ OFF RECORDS NOBODY TIED TO THE TABLE. Taking the scan's persisted
        recall verdict on trust lets a record that measured nothing but says "pass", a record
        left over from an earlier scan, and a table the scan's own acceptance check refused
        each license the strong code. `scan_evidence_status`
        requires every record to carry the digest of the table actually read, no acceptance
        verdict to be a miss, and the recall verdict to be re-derived from the record's own
        counts and to agree with what the record says.

    Two more, which remain even with every record tied to the table: its other numbers taken
    as written, and the recall bar read against the wrong quantity.

      * A BAR MET BY A POINT ESTIMATE. A recall sample that recovers every taxon it flags has
        seen no miss. It has not shown that misses are rarer than the bar allows, and a small
        sample cannot show that however clean it is. The licence therefore takes the recall
        verdict on the exact one-sided lower confidence bound, `recall_lower_bound`, which a
        sample clears only once it is large enough for "no miss seen" to mean "misses are
        rare". The scan's own acceptance gate stays on the point estimate, because it asks a
        different question — whether its exclusion list is fit for finding carriers — and a
        bound below the bar does not make that list wrong.

      * NUMBERS NOBODY RECOMPUTED. With every digest intact, a recall record whose rate
        disagrees with its own counts, a bar lowered in the very record whose counts it
        judges, more flagged taxa than the sample holds, and an acceptance table with its
        gated rows cut out would each license the strong code. Each of those numbers is
        compared with the counts it summarises, with the other records that carry it, or
        with the table that was read.

    EVERYTHING HERE IS PURE: plain numbers and labels in, NamedTuples out. No I/O, no
    files, no globals, no DataFrames. Two consumers need these predicates — the enrichment
    run and the sensitivity grid that varies one axis at a time — and a rule that exists in
    two copies is a rule that will eventually disagree with itself. The digest of a file is
    therefore computed by the caller, which has the file; this module only compares digests.

    NEUTRAL CODES, NEVER SENTENCES. Every function returns Symbols and label lists, so the
    step that writes prose maps a code to a phrase somewhere else. That keeps wording out
    of this repository, and it also makes a rule's output diffable across runs, which a
    sentence is not.

    UNKNOWN IS NOT FALSE. `missing`, `nothing` and NaN reach these functions whenever an
    estimate was not computed or a column was absent. None of them is read as a number and
    none of them licenses a claim: an unknown input either yields :undecidable or, where the
    contract is a Bool, yields the answer that withholds the claim together with the reason.

    The grade boundaries and the conservative band around each are the recorded rule's own
    parameters, and not quantities measured from data. SHARE_MAX is not one of them: the
    recorded rule sets no upper bound on the fraction, and the domain it guards is this
    module's own addition, made in the conservative direction.
=#

const SHARE_LARGE   = 0.75   # at or above this fraction, the strong grade
const SHARE_PARTIAL = 0.25   # at or above this fraction, the middle grade
const SHARE_BAND    = 0.05   # half-width of the conservative band around each boundary
const SHARE_MAX     = 1.0    # domain guard added here: a comparator effect above the
                             # reference effect is not a share of it

# The verdict vocabulary of the scan's acceptance table. Anything else is not a verdict.
const ACCEPTANCE_VERDICTS = ("PASS", "MISS", "DISCLOSED")

# The error rate of the recall bound: the one-sided 95 % bound. It is the level the bound
# is persisted at, so the bound re-derived here is the same quantity the record carries.
const RECALL_ALPHA = 0.05

# How closely a number a record writes has to match the same number derived here. A record
# writes a rate or a bar to six decimals and a bit-score cutoff to four, so a faithful record
# sits within half a unit of its last decimal. The bound is written at full precision, and
# both sides bisect to far below its tolerance. PARSE_SLACK is what reading a decimal back
# into binary can add to a difference: an honest rate exactly half a unit from its written
# form lands a hair past the half unit without it.
const RATE_TOL    = 5e-7
const BITS_TOL    = 5e-5
const BOUND_TOL   = 1e-9
const PARSE_SLACK = 1e-12

# An estimate that does not exist must never be read as a number: `missing`, `nothing`,
# NaN and anything non-numeric all mean "not computed". A Bool is rejected too, because a
# flag arriving where a bound was expected is a mix-up, not a bound of value 1.
_num(::Any) = nothing
_num(x::Real) = isnan(x) ? nothing : float(x)
_num(::Bool) = nothing

# Same, for a count of units: a negative count is not a count either.
_count(::Any) = nothing
_count(x::Real) = (isnan(x) || x < 0) ? nothing : float(x)
_count(::Bool) = nothing

# Same, for a count of taxa that a rate is taken over: a fraction of a taxon, or an infinite
# number of them, is a column mix-up rather than a count. So is a whole number past the range
# in which double-precision arithmetic holds every whole number exactly, which covers every
# value beyond Int range: the conversion would otherwise throw there, and an input this
# module cannot use has to come back as unknown with a reason, never as a crash. The rate and
# the bound are computed in that arithmetic, so a count past the range could not be used
# exactly even where it converts.
_whole(::Any) = nothing
_whole(x::Real) = (isfinite(x) && 0 <= x <= maxintfloat(Float64) && isinteger(x)) ? Int(x) : nothing
_whole(::Bool) = nothing

# A rate or bound as a record wrote it. NaN is kept, because NaN is what a record writes for
# a rate over an empty sample, and whether it belongs there is decided against the counts.
_written(::Any) = nothing
_written(x::Real) = Float64(x)
_written(::Bool) = nothing

# A text field as a record wrote it; blank is not a value.
_text(::Any) = nothing
_text(x::AbstractString) = (s = strip(x); isempty(s) ? nothing : String(s))

# Same, for a persisted verdict. Only a real Bool is a verdict; a 0/1 or a string is not.
_flag(x::Bool) = x
_flag(::Any) = nothing

_labels(::Nothing) = String[]
_labels(::Missing) = String[]
_labels(x::AbstractString) = isempty(strip(x)) ? String[] : String[strip(x)]
_labels(xs) = String[s for s in (strip(string(x)) for x in xs) if !isempty(s)]

# Sets leave this module sorted. A set compared in hash order is a test that passes on one
# Julia version and fails on the next, and a key written in hash order is not diffable.
_srt(xs) = sort!(unique(String.(xs)))

"""
    interval_excludes_zero(lo, hi) -> Symbol

`:excludes_zero`, `:includes_zero`, or `:undecidable` — the three states a bound pair can
be in, never two.

`:undecidable` covers every pair that is not an interval over the reals: either bound NaN,
`missing`, `nothing` or non-numeric (the estimator returns NaN bounds when an arm has no
hits at all, and a CSV column can arrive as `missing`), and also `lo > hi`, which means the
caller handed over its columns the wrong way round.

WHY THREE VALUES. A two-valued predicate maps a NaN interval to `false`, and a consumer then
prints "the interval includes zero" — a statement of fact about an estimate that was never
computed. An undecidable interval must be reportable as undecidable, or the
absence of an estimate becomes indistinguishable from evidence of a null result.

A bound of exactly zero counts as including zero: the recorded rules test the sign of a
bound strictly, so an interval that reaches the null value has not excluded it.
"""
function interval_excludes_zero(lo, hi)
    l = _num(lo); h = _num(hi)
    (l === nothing || h === nothing) && return :undecidable
    l > h && return :undecidable
    return (h < 0 || l > 0) ? :excludes_zero : :includes_zero
end

"""
    interval_code(lo, hi; excludes, includes) -> Symbol

The caller's code for a decision keyed on one interval: `excludes` when `(lo, hi)` excludes
zero, `includes` when it includes zero, and `:undecidable` whenever `interval_excludes_zero`
answers so — a NaN, missing, non-numeric or swapped bound pair.

A rule of the form "do this if the interval includes zero, that if it excludes zero" has two
decided outcomes and one that is not a decision. Written as a ternary next to each rule, the
mapping is one edit away from sending the uncomputed case to one of the decided codes, which
is how a Bool interval test comes to report estimates that do not exist. Here the third case
is fixed and only the two decided codes are the caller's, so every rule of this shape is
decided by one tested function and no second mapping has to live in a script.

THE CODES ARE CHECKED BEFORE THE INTERVAL IS. Both keywords are required and must be
Symbols. They must differ, because a rule whose two outcomes coincide does not depend on the
interval at all, which is a mistake at the call site rather than a rule. Neither may be
`:undecidable`, or an interval that was never computed would read the same as one that was.
Each of these is an `ArgumentError` whatever the bounds are, so a bad call fails on the
first run rather than on the first run whose interval happens to be decidable.
"""
function interval_code(lo, hi; excludes, includes)
    (excludes isa Symbol && includes isa Symbol) ||
        throw(ArgumentError("interval_code needs Symbol codes, got a $(typeof(excludes)) " *
                            "and a $(typeof(includes))"))
    excludes === includes &&
        throw(ArgumentError("interval_code: `excludes` and `includes` are both :$excludes, " *
                            "so the interval would decide nothing"))
    (excludes === :undecidable || includes === :undecidable) &&
        throw(ArgumentError("interval_code: :undecidable is reserved for an interval that " *
                            "was not computed and cannot be a decided code"))
    s = interval_excludes_zero(lo, hi)
    return s === :excludes_zero ? excludes : s === :includes_zero ? includes : :undecidable
end

_share_grade(x) = x >= SHARE_LARGE ? :large : x >= SHARE_PARTIAL ? :partial : :none

# The conservative band: a fraction within ±SHARE_BAND of a grade boundary has not resolved
# that boundary, so the weaker of the two grades it sits between applies. The two bands
# cannot overlap, and `nothing` means the fraction is clear of both.
function _band_grade(x)
    SHARE_LARGE   - SHARE_BAND <= x <= SHARE_LARGE   + SHARE_BAND && return :partial
    SHARE_PARTIAL - SHARE_BAND <= x <= SHARE_PARTIAL + SHARE_BAND && return :none
    return nothing
end

_in_share_domain(x) = 0 < x <= SHARE_MAX

"""
    sharing_wording_code(ratio, lo, hi; comparator_excludes_zero) -> (code, rule_used)

Grade how much of a reference arm's effect a comparator arm reproduces, as
`code ∈ (:none, :partial, :large)`, together with `rule_used` — the form of the rule that
decided, so the grade is auditable rather than merely asserted.

`ratio` is the comparator estimate over the reference estimate; `lo`/`hi` are the bounds of
an interval on that ratio, either of which may be NaN/`missing`/`nothing` when no interval
was computed. `comparator_excludes_zero` takes a Bool or the Symbol from
`interval_excludes_zero`, and describes the COMPARATOR's own estimate, not the ratio's.

THE FORMS, IN THE ORDER THEY GOVERN.

  `:comparator_null`  The comparator's own interval includes zero: there is no comparator
                      effect to be a fraction of, so no grade is licensed whatever the
                      ratio says. This overrides every form below.
  `:undecidable`      The comparator's status is itself unknown, or neither a ratio nor an
                      interval exists. Withholding the grade here is the whole point: an
                      expression that folded this state into `:comparator_null` would print
                      "the comparator interval includes zero" for an estimate that has never
                      been computed.
  `:out_of_domain`    A known ratio, or a known bound of the interval, lies outside
                      (0, SHARE_MAX]. See below.
  `:interval`         An interval on the ratio exists, and it supersedes the band — the
                      recorded rules say so explicitly. The grade is taken at the lower
                      bound, so an interval straddling a boundary takes the weaker grade.
  `:band`             No interval; the point fraction sits inside a conservative band.
  `:point`            No interval; the point fraction is clear of both bands.

OUTSIDE THE DOMAIN THE ANSWER IS :none, NOT THE STRONGEST GRADE, AND IT SAYS SO. A share is a
fraction of the reference effect, so this module grades only on (0, SHARE_MAX]. A ratio above
one (the comparator is the more enriched arm), a ratio at or below zero (the arms move in
opposite directions, or the comparator does not move), the huge ratio a near-zero
denominator produces, and an interval reaching past either end are all outside it. A bare
threshold comparison returns the strongest grade for exactly these inputs, including for a
reference arm with no effect at all. `code` is then :none and `rule_used` is
`:out_of_domain`, which is deliberately not the label of any grading form: a :none reached by
grading a small fraction and a :none that refused to grade are different findings, and a
table that gave both the same label could not be re-decided by its reader.

Every known value is checked, including a bound whose partner is missing and a point ratio
sitting beside an interval that is itself inside the domain. One value outside the domain
means the estimate is not a share, whichever form would otherwise have decided.
"""
function sharing_wording_code(ratio, lo, hi; comparator_excludes_zero)
    cmp = _zero_status(comparator_excludes_zero)
    cmp === :includes_zero && return (code = :none, rule_used = :comparator_null)
    cmp === :excludes_zero || return (code = :none, rule_used = :undecidable)

    r = _num(ratio); l = _num(lo); h = _num(hi)
    any(x -> x !== nothing && !_in_share_domain(x), (r, l, h)) &&
        return (code = :none, rule_used = :out_of_domain)

    if l !== nothing && h !== nothing && l <= h
        return (code = _share_grade(l), rule_used = :interval)
    end

    r === nothing && return (code = :none, rule_used = :undecidable)
    b = _band_grade(r)
    return b === nothing ? (code = _share_grade(r), rule_used = :point) :
                           (code = b, rule_used = :band)
end

_zero_status(x::Bool) = x ? :excludes_zero : :includes_zero
_zero_status(x::Symbol) = x in (:excludes_zero, :includes_zero) ? x : :undecidable
_zero_status(::Any) = :undecidable

"""
    depletion_tier(primary_lo, primary_hi, variant_intervals;
                 expected_variants = nothing) -> (tier, breakers)

Which robustness tier a depletion result reaches: `:t1` (every variant reproduces it),
`:t2` (the primary excludes zero but at least one variant does not — `breakers` names
them), `:t3` (the primary does not support the claim), or `:undecidable`.

`variant_intervals` is the one-at-a-time sensitivity grid as `name => (lo, hi)` pairs — a
Dict, a Vector of Pairs, a NamedTuple keyed by variant name, or `nothing` when the grid has
not been computed. A variant's value may also be a NamedTuple carrying `lo` and `hi`.
Anything else is an error rather than a guess.

WHAT MAY BE DECIDED WITHOUT THE GRID, AND WHAT MAY NOT. `:t3` needs no grid: if the primary
interval does not exclude zero on the depleted side there is no claim left for a variant to
break, and the recorded rules drop it there. `:t1` and `:t2` are *defined* by the grid, so
with no variants supplied this returns `:undecidable` rather than the tier the primary row
happens to look like. Returning a tier the data cannot support is the defect this function
exists to prevent, and a grid that does not exist yet supports neither strong tier.

An undecidable primary interval is `:undecidable` too, with no variants examined. The
caller can tell the two undecidable cases apart — no primary estimate versus no grid — by
asking `interval_excludes_zero` about the primary bounds, which is the same question this
function asked first.

A VARIANT WITH NO INTERVAL IS A BREAKER, not a variant that agrees. The recorded rules
require every variant to reproduce the depletion; one whose interval could not be computed
has not reproduced it, so it is named in `breakers` and the tier is `:t2`. `breakers` is
sorted and free of duplicates.

A VARIANT THAT IS ABSENT IS NOT ONE THAT AGREES EITHER. A grid step that skips a variant
which failed to run hands over a shorter grid, and every variant left in it may reproduce
the depletion. The optional `expected_variants` names the grid the caller set out to run;
while any of those names is missing from `variant_intervals`, `:t1` cannot be claimed and
the answer is `:undecidable`, unless a variant that IS present already breaks the result,
which decides `:t2` on its own. Without `expected_variants` the grid supplied is taken to
be the whole grid, so a caller that can lose variants should always pass it.
"""
function depletion_tier(primary_lo, primary_hi, variant_intervals; expected_variants = nothing)
    status = interval_excludes_zero(primary_lo, primary_hi)
    status === :undecidable && return (tier = :undecidable, breakers = String[])
    # includes zero, or excludes it on the ENRICHED side: either way no depletion to break
    (status === :excludes_zero && _num(primary_hi) < 0) ||
        return (tier = :t3, breakers = String[])

    grid = _variant_pairs(variant_intervals)
    isempty(grid) && return (tier = :undecidable, breakers = String[])
    brk = String[]; seen = String[]
    for (name, bounds) in grid
        vlo, vhi = _bounds(bounds)
        label = String(strip(string(name)))
        push!(seen, label)
        reproduces = interval_excludes_zero(vlo, vhi) === :excludes_zero &&
                     (v = _num(vhi)) !== nothing && v < 0
        reproduces || push!(brk, label)
    end
    isempty(brk) || return (tier = :t2, breakers = _srt(brk))
    isempty(setdiff(_labels(expected_variants), seen)) ||
        return (tier = :undecidable, breakers = String[])
    return (tier = :t1, breakers = String[])
end

_variant_pairs(::Nothing) = Pair[]
_variant_pairs(::Missing) = Pair[]
_variant_pairs(x::AbstractDict) = collect(x)
_variant_pairs(x::NamedTuple) = collect(pairs(x))
_variant_pairs(xs) = Pair[p isa Pair ? p :
    error("a sensitivity variant must be given as name => (lo, hi), got a $(typeof(p))")
    for p in xs]

_bounds(v::Tuple{Any,Any}) = (v[1], v[2])
_bounds(v::AbstractVector) = length(v) == 2 ? (v[1], v[2]) :
    error("a variant interval needs exactly two bounds, got $(length(v))")
_bounds(v::NamedTuple) = (hasproperty(v, :lo) && hasproperty(v, :hi)) ? (v.lo, v.hi) :
    error("a variant NamedTuple needs `lo` and `hi` fields, got $(keys(v))")
_bounds(v) = error("a variant interval must be (lo, hi), got a $(typeof(v))")

"""
    category_sets(primary_sig, anchor_sig, governed; primary_undecidable)
        -> (newly, preexisting, no_longer, legend_only, anchor_available,
            undecidable_in_primary)

Partition the categories whose primary interval excludes zero against the categories that
already did so in the earlier ANCHOR run, which is what makes "new" a measurable property
instead of an impression.

`primary_sig` and `anchor_sig` are the label sets of the two runs; `governed` names the
categories the rules already handle individually, which therefore need no standing of their
own here; `primary_undecidable` names the categories whose primary interval could not be
computed. All four may be any iterable of labels, and every returned set is a sorted,
de-duplicated `Vector{String}`.

  `newly`                   primary minus anchor: excludes zero now and did not before.
  `legend_only`             `newly` minus `governed`: the ones no other rule speaks for.
  `preexisting`             primary and anchor, minus `governed`: carried over, not new.
  `no_longer`               anchor minus primary, minus `primary_undecidable`: excluded
                            zero before, and now has an interval that includes it.
  `undecidable_in_primary`  anchor minus primary, within `primary_undecidable`: excluded
                            zero before, and has no primary interval to say anything now.

THE DEFECT THIS FORECLOSES. An expression that compares the current run against a hard-coded
tuple of names and calls everything else new nominates a category even on the anchor itself,
where nothing can be new; it can quietly demote a category already standing in the earlier
report, or present a carried-over one as a discovery. No constant can answer this question;
only the anchor set can.

A LOSS IS A MEASUREMENT. A category missing from the primary set is missing either because
its primary interval includes zero or because it has no interval at all, and only the first
says the earlier result went away. Reporting both as `no_longer` is the uncomputed-estimate-
read-as-a-null defect again, one level up, so the second are returned apart. A category the
primary run has no row for has no computed interval either, and belongs in
`primary_undecidable`. A category that is in `primary_sig` stays there even when it is also
undecidable, since it excluded zero on whichever unit was measured.

THE KEYWORD HAS NO DEFAULT, so the distinction cannot be forgotten by leaving it out. A caller
that cannot say which primary intervals were computed passes `nothing` (or `missing`), and
then no absence from the primary set is evidence: `no_longer` is empty and every anchor
category the primary set lacks is returned in `undecidable_in_primary`.

WHEN THE ANCHOR IS UNAVAILABLE (`nothing` or `missing`) every set is empty and
`anchor_available` is false. Empty is the only safe answer: treating an absent anchor as an
empty one would declare the entire primary set new. An EMPTY anchor set is not the same
input — it says the earlier run had nothing excluding zero, and then everything in the
primary set really is new.
"""
function category_sets(primary_sig, anchor_sig, governed; primary_undecidable)
    prim = _labels(primary_sig)
    gov  = _labels(governed)
    if anchor_sig === nothing || anchor_sig === missing
        return (newly = String[], preexisting = String[], no_longer = String[],
                legend_only = String[], anchor_available = false,
                undecidable_in_primary = String[])
    end
    anch   = _labels(anchor_sig)
    newly  = setdiff(prim, anch)
    absent = setdiff(anch, prim)
    unmeasured = (primary_undecidable === nothing || primary_undecidable === missing) ?
                 absent : intersect(absent, _labels(primary_undecidable))
    return (newly       = _srt(newly),
            preexisting = _srt(setdiff(intersect(prim, anch), gov)),
            no_longer   = _srt(setdiff(absent, unmeasured)),
            legend_only = _srt(setdiff(newly, gov)),
            anchor_available = true,
            undecidable_in_primary = _srt(unmeasured))
end

# A SHA-256 digest as a record carries it: 64 hexadecimal characters. Case and surrounding
# whitespace (a CRLF line end included) are not part of a digest; a truncated or padded one
# identifies nothing and is refused rather than compared by prefix.
_sha(::Any) = nothing
function _sha(x::AbstractString)
    s = lowercase(strip(x))
    return occursin(r"^[0-9a-f]{64}$", s) ? String(s) : nothing
end
_short(s::AbstractString) = first(s, 12) * "…"

# The persisted records' digests, as (record name, value) pairs. A value may itself be a list,
# for a record whose rows each carry the digest. `nothing` means none were supplied.
_record_pairs(::Nothing) = nothing
_record_pairs(::Missing) = nothing
_record_pairs(x::AbstractDict) = sort!([string(k) => v for (k, v) in x], by = first)
_record_pairs(x::NamedTuple) = [string(k) => v for (k, v) in pairs(x)]
_record_pairs(xs::Union{AbstractVector,Tuple}) =
    [p isa Pair ? (string(first(p)) => last(p)) : ("record $i" => p) for (i, p) in enumerate(xs)]
_record_pairs(x) = ["record 1" => x]

_digest_list(::Nothing) = Any[]
_digest_list(::Missing) = Any[]
_digest_list(xs::Union{AbstractVector,Tuple}) = collect(Any, xs)
_digest_list(x) = Any[x]

# Acceptance verdicts as strings. A blank cell is kept, as a verdict outside the vocabulary:
# an empty verdict is not a pass, and dropping it would turn one into nothing to object to.
_verdicts(::Nothing) = nothing
_verdicts(::Missing) = nothing
_verdicts(x::Union{AbstractString,Symbol}) = [String(strip(string(x)))]
_verdicts(xs::Union{AbstractVector,Tuple}) =
    [x === missing || x === nothing ? "" : String(strip(string(x))) for x in xs]
_verdicts(x) = [String(strip(string(x)))]

"""
    recall_lower_bound(k, n; alpha = 0.05) -> Float64

The exact one-sided Clopper-Pearson lower confidence bound on a recall observed as `k`
recovered out of `n`: the smallest true recall at which recovering `k` or more of `n` still
has probability `alpha`. At the default it is the one-sided 95 % lower bound.

  `n == 0`   NaN. An empty sample bounds nothing, and a number invented for it would read
             as evidence.
  `k == 0`   0.0. Recovering nothing is consistent with a recall of zero.
  `k == n`   `alpha^(1/n)`, in closed form: recovering all `n` has probability `p^n`.
  otherwise  the root in `p` of P(X >= k | n, p) = alpha, found by bisection on the binomial
             upper tail until the bracket is narrower than 1e-12. The tail rises
             monotonically in `p`, so the bracket always holds the root.

WHY A BOUND AND NOT THE RATE. A bar on the recall is a claim about the prefilter, not about
the sample: that it misses fewer carriers than the bar allows. A sample that recovered every
taxon it flagged has a point estimate of one whatever its size, so the point estimate clears
every bar on a sample of three as readily as on a sample of three thousand. The lower bound is
what the sample can vouch for, and it reaches a bar only once the sample is large enough for
the absence of a miss to say something.

A negative count, `k > n`, a Bool in place of a count, a count beyond what double-precision
arithmetic holds exactly, and an `alpha` outside (0, 1) are each an `ArgumentError`. Every one
is a mistake at the call site, and a number returned for it would bound a sample that does
not exist.

THE TAIL IS SUMMED WHERE ITS MASS IS. Each term starts in log space, so no binomial
coefficient overflows, and the sum runs from `k` away from the binomial's mode, where the
terms only fall, stopping once what is left cannot move the sum at double precision. The
cost grows with the spread of the binomial rather than with `n`, so an absurdly large count
gives a slow answer rather than an allocation the size of the count.
"""
function recall_lower_bound(k::Integer, n::Integer; alpha::Real = RECALL_ALPHA)
    (k isa Bool || n isa Bool) &&
        throw(ArgumentError("recall_lower_bound needs two counts, got a Bool"))
    (k < 0 || n < 0) &&
        throw(ArgumentError("recall_lower_bound: a count cannot be negative, got $k of $n"))
    k > n &&
        throw(ArgumentError("recall_lower_bound: $k recovered out of $n, which no sample can give"))
    n > maxintfloat(Float64) &&
        throw(ArgumentError("recall_lower_bound: a sample of $n is beyond what " *
                            "double-precision arithmetic holds exactly"))
    (alpha isa Bool || !(0 < alpha < 1)) &&
        throw(ArgumentError("recall_lower_bound: alpha must lie in (0, 1), got $alpha"))
    n == 0 && return NaN
    k == 0 && return 0.0
    a = Float64(alpha)
    k == n && return a^(1 / Float64(n))
    return _bisect_lower_bound(Int(k), Int(n), a)
end

# The bisection itself, apart from the closed forms, so a test can check the two against each
# other where both apply.
function _bisect_lower_bound(k::Int, n::Int, alpha::Float64)
    lo, hi = 0.0, 1.0
    while hi - lo > 1e-12
        mid = (lo + hi) / 2
        _binomial_upper_tail(k, n, mid) < alpha ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

# log(m!) for a whole m >= 0: summed directly while m is small, and from Stirling's series
# beyond, where the first term it leaves out is already below 1e-15 and falls from there.
const _LOGFACT_DIRECT = 20
function _logfact(m::Int)
    m <= _LOGFACT_DIRECT && return sum(log, 2:m; init = 0.0)
    x = Float64(m)
    return x * log(x) - x + 0.5 * log(2π * x) +
           1 / (12x) - 1 / (360x^3) + 1 / (1260x^5) - 1 / (1680x^7)
end

# P(X >= k) for X ~ Binomial(n, p), with 1 <= k <= n and p strictly inside (0, 1).
#
# Above the mode the terms fall from i = k upward, so the tail is summed that way. At or below
# it, the terms of the LOWER tail fall from i = k - 1 downward, and the upper tail is its
# complement; that complement is at least one half there, so taking it loses nothing. In both
# directions the ratio between neighbouring terms keeps shrinking, which bounds everything not
# yet added by t·r/(1 - r) and lets the sum stop as soon as that is negligible. A ratio that
# rounding puts at or above one is never taken for a small one.
function _binomial_upper_tail(k::Int, n::Int, p::Float64)
    lp, lq = log(p), log1p(-p)
    logterm(i) = _logfact(n) - _logfact(i) - _logfact(n - i) + i * lp + (n - i) * lq
    odds = p / (1 - p)
    if k + 1 > (n + 1.0) * p
        t = exp(logterm(k)); s = t; i = k
        while i < n
            r = (n - i) / (i + 1) * odds
            t *= r; s += t; i += 1
            r < 1 && t * r / (1 - r) <= 1e-17 * s && break
        end
        return s
    else
        t = exp(logterm(k - 1)); s = t; i = k - 1
        while i > 0
            r = i / (n - i + 1) / odds
            t *= r; s += t; i -= 1
            r < 1 && t * r / (1 - r) <= 1e-17 * s && break
        end
        return 1 - s
    end
end

_pct(alpha) = round(Int, 100 * (1 - alpha))

# A rate or bound a record wrote, against the one its counts give: NaN exactly where the
# derived value is NaN, and within `tol` everywhere else. With unusable counts there is
# nothing to compare against, and the counts already carry their own reason.
function _written_agrees!(reasons, key, what, raw, derived, tol)
    w = _written(raw)
    if w === nothing
        push!(reasons, "$key: unknown, or not a number")
    elseif derived !== nothing
        agrees = isnan(derived) ? isnan(w) : (!isnan(w) && abs(w - derived) <= tol)
        agrees || push!(reasons, isnan(derived) ?
            "$key: the record gives a $what of $w, but it flagged nothing, so its counts give none" :
            "$key: the record gives a $what of $w, but its own counts give $derived")
    end
    return reasons
end

# One field the recall record carries, against the same field in the other records. Every
# other record named has to carry a usable value, and every value has to agree with the
# recall record's own. Supplying no other record is a reason too: a field compared with
# nothing has not been checked.
function _others_agree!(reasons, key, what, own, others, parse, same)
    recs = _record_pairs(others)
    if recs === nothing || isempty(recs)
        push!(reasons, "$key: no other record's $what was supplied, so the recall record's " *
                       "is checked against nothing")
        return reasons
    end
    for (name, value) in recs
        vals = map(parse, _digest_list(value))
        if isempty(vals) || any(v -> v === nothing, vals)
            push!(reasons, "$key: $name carries no usable $what")
        elseif own !== nothing
            other = findfirst(v -> !same(v, own), vals)
            other === nothing || push!(reasons,
                "$key: $name gives a $what of $(vals[other]), but the recall record gives $own")
        end
    end
    return reasons
end

"""
    scan_evidence_status(; scan_sha256, record_sha256s,
                           acceptance_checks, acceptance_verdicts, required_checks,
                           flagged_row_value, table_n_flagged,
                           family_row_denominator, table_n_family,
                           recall_n_sampled, recall_n_full_flagged, recall_n_recovered,
                           recall_rate_field, recall_lower_bound_field, recall_bar,
                           recall_pass_field, recall_threshold_bits, recall_uniprot_release,
                           bar_in_other_records, threshold_bits_in_other_records,
                           uniprot_release_in_other_records)
        -> (ok, recall_pass_derived, reasons, recall_rate, recall_lower_bound)

Whether the records the background carrier scan persisted beside its table vouch for the
table that was actually read, and what the recall verdict is once it is re-derived from the
record's own counts, on the lower confidence bound, instead of read from its pass field.

THE TABLE AND ITS RECORDS.

  `scan_sha256`              the SHA-256 hex digest of the scan table, computed by the caller
                             from the very bytes it parsed — never taken from a record.
  `record_sha256s`           the digest each persisted record says it describes: a Dict or
                             NamedTuple keyed by record name, a Vector of `name => digest`
                             pairs, or a plain list of digests. A record whose rows each
                             carry a digest may give all of them as a Vector.
  `acceptance_checks`        the check column of the scan's acceptance table, row by row.
  `acceptance_verdicts`      its verdict column, in the same order.
  `required_checks`          the names of the gated checks the scan writes. The caller names
                             them because a check's wording belongs to the scan, and a module
                             that spelled it out would silently stop matching when it changed.
  `flagged_row_value`        the value of the acceptance row that counts the flagged taxa.
  `table_n_flagged`          the taxa with at least one hit, counted by the caller in the
                             table it parsed.
  `family_row_denominator`   the denominator of the acceptance row on family-signature
                             members at the length cut-off.
  `table_n_family`           the taxa with such a member, counted in the parsed table.

THE RECALL RECORD, AND WHAT THE OTHER RECORDS SAY OF THE SAME FIELDS.

  `recall_n_sampled`         taxa in the recall sample.
  `recall_n_full_flagged`    taxa of the sample that the full-proteome search flagged.
  `recall_n_recovered`       of those, the taxa the prefilter also recovered.
  `recall_rate_field`        the rate the record writes, NaN where nothing was flagged.
  `recall_lower_bound_field` the one-sided 95 % lower bound the record writes, NaN likewise.
  `recall_bar`               the recovery rate the check had to reach, as a fraction.
  `recall_pass_field`        the record's own pass verdict, as a Bool.
  `recall_threshold_bits`    the bit-score cutoff the record says it scored at.
  `recall_uniprot_release`   the release the record says its search space came from.
  `bar_in_other_records`, `threshold_bits_in_other_records`, `uniprot_release_in_other_records`
                             the same three fields as the other records carry them, in any
                             shape `record_sha256s` takes, with at least one record each.

`ok` is true only when all five of these hold.

  1. EVERY RECORD DESCRIBES THIS TABLE. Each digest equals `scan_sha256`. Records written by a
     run whose acceptance check failed, sitting beside a table left from an earlier run, carry
     a different digest; so does a record copied in from another scan. Records that agree
     with one another but not with the table read are foreign all the same.
  2. NOTHING WAS REFUSED, AND NO GATE WAS LEFT OUT. The verdicts are non-empty, every one is
     PASS, MISS or DISCLOSED, and none is MISS. A table that covers the frame but found no
     carriers anywhere passes every count a consumer can compute from it, and the scan's own
     acceptance check is the one place that can refuse it. But a table with its gated rows
     removed has no MISS left to find, so refusal alone is not enough: the check names must
     line up with the verdicts, at least one gated check must be required, and every
     required check must appear exactly once and not as a disclosure, which would mean its
     bar was never applied.
  3. THE ACCEPTANCE TABLE COUNTS THE TABLE READ. Its flagged-taxa value and its family
     denominator equal the same counts taken from the parsed table. A digest ties a record to
     bytes; this ties the record's numbers to what those bytes contain.
  4. THE RECALL RECORD AGREES WITH ITSELF AND WITH THE OTHER RECORDS. Its rate equals
     `n_recovered / n_full_flagged` to the six decimals it is written in, and its lower bound
     equals `recall_lower_bound(n_recovered, n_full_flagged)` to 1e-9, each NaN exactly when
     nothing was flagged. It flags no more taxa than it sampled. Its bar, cutoff and release
     equal the other records'. A bar read only from the record whose counts it judges can be
     lowered together with them, which is why the bar is held against another record's.
  5. THE RECALL VERDICT HOLDS, ON THE LOWER BOUND. `recall_pass_derived` is
     `n_full_flagged > 0 && recall_lower_bound(n_recovered, n_full_flagged) >= bar`, and it must
     be true. A record whose full search flagged nothing has measured nothing, so it cannot
     pass. A sample that recovered every flagged taxon but is too small for its bound to reach
     the bar has not shown that the bar is met, whatever its point estimate.

THE PASS FIELD IS JUDGED ON THE QUANTITY IT IS TAKEN ON. The record's pass field is the scan's
own gate, and that gate is the point estimate, `n_recovered / n_full_flagged >= bar`. The field
must equal that verdict re-derived from the counts, so a field that contradicts its own counts
is still refused, while a record that is honest about a point estimate its bound cannot back
is refused for the bound and not reported as inconsistent.

`recall_pass_derived` is `nothing` when the verdict cannot be derived: a count that is
unknown, not a whole non-negative number, or past the range in which double-precision
arithmetic holds whole numbers exactly, a recovered count larger than the flagged count (no rate of a sample can exceed one), or a bar
that is not a fraction in (0, 1] — a bar of zero passes every record, and a bar above one
passes none. `recall_rate` and `recall_lower_bound` are the rate and bound the counts give,
NaN over an empty sample and `nothing` where `recall_pass_derived` is, so a caller can report
the numbers the verdict rests on rather than the record's copies of them.

UNKNOWN IS NOT CLEAN HERE EITHER. Every keyword defaults to `nothing`, so a caller that has
not read one of the records cannot obtain `ok` by leaving it out. A missing or malformed
digest is a reason of its own, and so is a verdict outside the vocabulary, a blank cell
included. `reasons` names every failure, not only the first, each prefixed by the keyword it
concerns.

The caller passes `ok` on to `background_removal_complete` as `scan_evidence_ok` and
`recall_pass_derived` as `recall_pass`. A bound below the bar is inside `ok` too, so a caller
that forwarded the record's pass field in place of the derived verdict is refused all the same.
"""
function scan_evidence_status(; scan_sha256 = nothing, record_sha256s = nothing,
                              acceptance_checks = nothing, acceptance_verdicts = nothing,
                              required_checks = nothing,
                              flagged_row_value = nothing, table_n_flagged = nothing,
                              family_row_denominator = nothing, table_n_family = nothing,
                              recall_n_sampled = nothing, recall_n_full_flagged = nothing,
                              recall_n_recovered = nothing, recall_rate_field = nothing,
                              recall_lower_bound_field = nothing, recall_bar = nothing,
                              recall_pass_field = nothing, recall_threshold_bits = nothing,
                              recall_uniprot_release = nothing,
                              bar_in_other_records = nothing,
                              threshold_bits_in_other_records = nothing,
                              uniprot_release_in_other_records = nothing)
    # Every keyword arrives as whatever a record reader produced — a count, a float, a
    # string, `nothing` — and a body compiled afresh for each combination of those types costs
    # seconds per combination while buying nothing, since each value is classified by hand
    # below. One compiled body serves them all.
    @nospecialize
    reasons = String[]

    # 1. every record describes the table that was read
    scan = _sha(scan_sha256)
    scan === nothing &&
        push!(reasons, "scan_sha256: unknown, or not a SHA-256 hex digest")
    recs = _record_pairs(record_sha256s)
    if recs === nothing || isempty(recs)
        push!(reasons, "record_sha256s: no persisted record digest was supplied")
    else
        for (name, value) in recs
            digests = map(_sha, _digest_list(value))
            if isempty(digests) || any(d -> d === nothing, digests)
                push!(reasons, "record_sha256s: $name carries no usable SHA-256 digest")
            elseif scan !== nothing
                other = findfirst(d -> d != scan, digests)
                other === nothing || push!(reasons,
                    "record_sha256s: $name describes scan table $(_short(digests[other])), " *
                    "not the table read ($(_short(scan)))")
            end
        end
    end

    # 2. nothing was refused, and no gated check was left out
    verdicts = _verdicts(acceptance_verdicts)
    if verdicts === nothing || isempty(verdicts)
        push!(reasons, "acceptance_verdicts: no acceptance verdict was supplied")
        verdicts = nothing
    else
        n_unknown = count(v -> !(v in ACCEPTANCE_VERDICTS), verdicts)
        n_miss    = count(==("MISS"), verdicts)
        n_unknown > 0 && push!(reasons,
            "acceptance_verdicts: $n_unknown of $(length(verdicts)) verdicts are not " *
            "one of $(join(ACCEPTANCE_VERDICTS, ", "))")
        n_miss > 0 && push!(reasons,
            "acceptance_verdicts: $n_miss of $(length(verdicts)) acceptance checks MISS")
    end
    checks = _verdicts(acceptance_checks)
    if checks === nothing || isempty(checks)
        push!(reasons, "acceptance_checks: no acceptance check name was supplied")
        checks = nothing
    elseif verdicts !== nothing && length(checks) != length(verdicts)
        push!(reasons, "acceptance_checks: $(length(checks)) check names beside " *
                       "$(length(verdicts)) verdicts, so no verdict can be tied to its check")
    end
    required = unique(_labels(required_checks))
    if isempty(required)
        push!(reasons, "required_checks: no gated check was named, so a table with its " *
                       "gated rows removed would pass")
    elseif checks !== nothing
        aligned = verdicts !== nothing && length(checks) == length(verdicts)
        for name in required
            at = findall(==(name), checks)
            if isempty(at)
                push!(reasons, "acceptance_checks: the gated check \"$name\" is absent")
            elseif length(at) > 1
                push!(reasons, "acceptance_checks: the gated check \"$name\" appears " *
                               "$(length(at)) times, so no single verdict stands for it")
            elseif aligned && verdicts[only(at)] == "DISCLOSED"
                push!(reasons, "acceptance_checks: the gated check \"$name\" is DISCLOSED, " *
                               "so its bar was never applied")
            end
        end
    end

    # 3. the acceptance table counts the table that was read
    for (vkey, vraw, tkey, traw, what) in
            (("flagged_row_value", flagged_row_value, "table_n_flagged", table_n_flagged,
              "flagged taxa"),
             ("family_row_denominator", family_row_denominator, "table_n_family",
              table_n_family, "taxa with a family-signature member at the length cut-off"))
        v = _whole(vraw); t = _whole(traw)
        v === nothing && push!(reasons, "$vkey: unknown, or not a whole count")
        t === nothing && push!(reasons, "$tkey: unknown, or not a whole count")
        (v === nothing || t === nothing || v == t) || push!(reasons,
            "$vkey: the acceptance table gives $v $what, but the table read has $t")
    end

    # 4. the recall record agrees with itself and with the other records
    ns  = _whole(recall_n_sampled)
    nf  = _whole(recall_n_full_flagged)
    nr  = _whole(recall_n_recovered)
    bar = _num(recall_bar)
    ns === nothing &&
        push!(reasons, "recall_n_sampled: unknown, or not a whole count")
    nf === nothing &&
        push!(reasons, "recall_n_full_flagged: unknown, or not a whole count")
    nr === nothing &&
        push!(reasons, "recall_n_recovered: unknown, or not a whole count")
    bar_ok = bar !== nothing && 0 < bar <= 1
    bar_ok || push!(reasons, "recall_bar: unknown, or not a fraction in (0, 1]")
    ns !== nothing && nf !== nothing && nf > ns && push!(reasons,
        "recall_n_sampled: $nf taxa flagged in a sample of $ns, which no sample can give")
    nf !== nothing && nr !== nothing && nr > nf && push!(reasons,
        "recall_n_recovered: $nr taxa recovered out of $nf flagged, which no sample can give")
    counts_ok = nf !== nothing && nr !== nothing && nr <= nf
    rate = counts_ok ? (nf == 0 ? NaN : nr / nf) : nothing
    lb   = counts_ok ? recall_lower_bound(nr, nf) : nothing
    _written_agrees!(reasons, "recall_rate_field", "rate", recall_rate_field, rate,
                     RATE_TOL + PARSE_SLACK)
    _written_agrees!(reasons, "recall_lower_bound_field", "lower bound",
                     recall_lower_bound_field, lb, BOUND_TOL)

    bits    = _num(recall_threshold_bits)
    release = _text(recall_uniprot_release)
    bits === nothing &&
        push!(reasons, "recall_threshold_bits: unknown, or not a number")
    release === nothing &&
        push!(reasons, "recall_uniprot_release: unknown, or blank")
    within(tol) = (a, b) -> abs(a - b) <= tol
    _others_agree!(reasons, "bar_in_other_records", "bar", bar, bar_in_other_records,
                   _num, within(RATE_TOL + PARSE_SLACK))
    _others_agree!(reasons, "threshold_bits_in_other_records", "cutoff", bits,
                   threshold_bits_in_other_records, _num, within(BITS_TOL + PARSE_SLACK))
    _others_agree!(reasons, "uniprot_release_in_other_records", "release", release,
                   uniprot_release_in_other_records, _text, ==)

    # 5. the recall verdict on the lower bound; the record's pass field on its point estimate
    derived = (counts_ok && bar_ok) ? (nf > 0 && lb >= bar) : nothing
    point   = (counts_ok && bar_ok) ? (nf > 0 && rate >= bar) : nothing
    if derived === false
        push!(reasons, nf == 0 ?
            "recall_n_full_flagged: the full search flagged no taxon in the sample, so the " *
            "check measured nothing" :
            point ?
            "recall_n_recovered: $nr of $nf flagged taxa recovered, a rate at or above the " *
            "bar of $bar, but its exact one-sided $(_pct(RECALL_ALPHA)) % lower bound $lb is " *
            "below it: the sample is too small to show that the bar is met" :
            "recall_n_recovered: $nr of $nf flagged taxa recovered, below the bar of $bar")
    end
    field = _flag(recall_pass_field)
    if field === nothing
        push!(reasons, "recall_pass_field: unknown, or not a Bool")
    elseif point !== nothing && field != point
        push!(reasons, "recall_pass_field: the record says $field, but its own counts give " *
                       "$point on the point estimate the field is taken on")
    end

    return (ok = isempty(reasons), recall_pass_derived = derived, reasons = reasons,
            recall_rate = rate, recall_lower_bound = lb)
end

"""
    background_removal_complete(; left_hmm_hit, left_ipr_ge500, recall_pass,
                                frame_matches, scan_evidence_ok) -> (complete, reasons)

Whether the comparator background really is free of the enzyme carriers it is supposed to
be free of, and if not, why not.

  `left_hmm_hit`      carriers with an HMM hit still present in the ANALYSED background arm.
  `left_ipr_ge500`    carriers found by the family signature at the length cut-off, likewise
                      counted over the analysed arm.
  `recall_pass`       the verdict of the prefilter recall check, as a Bool — the one
                      `scan_evidence_status` re-derives on the lower confidence bound, not
                      the record's own field, which is taken on the point estimate.
  `frame_matches`     whether the scanned taxon set is the analysed frame, as a Bool.
  `scan_evidence_ok`  the `ok` of `scan_evidence_status`: whether the scan's persisted
                      records vouch for the table the counts above were taken from.

`complete` is true only when all five are known and all five are clean. `reasons` names
every failure, so a caller emits the cause rather than just the verdict.

WHY FIVE INPUTS AND NOT A FLAG. A predicate that reduces to "the exclusion flag was set"
grants the licence on no measurement at all: one of its counts is structurally zero for any
input, the other is populated only for taxa that already have a hit, a recall verdict that is
printed and never persisted is a verdict nothing can check, and whether the scan covered the
frame is never asked. Such a predicate answers "complete" for a background that still holds
carriers, and for an input from which nothing has been removed.

WHY THE FIFTH INPUT. Three of the other four are computed from the same scan table, so a
table that found nothing satisfies all three at once: no carrier is left because none was
ever reported, and the taxon set still matches. A carrier whose row is missing from the
table reads as no carrier too, which is why the frame comparison stays. Only evidence
from outside the table can tell "measured clean" from "measured nothing", and that
evidence is the scan's own acceptance record, tied to the table by its digest.

Where the exclusion itself is "any HMM hit", `left_hmm_hit` is zero by construction over
the analysed arm, and the completeness of that half rests on the recall check and the scan
evidence rather than on the count.

SO EVERY INPUT IS REQUIRED, AND UNKNOWN IS NOT CLEAN. A count that is `missing`, `nothing`,
NaN, negative or non-numeric, and a verdict that is not a Bool, each make the answer false
with a reason of their own. The keywords default to `nothing` so that a caller which has
not measured one of them cannot obtain a licence by leaving it out.
"""
function background_removal_complete(; left_hmm_hit = nothing, left_ipr_ge500 = nothing,
                                     recall_pass = nothing, frame_matches = nothing,
                                     scan_evidence_ok = nothing)
    reasons = String[]
    for (key, value, what) in
            (("left_hmm_hit",   left_hmm_hit,   "still carry an HMM hit"),
             ("left_ipr_ge500", left_ipr_ge500, "carry a family-signature entry at the length cut-off"))
        n = _count(value)
        if n === nothing
            push!(reasons, "$key: unknown, or not a count")
        elseif n > 0
            push!(reasons, "$key: $(_fmt(n)) analysed background taxa $what")
        end
    end
    for (key, value, what) in
            (("recall_pass",      recall_pass,      "the prefilter recall check did not pass"),
             ("frame_matches",    frame_matches,    "the scanned taxon set is not the analysed frame"),
             ("scan_evidence_ok", scan_evidence_ok,
              "the scan's persisted records, or the recall verdict derived from them, do not " *
              "support the stronger statement"))
        ok = _flag(value)
        if ok === nothing
            push!(reasons, "$key: unknown, or not a Bool")
        elseif !ok
            push!(reasons, "$key: $what")
        end
    end
    return (complete = isempty(reasons), reasons = reasons)
end

# A count for a reason string. Only a whole number double-precision arithmetic holds exactly
# is printed as an integer; converting a larger one would throw while the reason is being
# written, and the caller would get a crash in place of the reason.
_fmt(n::Real) = (isinteger(n) && abs(n) <= maxintfloat(Float64)) ? string(Int(n)) : string(n)

"""
    background_licence_code(complete) -> Symbol

`:strong` when the background removal is complete, `:fallback` otherwise — the only two
codes the recorded rules allow for this claim.

Takes the Bool or the whole NamedTuple from `background_removal_complete`. Anything that is
not exactly `true` is `:fallback`, so an unknown or malformed verdict cannot license the
strong code by accident. This function holds no wording: the neutral code is the output, and
the mapping from code to phrase lives with the step that writes prose.
"""
background_licence_code(complete::Bool) = complete ? :strong : :fallback
background_licence_code(x::NamedTuple) =
    background_licence_code(hasproperty(x, :complete) && x.complete === true)
background_licence_code(::Any) = :fallback

end # module
