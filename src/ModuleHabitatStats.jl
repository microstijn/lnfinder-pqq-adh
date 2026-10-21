module HabitatStats

export crude, standardised, bh, contrast_rows, flag_weights, collapse_to_genus,
       well_supported_mask

using HypothesisTests, Statistics

#=
    HabitatStats — the enrichment estimators, the multiple-testing correction and the row
    builder that scripts/habitat_enrichment.jl calls.

    WHY THEY LIVE HERE. Estimators defined inside an analysis script cannot be unit-tested,
    and cannot be re-used by a sensitivity grid without copying them. Everything here is
    pure: plain vectors in, NamedTuples out, no DataFrames, no I/O, no global row buffer to
    push into. One numeric path serves both weight representations, so a Bool hit vector
    and the equivalent 0/1 Float vector give bitwise identical results; that equality is
    asserted in test/test_habitat_stats.jl.

    TWO ESTIMATORS, ALWAYS REPORTED TOGETHER.

      CRUDE         a/n1 / (b/n2), unadjusted, as a log2 risk ratio with a Katz Wald
                    interval and a Fisher exact p.
      STANDARDISED  the background's share of the category REWEIGHTED to the foreground's
                    phylum distribution: p0_std = sum_k w_k (c_k/n0_k), w_k = n1_k/n1.

    The gap between them is the phylogenetic contribution. Habitat annotation tracks
    phylum, so a crude ratio can report "enrichment" that is nothing but the phylogenetic
    range of the protein family; reweighting makes that contribution legible instead of
    leaving it to be argued about.

    FRACTIONAL HITS. A hit vector may be Bool (a unit is in the category or it is not) or
    Real (a unit is split across categories, e.g. one part per tied term). See the
    `standardised` docstring for what fractional weights do to the intervals.

    THE UNIT MUST BE THE SAME ON BOTH SIDES. Aggregating the foreground to unique
    genus-by-category PAIRS while the background stays one row per taxon lets a genus
    count once per category it occupies on one side and once in total on the other. That
    is not a like-for-like comparison, and such a contrast cannot be read as a
    conservative bound. Collapsing BOTH arms with `collapse_to_genus` keeps n1 and n2 in
    the same currency.
=#

const MIN_CELL = 5        # a+b below this is reported but never tested
const MIN_BG_K = 20       # collapse phyla thinner than this in the background ...
const MIN_FG_K = 10       # ... or thinner than this in the foreground

_sv(x) = ismissing(x) ? "" : String(strip(string(x)))

"""
    crude(a, n1, b, n2) -> NamedTuple

Unadjusted log2 risk ratio with a Katz Wald 95% interval and a Fisher exact p-value.

`a` and `b` are hit totals (sums of weights) and may be fractional; `n1`/`n2` are unit
counts. Returns `(log2_rr, ci_lo, ci_hi, fisher_p, continuity_corrected,
rounded_for_fisher)`.

The Haldane correction (+0.5 on the hit cells, +1 on the denominators) fires ONLY when a
cell is empty. Applying it unconditionally shrinks every estimate towards zero, which
would silently damp real contrasts; applying it never leaves a zero cell with an infinite
ratio and an undefined interval.

Fisher's exact test is defined on counts, so a fractional `a` or `b` is rounded to the
nearest integer for the p-value ONLY — the point estimate and the interval use the
unrounded sums. `rounded_for_fisher` is true exactly when that rounding had something to
do, i.e. when a hit total was not already integral, and the p-value is then approximate.
The complements are taken as `n1 - round(a)` so the margins of the 2x2 table stay
consistent with its cells.
"""
function crude(a::Real, n1::Integer, b::Real, n2::Integer)
    ai = clamp(round(Int, a), 0, n1)
    bi = clamp(round(Int, b), 0, n2)
    rounded = !(isinteger(a) && isinteger(b))
    p = pvalue(FisherExactTest(ai, n1 - ai, bi, n2 - bi))
    corrected = (a == 0 || b == 0)
    ac = corrected ? a + 0.5 : float(a); bc = corrected ? b + 0.5 : float(b)
    N1 = corrected ? n1 + 1 : n1;        N2 = corrected ? n2 + 1 : n2
    l2 = log2((ac / N1) / (bc / N2))
    se = sqrt(max(1/ac + 1/bc - 1/N1 - 1/N2, 0.0)) / log(2)
    return (log2_rr = l2, ci_lo = l2 - 1.96se, ci_hi = l2 + 1.96se, fisher_p = p,
            continuity_corrected = corrected, rounded_for_fisher = rounded)
end

"""
    standardised(fg_phy, fg_hit, bg_phy, bg_hit; min_fg_k, min_bg_k) -> NamedTuple

Phylum-standardised log2 risk ratio: the background's share of the category reweighted to
the foreground's phylum mix, p0_std = sum_k w_k (c_k/n0_k) with w_k = n1_k/n1.

`fg_phy`/`bg_phy` are per-unit phylum labels; `fg_hit`/`bg_hit` are per-unit weights,
either Bool or Real. Returns `(log2_rr, ci_lo, ci_hi, n_strata, p1, p0, var,
n_zero_stratum_cells, min_bg_stratum_cell)`, where `var` is the variance of the log ratio
on the natural-log scale; the interval is that variance carried to log2 and taken 1.96
standard errors either side.

STRATUM HYGIENE. A phylum thinner than `min_fg_k` in the foreground or `min_bg_k` in the
background is relabelled "Other phyla" in BOTH arms before any counting, so a single-taxon
stratum cannot carry a large weight against a near-empty background cell. A stratum with
foreground units but no background units is dropped entirely; `n_strata` reports how many
strata actually contributed.

VARIANCE, AND WHY FRACTIONAL WEIGHTS MAKE IT CONSERVATIVE. The delta-method variance
`(1/a - 1/n1) + v0/p0^2` treats the weights as fixed and both hit totals as binomial: the
first term is the variance of log(a/n1) for a ~ Bin(n1, p1), and `v0` sums the binomial
variances of the per-stratum background shares. For a fractional weight vector those
totals are sums of bounded values in [0,1] rather than counts of Bernoulli draws, and such
a sum has variance no larger than the binomial one with the same mean. The formula is
therefore an UPPER BOUND on the true variance, and the intervals it returns are
conservative — wider than they need to be, never narrower. It is deliberately left alone
rather than adjusted, because adjusting it would need a model for how the weights were
split, and an interval that is too wide is the safer error.

THIN AND EMPTY BACKGROUND CELLS. That bound concerns the form of the variance, not how well
its plug-in value is estimated, and it says nothing about a thin cell. A stratum whose
background cell c_k is empty adds nothing to p0 and nothing to v0, because the plug-in
binomial variance of a share estimated as zero is zero, so the stratum that knows least
about the background contributes no uncertainty at all. A cell of one or two units is
better only in degree. The row totals cannot reveal either case: when one stratum
carries most of the weight, a background total `b` well clear of any gating threshold can
sit next to that heavy stratum's own cell being empty. Two disclosures therefore travel with every estimate,
taken over the strata that contributed (after the collapse, so a merged "Other phyla" is one
cell):

  - `n_zero_stratum_cells`, how many contributing strata have c_k == 0;
  - `min_bg_stratum_cell`, the smallest c_k among them, NaN when none contributed.

`min_bg_stratum_cell` is always a Float64 and is never rounded. A cell built from
fractional weights may hold a third of a unit, and rounding it would report an empty cell
that is not empty, or a whole unit that does not exist. Keeping one type also keeps a Bool
hit vector and its 0/1 Float equivalent bitwise identical on this field, as on the others.
The disclosures change no arithmetic: p0, v0 and the interval are computed exactly as
before.
"""
function standardised(fg_phy::AbstractVector, fg_hit::AbstractVector,
                      bg_phy::AbstractVector, bg_hit::AbstractVector;
                      min_fg_k::Integer = MIN_FG_K, min_bg_k::Integer = MIN_BG_K)
    n1 = length(fg_hit); a = sum(fg_hit)
    nanres(used, nzero, cmin) = (log2_rr = NaN, ci_lo = NaN, ci_hi = NaN, n_strata = used,
                                 p1 = NaN, p0 = NaN, var = NaN,
                                 n_zero_stratum_cells = nzero, min_bg_stratum_cell = cmin)
    (a == 0 || n1 == 0) && return nanres(0, 0, NaN)
    # A phylum column read straight from a CSV is Union{Missing,String}, so normalise here
    # rather than trusting every caller to do it: a grid variant that calls standardised
    # directly would otherwise hit a MethodError instead of the "Other phyla" collapse.
    # This is a no-op for the Vector{String} contrast_rows already passes.
    fg_phy = _sv.(fg_phy); bg_phy = _sv.(bg_phy)
    # stratum hygiene: collapse thin phyla so a single-unit stratum cannot dominate
    cnt_f = Dict{String,Int}(); cnt_b = Dict{String,Int}()
    for p in fg_phy; cnt_f[p] = get(cnt_f, p, 0) + 1; end
    for p in bg_phy; cnt_b[p] = get(cnt_b, p, 0) + 1; end
    lab(p) = (get(cnt_f, p, 0) >= min_fg_k && get(cnt_b, p, 0) >= min_bg_k) ? p : "Other phyla"
    F = lab.(fg_phy); B = lab.(bg_phy)
    strata = sort(unique(F))
    p1 = a / n1
    p0 = 0.0; v0 = 0.0; used = 0
    cells = Float64[]                              # c_k of each contributing stratum
    for k in strata
        n1k = count(==(k), F); n1k == 0 && continue
        idx = findall(==(k), B); n0k = length(idx)
        n0k == 0 && continue                       # no background for this stratum -> drop
        ck = sum(bg_hit[i] for i in idx)
        w = n1k / n1
        p0 += w * (ck / n0k)
        v0 += w^2 * ck * (n0k - ck) / n0k^3
        used += 1
        push!(cells, float(ck))
    end
    nzero = count(iszero, cells)
    cmin = isempty(cells) ? NaN : minimum(cells)
    (p0 <= 0 || p1 <= 0) && return nanres(used, nzero, cmin)
    l2 = log2(p1 / p0)
    var = (1/a - 1/n1) + v0 / p0^2
    se = sqrt(max(var, 0.0)) / log(2)
    return (log2_rr = l2, ci_lo = l2 - 1.96se, ci_hi = l2 + 1.96se, n_strata = used,
            p1 = p1, p0 = p0, var = var,
            n_zero_stratum_cells = nzero, min_bg_stratum_cell = cmin)
end

"""
    bh(ps) -> Vector{Float64}

Benjamini-Hochberg q-values, returned in the order the p-values were given. The running
minimum taken from the largest p-value downwards keeps them monotone, so a q-value can
never exceed that of a less significant test.
"""
function bh(ps::AbstractVector{<:Real})
    n = length(ps); ord = sortperm(ps); q = zeros(n); prev = 1.0
    for k in n:-1:1
        i = ord[k]; prev = min(prev, ps[i] * n / k); q[i] = prev
    end
    return q
end

"""
    flag_weights(values; splitcombo=false) -> Dict{String,Vector{Bool}}

Turn a column of category labels into one 0/1 weight vector per category.

`values` holds one label per unit; blanks and `missing` contribute nothing. With
`splitcombo=true` a label is split on `+` and the unit counts in every part, which is how
non-exclusive flags ("Hypersaline+Thermal") are read: those categories are not mutually
exclusive, and a unit belongs to each one it carries.

This helper cannot produce fractional weights — it reads an already-decided label. Build
the Dict directly (per-unit weights summing to 1 across the tied categories, say) to hand
`contrast_rows` a fractional unit.
"""
function flag_weights(values::AbstractVector; splitcombo::Bool = false)
    n = length(values)
    w = Dict{String,Vector{Bool}}()
    for (i, v) in enumerate(values)
        s = _sv(v)
        isempty(s) && continue
        for c in (splitcombo ? String.(split(s, '+')) : [s])
            haskey(w, c) || (w[c] = falses(n))
            w[c][i] = true
        end
    end
    return w
end

"""
    contrast_rows(fg_phy, fg_w, bg_phy, bg_w; contrast, unit, min_cell, min_fg_k, min_bg_k)

One result row per category, carrying both estimators and a BH q-value taken across the
categories of this contrast.

`fg_phy`/`bg_phy` are per-unit phylum labels and `fg_w`/`bg_w` map a category name to its
per-unit weight vector (see `flag_weights`); a category missing from one arm scores zero
there. `contrast` and `unit` are copied into every row so that several contrasts can be
concatenated into one table. `unit` is what states the currency of n1 and n2 ("taxa",
"genera", "sequences"), and rows carrying different units must never be compared with each
other.

Categories with no hits on either side are skipped. `tested` marks the rows whose combined
hit total reaches `min_cell`; only those enter the BH correction, and the rest carry a NaN
q-value rather than a number that would invite a comparison the data cannot support.

`n_zero_stratum_cells` and `min_bg_stratum_cell` are copied from `standardised` and are
appended AFTER every other column, so a reader that addresses columns by position still
finds the earlier ones where they were. `tested` and `continuity_corrected` look only at
the row totals `a` and `b`; these two columns are what show a row whose totals are healthy
while the standardised interval rests on an empty or near-empty stratum cell.
"""
function contrast_rows(fg_phy::AbstractVector, fg_w::AbstractDict,
                       bg_phy::AbstractVector, bg_w::AbstractDict;
                       contrast::AbstractString, unit::AbstractString,
                       min_cell::Real = MIN_CELL,
                       min_fg_k::Integer = MIN_FG_K, min_bg_k::Integer = MIN_BG_K)
    n1 = length(fg_phy); n2 = length(bg_phy)
    fphy = _sv.(fg_phy); bphy = _sv.(bg_phy)
    cats = sort(collect(union(keys(fg_w), keys(bg_w))))
    buf = NamedTuple[]
    for c in cats
        fhit = get(fg_w, c, falses(n1)); bhit = get(bg_w, c, falses(n2))
        length(fhit) == n1 ||
            error("foreground weights for '$c' cover $(length(fhit)) of $n1 units")
        length(bhit) == n2 ||
            error("background weights for '$c' cover $(length(bhit)) of $n2 units")
        a = sum(fhit); b = sum(bhit)
        (a + b) == 0 && continue
        cr = crude(a, n1, b, n2)
        st = standardised(fphy, fhit, bphy, bhit; min_fg_k, min_bg_k)
        push!(buf, (contrast = contrast, unit = unit, habitat = c,
                    a = a, n1 = n1, b = b, n2 = n2, pct_fg = 100a/n1, pct_bg = 100b/n2,
                    log2_rr = cr.log2_rr, ci_lo = cr.ci_lo, ci_hi = cr.ci_hi,
                    log2_rr_std = st.log2_rr, ci_lo_std = st.ci_lo, ci_hi_std = st.ci_hi,
                    n_strata = st.n_strata,
                    fisher_p = cr.fisher_p, q_bh = NaN,
                    tested = (a + b) >= min_cell,
                    continuity_corrected = cr.continuity_corrected,
                    rounded_for_fisher = cr.rounded_for_fisher,
                    n_zero_stratum_cells = st.n_zero_stratum_cells,
                    min_bg_stratum_cell = st.min_bg_stratum_cell))
    end
    isempty(buf) && return NamedTuple[]
    tested = [r.tested for r in buf]
    q = fill(NaN, length(buf))
    any(tested) && (q[tested] = bh([r.fisher_p for r in buf[tested]]))
    return NamedTuple[merge(r, (q_bh = q[i],)) for (i, r) in enumerate(buf)]
end

"""
    collapse_to_genus(genus, phylum, weights) -> NamedTuple

Collapse per-taxon rows to one row per genus, so both arms of a contrast can be counted in
the same currency.

A genus's weight in a category is the MEAN of its taxa's weights there, i.e. the share of
the genus that carries the category — so a genus is worth one unit in total however many
of its taxa were sequenced. Its phylum is the modal label of its taxa, ties resolved
alphabetically so the result does not depend on row order. Rows with a blank or missing
genus are dropped and counted in `n_dropped`.

Returns `(genus, phylum, weights, n_taxa, n_dropped)`, genera in sorted order, `weights`
being a Dict of category -> per-genus weight ready for `contrast_rows`.

The alternative — taking unique genus-by-category pairs — inflates an arm by the number of
categories each genus occupies, and is comparable only against an arm built the same way.
"""
function collapse_to_genus(genus::AbstractVector, phylum::AbstractVector,
                           weights::AbstractDict)
    length(genus) == length(phylum) ||
        error("genus and phylum differ in length ($(length(genus)) vs $(length(phylum)))")
    members = Dict{String,Vector{Int}}()
    n_dropped = 0
    for i in eachindex(genus)
        g = _sv(genus[i])
        if isempty(g)
            n_dropped += 1
            continue
        end
        push!(get!(members, g, Int[]), i)
    end
    gens = sort(collect(keys(members)))
    function modal(labels)
        cnt = Dict{String,Int}()
        for l in labels; cnt[l] = get(cnt, l, 0) + 1; end
        best = ""; bestn = -1
        for k in sort(collect(keys(cnt)))
            cnt[k] > bestn && (best = k; bestn = cnt[k])
        end
        return best
    end
    phy = [modal([_sv(phylum[i]) for i in members[g]]) for g in gens]
    w = Dict{String,Vector{Float64}}()
    for (cat, v) in weights
        length(v) == length(genus) ||
            error("weights for '$cat' cover $(length(v)) of $(length(genus)) rows")
        w[String(cat)] = [mean(float(v[i]) for i in members[g]) for g in gens]
    end
    return (genus = gens, phylum = phy, weights = w,
            n_taxa = [length(members[g]) for g in gens], n_dropped = n_dropped)
end

"""
    well_supported_mask(n_distinct; min_distinct=3) -> BitVector

Keep only the units whose call rests on at least `min_distinct` usable distinct terms.

A call resting on one or two terms is a single literature mention rather than a
distribution over evidence. Re-running a contrast on this subset answers whether a result
depends on the thinly annotated tail. `missing` counts as unsupported.
"""
well_supported_mask(n_distinct::AbstractVector; min_distinct::Integer = 3) =
    BitVector(!ismissing(x) && x >= min_distinct for x in n_distinct)

end # module
