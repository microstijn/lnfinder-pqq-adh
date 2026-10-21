#=
    validate_classifier.jl  — Control-set validation for the D-x-D Ln/Ca classifier
    ---------------------------------------------------------------------------------
    Purpose
    -------
    The whole LnFinder dataset rests on one equivalence (Voutsinos et al. 2025,
    doi:10.1093/ismejo/wraf057):

        D-x-D     at the calibrated active-site column  =>  Lanthanide
        D-x-[ATS] at the calibrated active-site column  =>  Calcium

    Nothing in the pipeline currently *checks* that this call is correct on proteins
    whose cofactor is experimentally known. This harness does exactly that. It also
    exposes a structural fragility in `calibrate_active_site`: it returns the FIRST
    model column where ref-1 is D-x-D and ref-2 is D-x-[ATS]. If more than one such
    "collision" column exists, the first one is not necessarily the catalytic one,
    and every downstream call inherits that mistake. `diagnose_calibration` reports
    how many collision columns exist so you can see whether the choice is unique.

    What it proves / does NOT prove
    -------------------------------
    - PROVES: whether D-x-D, read at the calibrated column, recovers the known
      cofactor on a labelled control set (sensitivity / specificity / accuracy).
    - DOES NOT prove functional Ln-dependence of any *uncharacterised* hit — that
      still needs structure (coordination geometry) and ultimately wet-lab assays.

    Human-in-the-loop steps (cannot be done from code alone)
    --------------------------------------------------------
    1. Confirm the labels in dat/control_set/known_cofactor_mdhs.tsv (domain expert).
    2. Provide an `hmmalign` binary (HMMER not on PATH here) OR a pre-aligned AFA.
    3. Network access for the UniProt sequence fetch.

    Pure-logic functions (`model_coords`, `diagnose_calibration`, `evaluate_controls`)
    are unit-tested offline in test/test_validate_classifier.jl — run that first.
=#

module ValidateClassifier

using DataFrames
using CSV
using FASTX
using HTTP
using CodecZlib
using JSON3
using SHA
using Printf
using HypothesisTests: BinomialTest, confint

# Robust GET → (String, response). UniProt may gzip the body even when not asked;
# request identity encoding and gunzip defensively if a gzip magic header still
# appears. The response is handed back too, because the retrieval manifest records the
# `X-UniProt-Release` header, and because a 200 with an empty body is meaningful here.
function _http_fetch(url::AbstractString)
    resp = HTTP.get(url, ["Accept-Encoding" => "identity"]; status_exception = false, readtimeout = 30)
    resp.status == 200 || return ("", resp)
    body = resp.body
    if length(body) >= 2 && body[1] == 0x1f && body[2] == 0x8b
        return (String(transcode(GzipDecompressor, body)), resp)
    end
    return (String(body), resp)
end

_http_text(url::AbstractString) = first(_http_fetch(url))

# Reuse the *production* calibration/parsing code so we test the real thing.
include("ModuleHmmer.jl")
using .ModuleHmmer: calibrate_active_site, parse_hmmalign_results

export model_coords, diagnose_calibration, evaluate_controls
export fetch_control_sequences, run_hmmalign, reorder_references!, validate, validate_from_afa, win_to_wsl
export resolve_control_accessions
export is_heldout, heldout_mask, stratified_summary, domain_counts, switch_residue_counts
export identity_clusters

# ---------------------------------------------------------------------------
# Pure helpers (offline-testable)
# ---------------------------------------------------------------------------

"""
    model_coords(seq::AbstractString) -> String

Strip insert states (lowercase) exactly as the production pipeline does, leaving
the HMM match-state string (uppercase + '-'). Mirrors the filter used inside
`calibrate_active_site` / `parse_hmmalign_results` so diagnostics agree with calls.
"""
model_coords(seq::AbstractString) = filter(c -> isuppercase(c) || c == '-', seq)

"""
    collision_columns(ref_ln_model::AbstractString, ref_ca_model::AbstractString) -> Vector{Int}

Return EVERY model index `i` where ref-1 is D-x-D and ref-2 is D-x-[ATS] — i.e.
every column `calibrate_active_site` could legitimately have returned. The
production code returns `first(...)`; if this vector has length > 1 the calibration
is ambiguous and the chosen column should be confirmed against a crystal structure.
"""
function collision_columns(ref_ln_model::AbstractString, ref_ca_model::AbstractString)
    @assert length(ref_ln_model) == length(ref_ca_model) "reference model strings differ in length"
    hits = Int[]
    for i in 1:(length(ref_ln_model) - 2)
        is_ln = ref_ln_model[i] == 'D' && ref_ln_model[i+2] == 'D'
        is_ca = ref_ca_model[i] == 'D' && ref_ca_model[i+2] in ('A', 'T', 'S')
        (is_ln && is_ca) && push!(hits, i)
    end
    return hits
end

"""
    diagnose_calibration(afa_path) -> NamedTuple

Read the first two records of an AFA (seq1 = Ln ref, seq2 = Ca ref, the order
`calibrate_active_site` assumes) and report all collision columns and the one the
production code will pick. Prints a warning when the choice is not unique.
"""
function diagnose_calibration(afa_path::AbstractString)
    seqs = String[]
    FASTA.Reader(open(afa_path)) do reader
        for (i, rec) in enumerate(reader)
            push!(seqs, FASTA.sequence(String, rec))
            i >= 2 && break
        end
    end
    length(seqs) < 2 && error("AFA has < 2 sequences; cannot diagnose calibration.")

    ref_ln = model_coords(seqs[1])
    ref_ca = model_coords(seqs[2])
    cols = collision_columns(ref_ln, ref_ca)
    chosen = isempty(cols) ? nothing : first(cols)

    println("── Calibration diagnostic ──")
    if isempty(cols)
        println("  ❌ No D-x-D / D-x-[ATS] collision column found — calibration will ERROR.")
    elseif length(cols) == 1
        @printf("  ✅ Unique collision column at model index %d (Ln=%s | Ca=%s).\n",
                chosen, ref_ln[chosen:chosen+2], ref_ca[chosen:chosen+2])
    else
        @printf("  ⚠️  %d collision columns: %s — production picks the FIRST (%d).\n",
                length(cols), string(cols), chosen)
        println("     The chosen column is NOT guaranteed catalytic. Confirm against the")
        println("     K-Ln / K-Ca crystal structure before trusting downstream calls.")
    end
    return (chosen = chosen, collisions = cols, ref_ln = ref_ln, ref_ca = ref_ca)
end

# ---------------------------------------------------------------------------
# Coordinate guard: an active-site index is only meaningful in model coordinates
# ---------------------------------------------------------------------------

"""
    _hmm_length(hmm_path) -> Int

Number of match states of a HMMER3 profile, read from its `LENG` header line. That
number is also the number of columns any `hmmalign` alignment to this profile must
have, so it is the yardstick every match-state coordinate is measured against.

Deliberately standalone — a path in, an integer out, no other state from this file —
so a second pipeline aligning to the same profile can assert the same invariant on
its own alignment without pulling in the validation flow.

Only the header of the first model is read; the scan stops where the emission table
(or the model terminator) begins, so a stray `LENG`-like line inside the body cannot
be mistaken for the header field. A profile without a `LENG` line raises rather than
returning a silent zero.
"""
function _hmm_length(hmm_path::AbstractString)
    open(hmm_path) do io
        for line in eachline(io)
            m = match(r"^LENG\s+(\d+)\s*$", line)
            m === nothing || return parse(Int, m.captures[1])
            (startswith(line, "HMM ") || startswith(line, "//")) && break
        end
        error("no LENG line in the HMMER profile header of $hmm_path — " *
              "the model length cannot be read, so match coordinates cannot be checked.")
    end
end

"""
    _afa_match_length(afa_path) -> Int

Number of HMM match columns in an aligned FASTA: the length of its first record once
insert states (lowercase) are stripped, i.e. exactly what `model_coords` leaves. Every
record of an `hmmalign` AFA shares that count, so the first record speaks for the file.
"""
function _afa_match_length(afa_path::AbstractString)
    seq = ""
    FASTA.Reader(open(afa_path)) do reader
        for rec in reader
            seq = FASTA.sequence(String, rec)
            break
        end
    end
    isempty(seq) && error("$afa_path holds no sequence record; cannot count match columns.")
    return length(model_coords(seq))
end

"""
    _assert_match_columns(afa_path, expected_len) -> Int

Fail unless the alignment carries exactly `expected_len` match columns, and print the
two numbers when it does.

The active-site column is an INDEX into the match states of the model. An alignment
built against a different profile, or one whose columns were re-derived after the fact,
still answers every index — it simply answers at a different residue in every sequence,
and the classifier reads a position that is not the catalytic one. Nothing throws, no
count changes, and the calls are quietly wrong. Comparing the alignment's match count
with the model's `LENG` is the cheapest way to make that mismatch loud.
"""
function _assert_match_columns(afa_path::AbstractString, expected_len::Integer)
    got = _afa_match_length(afa_path)
    got == expected_len || error(
        "match-column mismatch: $afa_path has $got match columns, the model has " *
        "$expected_len (LENG). Column indices read from this alignment are NOT model " *
        "coordinates, so the calibrated active-site index would point at the wrong " *
        "residue. Re-align against the same profile, or pass the hmm/expected_len " *
        "that actually describes this alignment.")
    @printf("  ✅ match columns = %d (HMM LENG %d)\n", got, expected_len)
    return got
end

# Path of the profile this harness is built around; used only as the fallback yardstick
# when a caller supplies neither an HMM nor an expected length.
const _DEFAULT_HMM = normpath(joinpath(@__DIR__, "..", "dat", "models",
                                       "methanol_dehydrogenase_pqq_xoxF_mxaF.hmm"))

# Match-column count of the profile above. Used only when that file cannot be read, so
# that a missing model downgrades the guard's provenance, never the guard itself.
const _MODEL_LENG = 589

"""
    _resolve_expected_len(hmm, expected_len) -> Union{Int,Nothing}

Decide how many match columns an alignment must have: an explicit `expected_len` wins,
then the `LENG` of an explicitly named `hmm`. Given neither, fall back to the profile
above and say so in a warning.

There is deliberately no path that returns `nothing`. Scoring an alignment whose columns
were never checked is the failure this guard exists to close, so when even the bundled
profile is missing the harness falls back to `_MODEL_LENG`, the length that profile has
carried since the census was built. A caller aligning against a different model has to
say so, by passing `hmm` or `expected_len`; it cannot happen by omission.
"""
function _resolve_expected_len(hmm::Union{Nothing,AbstractString},
                               expected_len::Union{Nothing,Integer})
    expected_len === nothing || return Int(expected_len)
    hmm === nothing || return _hmm_length(hmm)
    if isfile(_DEFAULT_HMM)
        len = _hmm_length(_DEFAULT_HMM)
        @warn "no hmm/expected_len given — guarding against the default profile" model = _DEFAULT_HMM leng = len
        return len
    end
    @warn "no hmm/expected_len given and the default profile is missing — guarding against " *
          "the length it has carried since the census was built" model = _DEFAULT_HMM leng = _MODEL_LENG
    return _MODEL_LENG
end

"""
    evaluate_controls(pred_df, labels) -> NamedTuple

Join predicted classifications (`accession`, `classification`) to known labels
(`accession`, `cofactor` ∈ {"Lanthanide","Calcium"}) and compute a confusion matrix
plus sensitivity (Ln recall), specificity (Ca recall) and accuracy.

Predictions of :Gap / :Truncated / :Unknown are kept in the table (they are
informative failures) but counted as *misclassifications* against their known label.

A labelled accession with no prediction is the silent-n-shrink failure mode: the inner
join would simply drop the row, and every rate below would be computed on a smaller
denominator than the control set actually has. Those accessions are collected BEFORE
the join. `strict = true` (the default) throws an `ArgumentError` naming them;
`strict = false` prints them and returns them in the `missing` field.
"""
function evaluate_controls(pred_df::DataFrame, labels::DataFrame; strict::Bool = true)
    # Compared exactly as the join compares them (no trimming): a key that only looks
    # equal is a key the join would drop.
    predicted = Set(string.(pred_df.accession))
    missing_accs = unique([string(a) for a in labels.accession if !(string(a) in predicted)])
    if !isempty(missing_accs)
        detail = "labelled accession(s) with no prediction: " * join(missing_accs, ", ")
        strict && throw(ArgumentError(detail *
            " — they would shrink n silently. Retrieve or remove them, or pass " *
            "strict = false to score the remainder knowingly."))
        println("  ⚠ $(length(missing_accs)) $detail")
    end

    joined = innerjoin(labels, pred_df, on = :accession, makeunique = true)
    isempty(joined) && error("No overlap between control labels and predictions on :accession.")

    pred_sym(c) = c isa Symbol ? c : Symbol(c)
    joined.pred = pred_sym.(joined.classification)
    joined.known = Symbol.(joined.cofactor)
    joined.correct = joined.pred .== joined.known

    # Confusion counts on the two real classes.
    tp = count(r -> r.known == :Lanthanide && r.pred == :Lanthanide, eachrow(joined)) # Ln called Ln
    fn = count(r -> r.known == :Lanthanide && r.pred != :Lanthanide, eachrow(joined)) # Ln missed
    tn = count(r -> r.known == :Calcium    && r.pred == :Calcium,    eachrow(joined)) # Ca called Ca
    fp = count(r -> r.known == :Calcium    && r.pred == :Lanthanide, eachrow(joined)) # Ca called Ln

    sens = (tp + fn) == 0 ? NaN : tp / (tp + fn)
    spec = (tn + fp) == 0 ? NaN : tn / (tn + fp)
    acc  = count(joined.correct) / nrow(joined)

    println("── Control-set evaluation (n = $(nrow(joined))) ──")
    @printf("  Ln  recall (sensitivity): %.3f   (%d/%d)\n", sens, tp, tp + fn)
    @printf("  Ca  recall (specificity): %.3f   (%d/%d)\n", spec, tn, tn + fp)
    @printf("  Overall accuracy:         %.3f   (%d/%d)\n", acc, count(joined.correct), nrow(joined))
    bad = filter(r -> !r.correct, joined)
    if nrow(bad) > 0
        println("  Misclassified / non-call:")
        for r in eachrow(bad)
            @printf("    %-14s known=%-11s pred=%-11s site=%s\n",
                    r.accession, r.known, r.pred, get(r, :active_site, "?"))
        end
    end
    return (table = joined, sensitivity = sens, specificity = spec, accuracy = acc,
            tp = tp, fp = fp, tn = tn, fn = fn, missing = missing_accs)
end

# ---------------------------------------------------------------------------
# Tiered, stratified scoring
# ---------------------------------------------------------------------------

"""
    REFERENCE_ROLE

The one `role` value that takes a control row OUT of the held-out set. Every other value,
including no value at all, leaves the row held out.
"""
const REFERENCE_ROLE = "reference"

"""
    CONTESTED_TIER

The `tier` value of a control row whose published cofactor assignment is disputed. Such a
row is listed in the stratified table together with the call the classifier made for it,
and is never part of a numerator or a denominator.
"""
const CONTESTED_TIER = "X"

_norm(x) = lowercase(strip(string(x)))
_norm_tier(x) = uppercase(strip(string(x)))

"""
    is_heldout(role) -> Bool

Held out means exactly `role != "reference"`, compared case- and whitespace-insensitively.

This is a named function rather than an inline comparison because it is load-bearing. A
control row is given the `reference` role precisely when the classifier's calibration was
derived from it: the calibrated column is read off that row, so the row cannot also be
evidence that the column is right. Scoring it as held out would grade the rule on the
sequence that set the rule. Every held-out stratum below is therefore built from this one
predicate — a reference row cannot reach a held-out stratum through any tier, cofactor or
domain, because no stratum re-derives the test for itself.

A control table with no `role` column has nothing demoted, so every row is held out.
"""
is_heldout(role) = _norm(role) != REFERENCE_ROLE

# A control table may or may not carry role / tier / domain / organism columns. An absent
# column reads as empty strings instead of throwing, so a minimal accession+cofactor set
# still scores, and `missing` cells read as empty rather than as the text "missing".
function _text_col(df::AbstractDataFrame, col::Symbol)
    String(col) in names(df) || return fill("", nrow(df))
    return String[String(strip(string(coalesce(v, "")))) for v in df[!, col]]
end

"""
    heldout_mask(df) -> BitVector

`true` for every row of a labelled control table that is held out, per [`is_heldout`](@ref).
Every held-out stratum is an AND with this mask, which is what makes the role check
impossible to bypass.
"""
heldout_mask(df::AbstractDataFrame) = BitVector(is_heldout.(_text_col(df, :role)))

"""
    _label_table(pred_df, labels) -> DataFrame

Join predictions to labels for LISTING: the same join `evaluate_controls` scores on, plus
`pred`, `known` and `correct`, but no rates and no printing. Rows that no held-out rate counts —
the calibration references and the contested tier — still have to appear in the written
table, which is why this is separate from the scoring path. Labelled accessions with no
prediction are named in a warning; `evaluate_controls` is the one that refuses to score
without them.
"""
function _label_table(pred_df::DataFrame, labels::AbstractDataFrame)
    predicted = Set(string.(pred_df.accession))
    absent = unique([string(a) for a in labels.accession if !(string(a) in predicted)])
    isempty(absent) || @warn "label rows without a prediction are absent from the listing" accessions = absent
    joined = innerjoin(labels, pred_df, on = :accession, makeunique = true)
    isempty(joined) && error("No overlap between control labels and predictions on :accession.")
    joined.pred = [c isa Symbol ? c : Symbol(c) for c in joined.classification]
    joined.known = Symbol.(joined.cofactor)
    joined.correct = joined.pred .== joined.known
    return joined
end

"""
    _clopper_pearson(k, n; level = 0.95) -> (lo, hi)

Exact binomial (Clopper–Pearson) interval for `k` successes out of `n`, or
`(missing, missing)` for an empty stratum.

Exact rather than normal-approximate because the strata here are small and often all
correct, and the Wald interval of a k == n stratum has zero width — it would report a
certainty the data cannot carry. The exact lower bound is the number worth quoting.
"""
function _clopper_pearson(k::Integer, n::Integer; level::Real = 0.95)
    n == 0 && return (missing, missing)
    lo, hi = confint(BinomialTest(Int(k), Int(n)); level = level, method = :clopper_pearson)
    return (lo, hi)
end

const _EXPECTED_DOMAINS = ("Bacteria", "Archaea", "Eukaryota")

"""
    domain_counts(joined) -> DataFrame(domain, n)

Composition of a labelled table by taxonomic `domain`, with EVERY expected domain listed
even when it holds no rows. A domain that is simply missing from a printout reads as "not
looked at", whereas a printed zero is a statement about the set — and the scope of a
control set is exactly the kind of claim a reader has to be able to check. Values not in
the expected list are appended; rows with no domain are counted as `unassigned`.
"""
function domain_counts(joined::AbstractDataFrame)
    tally = Dict{String,Int}()
    for v in _text_col(joined, :domain)
        key = isempty(v) ? "unassigned" : v
        for d in _EXPECTED_DOMAINS
            if lowercase(key) == lowercase(d)
                key = d
                break
            end
        end
        tally[key] = get(tally, key, 0) + 1
    end
    out = DataFrame(domain = String[], n = Int[])
    for d in _EXPECTED_DOMAINS
        push!(out, (d, get(tally, d, 0)))
    end
    for d in sort(collect(keys(tally)))
        d in _EXPECTED_DOMAINS || push!(out, (d, tally[d]))
    end
    return out
end

const _SWITCH_RESIDUES = ("A", "S", "T")

"""
    switch_residue_counts(joined) -> DataFrame(residue, n)

Tally of the third ("switch") residue of the active-site triad over the rows whose KNOWN
cofactor is calcium. The Ca call is the disjunction D-x-[ATS], so this says which of its
three branches the control set actually exercises: a branch with no rows behind it has not
been tested, and is listed with a zero rather than left out.
"""
function switch_residue_counts(joined::AbstractDataFrame)
    known = "known" in names(joined) ? joined.known : Symbol.(joined.cofactor)
    sites = _text_col(joined, :active_site)
    tally = Dict{String,Int}()
    for (k, s) in zip(known, sites)
        k == :Calcium || continue
        res = length(s) >= 3 ? string(s[3]) : "?"
        tally[res] = get(tally, res, 0) + 1
    end
    out = DataFrame(residue = String[], n = Int[])
    for r in _SWITCH_RESIDUES
        push!(out, (r, get(tally, r, 0)))
    end
    for r in sort(collect(keys(tally)))
        r in _SWITCH_RESIDUES || push!(out, (r, tally[r]))
    end
    return out
end

# Identity over the columns where at least one of the two records carries a residue.
function _pair_identity(a::Vector{Char}, b::Vector{Char})
    same = 0
    comparable = 0
    for i in 1:min(length(a), length(b))
        (a[i] == '-' && b[i] == '-') && continue
        comparable += 1
        a[i] == b[i] && (same += 1)
    end
    return comparable == 0 ? 0.0 : same / comparable
end

# Key a record the way `parse_hmmalign_results` keys it, so the two agree on what an
# accession is (second '|' field, else the leading alphanumeric run).
function _accession_key(header::AbstractString)
    id = String(first(split(strip(String(header)))))
    if occursin('|', id)
        parts = split(id, '|')
        return length(parts) >= 2 ? String(parts[2]) : String(parts[1])
    end
    m = match(r"^([A-Z0-9]+)", id)
    return m === nothing ? id : String(m.captures[1])
end

"""
    identity_clusters(afa_path, accessions; threshold = 0.9) -> (clusters, n_eff)

Single-linkage clusters of the named records at `threshold` or more pairwise identity over
the match columns, plus `n_eff`, the number of clusters.

An interval computed from n rows assumes n independent draws. Near-identical proteins are
closer to one draw than to several, so a control set that contains a group of them has
fewer effective observations than rows, and the interval is narrower than the evidence
warrants. `n_eff` is the number that belongs next to any interval quoted from this set.
"""
function identity_clusters(afa_path::AbstractString, accessions; threshold::Real = 0.9)
    want = Set(String[String(strip(string(a))) for a in accessions])
    ids = String[]
    seqs = Vector{Vector{Char}}()
    FASTA.Reader(open(afa_path)) do reader
        for rec in reader
            key = _accession_key(FASTA.identifier(rec))
            (key in want && !(key in ids)) || continue
            push!(ids, key)
            push!(seqs, collect(model_coords(FASTA.sequence(String, rec))))
        end
    end
    n = length(ids)
    parent = collect(1:n)
    function root(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    for i in 1:n, j in (i + 1):n
        _pair_identity(seqs[i], seqs[j]) >= threshold || continue
        ri, rj = root(i), root(j)
        ri == rj || (parent[max(ri, rj)] = min(ri, rj))
    end
    groups = Dict{Int,Vector{String}}()
    for i in 1:n
        push!(get!(groups, root(i), String[]), ids[i])
    end
    clusters = [sort(g) for g in values(groups)]
    sort!(clusters, by = first)
    return (clusters = clusters, n_eff = length(clusters))
end

"""
    stratified_summary(joined; level = 0.95, quiet = false) -> DataFrame

One row per stratum — `stratum`, `n`, `correct`, `accuracy`, `ci_lo`, `ci_hi`, `note` —
with exact (Clopper–Pearson) intervals at `level`. `joined` is a labelled table carrying
`correct` (see `_label_table` or the `table` field of `evaluate_controls`).

The strata, in order:
  - `reference`   — the calibration rows (`role == "reference"`), reported apart because
                    they are not evidence about held-out proteins;
  - `heldout_<t>` — one per tier present among the held-out rows, plus the cofactor split
                    of tier A, so the strongest evidence can be read on its own;
  - `tier_X`      — the contested rows, listed with the calls they received and never
                    scored;
  - `heldout_all` — the pooled held-out rows, flagged as mixed evidence: tiers differ in
                    how directly the cofactor was established, so one number over all of
                    them is a summary, not the headline.

A single pooled accuracy hides all of that, which is the reason this function exists.
Also prints the composition by `domain` and, for the known-Ca rows, the switch-residue
tally.
"""
function stratified_summary(joined::AbstractDataFrame; level::Real = 0.95, quiet::Bool = false)
    "correct" in names(joined) ||
        error("stratified_summary needs a `correct` column; pass the table from " *
              "`evaluate_controls` or `_label_table`.")
    heldout = heldout_mask(joined)
    tiers = _norm_tier.(_text_col(joined, :tier))
    contested = tiers .== CONTESTED_TIER
    known = "known" in names(joined) ? Symbol.(joined.known) : Symbol.(joined.cofactor)
    correct = BitVector(joined.correct)

    rows = DataFrame(stratum = String[], n = Int[], correct = Union{Int,Missing}[],
                     accuracy = Union{Float64,Missing}[], ci_lo = Union{Float64,Missing}[],
                     ci_hi = Union{Float64,Missing}[], note = String[])
    function scored_row!(name, mask; note = "")
        n = count(mask)
        k = count(correct .& mask)
        lo, hi = _clopper_pearson(k, n; level = level)
        push!(rows, (name, n, n == 0 ? missing : k, n == 0 ? missing : k / n, lo, hi, note))
    end

    # Contested rows are subtracted from the reference stratum too: "never counted" has to
    # mean never, not "never, except in the one stratum nobody reads".
    scored_row!("reference", .!heldout .& .!contested;
                note = "calibration source — not held out")

    for t in sort(unique(tiers[heldout .& .!contested]))
        isempty(t) && continue
        mask = heldout .& (tiers .== t)
        scored_row!("heldout_$t", mask)
        if t == "A"
            scored_row!("heldout_A_Ca", mask .& (known .== :Calcium))
            scored_row!("heldout_A_Ln", mask .& (known .== :Lanthanide))
        end
    end

    # The contested rows are reported by the call they received, which is the only thing
    # that can be said about them without settling the dispute over their label.
    call_col = count(contested) == 0 ? nothing :
               "pred" in names(joined) ? joined.pred :
               "classification" in names(joined) ? joined.classification : nothing
    calls = call_col === nothing ? String[] : sort(unique(String.(string.(call_col[contested]))))
    push!(rows, ("tier_$CONTESTED_TIER", count(contested), missing, missing, missing, missing,
                 "contested label — listed, never scored" *
                 (isempty(calls) ? "" : "; calls: " * join(calls, ", "))))

    scored_row!("heldout_all", heldout .& .!contested; note = "mixed evidence (tiers pooled)")

    if !quiet
        pct = round(Int, 100 * level)
        println("── Stratified summary (exact $(pct)% CIs) ──")
        @printf("  %-14s %5s %8s %9s  %-16s %s\n", "stratum", "n", "correct", "accuracy",
                "$(pct)% CI", "note")
        for r in eachrow(rows)
            acc = r.accuracy === missing ? "        —" : @sprintf("%9.3f", r.accuracy)
            ci  = r.ci_lo === missing ? "—" : @sprintf("[%.3f, %.3f]", r.ci_lo, r.ci_hi)
            cor = r.correct === missing ? "—" : string(r.correct)
            @printf("  %-14s %5d %8s %s  %-16s %s\n", r.stratum, r.n, cor, acc, ci, r.note)
        end
        println("  Domain composition:")
        for r in eachrow(domain_counts(joined))
            @printf("    %-12s %d\n", r.domain, r.n)
        end
        println("  Ca switch residues (known-Ca rows):")
        for r in eachrow(switch_residue_counts(joined))
            @printf("    %-12s %d\n", r.residue, r.n)
        end
    end
    return rows
end

# The per-row audit trail written next to the summary. Every column appears whether or not
# the control table carried it, so the file has one shape across runs and a reader can tell
# an empty field from a column that was never there.
const _LABELLED_COLUMNS = ("accession", "cofactor", "role", "tier", "domain", "organism",
                           "active_site", "classification", "correct")

function _labelled_output(labelled::AbstractDataFrame)
    out = DataFrame()
    for col in _LABELLED_COLUMNS
        out[!, col] = col in names(labelled) ? collect(labelled[!, col]) : fill("", nrow(labelled))
    end
    return out
end

# ---------------------------------------------------------------------------
# Network / external-tool steps (guarded; run against your real data)
# ---------------------------------------------------------------------------

# Sequence retrieval: UniProtKB, with a UniParc fallback for withdrawn accessions and
# an optional NCBI protein route.

const _NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

"""
    _json_object(json) -> Union{JSON3.Object,Nothing}

Parse a JSON body into an object, or `nothing` for an empty body, a non-object, or
anything unparseable. Keeps the callers below free of try/catch noise.
"""
function _json_object(json::AbstractString)
    isempty(strip(json)) && return nothing
    obj = try
        JSON3.read(json)
    catch
        return nothing
    end
    return obj isa JSON3.Object ? obj : nothing
end

"""
    _inactive_uniparc_id(json) -> Union{String,Nothing}

Return the UniParc id of a UniProtKB entry that has been withdrawn, given that entry's
JSON view. An accession deleted from UniProtKB still answers `/uniprotkb/<acc>.fasta`
with HTTP 200 and a ZERO-BYTE body, so the JSON view is the only place the withdrawal
is visible: `entryType` reads "Inactive" and `extraAttributes.uniParcId` carries the
UPI whose sequence archive still holds the protein. Returns `nothing` for an active
entry, an empty body, or a body without that field.
"""
function _inactive_uniparc_id(json::AbstractString)
    obj = _json_object(json)
    obj === nothing && return nothing
    get(obj, :entryType, "") == "Inactive" || return nothing
    extra = get(obj, :extraAttributes, nothing)
    extra === nothing && return nothing
    upi = get(extra, :uniParcId, nothing)
    (upi === nothing || isempty(String(upi))) && return nothing
    return String(upi)
end

"""
    _inactive_reason(json) -> String

The `inactiveReason.inactiveReasonType` of a UniProtKB entry JSON ("DELETED",
"MERGED", ...), or "unknown". Recorded in the rewritten header so the provenance of a
fallback sequence stays readable in the alignment itself.
"""
function _inactive_reason(json::AbstractString)
    obj = _json_object(json)
    obj === nothing && return "unknown"
    reason = get(obj, :inactiveReason, nothing)
    reason === nothing && return "unknown"
    kind = get(reason, :inactiveReasonType, nothing)
    return kind === nothing ? "unknown" : String(kind)
end

"""
    _rewrite_header(acc, alt_id, tag; db = "up") -> String

Build a FASTA header keyed by the CONTROL-SET accession, e.g.
`>up|<acc>|<UPI> uniparc (UniProtKB deleted: DELETED)`.

`parse_hmmalign_results` keys a record on the second '|' field, so a sequence fetched
under a different identifier (a UniParc UPI, an NCBI protein id) MUST be re-keyed here:
left alone the row would be keyed by that foreign id, would never join to its label,
and would disappear from the scored set without anything failing.
"""
function _rewrite_header(acc::AbstractString, alt_id::AbstractString, tag::AbstractString;
                         db::AbstractString = "up")
    head = ">" * db * "|" * strip(acc) * "|" * strip(alt_id)
    return isempty(strip(tag)) ? head : head * " " * strip(tag)
end

"""
    _first_fasta_record(body) -> (header, sequence)

Split a FASTA body into its first record's header (without the leading '>') and its
unwrapped sequence. Returns `("", "")` when the body holds no record.
"""
function _first_fasta_record(body::AbstractString)
    lines = split(replace(String(body), "\r\n" => "\n"), '\n')
    i = findfirst(l -> startswith(l, ">"), lines)
    isnothing(i) && return ("", "")
    header = String(lines[i][2:end])
    buf = IOBuffer()
    for line in lines[(i + 1):end]
        startswith(line, ">") && break
        write(buf, strip(line))
    end
    return (header, String(take!(buf)))
end

"""
    _write_fasta_record(io, header_line, seq; width = 60)

Write one wrapped FASTA record. `header_line` may or may not carry its own '>'.
"""
function _write_fasta_record(io::IO, header_line::AbstractString, seq::AbstractString;
                             width::Int = 60)
    println(io, startswith(header_line, ">") ? header_line : ">" * header_line)
    for i in 1:width:length(seq)
        println(io, seq[i:min(i + width - 1, lastindex(seq))])
    end
    return nothing
end

"""
    _source_map(accessions, sources) -> Dict{String,String}

Normalise the `sources` argument of `fetch_control_sequences` — `nothing`, a Dict-like
accession => database, or a vector parallel to `accessions` — into a lookup.
Recognised databases: "uniprot" (the default for an unlisted accession) and "ncbi".
"""
function _source_map(accessions, sources)
    out = Dict{String,String}()
    sources === nothing && return out
    if sources isa AbstractDict
        for (acc, db) in sources
            out[String(strip(string(acc)))] = lowercase(strip(string(db)))
        end
    else
        length(sources) == length(accessions) ||
            error("`sources` must be a Dict, or a vector as long as `accessions`.")
        for (acc, db) in zip(accessions, sources)
            out[String(strip(string(acc)))] = lowercase(strip(string(db)))
        end
    end
    return out
end

# UniProtKB first; on the zero-byte body of a withdrawn entry, fall back to UniParc.
function _fetch_uniprot(acc::AbstractString)
    text, resp = _http_fetch("https://rest.uniprot.org/uniprotkb/$(acc).fasta")
    release = HTTP.header(resp, "X-UniProt-Release", "")
    if !isempty(strip(text))
        header, seq = _first_fasta_record(text)
        isempty(seq) && return nothing
        return (db = "uniprotkb", upi = "", header = ">" * header, seq = seq, release = release)
    end

    json, jresp = _http_fetch("https://rest.uniprot.org/uniprotkb/$(acc)?format=json")
    isempty(release) && (release = HTTP.header(jresp, "X-UniProt-Release", ""))
    upi = _inactive_uniparc_id(json)
    upi === nothing && return nothing

    reason = _inactive_reason(json)
    archived, uresp = _http_fetch("https://rest.uniprot.org/uniparc/$(upi).fasta")
    _, seq = _first_fasta_record(archived)
    isempty(seq) && return nothing
    isempty(release) && (release = HTTP.header(uresp, "X-UniProt-Release", ""))
    return (db = "uniparc", upi = upi,
            header = _rewrite_header(acc, upi, "uniparc (UniProtKB deleted: $reason)"),
            seq = seq, release = release)
end

# NCBI protein, for control rows whose sequence has no UniProtKB entry at all.
function _fetch_ncbi(acc::AbstractString)
    url = "$(_NCBI_EFETCH)?db=protein&id=$(HTTP.escapeuri(acc))&rettype=fasta&retmode=text"
    text, _ = _http_fetch(url)
    header, seq = _first_fasta_record(text)
    isempty(seq) && return nothing
    parts = split(header, ' ', limit = 2)
    ncbi_id = String(parts[1])
    description = length(parts) > 1 ? String(parts[2]) : ""
    return (db = "ncbi", upi = "",
            header = _rewrite_header(acc, ncbi_id, "ncbi ($description)"; db = "ncbi"),
            seq = seq, release = "")
end

"""
    fetch_control_sequences(accessions; out_fasta, sources = nothing) -> DataFrame

Download protein sequences for the labelled control set into `out_fasta` and return the
retrieval manifest: one row per accession with `accession`, `db`, `upi`, `header`,
`length`, `seq_sha256` and `uniprot_release` (the `X-UniProt-Release` response header),
so a later re-run can be compared against exactly what was scored. The digest is
SHA-256 because that is the hash the dependency set already provides.

UniProtKB is the default source. A withdrawn accession answers with HTTP 200 and an
empty body, so an empty body is followed up with the JSON view; an inactive entry is
re-fetched from UniParc and its header rewritten to carry the control-set accession
(see `_rewrite_header`). Accessions marked "ncbi" in `sources` (a Dict
accession => database, or a vector parallel to `accessions`) are fetched from NCBI
efetch instead.

Errors if ANY accession yields no sequence: a labelled row that cannot be retrieved has
to stop the run, because skipping it would quietly shrink the denominator of every
number computed downstream. Prints the header returned for each one so a domain expert
can confirm the label actually matches the intended (gene, organism). Requires network.
"""
function fetch_control_sequences(accessions::Vector{<:AbstractString};
                                 out_fasta::AbstractString,
                                 sources = nothing)
    source_of = _source_map(accessions, sources)
    manifest = DataFrame(accession = String[], db = String[], upi = String[],
                         header = String[], length = Int[], seq_sha256 = String[],
                         uniprot_release = String[])
    unretrieved = String[]
    open(out_fasta, "w") do io
        for acc in unique(String.(strip.(accessions)))
            isempty(acc) && continue
            db = get(source_of, acc, "uniprot")
            record = db == "ncbi" ? _fetch_ncbi(acc) : _fetch_uniprot(acc)
            if record === nothing
                @warn "no sequence retrieved" accession = acc source = db
                push!(unretrieved, acc)
                continue
            end
            _write_fasta_record(io, record.header, record.seq)
            push!(manifest, (acc, record.db, record.upi, record.header,
                             length(record.seq), bytes2hex(sha256(record.seq)),
                             record.release))
            println("  fetched $acc  [$(record.db)]  →  $(record.header)")
        end
    end
    isempty(unretrieved) || error(
        "No sequence for $(length(unretrieved)) labelled accession(s): " *
        join(unretrieved, ", ") * ". Fix or remove those rows — dropping them here " *
        "would shrink the control set without anything failing.")
    println("✅ control sequences written to $out_fasta  (n = $(nrow(manifest)))")
    return manifest
end

"""
    win_to_wsl(path) -> String

Translate a Windows path (`C:\\Users\\...`) to its WSL mount form (`/mnt/c/Users/...`).
Paths with spaces are fine: they are passed as single argv elements by the Cmd object,
so no manual quoting is needed. Non-Windows paths are returned unchanged.
"""
function win_to_wsl(path::AbstractString)
    q = replace(String(path), '\\' => '/')
    m = match(r"^([A-Za-z]):/(.*)$", q)
    return m === nothing ? q : "/mnt/" * lowercase(m.captures[1]) * "/" * m.captures[2]
end

# ---------------------------------------------------------------------------
# Accession resolver: turn (gene, organism) into a proposed UniProt accession
# for a human to confirm. Removes the "type 9 accessions by hand" step.
# ---------------------------------------------------------------------------

const _UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# Genus + species only (drop strain) — used as a looser fallback organism query.
_genus_species(org::AbstractString) = join(split(strip(org))[1:min(2, length(split(strip(org))))], " ")

# Map a gene symbol to a protein-name keyword for the looser fallback search.
function _mdh_keyword(gene::AbstractString)
    g = lowercase(gene)
    (occursin("xox", g) || occursin("mxa", g) || occursin("mdh", g) || occursin("mox", g)) && return "methanol dehydrogenase"
    (occursin("adh", g) || occursin("exa", g)) && return "alcohol dehydrogenase"
    occursin("gcd", g) && return "glucose dehydrogenase"
    return nothing
end

# The classifier needs the catalytic LARGE subunit (carries the active-site D-x-D),
# not the small subunit (mxaI / "subunit 2 / beta"). Bias selection accordingly.
function _subunit_bonus(protein_name::AbstractString)
    n = lowercase(protein_name)
    (occursin("small", n) || occursin("beta", n) || occursin("subunit 2", n) ||
        occursin("associated peptide", n)) && return -3
    (occursin("large", n) || occursin("alpha", n) || occursin("subunit 1", n)) && return 1
    return 0
end

# An intact PQQ-dependent dehydrogenase of this family is far longer than this. An entry
# below the bar is a fragment, an accessory peptide, or a different protein, and none of
# those can be aligned across the profile that the active-site column indexes into.
const _MIN_INTACT_LENGTH = 400

# UniProt's `fragment` field is empty for a complete entry and carries a word otherwise.
# Anything neither empty nor explicitly negative counts as a fragment, so a change of
# wording upstream fails towards caution rather than towards a silent zero penalty.
const _NOT_FRAGMENT = ("", "no", "false", "0", "n", "-")
_is_fragment(v) = !(lowercase(strip(string(coalesce(v, "")))) in _NOT_FRAGMENT)

# Sequence length as an integer, or `nothing` when the field is absent or unparseable: an
# unknown length must not be scored as though it were a short one.
function _candidate_length(v)
    v === missing && return nothing
    v isa Integer && return Int(v)
    return tryparse(Int, strip(string(v)))
end

_candidate_lengths(df::AbstractDataFrame) =
    "length" in names(df) ? Any[_candidate_length(v) for v in df[!, :length]] :
                            Any[nothing for _ in 1:nrow(df)]

"""
    _rank_candidates(df, organism) -> DataFrame

Rank the UniProt hits for one control row: the same columns plus `score`, best first.

Pure — a table in, a table out, no network — so the ranking can be exercised offline on a
fixture of the exact hits that have to come out in a given order.

The score sums independent statements about a hit:
  * `+2` reviewed (Swiss-Prot): a curated entry over an automatic one;
  * `+1` the hit organism contains the target genus and species;
  * `±`  the subunit bias of [`_subunit_bonus`](@ref) — the catalytic large subunit carries
    the active-site triad, the small subunit does not;
  * `-5` the entry is a FRAGMENT: a partial sequence need not contain the active-site
    column at all, so it cannot answer the question the control set asks;
  * `-2` fewer than `_MIN_INTACT_LENGTH` residues.

The two penalties are what reviewed-plus-organism alone cannot express: to that score a
short, unreviewed fragment from the right organism ranks exactly as well as a complete
entry, and a fragment that reaches the control set carries its label onto a sequence that
need never have been aligned over the catalytic column.

Ties keep the order UniProt returned, so equally good hits are never reshuffled.
"""
function _rank_candidates(df::AbstractDataFrame, organism::AbstractString)
    gs = lowercase(_genus_species(organism))
    reviewed  = _text_col(df, :reviewed)
    organisms = _text_col(df, :organism)
    proteins  = _text_col(df, :protein_name)
    fragments = _text_col(df, :fragment)
    lengths   = _candidate_lengths(df)
    score = Int[(lowercase(reviewed[i]) == "reviewed" ? 2 : 0) +
                (occursin(gs, lowercase(organisms[i])) ? 1 : 0) +
                _subunit_bonus(proteins[i]) +
                (_is_fragment(fragments[i]) ? -5 : 0) +
                (lengths[i] !== nothing && lengths[i] < _MIN_INTACT_LENGTH ? -2 : 0)
                for i in 1:nrow(df)]
    perm = sortperm(score; rev = true)
    ranked = df[perm, :]
    ranked.score = score[perm]
    return ranked
end

"""
    _search_tiers(gene, organism, protein_keyword) -> Vector{Tuple{String,String}}

The ordered `(mode, query)` pairs `_uniprot_search` tries, stopping at the first that
returns anything. The keyword tiers appear only when a protein-name keyword is known:

  1. `name_reviewed`       — reviewed entries whose PROTEIN NAME matches the keyword, in
     the full organism;
  2. `gene`                — gene symbol + full organism;
  3. `name_reviewed+genus` — as 1, on genus/species (drops the strain);
  4. `gene+genus`          — gene symbol + genus/species;
  5. `name`                — protein-name keyword + full organism, reviewed or not;
  6. `name+genus`          — protein-name keyword + genus/species.

Two orderings are at work, and they are independent. Within one organism specificity the
reviewed protein-name question is asked BEFORE the gene one, because the gene tiers cannot
see an entry that carries no gene name at all: a curated, full-length entry with an empty
`gene_names` field never enters the candidate list, the search stops at the first tier that
answered, and the row is resolved from whatever unreviewed material happens to carry the
symbol. Asking Swiss-Prot by protein name first puts the curated entry in front of
[`_rank_candidates`](@ref), which is the only place it can be weighed at all.

Then the exact organism is tried before genus/species, and that is why the reviewed tier
needs both forms. UniProt writes a strain as a parenthesised qualifier rather than as part
of the organism name, so a control row naming its strain matches NOTHING under the exact
form and everything under the looser one — a reviewed-only tier that existed solely in the
exact form would silently never fire for exactly the rows it was added for.

Name-based matches are looser than gene-based ones in one specific way: a protein-name
keyword is shared by the whole enzyme family, so a hit can be the right family and the
wrong paralogue. `_uniprot_search` reports the mode for that reason, and every proposal is
flagged for human confirmation.
"""
function _search_tiers(gene::AbstractString, organism::AbstractString,
                       protein_keyword::Union{Nothing,AbstractString})
    gs = _genus_species(organism)
    reviewed_name(org) = "(reviewed:true) AND (protein_name:\"$protein_keyword\") " *
                         "AND (organism_name:\"$org\")"
    tiers = Tuple{String,String}[]
    protein_keyword === nothing || push!(tiers, ("name_reviewed", reviewed_name(organism)))
    push!(tiers, ("gene", "(gene:$gene) AND (organism_name:\"$organism\")"))
    protein_keyword === nothing || push!(tiers, ("name_reviewed+genus", reviewed_name(gs)))
    push!(tiers, ("gene+genus", "(gene:$gene) AND (organism_name:\"$gs\")"))
    if protein_keyword !== nothing
        push!(tiers, ("name",       "(protein_name:\"$protein_keyword\") AND (organism_name:\"$organism\")"))
        push!(tiers, ("name+genus", "(protein_name:\"$protein_keyword\") AND (organism_name:\"$gs\")"))
    end
    return tiers
end

"""
    _uniprot_search(gene, organism; size, protein_keyword) -> (df, mode)

Search UniProt through the tiers of [`_search_tiers`](@ref) and return the first tier that
hits, tagged with the `mode` used, so callers can flag the looser, name-based matches for
extra scrutiny.

`fragment` and `length` are requested alongside the descriptive fields because
[`_rank_candidates`](@ref) scores on them; without those two columns a partial entry is
indistinguishable from a complete one.
"""
function _uniprot_search(gene::AbstractString, organism::AbstractString;
                         size::Int = 5, protein_keyword::Union{Nothing,AbstractString} = nothing)
    fields = "accession,protein_name,gene_names,organism_name,reviewed,length,fragment"
    function query(qstr)
        url = "$(_UNIPROT_SEARCH)?query=$(HTTP.escapeuri(qstr))&format=tsv&fields=$fields&size=$size"
        body = _http_text(url)
        isempty(strip(body)) && return DataFrame()
        df = CSV.read(IOBuffer(body), DataFrame; delim = '\t')
        nrow(df) == 0 && return DataFrame()
        # Returned columns are human-readable; we requested a fixed order, so rename by position.
        target = [:accession, :protein_name, :gene_names, :organism, :reviewed, :length, :fragment]
        rename!(df, names(df) .=> target[1:ncol(df)])
        return df
    end
    for (mode, q) in _search_tiers(gene, organism, protein_keyword)
        df = query(q)
        nrow(df) > 0 && return (df = df, mode = mode)
    end
    return (df = DataFrame(), mode = "")
end

"""
    resolve_control_accessions(control_tsv; out_tsv=nothing, per_query=10) -> DataFrame

For every row in the control TSV whose `accession` is "TODO"/empty, search UniProt by
(`gene`, `organism`) through the tiers of [`_search_tiers`](@ref) and propose the best
candidate under [`_rank_candidates`](@ref): reviewed entries and organism-name matches up,
fragments and short entries down. Prints each proposal with its UniProt protein name +
organism so a domain expert can CONFIRM (or pick from the listed alternatives) rather than
type accessions by hand. Rows that already have an accession (e.g. P16027) are left as-is.

Every proposal is a proposal. The search can only say which entry best matches a name and
an organism; whether that protein is the one whose cofactor was measured is a judgement no
query makes, which is why nothing here writes a label.

If `out_tsv` is given, writes a resolved copy with accessions filled in (still verify it).
Returns the proposals table. Requires network.
"""
function resolve_control_accessions(control_tsv::AbstractString;
                                    out_tsv::Union{Nothing,AbstractString} = nothing,
                                    per_query::Int = 10)
    labels = CSV.read(control_tsv, DataFrame; delim = '\t', comment = "#")
    @assert all(in(names(labels)), ["accession", "gene", "organism"]) "TSV needs accession,gene,organism columns"
    resolved = copy(labels)
    resolved.accession = String.(string.(resolved.accession))  # widen from inferred String7
    proposals = DataFrame(row = Int[], gene = String[], target_organism = String[],
                          cofactor = String[], role = String[], accession = String[],
                          hit_protein = String[], hit_organism = String[],
                          reviewed = String[], match_mode = String[],
                          n_alternatives = Int[], alternatives = String[])

    for (i, r) in enumerate(eachrow(labels))
        acc = strip(String(r.accession))
        (acc != "TODO" && !isempty(acc)) && continue          # already known — keep it
        gene = strip(String(r.gene)); org = strip(String(r.organism))
        res = _uniprot_search(gene, org; size = (per_query == 0 ? 5 : per_query),
                              protein_keyword = _mdh_keyword(gene))
        df, mode = res.df, res.mode

        cof  = hasproperty(r, :cofactor) ? String(r.cofactor) : ""
        role = hasproperty(r, :role) ? String(r.role) : ""
        if nrow(df) == 0
            push!(proposals, (i, gene, org, cof, role, "NOT_FOUND", "", "", "", "", 0, ""))
            continue
        end
        ranked = _rank_candidates(df, org)
        best = ranked[1, :]
        alts = join(["$(row.accession)($(row.reviewed))" for row in eachrow(ranked[2:min(end, 4), :])], "; ")
        push!(proposals, (i, gene, org, cof, role, String(best.accession),
                          String(best.protein_name), String(best.organism),
                          String(best.reviewed), mode, nrow(ranked) - 1, alts))
        resolved[i, :accession] = String(best.accession)
    end

    println("\n── Proposed control accessions (CONFIRM before use; ⚠ = looser name-based match) ──")
    for p in eachrow(proposals)
        flag = startswith(p.match_mode, "name") ? "⚠" : " "
        @printf("  %s[%2d] %-7s %-30s → %-12s %-10s %-11s %s\n",
                flag, p.row, p.gene, first(p.target_organism, 30), p.accession,
                p.reviewed, p.match_mode, p.hit_protein)
        isempty(p.alternatives) || println("            alternatives: ", p.alternatives)
    end

    if out_tsv !== nothing
        open(out_tsv, "w") do io
            println(io, "# Auto-resolved control set — VERIFY every accession against the printed name/organism.")
            println(io, join(names(resolved), '\t'))
            for row in eachrow(resolved)
                println(io, join([string(row[c]) for c in names(resolved)], '\t'))
            end
        end
        println("\n✅ wrote resolved copy to $out_tsv")
    end
    return proposals
end

"""
    run_hmmalign(hmm, fasta; hmmalign_bin="hmmalign", wsl=false, out_afa) -> out_afa

Align control sequences to the model so match-state coordinates are meaningful
(a generic MSA would NOT reproduce the HMM column numbering the pipeline relies on).

HMMER has no native Windows build. Set `wsl=true` to invoke a HMMER installed in WSL
(Ubuntu `apt install hmmer`); the three paths are auto-translated to `/mnt/<drive>`.
Otherwise `hmmalign_bin` is run directly (e.g. on Linux/macOS or a server).

`--amino` asserts that both files are protein, so alphabet autodetection cannot reach a
different answer for a short or unusual input. There is no option asking hmmalign for
the full set of model columns, so the column count is checked after the fact instead:
`_assert_match_columns` compares it with the model's `LENG`.
"""
function run_hmmalign(hmm::AbstractString, fasta::AbstractString;
                      hmmalign_bin::AbstractString = "hmmalign",
                      wsl::Bool = false,
                      out_afa::AbstractString)
    cmd = if wsl
        h, f, o = win_to_wsl(hmm), win_to_wsl(fasta), win_to_wsl(out_afa)
        `wsl $hmmalign_bin --amino --trim --outformat afa -o $o $h $f`
    else
        `$hmmalign_bin --amino --trim --outformat afa -o $out_afa $hmm $fasta`
    end
    println("▶ ", cmd)
    run(cmd)
    return out_afa
end

"""
    reorder_references!(afa_path, ln_ref_acc, ca_ref_acc)

`calibrate_active_site` assumes record 1 = Ln ref and record 2 = Ca ref. hmmalign
preserves input order, but this rewrites the AFA to force the two references to the
top regardless, so calibration is deterministic.
"""
function reorder_references!(afa_path::AbstractString, ln_ref_acc::AbstractString, ca_ref_acc::AbstractString)
    recs = NamedTuple{(:id, :seq), Tuple{String, String}}[]
    FASTA.Reader(open(afa_path)) do reader
        for rec in reader
            push!(recs, (id = FASTA.identifier(rec), seq = FASTA.sequence(String, rec)))
        end
    end
    keyacc(id) = occursin('|', id) ? split(id, '|')[2] : split(id, ' ')[1]
    ln = findfirst(r -> occursin(ln_ref_acc, keyacc(r.id)), recs)
    ca = findfirst(r -> occursin(ca_ref_acc, keyacc(r.id)), recs)
    (isnothing(ln) || isnothing(ca)) && error("Reference accessions not found in AFA: $ln_ref_acc / $ca_ref_acc")
    order = vcat(ln, ca, setdiff(1:length(recs), [ln, ca]))
    open(afa_path, "w") do io
        for i in order
            write(io, ">", recs[i].id, "\n", recs[i].seq, "\n")
        end
    end
    return afa_path
end

"""
    validate(; control_tsv, hmm, work_dir, hmmalign_bin, ln_ref_acc, ca_ref_acc) -> NamedTuple

End-to-end driver: load labels → fetch sequences → hmmalign → check match columns
against the model's `LENG` → reorder refs → diagnose calibration → calibrate → parse →
evaluate (all rows, then the held-out rows) → stratify. A `source` column in the control
TSV (e.g. "ncbi") routes individual rows to another database; `include_tiers` keeps only
the listed evidence tiers.

Writes <work_dir>/control_fetch_manifest.csv, control_predictions.csv,
control_predictions_labelled.csv (one row per labelled protein) and control_summary.csv
(one row per stratum, with exact CIs), and returns the manifest as `fetch_manifest`.
"""
function validate(; control_tsv::AbstractString,
                    hmm::AbstractString,
                    work_dir::AbstractString = mktempdir(),
                    hmmalign_bin::AbstractString = "hmmalign",
                    wsl::Bool = false,
                    ln_ref_acc::AbstractString,
                    ca_ref_acc::AbstractString,
                    exclude_roles::Vector{<:AbstractString} = String[],
                    include_tiers::Union{Nothing,AbstractVector} = nothing)
    mkpath(work_dir)
    labels = CSV.read(control_tsv, DataFrame; delim = '\t', comment = "#")
    @assert all(in(names(labels)), ["accession", "cofactor"]) "control TSV needs accession,cofactor columns"
    accs = String.(strip.(labels.accession))
    accs = filter(a -> !isempty(a) && a != "TODO", accs)

    fasta = joinpath(work_dir, "control.fasta")
    afa   = joinpath(work_dir, "control_aligned.afa")
    sources = "source" in names(labels) ?
        Dict(String(strip(string(r.accession))) => String(strip(string(r.source))) for r in eachrow(labels)) :
        nothing
    manifest = fetch_control_sequences(accs; out_fasta = fasta, sources = sources)
    CSV.write(joinpath(work_dir, "control_fetch_manifest.csv"), manifest)
    run_hmmalign(hmm, fasta; hmmalign_bin = hmmalign_bin, wsl = wsl, out_afa = afa)
    scored = _score_afa(afa, labels, work_dir, ln_ref_acc, ca_ref_acc;
                        exclude_roles = exclude_roles, include_tiers = include_tiers,
                        expected_len = _hmm_length(hmm))
    return merge(scored, (fetch_manifest = manifest,))
end

"""
    validate_from_afa(; control_afa, control_tsv, ln_ref_acc, ca_ref_acc, work_dir,
                        hmm = nothing, expected_len = nothing) -> NamedTuple

Score an ALREADY-ALIGNED control AFA — decouples HMMER entirely from the critical
path. Align your control sequences wherever HMMER lives (WSL, a server, anywhere),
then hand the resulting `.afa` here. Checks the match-column count, reorders the two
references to the top, diagnoses calibration, calibrates, parses and scores against
`control_tsv` — as a whole, on the held-out rows alone, and stratum by stratum.

Because the alignment arrives from elsewhere, the profile it was built against is the
one thing this entry point cannot see. Name it with `hmm` (its `LENG` is then the
yardstick) or give `expected_len` outright; with neither, the harness falls back to its
default profile and warns, so an alignment in foreign coordinates still cannot be
scored in silence.
"""
function validate_from_afa(; control_afa::AbstractString,
                             control_tsv::AbstractString,
                             ln_ref_acc::AbstractString,
                             ca_ref_acc::AbstractString,
                             work_dir::AbstractString = mktempdir(),
                             exclude_roles::Vector{<:AbstractString} = String[],
                             include_tiers::Union{Nothing,AbstractVector} = nothing,
                             hmm::Union{Nothing,AbstractString} = nothing,
                             expected_len::Union{Nothing,Integer} = nothing)
    mkpath(work_dir)
    labels = CSV.read(control_tsv, DataFrame; delim = '\t', comment = "#")
    afa = joinpath(work_dir, "control_aligned.afa")
    cp(control_afa, afa; force = true)
    return _score_afa(afa, labels, work_dir, ln_ref_acc, ca_ref_acc;
                      exclude_roles = exclude_roles, include_tiers = include_tiers,
                      expected_len = _resolve_expected_len(hmm, expected_len))
end

# Shared scoring tail: check columns → reorder refs → diagnose → calibrate → parse →
# label → evaluate twice → stratify → write the labelled predictions and the summary.
#
# Two evaluations are always computed, because they answer different questions:
#   `all`     — every scorable labelled row, calibration references included;
#   `heldout` — `role != reference` only (see `is_heldout`), the rows the classifier was
#               not built from, and therefore the only ones that test it.
# Contested rows (tier X) are outside both, and are listed by `stratified_summary` with
# the calls they received.
#
# `include_tiers` keeps only the listed tiers (default: every tier). `exclude_roles` is
# the older, coarser switch and still works; it can only remove rows, never promote one
# into a rate. `expected_len` is the model's match-state count; that check runs before
# anything reads a column index, and before the file is rewritten.
function _score_afa(afa, labels, work_dir, ln_ref_acc, ca_ref_acc;
                    exclude_roles = String[],
                    include_tiers::Union{Nothing,AbstractVector} = nothing,
                    level::Real = 0.95,
                    identity_threshold::Real = 0.9,
                    expected_len::Union{Nothing,Integer} = nothing)
    if expected_len === nothing
        @warn "no expected match-column count given — coordinate guard SKIPPED for $afa"
    else
        _assert_match_columns(afa, expected_len)
    end
    reorder_references!(afa, ln_ref_acc, ca_ref_acc)
    diag = diagnose_calibration(afa)
    idx  = calibrate_active_site(afa)
    preds = parse_hmmalign_results(afa, idx)
    CSV.write(joinpath(work_dir, "control_predictions.csv"), preds)

    kept = labels
    if include_tiers !== nothing
        if "tier" in names(kept)
            want = Set(_norm_tier.(include_tiers))
            n_before = nrow(kept)
            kept = kept[in.(_norm_tier.(_text_col(kept, :tier)), Ref(want)), :]
            println("(keeping tiers $(join(sort(collect(want)), ", ")): $n_before → $(nrow(kept)) labelled rows)")
        else
            @warn "include_tiers given but the control table has no `tier` column — every row kept"
        end
    end
    if !isempty(exclude_roles) && "role" in names(kept)
        drop = Set(lowercase.(exclude_roles))
        n_before = nrow(kept)
        kept = kept[.!in.(_norm.(_text_col(kept, :role)), Ref(drop)), :]
        println("(excluding roles $(collect(drop)): $n_before → $(nrow(kept)) labelled rows)")
    end

    scorable = kept[_norm_tier.(_text_col(kept, :tier)) .!= CONTESTED_TIER, :]
    n_contested = nrow(kept) - nrow(scorable)
    n_contested > 0 && println("($n_contested contested tier-$CONTESTED_TIER row(s): listed, never scored)")

    println("── all labelled rows (calibration references included) ──")
    ev_all = evaluate_controls(preds, scorable)

    heldout_labels = scorable[heldout_mask(scorable), :]
    ev_heldout = if nrow(heldout_labels) == 0
        @warn "no row has a role other than \"$REFERENCE_ROLE\" — the held-out evaluation is empty"
        nothing
    elseif nrow(heldout_labels) == nrow(scorable)
        println("── held-out rows (role != $REFERENCE_ROLE): the same rows as above ──")
        ev_all
    else
        println("── held-out rows (role != $REFERENCE_ROLE) ──")
        evaluate_controls(preds, heldout_labels)
    end

    labelled = _label_table(preds, kept)
    summary = stratified_summary(labelled; level = level)

    redundancy = identity_clusters(afa, heldout_labels.accession; threshold = identity_threshold)
    @printf("  Held-out records: n = %d, clusters at ≥ %d%% identity: n_eff = %d\n",
            nrow(heldout_labels), round(Int, 100 * identity_threshold), redundancy.n_eff)
    for c in redundancy.clusters
        length(c) > 1 && println("    near-identical: ", join(c, ", "))
    end

    CSV.write(joinpath(work_dir, "control_predictions_labelled.csv"), _labelled_output(labelled))
    CSV.write(joinpath(work_dir, "control_summary.csv"), summary)

    return (diagnostic = diag, calibrated_index = idx, predictions = preds,
            evaluation = ev_all, evaluation_all = ev_all, evaluation_heldout = ev_heldout,
            labelled = labelled, summary = summary, redundancy = redundancy)
end

end # module
