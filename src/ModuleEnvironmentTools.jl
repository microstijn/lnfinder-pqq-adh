module EnvironmentTools

export fetch_omnicrobe_env!, agg_habitat, calculate_enrichment_stats, compare_archaea_bacteria
export aggregate_habitats, category_from_path, habitat_type, habitat_chemistry, coarse_pillar
export is_isolation_source, is_generic_term

using HTTP, JSON3, CSV, DataFrames, ProgressMeter, HypothesisTests

# We cache OBT paths to avoid redundant API calls for the same habitat
const OBT_CACHE = Dict{String, String}()

#=
    OBT branch ids → habitat TYPE.

    Derived empirically from the OBT ontology paths of the habitat terms the input data
    returns, not assumed. Four branch ids are deliberately NOT in the table, because keying
    it at their level mis-assigns wholesale:

        OBT:000013 is not "Aquatic"     — it is the NATURAL-ENVIRONMENT branch and
                                          contains soil, rock, forest and air as
                                          well as lake and sea, so `soil` and even
                                          `terrestrial habitat` would score Aquatic.
        OBT:000006 is not "Terrestrial" — it is the MAN-MADE branch (saltern, mine,
                                          sewage, bioreactor, hospital, landfill).
        OBT:001660 / OBT:002027 / OBT:000987 / OBT:000763 match NOTHING in the returned
                                          paths at all, so a table keyed on them leaves
                                          Host-associated and Anthropogenic unfired.

    The real pillars sit one level down, e.g. soil = /OBT:000001/OBT:000013/OBT:000158/…
    and lake = /OBT:000001/OBT:000013/OBT:000047/…

    Checked in order; first match wins, so the specific anatomy branches are tested before
    the general "living organism" one.
=#
    # Official ontology names, fetched from /api/get/obt/<id> — quoted here so a future
    # reader can check the label against the source without a network call.
const OBT_TYPE = [
    # Aquatic must include the SUBSTRATE nodes below, not just "aquatic environment".
    # OBT:000079 "environmental matter" is not labelled wholesale: only `rock` is geological,
    # while the branch also holds aquatic sediments and the bare term `water`, which
    # "Geological" would classify as rock. Same failure shape as labelling the
    # natural-environment parent Aquatic, so the branch is split by child instead of
    # labelled at the parent.
    ("Aquatic",        ["OBT:000047",                   # "aquatic environment"
                        "OBT:000468",                   # "water"
                        "OBT:001430",                   # "aquatic sediment"
                        "OBT:000050"]),                 # "artificial water environment" —
                                                        #   salterns; JUDGEMENT CALL, see below
    ("Terrestrial",    ["OBT:000158",                   # "terrestrial habitat"
                        "OBT:000428"]),                 # "soil part"
    ("Geological",     ["OBT:000350"]),                 # "mineral matter" — rock proper, the
                                                        #   only genuinely geological child
    ("Atmospheric",    ["OBT:000052"]),                 # "atmosphere part" — air, aerosol
    ("Host-plant",     ["OBT:000140",                   # "plant habitat" (root, rhizosphere)
                        "OBT:000393"]),                 # "plant" — the whole organism
    ("Host-animal",    ["OBT:000039",                   # "animal part" (gut, rumen)
                        "OBT:000193"]),                 # "animal" — the whole organism
    ("Host-other",     ["OBT:000010",                   # "living organism"
                        "OBT:000014"]),                 # "part of living organism" (PARENT of
                                                        #   000039/000140 — hence tested last)
    ("Agricultural",   ["OBT:000004"]),                 # "animal husbandry and agricultural habitat"
    ("Lab-medium",     ["OBT:000007"]),                 # "experimental medium" — nutrient
                                                        #   broth, yeast extract, cell culture
    ("Anthropogenic",  ["OBT:000006",                   # "artificial environment"
                        "OBT:000011"]),                 # "medical environment"
    ("Food/substance", ["OBT:000008"]),                 # "food" (also carries `salt` itself)
]

#=
    ISOLATION SOURCES, not habitats. "Food/substance" and "Lab-medium" record where a strain
    was OBTAINED — fermented seafood, nutrient broth, a culture flask — not where the
    organism lives. Counting them as habitats inflates the anthropogenic pillar with what is
    really sampling and cultivation history, so ecological analyses should exclude them
    (and say how many they dropped). They stay as first-class types rather than being folded
    away, so the exclusion is visible and reversible rather than silent.
=#
const ISOLATION_SOURCE = Set(["Food/substance", "Lab-medium"])

"""
    is_isolation_source(habitat_type) -> Bool

True for habitat types that record where a strain was obtained rather than where the
organism lives. Filter these out of ecological enrichment; report the count dropped.
"""
is_isolation_source(t) = !ismissing(t) && String(t) in ISOLATION_SOURCE
#=
    On salterns being Aquatic rather than Anthropogenic: this looked like a judgement call
    but the ontology settles it. OBT:003164 "saltern" carries THREE lineages —
        /OBT:000001/OBT:000006/OBT:000050/OBT:003164                       (artificial)
        /OBT:000001/OBT:000009/OBT:000090/OBT:000302/OBT:000667/OBT:003164 (extreme/haline)
        /OBT:000001/OBT:000013/OBT:000047/OBT:000326/OBT:000705/OBT:003164 (AQUATIC)
    — so OBT already classes a saltern as an aquatic environment, and the multi-lineage
    lookup finds it via OBT:000047 regardless of the OBT:000050 entry above. Nothing is
    lost either way: on the two-axis scheme it reads Aquatic × Hypersaline, and the
    constructed-ness is recoverable from the artificial lineage if ever wanted.
=#

#=
    OBT branch ids → habitat CHEMISTRY (the second, orthogonal axis).

    OBT:000009 is an environmental-CONDITION branch, which is what makes two axes possible
    without inventing a scheme: salinity, alkalinity and oxygen state are already modelled
    there, separately from habitat type. A solar saltern is Anthropogenic × Hypersaline; a
    hypersaline lake is Aquatic × Hypersaline; saline soil is Terrestrial × Hypersaline.
    No single-axis pillar list can express that, because salinity cuts across every
    habitat type.
=#
    # These ids are the SPECIFIC condition nodes, not their parents. Three obvious-looking
    # parents were checked and rejected, because each would have mislabelled wholesale:
    #   OBT:002018 "environmental water with chemical property" — any chemistry, not saline
    #   OBT:000141 "polluted environment"   — would flag every polluted site as hypersaline
    #   OBT:000090 "extreme environment"    — thermal/acidic/etc., not specifically alkaline
    #   OBT:000078 "environment wrt oxygen level" — the parent of BOTH aerobic and anaerobic
    # NB the parent OBT:000009 is "habitat wrt chemico-PHYSICAL property", so this axis
    # legitimately carries physical conditions (temperature, pressure, osmolarity) as well
    # as chemical ones. The field is called `chemistry` for brevity; read it as "condition".
    # Extensible: acidity, pressure (OBT:000307) and osmolarity (OBT:000306) have nodes too
    # and can be added the same way when a dataset needs them.
const OBT_CHEMISTRY = [
    ("Hypersaline", ["OBT:002592",                   # "saline water"
                     "OBT:003115",                   # "salt contaminated soil"
                     "OBT:001751",                   # "saline sediment"
                     "OBT:003164"]),                 # "saltern" (under artificial water env.)
    ("Alkaline",    ["OBT:000187"]),                 # "alkaline environment"
    ("Anoxic",      ["OBT:000191"]),                 # "anaerobic environment"
    ("Thermal",     ["OBT:000309",                   # "high temperature environment"
                     "OBT:000625"]),                 # "extreme high temperature environment"
    ("Cold",        ["OBT:000238"]),                 # "cold temperature environment"
]

#=
    Text fallback for chemistry. Necessary because OmniCrobe gives a term ONE path: many
    saline habitats are filed under their type branch instead of the condition branch
    ("lake, Salt Lake" sits under Aquatic; "solar saltern" under man-made), so the
    condition would be lost. Ontology first, text only where the ontology is silent.
=#
const CHEM_KEYWORDS = [
    ("Hypersaline", ["hypersalin", "salin", "saltern", "brine", "halite", "salt lake",
                     "salt mine", "salt soil", "salt marsh", "salt pan", "solar salt",
                     "salt-fermented", "salt contaminated"]),
    ("Alkaline",    ["soda lake", "alkaline", "alkalin"]),
    ("Thermal",     ["hot spring", "hydrothermal", "geothermal", "thermal spring"]),
    ("Anoxic",      ["anaerobic", "anoxic"]),
]

#=
    Word-boundary chemistry patterns. Bare `salt` is sometimes the entire habitat statement
    OmniCrobe returns for a halophile, including *salifodinae* ("of the salt mine") strains.
    That is a chemistry statement, not a habitat, and OBT files it under "food"
    (OBT:000008), so without this a salt-mine archaeon reads as an anthropogenic food
    habitat.

    Anchored because a plain substring would fire on `basalt`, which is a rock.
=#
const CHEM_REGEX = [
    ("Hypersaline", [r"\bsalts?\b", r"\bhalite\b", r"\bbrine\b"]),
]

# Collapse the fine types onto the original five pillars, so existing analyses stay
# comparable with the five-pillar scheme.
const COARSE = Dict(
    "Aquatic" => "Aquatic",
    "Terrestrial" => "Terrestrial", "Geological" => "Terrestrial",
    # Air belongs to none of the original four; mapping it to Other is deliberate, not a
    # gap. Give it its own pillar if aerosol/dust habitats ever matter to the analysis.
    "Atmospheric" => "Other",
    "Host-plant" => "Host-associated", "Host-animal" => "Host-associated",
    "Host-other" => "Host-associated",
    "Agricultural" => "Anthropogenic", "Anthropogenic" => "Anthropogenic",
    "Food/substance" => "Anthropogenic", "Lab-medium" => "Anthropogenic",
    "Other" => "Other",
)

"""
    habitat_type(path) -> String

Fine-grained habitat type from an OBT ontology path; "Other" if it matches no known
branch. Pure — the ontology climb happens elsewhere.
"""
function habitat_type(path::AbstractString)
    for (cat, ids) in OBT_TYPE
        any(r -> occursin(r, path), ids) && return cat
    end
    return "Other"
end

"""
    habitat_chemistry(path, name="") -> String

Chemistry/condition flags for a habitat, from the OBT condition branch where present and
from the term text otherwise. Returns "" when none apply, or several joined by `+`
(a soda lake is `Alkaline+Hypersaline`).
"""
function habitat_chemistry(path::AbstractString, name::AbstractString = "")
    flags = String[]
    for (chem, ids) in OBT_CHEMISTRY
        any(r -> occursin(r, path), ids) && push!(flags, chem)
    end
    n = lowercase(name)
    if !isempty(n)
        for (chem, kws) in CHEM_KEYWORDS
            (chem in flags) && continue
            any(k -> occursin(k, n), kws) && push!(flags, chem)
        end
        for (chem, res) in CHEM_REGEX
            (chem in flags) && continue
            any(r -> occursin(r, n), res) && push!(flags, chem)
        end
    end
    return join(sort(unique(flags)), "+")
end

"""
    coarse_pillar(fine) -> String

Collapse a fine habitat type onto the original five pillars (Aquatic / Terrestrial /
Host-associated / Anthropogenic / Other).
"""
coarse_pillar(fine::AbstractString) = get(COARSE, fine, "Other")

"""
    category_from_path(path) -> String

Backward-compatible five-pillar call: `coarse_pillar(habitat_type(path))`. It resolves
through the branch table above, so `soil` reaches Terrestrial and an engineered site reaches
Anthropogenic; a map keyed at OBT:000013 / OBT:000006 instead would score `soil` Aquatic and
the engineered site Terrestrial — two different mis-assignments, not an exchange of the two
labels. The branch-table header above says why.
"""
category_from_path(path::AbstractString) = coarse_pillar(habitat_type(path))

#=
    GENERIC TERMS: ontology nodes so high in the tree that they carry no ecological
    information. OBT:000010 "living organism" and OBT:000061 "cell" are the two the habitat
    data returns as terms in their OWN right — a relation typed to either says "this
    organism was found in or on some living thing", which is true of almost every
    host-associated record and specific to nothing. Left in, they hand the Host-other
    bucket to taxa whose only informative habitat statements are elsewhere.

    Matched on the TERM'S OWN id, never on its lineage. OBT:000010 is the ancestor of every
    plant and animal host term, so a substring test against the path would delete the whole
    host-associated pillar instead of the two uninformative nodes. The term's own id is the
    leaf of each lineage, which is why the check reads the last segment.
=#
const GENERIC_TERMS = Set(["OBT:000010",      # "living organism"
                           "OBT:000061"])     # "cell"

"""
    is_generic_term(obtid, path="") -> Bool

True for ontology nodes that name life in general rather than a habitat (see
`GENERIC_TERMS`). Both the id and the leaf of every lineage in `path` are checked, so the
predicate works whether the caller holds the OBT id or only the term's resolved path. The
lineage INTERIOR is deliberately not checked: the generic ids sit above the host-associated
terms, and matching them there would drop those terms too.
"""
function is_generic_term(obtid::AbstractString, path::AbstractString = "")
    String(obtid) in GENERIC_TERMS && return true
    for lineage in split(path, ' '; keepempty = false)
        i = findlast('/', lineage)
        leaf = i === nothing ? lineage : lineage[nextind(lineage, i):end]
        String(leaf) in GENERIC_TERMS && return true
    end
    return false
end

#=
    THE RULE NAMES, one tuple per knob. The check that refuses an unknown name and the
    message that lists the valid ones are both read off these tuples, so a rule added to the
    code cannot be missing from the message. Without that check, any symbol `resolve` does not
    recognise falls through to the fractional split, and a misspelt rule then produces a
    plausible table under a rule nobody asked for.
=#
const TIE_RULES = (:mentions_fractional, :fractional, :alphabetical, :legacy_dict,
                   :unresolved, :share, :any_type)
const CHEMISTRY_RULES = (:any_term, :share, :winner_fractional)

# The three tie rules that decide the TYPE axis and nothing else (see `aggregate_habitats`).
const TYPE_ONLY_RULES = (:unresolved, :share, :any_type)

# The share a category or a flag has to reach under either `:share` rule. A parameter of the
# rule, not a property of any dataset; it lives here so a caller names it instead of retyping
# the number, and two call sites cannot quietly disagree about it.
const SHARE_THRESHOLD = 0.25

# The label of a taxon a rule declines to call: a tie that survives the mention step under
# `:unresolved`, or evidence spread so thinly under `:share` that no category reaches the
# share. No branch table produces it, so it can never be scored as a category; and it is not
# blank, so a caller that drops blank calls keeps the taxon in its arm, where it scores zero.
# Dropping it instead would let a rule that is meant to change only the assignment shrink
# the denominator as well.
const UNRESOLVED = "Unresolved"

_rule_names(rules) = join((":" * String(r) for r in rules), ", ")

function _check_rule(knob::AbstractString, rule::Symbol, rules)
    rule in rules || error("$knob must be one of $(_rule_names(rules)); got :$rule")
    return rule
end

# A share is a fraction. 25 (a percentage) or 0 would each turn the share rule into a
# different rule without saying so — nothing clears 25, everything clears 0 — so both are
# refused rather than run.
function _check_share(t::Real)
    (isfinite(t) && 0 < t <= 1) || error("share_threshold must lie in (0, 1]; got $t")
    return t
end

"""
    resolve(primary, mentions, rule; share_threshold=SHARE_THRESHOLD)
        -> (label, weights, n_tied, resolved_by)

Pick the winning category from a tally, with any tie made explicit instead of being left to
`Dict` iteration order.

`primary` is the tally the winner is ranked on (distinct terms per category by default);
`mentions` is the number of asserted statements per category, used only to break a tie in
`primary`. Either may equally be a tally read back from a written call table — that is how a
variant rule is recomputed without refetching, and it is why the rules live here rather than
being rewritten by each consumer. Returns the winning `label`, the `weights` each winner
contributes downstream
(`1/k` each for a `k`-way residual tie, so a split taxon is counted as split rather than
rounded onto one pillar), how many categories were still tied when the rule had to step in
(`n_tied`), and which clause settled it (`resolved_by`).

Rules (`TIE_RULES`), each with what ONE taxon's `weights` sum to. That sum is part of the rule,
not a detail of it: the estimators take Real weights and score every category in a table of
its own, so a rule that quietly handed a taxon extra weight would read as extra hits with
nothing downstream to notice.
  • `:mentions_fractional` (the primary rule) — most distinct terms, then most statements,
    then split `1/k`. Sums to 1.
  • `:fractional` — most distinct terms, then split `1/k`; no mention tie-break. Sums to 1.
  • `:alphabetical` — as `:mentions_fractional`, but a residual tie goes to the first label
    in sort order instead of being split. `n_tied` still reports the tie. Sums to 1.
  • `:unresolved` — as `:mentions_fractional`, but a tie that survives the mention step is
    left uncalled: label `UNRESOLVED`, no weight anywhere. Sums to 1 for a decided call and
    to 0 for a tied one — the taxon is still one taxon of its arm, scoring zero everywhere.
  • `:share` — MULTI-LABEL: weight 1 to every category holding at least `share_threshold` of
    `primary`. Sums to a whole number from 0 to `floor(1 / share_threshold)`. It is 0, and the
    label `UNRESOLVED`, when the evidence is spread so thinly that no category reaches the
    share (five categories of one term each, at the default share).
  • `:any_type` — MULTI-LABEL: weight 1 to every category `primary` holds at all. Sums to the
    number of categories present, at least 1. Not the tied set: a category trailing by any
    margin still counts.
  • `:legacy_dict` — REPLICATION ONLY. Computes
    `first(sort(collect(t), by = kv -> -kv[2]))`: a stable sort over a `Dict`, so among equal
    tallies the winner is whichever key the hash table happens to yield first. `n_tied > 1`
    here means the answer was decided by hash order. Kept so the committed numbers can be
    regenerated and diffed, never for new analysis. Sums to 1.

Under every rule each single weight lies in (0, 1]. The multi-label rules give a taxon more
than one category, never more than one unit OF a category, so a category's hit total can never
exceed the number of taxa it is scored over. An empty tally returns `("", empty, 0, :none)`
under every rule. An unknown rule raises, naming the valid ones.

A residual tie is labelled `join(sort(tied), "|")` — a value no single category can take, so
downstream code cannot mistake a tie for a call. A multi-label call is labelled
`join(sort(labelled), "+")`, the separator the chemistry axis already uses for conditions that
hold together, so two categories held together cannot be read as a tie that was split. The
multi-label rules report the tie at the top of `primary` in `n_tied`, and `resolved_by` is
`:none` when they label exactly the one outright winner, their own name otherwise.
"""
function resolve(primary::AbstractDict, mentions::AbstractDict, rule::Symbol;
                 share_threshold::Real = SHARE_THRESHOLD)
    _check_rule("tie_rule", rule, TIE_RULES)
    _check_share(share_threshold)
    isempty(primary) && return ("", Dict{String,Float64}(), 0, :none)

    top  = maximum(values(primary))
    tied = sort!(String[k for (k, v) in primary if v == top])

    if rule === :legacy_dict
        lab = first(sort(collect(primary), by = kv -> -kv[2]))[1]
        return (lab, Dict{String,Float64}(lab => 1.0), length(tied), :none)
    end

    # The multi-label rules never choose among the tied — they label every category meeting
    # their condition — so neither the mention step nor a split has anything to decide. The
    # share is compared as a quotient against the threshold, not as `v >= t * total`: the
    # quotient of two integers is correctly rounded, so a category sitting exactly on a
    # decimal threshold such as 0.3 stays on it, where the product can round past it.
    if rule === :share || rule === :any_type
        total = sum(values(primary))
        labelled = sort!(String[k for (k, v) in primary
                                if v > 0 && (rule === :any_type || v / total >= share_threshold)])
        lab = isempty(labelled) ? UNRESOLVED : join(labelled, "+")
        return (lab, Dict{String,Float64}(c => 1.0 for c in labelled), length(tied),
                length(labelled) == 1 ? :none : rule)
    end

    by = :none
    if length(tied) > 1 && rule in (:mentions_fractional, :alphabetical, :unresolved)
        m = maximum(get(mentions, k, 0) for k in tied)
        narrowed = String[k for k in tied if get(mentions, k, 0) == m]
        length(narrowed) < length(tied) && (by = :mentions)
        tied = narrowed
    end

    n_tied = length(tied)
    winners = tied
    if n_tied > 1
        if rule === :alphabetical
            winners = [tied[1]]                      # `tied` is sorted
            by = :alphabetical
        elseif rule === :unresolved
            # Neither split nor awarded: no weight, and a label no category can take.
            return (UNRESOLVED, Dict{String,Float64}(), n_tied, :unresolved)
        else
            by = :fractional
        end
    end

    k = length(winners)
    return (k == 1 ? winners[1] : join(winners, "|"),
            Dict{String,Float64}(c => 1 / k for c in winners), n_tied, by)
end

"""
    resolve_chemistry(combos, rule; share_threshold=SHARE_THRESHOLD) -> Dict{String,Float64}

Per-flag weights on the chemistry axis: which conditions a taxon is scored for, and with how
much weight.

`combos` is the tally of chemistry-bearing terms per flag COMBINATION, keyed as
`habitat_chemistry` writes them (`"Alkaline+Hypersaline" => 1`). Pass counts, not rounded
shares: a flag sitting exactly on the threshold lands on either side of it depending on how
its shares were rounded and summed.

The axis is non-exclusive — the conditions are not alternatives to each other — so every flag
is scored in a table of its own and the weights are NOT a distribution over flags. Rules
(`CHEMISTRY_RULES`), each with what one taxon's weights sum to:
  • `:any_term` (the primary rule) — weight 1 to every flag that any term carries, i.e.
    exactly the flags of `chemistry_flags`. Sums to the number of distinct flags, at least 1.
  • `:share` — weight 1 to every flag named by at least `share_threshold` of the
    chemistry-bearing terms; a combination counts toward EACH flag it names. Sums to a whole
    number from 0 to the number of flags present; 0 when no flag reaches the share.
  • `:winner_fractional` — the exclusive winner: the combinations tied at the top of the tally,
    each flag weighted by the fraction of those `k` combinations that name it (tied between
    `Thermal` and `Alkaline+Thermal`, Thermal carries 1 and Alkaline 1/2). Per combination the
    weights sum to 1; per flag they sum to the mean number of flags in the tied combinations,
    at least 1. There is no mention step, unlike the back-compatible `chemistry` winner: the
    tie is settled on the tally alone, which is all `chemistry_profile` records, so the rule
    can be recomputed from the written table.

Every single weight lies in (0, 1], for the reason `resolve` gives. A tally with no
chemistry-bearing term returns an empty map under every rule; an unknown rule raises, naming
the valid ones.
"""
function resolve_chemistry(combos::AbstractDict, rule::Symbol;
                           share_threshold::Real = SHARE_THRESHOLD)
    _check_rule("chemistry_rule", rule, CHEMISTRY_RULES)
    _check_share(share_threshold)
    held = [(String(c), v) for (c, v) in combos if !isempty(c) && v > 0]
    out = Dict{String,Float64}()
    isempty(held) && return out
    flags_of(c) = unique!(String.(split(c, '+'; keepempty = false)))

    if rule === :any_term
        for (c, _) in held, f in flags_of(c)
            out[f] = 1.0
        end
    elseif rule === :share
        total = sum(last, held)
        mass = Dict{String,Float64}()
        for (c, v) in held, f in flags_of(c)
            mass[f] = get(mass, f, 0.0) + v
        end
        for (f, m) in mass
            m / total >= share_threshold && (out[f] = 1.0)
        end
    else
        top = maximum(last, held)
        tied = [c for (c, v) in held if v == top]
        # Counted first and divided once. Adding 1/k per combination drifts: seven additions
        # of 1/7 come to 0.9999999999999998, so a flag every tied combination names would stop
        # being a whole hit, and the 0/1 fast path downstream would never see it as one.
        naming = Dict{String,Int}()
        for c in tied, f in flags_of(c)
            naming[f] = get(naming, f, 0) + 1
        end
        for (f, n) in naming
            out[f] = n / length(tied)
        end
    end
    return out
end

"""
    aggregate_habitats(relations, path_of; weighting=:distinct, max_distinct=50,
                       isolation=:term, drop_generic=true,
                       tie_rule=:mentions_fractional, chem_text=:forms_union,
                       chemistry_rule=:any_term,
                       share_threshold=SHARE_THRESHOLD) -> NamedTuple

Collapse ALL of a taxon's OmniCrobe habitat relations into one call, deterministically.

`relations` is one entry per returned relation (`obtid`, `name`, and — from
`parse_relations` — `forms`); `path_of(obtid) -> String` resolves the OBT path (injected so
this is testable offline, and so the caller controls caching).

Why this exists: taking `json[1]` — an arbitrary single relation out of up to thousands —
would assign *Methylorubrum extorquens*, a model phyllosphere methylotroph, the habitat
"blood"; everything downstream inherits that draw.

THE RULE, in the order it is applied:

1. **Rank** the distinct terms by mentions, ties broken by OBT id. The ranking decides what a
   `max_distinct` cap drops, so leaving equal-mention terms in `Dict` order would make the
   cap's victims depend on hash layout.
2. **Filter, then truncate.** Walk that ranking, resolving paths lazily, and keep terms until
   `max_distinct` usable ones are held; terms never looked at are counted in `n_truncated`.
   With `drop_generic = true` the generic nodes (`GENERIC_TERMS`) go, and with
   `isolation = :term` so do isolation sources — off the TYPE axis, because they record where
   a strain was OBTAINED, not where it lives. `Lab-medium` leaves both axes; `Food/substance`
   stays on the chemistry axis, since OBT files the bare term `salt` there and that is a
   salinity statement. `isolation = :taxon` is what the `legacy` configuration selects:
   nothing is filtered and the caller drops whole taxa afterwards on the winning type.
3. **Tally** the kept terms per category, weighting each term by 1 (`weighting = :distinct`,
   the default — robust to literature attention, so a taxon with very many relations about one
   well-studied context does not win on volume alone) or by its mentions
   (`weighting = :mentions`).
4. **Resolve** the winner with `tie_rule` (see `resolve`) and the chemistry weights with
   `chemistry_rule` (see `resolve_chemistry`). The four rules that pick ONE category —
   `:mentions_fractional`, `:fractional`, `:alphabetical`, `:legacy_dict` — settle the fine
   type, the coarse pillar and the chemistry winner alike, which is the behaviour the
   replication configuration reproduces. The three `TYPE_ONLY_RULES` added for the one-axis-at-
   a-time sensitivity runs settle the TYPE axis only: under them the chemistry winner is taken
   on the primary rule, so that moving the type rule cannot move a chemistry field and a
   variant stays one move off the specification it varies.

Returns the winning `category` and fine `habitat_type` (a residual tie reads `"A|B"`, a
multi-label call `"A+B"`, an uncalled taxon `UNRESOLVED`), a representative `name` (the
most-mentioned term *within* the winning type, `missing` when no category was called), the
`confidence` values, the non-exclusive `chemistry_flags`, and the raw counts every downstream
sensitivity rule needs — `type_counts`, `type_mentions`, `type_weights`, `chemistry_weights`,
`chemistry_flag_counts`, `n_isolation_terms`, `n_generic_terms`, `isolation_only`,
`generic_only` — so a variant rule can be recomputed without refetching anything.

Chemistry is NON-EXCLUSIVE: `chemistry_flags` is the union of every flag over the kept terms
(a taxon in a hot brine reads `Hypersaline+Thermal`), while the back-compatible `chemistry`
field still reports the single winning flag combination. `chemistry_weights` is what
`chemistry_rule` decides — how much weight each flag carries — and it is the only field that
rule touches: the flags, the winner, the counts and the profile are the raw evidence and stay
as they are under every chemistry rule. `chem_text` chooses which text the keyword fallback
reads: `:forms_union` (every surface form, the default), `:object_name` (the ontology term
alone) or `:first_form` (the first surface form — the name the committed tables carry, so
replication only); it decides which terms carry a flag, which is a different question from how
the flags a term carries are weighted, so the two knobs stay apart.
"""
function aggregate_habitats(relations, path_of; weighting::Symbol = :distinct,
                            max_distinct::Int = 50, isolation::Symbol = :term,
                            drop_generic::Bool = true,
                            tie_rule::Symbol = :mentions_fractional,
                            chem_text::Symbol = :forms_union,
                            chemistry_rule::Symbol = :any_term,
                            share_threshold::Real = SHARE_THRESHOLD)
    weighting in (:distinct, :mentions) || error("weighting must be :distinct or :mentions")
    isolation in (:term, :taxon) || error("isolation must be :term or :taxon")
    _check_rule("tie_rule", tie_rule, TIE_RULES)
    _check_rule("chemistry_rule", chemistry_rule, CHEMISTRY_RULES)
    # Checked even when neither rule is a share rule: a threshold given in percent is a
    # mistake worth hearing about at the call that made it, not two configurations later.
    _check_share(share_threshold)
    chem_text in (:forms_union, :object_name, :first_form) ||
        error("chem_text must be :forms_union, :object_name or :first_form")
    # `nothing` means the FETCH FAILED (see `fetch_relations`); it does not mean "this taxon
    # has no habitat". Silently aggregating it as an empty vector is exactly the confusion
    # the fetch signature exists to prevent, so refuse it here rather than return a
    # blank call that is indistinguishable from a genuinely unannotated taxon.
    relations === nothing &&
        throw(ArgumentError("aggregate_habitats received `nothing`: the relation fetch " *
                            "failed. Handle the failure (retry it, or record the taxon as " *
                            "unfetched) instead of aggregating it as 'no habitat'."))
    empty_result = (name = missing, category = missing, confidence = 0.0,
                    habitat_type = missing, type_confidence = 0.0,
                    chemistry = "", chemistry_share = 0.0, chemistry_coverage = 0.0,
                    n_relations = 0, n_distinct = 0, n_truncated = 0,
                    profile = "", type_profile = "", chemistry_profile = "",
                    type_counts = "", type_mentions = "",
                    type_weights = Dict{String,Float64}(), n_tied = 0, resolved_by = :none,
                    chemistry_flags = "", chemistry_weights = Dict{String,Float64}(),
                    chemistry_flag_counts = "", n_chem_terms = 0,
                    n_isolation_terms = 0, n_generic_terms = 0,
                    isolation_only = false, generic_only = false)
    isempty(relations) && return empty_result

    mentions = Dict{String,Int}(); label = Dict{String,String}()
    forms = Dict{String,Vector{String}}(); first_form = Dict{String,String}()
    for r in relations
        oid = String(r.obtid)
        mentions[oid] = get(mentions, oid, 0) + 1
        get!(label, oid, String(r.name))
        fs = hasproperty(r, :forms) ? String.(r.forms) : String[]
        get!(first_form, oid, isempty(fs) ? String(r.name) : fs[1])
        if !isempty(fs)
            v = get!(forms, oid, String[])
            for f in fs
                f in v || push!(v, f)
            end
        end
    end
    n_rel = length(relations); n_dist = length(mentions)

    # Which text the chemistry keyword fallback reads. The forms are joined with a separator
    # a multi-word keyword cannot cross, or "…salt" next to "lake…" would manufacture a
    # "salt lake" that neither form states.
    function chem_source(oid)
        chem_text === :object_name && return label[oid]
        chem_text === :first_form && return first_form[oid]
        fs = get(forms, oid, String[])
        return isempty(fs) ? label[oid] : join(sort(fs), " | ")
    end

    # Rank by mentions, ties by id: a cap must drop the long tail reproducibly, and the second
    # key is what makes "which of the equal-mention terms survives" independent of hash order.
    # Then FILTER BEFORE TRUNCATING — walk the ranking, resolving paths lazily, until
    # `max_distinct` usable terms are held, so the cap spends its budget on terms that can
    # contribute rather than on generic nodes and culture media.
    #
    # `:legacy_dict` selects the `Dict`-order ranking along with the matching tie-break,
    # because the two are one effect seen twice: sorting on mentions alone is stable, so among
    # equal-mention terms the surviving 50 are whichever ones the hash table happens to
    # yield first. Replicating the committed numbers means replicating that too. Nothing but
    # replication should use this branch (see `resolve`).
    legacy = tie_rule === :legacy_dict
    ranked = legacy ? sort(collect(keys(mentions)), by = o -> -mentions[o]) :
                      sort(collect(keys(mentions)), by = o -> (-mentions[o], o))
    paths = Dict{String,String}(); ftype = Dict{String,String}()
    kept = String[]          # usable on at least one axis
    type_terms = String[]    # usable on the TYPE axis
    n_examined = 0; n_generic = 0; n_isolation = 0
    for oid in ranked
        length(kept) >= max_distinct && break
        n_examined += 1
        p = path_of(oid); paths[oid] = p
        ft = habitat_type(p); ftype[oid] = ft
        generic = is_generic_term(oid, p); generic && (n_generic += 1)
        isolat = is_isolation_source(ft); isolat && (n_isolation += 1)
        drop_generic && generic && continue
        if isolation === :term
            ft == "Lab-medium" && continue        # a culture flask is neither a habitat nor
            push!(kept, oid)                      #   a condition: it leaves both axes
            isolat || push!(type_terms, oid)      # `salt` still speaks for the chemistry axis
        else
            push!(kept, oid); push!(type_terms, oid)
        end
    end
    truncated = n_dist - n_examined
    chem_of = Dict(oid => habitat_chemistry(paths[oid], chem_source(oid)) for oid in kept)

    w(oid) = weighting === :mentions ? mentions[oid] : 1

    # Three tallies per axis: the one the winner is ranked on, plus the raw distinct-term and
    # mention counts, so a downstream variant rule can be recomputed without refetching.
    function tally(terms, f)
        prim = Dict{String,Int}(); dist = Dict{String,Int}(); ment = Dict{String,Int}()
        for oid in terms
            v = f(oid)
            isempty(v) && continue
            prim[v] = get(prim, v, 0) + w(oid)
            dist[v] = get(dist, v, 0) + 1
            ment[v] = get(ment, v, 0) + mentions[oid]
        end
        return prim, dist, ment
    end
    # Sorted by weight and then by name: two runs over the same evidence must write the same
    # string, which a `-value`-only sort does not guarantee for equal values.
    ord(t) = sort(collect(t), by = kv -> (-kv[2], kv[1]))
    counts_str(t) = join([string(k, ":", v) for (k, v) in ord(t)], ";")
    # The three `legacy` profile strings keep the hash-order sort under the replication
    # rule, so the committed output can be diffed field by field and not only on the winner.
    # The count fields above have no `legacy` counterpart and stay deterministic.
    prof(t) = (s = sum(values(t)); s == 0 ? "" :
               join([string(k, ":", round(v / s, digits = 3))
                     for (k, v) in (legacy ? sort(collect(t), by = kv -> -kv[2]) : ord(t))], ";"))
    # Share of ONE winner, not of the tied set: for a k-way tie this stays the per-category
    # share, so `confidence` keeps meaning "how much of the evidence backs this call".
    share1(t, weights) = (s = sum(values(t));
                          (isempty(weights) || s == 0) ? 0.0 :
                          t[first(sort(collect(keys(weights))))] / s)

    fprim, fdist, fment = tally(type_terms, oid -> ftype[oid])                  # fine type
    cprim, _,     cment = tally(type_terms, oid -> coarse_pillar(ftype[oid]))   # 5-pillar
    xprim, _,     xment = tally(kept,       oid -> chem_of[oid])                # chemistry

    # A type-axis rule settles the type axis. The chemistry WINNER under one of those three
    # falls back to the primary rule, because they have no single-winner reading: a multi-label
    # rule would join flag combinations into a string that no longer parses as combinations,
    # and `:unresolved` would blank a chemistry call over a tie in the type tally's rule.
    # Keeping the fallback also keeps a type variant from moving a chemistry row.
    chem_rule = tie_rule in TYPE_ONLY_RULES ? :mentions_fractional : tie_rule
    fine,   fweights, n_tied, resolved_by = resolve(fprim, fment, tie_rule; share_threshold)
    coarse, cweights, _, _                = resolve(cprim, cment, tie_rule; share_threshold)
    chem,   xweights, _, _                = resolve(xprim, xment, chem_rule; share_threshold)
    chem_weights = resolve_chemistry(xprim, chemistry_rule; share_threshold)

    # NON-EXCLUSIVE chemistry: a term that is both alkaline and hypersaline votes for both
    # flags, which is what a condition axis means - the conditions are not alternatives to
    # each other. The winning COMBINATION stays in `chemistry` for back-compatibility.
    flag_counts = Dict{String,Int}(); n_chem_terms = 0
    for oid in kept
        isempty(chem_of[oid]) && continue
        n_chem_terms += 1
        for f in split(chem_of[oid], '+')
            flag_counts[String(f)] = get(flag_counts, String(f), 0) + 1
        end
    end

    # Representative name: most-mentioned term that actually sits in the winning fine type,
    # so the name never comes from a different pillar than the one reported.
    in_best = [o for o in type_terms if haskey(fweights, ftype[o])]
    rep_order = legacy ? (o -> -mentions[o]) : (o -> (-mentions[o], o))
    rep = isempty(in_best) ? missing : label[first(sort(in_best, by = rep_order))]
    # How much of the evidence carries ANY chemistry flag — distinguishes "not saline" from
    # "we only had two terms and neither said anything about chemistry".
    chem_any = sum(values(xprim)) / max(sum((w(o) for o in kept), init = 0), 1)

    return (name = rep,
            category = isempty(coarse) ? missing : coarse,
            confidence = share1(cprim, cweights),
            habitat_type = isempty(fine) ? missing : fine,
            type_confidence = share1(fprim, fweights),
            chemistry = chem,
            chemistry_share = isempty(chem) ? 0.0 : share1(xprim, xweights),
            chemistry_coverage = chem_any,
            n_relations = n_rel, n_distinct = n_dist, n_truncated = truncated,
            profile = prof(cprim), type_profile = prof(fprim),
            chemistry_profile = prof(xprim),
            type_counts = counts_str(fdist), type_mentions = counts_str(fment),
            type_weights = fweights, n_tied = n_tied, resolved_by = resolved_by,
            chemistry_flags = join(sort(collect(keys(flag_counts))), "+"),
            chemistry_weights = chem_weights,
            chemistry_flag_counts = counts_str(flag_counts), n_chem_terms = n_chem_terms,
            n_isolation_terms = n_isolation, n_generic_terms = n_generic,
            # "every term this taxon has is one of these" — how a caller tells a call blanked
            # by filtering from a taxon that was never annotated at all.
            isolation_only = n_isolation > 0 && n_isolation == n_examined,
            generic_only = n_generic > 0 && n_generic == n_examined)
end

"""
    get_obt_path(obtid; strict=false) -> String

Every OBT lineage of a term, joined by a space, memoised in `OBT_CACHE`.

`strict = true` is the OFFLINE mode: an id that is not already in the cache raises instead
of reaching for the network. Use it for replayable re-runs. With `strict = false` an
unreachable API and a genuinely path-less term both return `""`, and an empty path types
the term as "Other" without anything in the output saying that a lookup was lost.
"""
function get_obt_path(obtid::AbstractString; strict::Bool = false)
    oid = String(obtid)
    haskey(OBT_CACHE, oid) && return OBT_CACHE[oid]
    strict && error("OBT path for $oid is not cached and strict = true forbids a network " *
                    "lookup. Prime the cache with load_obt_cache!(csv), or pass strict = false.")

    url = "https://omnicrobe.migale.inrae.fr/api/get/obt/$oid"
    try
        r = HTTP.get(url, require_ssl_verification=false, status_exception=false)
        if r.status == 200
            json = JSON3.read(r.body)
            # OBT is a DAG, not a tree: a term can hang under SEVERAL parents and the API
            # returns every lineage. Taking path[1] silently discards the rest — `saline
            # lake` returns its chemico-physical lineage first, so its aquatic lineage would
            # be invisible and it would type as "Other" despite obviously being a lake. Join
            # them all; both classifiers match by substring, so a term is assigned from the
            # union of its lineages (type precedence in OBT_TYPE breaks any tie).
            path_str = isempty(json.path) ? "" : join(String.(json.path), " ")
            OBT_CACHE[oid] = path_str
            return path_str
        end
    catch
    end
    return ""
end

"""
    load_obt_cache!(csv) -> Int

Prime `OBT_CACHE` from a two-column CSV (`obtid`, `path`) and return the cache size. A
missing file is not an error — it returns 0, so a cold first run simply starts empty.

The ontology climb is the slow part of a habitat run and its answers never change, so every
script shares one on-disk cache. The loader and the writer live here rather than being
pasted into each script: two copies that drift are two runs of "the same" pipeline that are
no longer comparable.
"""
function load_obt_cache!(csv::AbstractString)
    isfile(csv) || return 0
    for r in eachrow(CSV.read(csv, DataFrame))
        OBT_CACHE[String(r.obtid)] = ismissing(r.path) ? "" : String(r.path)
    end
    return length(OBT_CACHE)
end

"""
    save_obt_cache(csv) -> Int

Write `OBT_CACHE` to `csv` as (`obtid`, `path`) and return the number of entries written.
Rows are sorted by id, so two runs that reached the same terms produce the same file and a
diff shows what a run actually added.
"""
function save_obt_cache(csv::AbstractString)
    dir = dirname(csv)
    isempty(dir) || mkpath(dir)
    ids = sort!(collect(keys(OBT_CACHE)))
    CSV.write(csv, DataFrame(obtid = ids, path = [OBT_CACHE[i] for i in ids]))
    return length(ids)
end

"""
    Relation

One parsed OmniCrobe habitat relation: the OBT id, the ontology term `name`, every surface
`forms` string the statement used, and the `source` that asserted it — a literature id
(PubMed, EPMC) or a culture-collection / sequence record (GenBank, DSMZ, CIRM). `source` is
carried because a text-mined sentence and a curated strain record are not the same evidence,
and a sensitivity analysis has to be able to separate them afterwards.
"""
const Relation = NamedTuple{(:obtid, :name, :forms, :source),
                            Tuple{String, String, Vector{String}, String}}

# JSON values arrive loosely typed; never let a stray number or null throw mid-parse.
_str(x) = x === nothing ? "" : (x isa AbstractString ? String(x) : string(x))

# `obt_object` is a single object carrying `name`. The legacy alias is accepted in both the
# object and the list-of-objects shape, since it was never observed and its shape is a guess.
function _term_name(v)
    v === nothing && return ""
    v isa AbstractString && return String(v)
    v isa AbstractDict && return _str(get(v, :name, nothing))
    v isa AbstractVector && return isempty(v) ? "" : _term_name(first(v))
    return ""
end

_as_string(b::AbstractString) = String(b)
_as_string(b::AbstractVector{UInt8}) = String(copy(b))

"""
    parse_relations(body) -> Vector{Relation}

Parse an OmniCrobe `/api/search/relations` response body. Pure: no network and no cache, so
it can be tested against a fixture of the real payload and replayed over stored bodies.

`name` comes from `obt_object.name` — the SINGULAR key the API really returns. A parser that
tests for `obt_objects` instead never fires on that branch, because the API does not return
that key, and every name then falls through to `obt_forms[1]`: not the ontology term but the
first surface form of the statement ("spring, mesophilic saline spring" where the term is
"spring").
`obt_objects` is kept as a legacy alias and `obt_forms[1]` stays the last-resort fallback.

Items with no `obtid`, or with no usable name at all, are skipped. A body that is not a JSON
array of relations raises, rather than quietly parsing as zero relations.
"""
function parse_relations(body::AbstractString)
    json = JSON3.read(body)
    json isa AbstractVector ||
        throw(ArgumentError("expected a JSON array of relations, got $(typeof(json))"))

    out = Relation[]
    for item in json
        item isa AbstractDict || continue
        haskey(item, :obtid) || continue
        oid = _str(item.obtid)
        isempty(oid) && continue

        forms = String[]
        fs = get(item, :obt_forms, nothing)
        if fs isa AbstractVector
            for f in fs
                s = _str(f)
                isempty(s) || push!(forms, s)
            end
        end

        nm = _term_name(get(item, :obt_object, nothing))
        isempty(nm) && (nm = _term_name(get(item, :obt_objects, nothing)))   # legacy alias
        isempty(nm) && !isempty(forms) && (nm = forms[1])
        isempty(nm) && continue

        push!(out, (obtid = oid, name = nm, forms = forms,
                    source = _str(get(item, :source, nothing))))
    end
    return out
end

"""
    fetch_relations(taxid; cache_dir=nothing, getter=HTTP.get, retries=0, retry_delay=1.0,
                    backoff=2.0, pause=0, user_agent) -> Vector{Relation} | Nothing

Every habitat relation OmniCrobe holds for `taxid`.

Returns `nothing` when the FETCH FAILED, and an empty vector when the taxon genuinely has no
habitat relation. Collapsing the two into one value — a bare `catch` returning an empty
vector — lets a network blip enter the analysis as the ecological statement "this organism
has no known habitat". Callers must handle `nothing`; where the collapsed reading is really
wanted, it is `something(fetch_relations(t), EnvironmentTools.Relation[])`.

`cache_dir` makes a run idempotent and replayable offline: `<cache_dir>/<taxid>.json` is
read when it exists, and the RAW response body is written there after a successful fetch
(status 200 only). Storing the raw body rather than the parsed rows means a later parser fix
can be replayed over the same bytes without refetching. `getter` is injected for tests.

`retries` (default 0: one request, nothing asked again unless the caller asks for it)
re-asks only on a transient failure — status 429, 500, 502, 503 or 504, or a request that
raised before any status came back — waiting `retry_delay * backoff^(attempt-1)` seconds
before each new attempt; any other non-200 is returned as a failure at once. `pause` sleeps
after a fetch that actually went out and succeeded (never after a cache hit), so a long run
stays polite, and `user_agent` names this client to the service.
"""
function fetch_relations(taxid::Int; cache_dir::Union{AbstractString,Nothing} = nothing,
                         getter = HTTP.get, retries::Int = 0, retry_delay::Real = 1.0,
                         backoff::Real = 2.0, pause::Real = 0,
                         user_agent::AbstractString = "LnFinder/1.0 (research; serial)")
    cached = cache_dir === nothing ? nothing : joinpath(cache_dir, "$taxid.json")
    if cached !== nothing && isfile(cached)
        try
            return parse_relations(read(cached, String))
        catch e
            @warn "Unreadable relation cache entry; refetching" file = cached exception = e
        end
    end

    url = "https://omnicrobe.migale.inrae.fr/api/search/relations?taxid=ncbi%3A$taxid&type=habitat"
    # A TRANSIENT FAILURE AND A REFUSAL NEED OPPOSITE RESPONSES. 429 and the 5xx statuses
    # below are the service under load, and a request that raised is a connection that
    # dropped; both are worth asking again. Any other non-200 is an answer about this
    # request, and repeating it only repeats the mistake. `retries` defaults to 0, so a caller
    # that does not ask for retries, and every injected test getter, makes exactly one
    # request per uncached taxon -- the long network runs, where one blip otherwise discards
    # the whole fetch, opt in explicitly.
    body = nothing
    for attempt in 1:(retries + 1)
        transient = false
        body = try
            r = getter(url; require_ssl_verification = false, status_exception = false,
                       headers = ["User-Agent" => user_agent])
            if r.status == 200
                _as_string(r.body)
            elseif r.status in (429, 500, 502, 503, 504)
                transient = true
                @warn "OmniCrobe relations: retryable status" taxid status = r.status attempt
                nothing
            else
                @warn "OmniCrobe relations: non-200 status" taxid status = r.status
                return nothing
            end
        catch e
            transient = true
            @warn "OmniCrobe relations: fetch failed" taxid exception = e attempt
            nothing
        end
        body === nothing || break
        (!transient || attempt == retries + 1) && return nothing
        sleep(float(retry_delay) * float(backoff)^(attempt - 1))
    end
    body === nothing && return nothing
    pause > 0 && sleep(pause)

    rels = try
        parse_relations(body)
    catch e
        @warn "OmniCrobe relations: unparseable body" taxid exception = e
        return nothing
    end

    if cached !== nothing
        try
            mkpath(dirname(cached))
            tmp = cached * ".tmp"          # write-then-rename: a killed run leaves no
            write(tmp, body)               # half-written file that later parses as empty
            mv(tmp, cached; force = true)
        catch e
            @warn "Could not write the relation cache entry" file = cached exception = e
        end
    end
    return rels
end

# Aggregate over ALL relations rather than trusting the first one. See aggregate_habitats.
# The third element reports whether the FETCH succeeded: a failed fetch also yields a blank
# call, and upstream has to be able to tell the two apart.
function query_habitat_with_hierarchy(taxid::Int; weighting::Symbol = :distinct,
                                      max_distinct::Int = 50,
                                      isolation::Symbol = :term, drop_generic::Bool = true,
                                      tie_rule::Symbol = :mentions_fractional,
                                      chem_text::Symbol = :forms_union,
                                      cache_dir::Union{AbstractString,Nothing} = nothing,
                                      path_of = get_obt_path)
    rels = fetch_relations(taxid; cache_dir = cache_dir)
    rels === nothing && return (taxid, aggregate_habitats(Relation[], path_of), false)
    agg = aggregate_habitats(rels, path_of; weighting = weighting, max_distinct = max_distinct,
                             isolation = isolation, drop_generic = drop_generic,
                             tie_rule = tie_rule, chem_text = chem_text)
    return (taxid, agg, true)
end

# A tie's 1/k weights, flattened for a CSV column. Same "Key:value;…" shape and same
# (-weight, name) ordering as the count and profile strings, so one parser reads them all.
_weights_str(d::AbstractDict) =
    join([string(k, ":", round(v, digits = 3))
          for (k, v) in sort(collect(d), by = kv -> (-kv[2], kv[1]))], ";")

"""
    fetch_omnicrobe_env!(df; taxid_col=:TaxID, weighting=:distinct, max_distinct=50,
                         isolation=:term, drop_generic=true,
                         tie_rule=:mentions_fractional, chem_text=:forms_union,
                         cache_dir=nothing)

Annotate `df` with OmniCrobe habitat, aggregated over every relation a taxon has (see
`aggregate_habitats`, whose keyword arguments are passed straight through). Adds, besides the
original two columns:

  • `habitat_confidence`  — winning pillar's share of the tallied weight (0–1). **Filter on
    this before drawing ecological conclusions**; a 0.35 call is a coin-toss between pillars.
  • `habitat_n_relations` / `habitat_n_distinct` — how much evidence the call rests on.
  • `habitat_profile`     — the full pillar distribution, e.g. "Aquatic:0.5;Terrestrial:0.3".
  • `habitat_type_counts` / `habitat_type_weights` — the raw distinct-term count per fine
    type, and the 1/k weight each winner carries when the tie rule had to split one. Together
    they let a downstream table be recomputed under a different rule without refetching.
  • `habitat_chemistry_flags` / `habitat_chemistry_flag_counts` — the NON-EXCLUSIVE chemistry
    union and how many terms voted for each flag, next to the single-winner
    `habitat_chemistry`.

`habitat_specific` is the most-mentioned term *within* the winning pillar (not the first
relation returned, which ignores the pillar), so `agg_habitat` consumes it directly.

`cache_dir` is handed straight to `fetch_relations`, so a re-run reads the stored response
bodies instead of the API. TaxIDs whose fetch FAILED are counted and reported: their rows
are blank, and a blank row from a failed fetch must not be read as "no habitat".
"""
function fetch_omnicrobe_env!(df::DataFrame; taxid_col=:TaxID, weighting::Symbol = :distinct,
                              max_distinct::Int = 50,
                              isolation::Symbol = :term, drop_generic::Bool = true,
                              tie_rule::Symbol = :mentions_fractional,
                              chem_text::Symbol = :forms_union,
                              cache_dir::Union{AbstractString,Nothing} = nothing)
    unique_tids = unique(collect(skipmissing(df[!, taxid_col])))
    println("🌐 Fetching habitats and OBT hierarchies for $(length(unique_tids)) TaxIDs...")

    results = @showprogress asyncmap(t -> query_habitat_with_hierarchy(t; weighting,
                                                                      max_distinct, isolation,
                                                                      drop_generic, tie_rule,
                                                                      chem_text, cache_dir),
                                     unique_tids)

    tid_map = Dict(tid => agg for (tid, agg, _) in results)
    n_failed = count(r -> !r[3], results)
    blank = aggregate_habitats(Relation[], identity)
    getagg(t) = ismissing(t) ? blank : get(tid_map, Int(t), blank)

    df.habitat_specific    = [getagg(t).name              for t in df[!, taxid_col]]
    df.habitat_category    = [getagg(t).category          for t in df[!, taxid_col]]
    df.habitat_confidence  = [getagg(t).confidence        for t in df[!, taxid_col]]
    df.habitat_type        = [getagg(t).habitat_type      for t in df[!, taxid_col]]
    df.habitat_chemistry   = [getagg(t).chemistry         for t in df[!, taxid_col]]
    df.chemistry_share     = [getagg(t).chemistry_share   for t in df[!, taxid_col]]
    df.habitat_n_relations = [getagg(t).n_relations       for t in df[!, taxid_col]]
    df.habitat_n_distinct  = [getagg(t).n_distinct        for t in df[!, taxid_col]]
    df.habitat_profile     = [getagg(t).profile           for t in df[!, taxid_col]]
    df.habitat_type_profile= [getagg(t).type_profile      for t in df[!, taxid_col]]
    df.habitat_type_counts = [getagg(t).type_counts       for t in df[!, taxid_col]]
    df.habitat_type_weights= [_weights_str(getagg(t).type_weights) for t in df[!, taxid_col]]
    df.habitat_chemistry_flags       = [getagg(t).chemistry_flags for t in df[!, taxid_col]]
    df.habitat_chemistry_flag_counts = [getagg(t).chemistry_flag_counts
                                        for t in df[!, taxid_col]]

    mapped = count(!ismissing, df.habitat_category)
    conf = df.habitat_confidence[.!ismissing.(df.habitat_category)]
    ntrunc = count(t -> !ismissing(t) && getagg(t).n_truncated > 0, df[!, taxid_col])
    println("✅ $mapped/$(nrow(df)) rows mapped to a pillar; median confidence " *
            (isempty(conf) ? "n/a" : string(round(sort(conf)[max(1, end ÷ 2)], digits = 2))))
    ntrunc > 0 && println("⚠ $ntrunc rows had >$(max_distinct) distinct habitat terms; " *
                          "the long tail was dropped (most-mentioned kept).")
    n_failed > 0 && println("⚠ $n_failed/$(length(unique_tids)) TaxIDs could not be fetched. " *
                            "Their rows are blank, which is NOT evidence of 'no habitat'; " *
                            "re-run to fill them before using the ecology columns.")
    return df
end

const _KW_ANTHRO = ["contaminated", "food", "broth", "agar", "industrial", "waste", "beer",
                    "wine", "storage", "medium", "sewage", "polluted", "mine", "arsenic",
                    "oil-", "hch-", "experimental"]
# Host-associated, split by which kind of host — see `split_host` below.
const _KW_PLANT  = ["plant", "root", "leaf", "rhizosphere", "nodule", "seed", "fruit", "tree"]
const _KW_ANIMAL = ["human", "animal", "feces", "gut", "blood", "insect", "patient",
                    "sputum", "lung", "skin", "intestine", "respiratory"]
const _KW_HOST_X = ["host", "biota", "organism", "cell"]          # host unspecified
const _KW_AQUA   = ["water", "marine", "lake", "sea", "spring", "ocean", "ice", "glacier",
                    "lagoon", "freshwater", "brine", "hydrothermal", "aquatic", "hotspring",
                    "trench", "coastal"]
const _KW_TERRA  = ["soil", "paddy", "field", "permafrost", "sediment", "mud", "rhizoplane",
                    "terrestrial"]
# Plant SURFACES. Under the original five-pillar convention "rhizoplane" fell to
# Terrestrial and "phyllosphere" to Other; both are left exactly there by default and are
# only reclassified when the caller explicitly asks for the host split.
const _KW_PLANT_X = ["phyllosphere", "rhizoplane"]

"""
    agg_habitat(name; split_host=false)

Heuristic roll-up of a specific OBT leaf name into ecological pillars.

`split_host=false` (default) reproduces the original five-pillar behaviour exactly.

`split_host=true` splits **Host-associated** into `Plant-associated`, `Animal-associated`
and `Host-other`. The lumped pillar mixes rhizosphere/root/nodule with gut/blood/patient;
the split keeps plant- and animal-associated hosts apart.
"""
function agg_habitat(name::Union{String, Missing}; split_host::Bool = false)
    if ismissing(name) || name == "Unknown"
        return "_Not Mapped_"
    end

    n = lowercase(name)
    hit(kws) = any(occursin.(kws, n))

    # Priority 1: Anthropogenic
    if hit(_KW_ANTHRO)
        return "Anthropogenic"
    # Plant surfaces — only reclassified when the host split is requested (see _KW_PLANT_X).
    elseif split_host && hit(_KW_PLANT_X)
        return "Plant-associated"
    # Priority 2: Host-associated (optionally resolved by host kind; plant before animal)
    elseif hit(_KW_PLANT) || hit(_KW_ANIMAL) || hit(_KW_HOST_X)
        split_host || return "Host-associated"
        hit(_KW_PLANT)  && return "Plant-associated"
        hit(_KW_ANIMAL) && return "Animal-associated"
        return "Host-other"
    # Priority 3: Aquatic
    elseif hit(_KW_AQUA)
        return "Aquatic"
    # Priority 4: Terrestrial
    elseif hit(_KW_TERRA)
        return "Terrestrial"
    else
        return "Other"
    end
end

function calculate_enrichment_stats(df; top_n::Int = 5)
    # Filter for hits that actually have a habitat and a valid position-2 residue. WHICH
    # residues the table covers is not fixed here: it is the `top_n` commonest position-2
    # residues OF THE TABLE PASSED IN, so the selection is a property of the caller's data
    # and not a constant of this repository, and ties at the cut are broken alphabetically
    # so the choice is reproducible. The residues kept are the row labels of the table
    # returned, so it states what it was built on. Valid means a standard one-letter code:
    # the placeholder written for a sequence with no second motif residue is not one.
    aa = Set(string.(collect("ACDEFGHIKLMNPQRSTVWY")))
    with_habitat = filter(row -> row.report_group ∉ ["_Not Mapped_", "Other"] &&
                                 row.pos2_aa ∈ aa, df)
    counts = combine(groupby(with_habitat, :pos2_aa), nrow => :n)
    sort!(counts, [order(:n, rev = true), :pos2_aa])
    top_aas = collect(first(counts.pos2_aa, min(top_n, nrow(counts))))
    sub = filter(row -> row.pos2_aa ∈ top_aas, with_habitat)
    
    # Create contingency table
    ct = combine(groupby(sub, [:report_group, :pos2_aa]), nrow => :Count)
    matrix_df = unstack(ct, :pos2_aa, :report_group, :Count, fill=0)
    
    # Math for Z-scores
    obs = Array(matrix_df[:, 2:end])
    row_sums = sum(obs, dims=2)
    col_sums = sum(obs, dims=1)
    total = sum(obs)

    # Calculate Z-score matrix
    expected = (row_sums * col_sums) ./ total
    z_scores = (obs .- expected) ./ sqrt.(expected)
    
    # Pearson Chi-Squared Test
    pval = pvalue(ChisqTest(obs))
    
    return matrix_df, obs, z_scores, pval
end

"""
    compare_archaea_bacteria(df::DataFrame) -> (table, p)

Among Lanthanide-classified centroids in Bacteria and Archaea, test whether the D-x-D
spacer residue (`pos2_aa`) differs between the two domains: a Pearson chi-squared test on
the kingdom x position-2 count table.

Exploratory only. Where the archaeal counts are small, expected cell counts fall below 5,
and the asymptotic p-value is not reliable there; use an exact or Monte Carlo test before
reporting one.
"""
function compare_archaea_bacteria(df::DataFrame)
    # Lanthanide-classified centroids in Bacteria and Archaea only
    sub = filter(row -> row.classification == :Lanthanide && 
                        row.kingdom in ["Bacteria", "Archaea"], df)
    
    # Contingency table: kingdom x position-2 residue
    ct = combine(groupby(sub, [:kingdom, :pos2_aa]), nrow => :Count)
    matrix_df = unstack(ct, :pos2_aa, :kingdom, :Count, fill=0)
    
    # Pearson chi-squared test on the counts
    obs = Matrix{Int}(matrix_df[:, [:Bacteria, :Archaea]])
    p_val = pvalue(ChisqTest(obs))
    
    return matrix_df, p_val
end

end