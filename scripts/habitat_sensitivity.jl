#=
    habitat_sensitivity.jl — the one-at-a-time sensitivity grid around the registered
    habitat specification, the robustness flag it decides, and the rendered panel.

    WHY THIS STEP EXISTS. Every estimate the enrichment step reports rests on rule choices a
    reader cannot see: which habitat wins when a taxon's terms tie, whether a chemistry flag
    needs a share of the evidence or a single mention, which text a keyword fallback reads,
    and which background counts as "carries no such enzyme". None of those choices is
    visible in an interval, and each of them can move one. A result that holds only under its
    own tie rule is not a result. This step varies each choice ONE AT A TIME off the
    registered specification, reports what survives, and names — as data, per category — the
    variants that break what does not.

    THE GRID IS NOT COMPUTED HERE. Every level reaches the estimators through code that
    already exists and is already tested: the aggregation rules through
    `EnvironmentTools.resolve` / `resolve_chemistry`, the arms and the rows through
    `analyse` and `contrast_rows` in scripts/habitat_enrichment.jl and
    src/ModuleHabitatStats.jl, the interval predicates and the tier through
    src/ModuleOutcomeRules.jl. A second implementation of any of those would be the one that
    drifts, and a grid computed by drifted code measures the drift rather than the rule.

    THE TWO SEAMS THAT ARE NOT ARGUMENTS OF `analyse`, and how they are crossed.

      TYPE RULE.  `analyse` reads the taxon-call table's `type_weights` column verbatim; it
      has no tie-rule argument. A type variant is therefore computed by rewriting that column
      from the same table's own `type_counts` and `type_mentions` tallies through
      `resolve`, and calling `analyse` unchanged. The reconstruction is checked before it is
      used: recomputing the REGISTERED rule from those tallies must reproduce the committed
      column exactly, or the rewrite is not reading the tallies the table was built from.

      CHEMISTRY RULE. `analyse` reads a chemistry LABEL column through `flag_weights`, which
      cannot express a fractional weight — and silently mis-parses a weight string, since
      "Flag:0.5" is a label to it. A chemistry variant is therefore carried by the weight
      column instead, which `weight_map` parses, and the chemistry rows of that run are the
      run's `taxa` rows read under their flag names. The route is checked against the one it
      replaces: under the registered chemistry rule the two must agree cell for cell.

    WEIGHTS ARE NOT ROUNDED ON THE WAY IN. The published call table rounds its weight column
    for reading; this grid compares hit totals, and a three-decimal round trip changes one.
    A rewritten column is therefore written at full precision, and only the CHECK against the
    published column uses the published format. The registered specification itself is run on
    the committed column verbatim, so the primary rows of this grid are the enrichment step's
    own rows and not a recomputation of them.

    WHAT DECIDES A FLAG, AND WHAT ONLY DESCRIBES. The flag reads the primary interval and the
    one-at-a-time variants, on the Ln arm against the background, on the two units whose
    denominator is a taxon. The unit checks — the symmetric genus unit, the per-sequence
    unit, and the subset of taxa whose call rests on several distinct terms — are reported
    for the same grid and marked as NOT feeding the flag, because a unit check answers a
    different question from a rule variant: it changes what is being counted, not which rule
    counted it. The anchors and the disclosed background lines sit outside the variant set
    for the same reason and are marked the same way.

    A UNIT CHECK IS A UNIT, NOT A FIFTEENTH RULE. All three are reported for the registered
    specification AND for every one-at-a-time variant, under the variant's own name, and what
    separates their rows is the unit they carry. The genus and per-sequence units cost nothing
    — the arm builder emits them beside the taxon unit — while the thin-evidence subset is a
    second run of the same specification over fewer units, and its rows say so in the only
    place that cannot be lost: the name of the unit whose n1 and n2 they report.

    THE ANCHORS ARE THE WIRING CHECK. Two specifications outside the grid reproduce earlier
    runs of this analysis: the all-legacy specification, compared row by row against an
    earlier run's rows table, and an earlier sensitivity run's own chemistry rules on the
    legacy arms, compared against transcribed targets. Targets are also taken on the REGISTERED
    specification's own rows, because an anchor that only reproduces the earlier pipeline says
    nothing about the rows that become colours: the anchors run a different configuration, a
    different chemistry column and no background exclusion, so a wiring error confined to the
    registered specification's own inputs would meet every one of them. Neither target is
    written in this file — a target typed into a script is a number about this dataset, and
    the script must carry none — so both arrive by option, and a run given neither says so in
    its own output rather than passing quietly.

    WHETHER THE ANCHORS ARE BUILT is set by ANCHOR_LEVELS_RETIRED. Both anchors reproduce
    rows of a census other than the one a run reads, so while that constant holds they are
    not built: the row-by-row comparison covers zero rows and records itself as not
    applicable, --anchor-reference is not required and not opened, and the targets on the
    registered specification and its levels are the wiring check that remains.

    AND THE SAME FILE SAYS WHAT THE SET IS. Everything in this script that could contradict a
    deleted or substituted level is enumerated from the tuples the level would be deleted
    from, so an edit there moves the wiring and the vocabulary together and nothing is left to
    disagree. The fourteen level NAMES — rule names, not measurements — therefore arrive with
    the targets, and a run whose enumerated set differs from the declared one in either
    direction stops before it computes.

    Run:  julia --project=. scripts/habitat_sensitivity.jl [options]
          LNFINDER_FIGDIR points every output at the figure directory; see --figdir.
=#

const REPO = normpath(joinpath(@__DIR__, ".."))

const USAGE = """
Usage: julia --project=. scripts/habitat_sensitivity.jl [options]

Options:
  --calls <path>             taxon-call table (default: <data>/habitat/habitat_taxon_calls.csv)
  --scan <path>              background carrier scan table
  --census <path>            census table, for the disclosed calcium comparator arm
  --figdir <path>            where every output goes (default: \$LNFINDER_FIGDIR, else figs/)
  --anchor-reference <path>  an earlier run's rows CSV; the all-legacy anchor is compared
                             against it row by row. While the anchor levels are evaluated a
                             run without it is refused, for the same reason as one without
                             --anchor-targets: a check one dropped flag removes is a check
                             nobody misses. While they are retired (ANCHOR_LEVELS_RETIRED)
                             there are no anchor rows, so it is not required: the run records
                             the comparison as not applicable, and a path given anyway is
                             recorded by name and not opened
  --anchor-targets <path>    the transcribed anchor targets, their named exceptions, the NAMES
                             of the registered one-at-a-time levels, the file the targets were
                             transcribed FROM with its digest, and the three counts the file
                             declares about itself, as kind,anchor,unit,category,quantity,
                             value,tol,note,digest with kind in
                             target|exception|variant|expect|source
  --no-anchors               run with no target file and no reference at all, and say so in
                             the output. Neither half of the wiring check is then made and the
                             registered set is compared against nothing; this has to be typed,
                             because a check that is skipped by default is a check nobody runs
  --no-figure                write the tables and skip the rendered panel. Refused if the
                             figure directory already holds a panel of the current name, which
                             this run would leave standing beside its own tables
  -h, --help                 print this text and exit

Exit status: 0 on a completed run; 1 on a command-line error, a missing input, a failed
guard, or an anchor target that was not reproduced.
"""

# --help must cost nothing: the analysis module below takes seconds to load.
if abspath(PROGRAM_FILE) == (@__FILE__) && any(a -> a in ("-h", "--help"), ARGS)
    print(USAGE)
    exit(0)
end

using Pkg; Pkg.activate(REPO)
using DataFrames, CSV, Printf, Dates, SHA, Statistics

# Where every input of this step lives, by NAME; the layout is declared once, in
# src/ModuleDataInputs.jl, and documented once, in docs/DATA.md.
include(joinpath(REPO, "src", "ModuleDataInputs.jl"))

# The enrichment step goes into a module of its own, behind its own program-file guard, so
# including it defines its functions and runs nothing. It is the ONLY arm builder in this
# pipeline, and it brings its libraries with it: a second copy of `analyse`, of the
# estimators or of the exclusion rules would be the copy that drifts.
module Enrichment end
Base.include(Enrichment, joinpath(REPO, "scripts", "habitat_enrichment.jl"))
const HE = Enrichment
const ET = Enrichment.EnvironmentTools
const HS = Enrichment.HabitatStats
const OR = Enrichment.OutcomeRules

sv(x) = ismissing(x) ? "" : String(strip(string(x)))
hascol(df::AbstractDataFrame, name::Symbol) = name in propertynames(df)
fnum(x) = x === missing || x === nothing ? NaN :
          x isa Real ? float(x) : something(tryparse(Float64, sv(x)), NaN)
# WHAT A DIGEST LOOKS LIKE, asked in one place. Every provenance row that carries a sha256 is
# a claim about which bytes an input was, and the empty string satisfies "a digest was
# recorded" while saying nothing at all — which is exactly what a run over a file that was not
# there writes. A digest field is therefore checked for being a digest wherever one is
# required, rather than only for being present.
is_sha256(s::AbstractString) = length(s) == 64 && all(isxdigit, s)

# ---------------------------------------------------------------------------
# The specifications. Nothing below is a quantity measured on any data: these are rule
# names, configuration names and the names of the columns each level reaches through.
# ---------------------------------------------------------------------------

# The configuration column of the taxon-call table that carries the registered levels of the
# text, generic-term and isolation axes. The other five configurations ARE the variants of
# those three axes, which is why those axes need no rule argument here.
const PRIMARY_CONFIG = "primary"

# The type-axis levels, in the order the registered set names them. Each is a rule of
# `EnvironmentTools.resolve`; `nothing` means "the configuration's own committed column",
# which is the registered level.
const TYPE_LEVELS = (:fractional, :alphabetical, :unresolved, :share, :any_type)
# The chemistry-axis levels, each a rule of `EnvironmentTools.resolve_chemistry`.
const CHEMISTRY_LEVELS = (:share, :winner_fractional)
# The chemistry-text and term-filter axes, as the configurations that carry them.
const TEXT_CONFIGS      = ("text_first", "text_object")
const GENERIC_CONFIG    = "generic_none"
const ISOLATION_CONFIG  = "iso_taxon"
# The background levels, as modes of the enrichment step's own exclusion function. The
# registered primary is that function's first mode and is NOT a variant.
const BACKGROUND_LEVELS = ("ipr017512", "none", HE.MOVEIN_MODE)
# Background rules reported beside the grid and deliberately outside it: they are disclosed
# lines, not one-at-a-time variants of the registered specification.
const DISCLOSED_BACKGROUND_LEVELS = ("lnca", "lntype")

# An earlier sensitivity run's own chemistry rules on the legacy arms. Three are rules of
# `resolve_chemistry`; the fourth settles the tie among chemistry COMBINATIONS
# alphabetically, which is `resolve`'s alphabetical rule applied to the combination tally —
# see `chemistry_weights`.
const AUDIT_CHEMISTRY_RULES = (:alphabetical_combination, :winner_fractional, :any_term, :share)
const COMBINATION_ALPHABETICAL = :alphabetical_combination

# THE ANCHOR LEVELS ARE NOT EVALUATED AGAINST THE CENSUS THIS RUN READS.
#
# `anchor:all_legacy` reproduces an earlier run's rows, and the four `audit:*` levels reproduce
# an earlier sensitivity run's chemistry rules on the earlier arms. Both reproduce rows of a
# census other than the one this run reads. Asked of a run reading a different census they are
# a DIFFERENT test, and a pass would be misread as evidence that the two censuses give the same
# rows - so they are not evaluated.
#
# The definitions below are KEPT, not deleted, so that a run against the census they describe
# can still produce them by flipping this one constant. Retiring them empties `anchor_rows`, so
# the `--anchor-reference` comparison covers zero rows; the targets file states that as
# `reference_rows,0` rather than leaving it to be inferred, and the run records the retirement
# in its provenance (`anchor_levels_retired`, and an `anchor_reference_status` that says why)
# so a reader of the output sees a decision rather than an absence.
const ANCHOR_LEVELS_RETIRED = true

# How many distinct usable terms a call must rest on for the thinly-annotated-tail unit
# check. A parameter of that check, named here and passed to the library's own mask so the
# label and the mask cannot drift apart.
const MIN_DISTINCT_TERMS = 3
# The subset is reported as a UNIT of every specification rather than as a specification of
# its own, because that is what it is: the same rule counted over fewer units. Suffixing the
# unit name keeps it in the same row key as its parent, so the grid carries one series per
# specification and a reader cannot line a subset row up against a whole-arm row by accident.
const MIN_DISTINCT_SUFFIX = " [>=$MIN_DISTINCT_TERMS terms]"

# The two chemistry routes are two computations of ONE quantity, not two measurements, so
# what may separate them is the order the floats were summed and nothing else.
const CHEMISTRY_ROUTE_TOL = 1e-9

# The names of the two derived units, built from the unit the arm builder emits so a rename
# there cannot leave this file describing a unit that no longer exists.
chemistry_unit(unit::AbstractString) = unit * " (chemistry)"
is_chemistry_unit(unit::AbstractString) = endswith(unit, " (chemistry)")
# The unit a written row carries may wear the subset suffix, so anything that asks WHICH KIND
# of unit a row is on reads the stem rather than the whole name: "taxa (chemistry) [>=3
# terms]" is a chemistry unit, and a test that only looked at the end of the name would say
# it is not.
unit_stem(unit::AbstractString) = replace(String(unit), MIN_DISTINCT_SUFFIX => "")

const ROLE_PRIMARY   = "primary"      # the registered specification itself
const ROLE_VARIANT   = "variant"      # a one-at-a-time level of the registered set
const ROLE_UNIT      = "unit_check"   # a change of what is counted, not of the rule
const ROLE_ANCHOR    = "anchor"       # an earlier specification, reproduced
const ROLE_AUDIT     = "audit"        # an earlier sensitivity run, reproduced
const ROLE_DISCLOSED = "disclosed"    # reported beside the grid, outside it
# ROLE_AUDIT's value and the `audit:` prefix of the specification names it labels are the
# names the recorded target file spells in its own `anchor` column, so they are an interface:
# renaming either invalidates every target written against them. What the role means is its
# comment above, not the word.

"""
One specification of the whole analysis: which configuration of the call table, which rule
each of the two aggregation axes runs under, which background rule, which subset of the
taxa, and what the run is FOR.

`tie_rule === nothing` and `chem_rule === nothing` each mean "the level the call table
already carries", which is how the registered specification is expressed — it is run on the
committed columns rather than on a recomputation of them.
"""
struct Spec
    name::String
    axis::String
    level::String
    config::String
    tie_rule::Union{Symbol,Nothing}
    chem_col::Symbol
    chem_rule::Union{Symbol,Nothing}
    bg_mode::String
    subset::Symbol                     # :all | :min_distinct
    role::String
    emit::Symbol                       # :both | :chemistry_only
    unit_suffix::String                # what the subset does to the name of every unit
    note::String
end

# A specification's identity in the run cache: the name alone is not unique, because the
# thin-evidence subset of a specification carries the SAME name (it is the same rule) and
# differs only in which units it was counted over.
run_key(spec::Spec) = (spec.name, spec.unit_suffix)

primary_spec() = Spec("primary", "primary", "registered", PRIMARY_CONFIG, nothing,
                      :chemistry_flags, nothing, HE.PRIMARY_MODE, :all, ROLE_PRIMARY, :both,
                      "", "the registered specification, run on the committed call columns")

"""
    sensitivity_specs() -> Vector{Spec}

The registered sensitivity set, one level at a time off the primary: the type rule, the
chemistry rule, the chemistry text, the generic-term filter, the isolation filter and the
background exclusion. The order is the order the set is registered in, and the count is a
property of these tuples rather than a number written down beside them.
"""
function sensitivity_specs()
    p = primary_spec()
    out = Spec[]
    for r in TYPE_LEVELS
        push!(out, Spec("T:$r", "T", String(r), p.config, r, p.chem_col, nothing, p.bg_mode,
                        :all, ROLE_VARIANT, :both, "",
                        "type rule rewritten from the table's tallies"))
    end
    for r in CHEMISTRY_LEVELS
        push!(out, Spec("C:$r", "C", String(r), p.config, nothing, p.chem_col, r, p.bg_mode,
                        :all, ROLE_VARIANT, :both, "",
                        "chemistry rule carried by the weight column"))
    end
    for c in TEXT_CONFIGS
        push!(out, Spec("X:$c", "X", c, c, nothing, p.chem_col, nothing, p.bg_mode, :all,
                        ROLE_VARIANT, :both, "",
                        "chemistry text, as a configuration of the call table"))
    end
    push!(out, Spec("G:$GENERIC_CONFIG", "G", GENERIC_CONFIG, GENERIC_CONFIG, nothing,
                    p.chem_col, nothing, p.bg_mode, :all, ROLE_VARIANT, :both, "",
                    "generic terms kept, as a configuration of the call table"))
    push!(out, Spec("I:$ISOLATION_CONFIG", "I", ISOLATION_CONFIG, ISOLATION_CONFIG, nothing,
                    p.chem_col, nothing, p.bg_mode, :all, ROLE_VARIANT, :both, "",
                    "isolation filtered at taxon level, as a configuration of the call table"))
    for m in BACKGROUND_LEVELS
        push!(out, Spec("B:$m", "B", m, p.config, nothing, p.chem_col, nothing, m, :all,
                        ROLE_VARIANT, :both, "",
                        "background exclusion mode of the enrichment step"))
    end
    return out
end

"""
    unit_check_specs(base) -> Vector{Spec}

The unit checks that need a run of their own, one per specification of `base`.

The symmetric genus unit and the per-sequence unit need no run at all — the arm builder emits
them beside the taxon unit for EVERY specification, so they are reported across the whole grid
at no extra cost. The thin-evidence subset does need one: it changes which taxa are in the
arms, and it has to change BOTH arms, so it is a filter on the call table applied before the
arms are built.

It is nevertheless a UNIT and not a specification. It is run for the registered specification
and for every one-at-a-time variant, under the SAME name, and what distinguishes its rows is
the suffix on every unit they carry. That is the honest shape: "how far does the rule choice
move the estimate" and "how far does dropping the thinly annotated tail move it" are two
questions about one specification, and a subset row that wore a variant's name could be read
as a fifteenth rule.
"""
unit_check_specs(base::AbstractVector{Spec}) =
    [Spec(s.name, s.axis, s.level, s.config, s.tie_rule, s.chem_col, s.chem_rule, s.bg_mode,
          :min_distinct, ROLE_UNIT, s.emit, MIN_DISTINCT_SUFFIX,
          "both arms filtered to calls resting on several distinct terms") for s in base]

"""
    anchor_specs(; retired = ANCHOR_LEVELS_RETIRED) -> Vector{Spec}

The two anchors. The all-legacy specification is the earlier pipeline in full — the legacy
configuration, its exclusive winner-chemistry column and no background exclusion at all —
and it is compared against an earlier run's rows. The audit rows are that specification's
arms with the chemistry rules of an earlier sensitivity run, which are the rules this
grid's chemistry axis
generalises; they are emitted on the chemistry unit alone, because a chemistry rule moves
nothing else and a second copy of the legacy type rows under four more names would only
invite them to be read as four measurements. With `retired` true, which is what
`ANCHOR_LEVELS_RETIRED` gives a run by default, it returns no levels at all.
"""
function anchor_specs(; retired::Bool = ANCHOR_LEVELS_RETIRED)
    retired && return Spec[]
    out = Spec[Spec("anchor:all_legacy", "anchor", "all_legacy", HE.ANCHOR_CONFIG, nothing,
                    HE.ANCHOR_CHEM_COL, nothing, "none", :all, ROLE_ANCHOR, :both, "",
                    "the earlier specification in full, compared against an earlier run")]
    for r in AUDIT_CHEMISTRY_RULES
        push!(out, Spec("audit:$r", "audit", String(r), HE.ANCHOR_CONFIG, nothing,
                        HE.ANCHOR_CHEM_COL, r, "none", :all, ROLE_AUDIT, :chemistry_only, "",
                        "an earlier sensitivity rule on the earlier arms"))
    end
    return out
end

"""
    disclosed_specs() -> Vector{Spec}

Background rules reported beside the grid and outside it. They answer "what counts as a
carrier" with a different rule from the registered one, and the recorded decision keeps the
registered rule as the primary and these as disclosed lines; emitting them inside the
variant set would let a disclosed line decide a flag.
"""
disclosed_specs() = [Spec("disclosed:B:$m", "disclosed", m, PRIMARY_CONFIG, nothing,
                          :chemistry_flags, nothing, m, :all, ROLE_DISCLOSED, :both, "",
                          "a disclosed background rule, reported beside the grid")
                     for m in DISCLOSED_BACKGROUND_LEVELS]

function all_specs(; retired::Bool = ANCHOR_LEVELS_RETIRED)
    graded = vcat([primary_spec()], sensitivity_specs())
    return vcat(graded, unit_check_specs(graded), anchor_specs(; retired), disclosed_specs())
end

variant_names(specs) = [s.name for s in specs if s.role == ROLE_VARIANT]

"""
    grid_configurations(specs) -> Vector{String}
    registered_type_rule_configurations(specs) -> Vector{String}

Which configurations of the call table this grid reads, and which of them carry a committed
`type_weights` column the REGISTERED type rule wrote.

Both are derived from the specifications rather than listed, so a configuration that joins or
leaves the grid moves these lists with it. The second is the smaller: the earlier pipeline's
configuration carries its OWN type rule in that column on purpose — that is what the anchor
reproduces — so recomputing the registered rule over it would disagree by design. Every other
configuration the grid reads is one whose committed column the registered rule produced, and
a variant of the text, generic-term or isolation axis IS that column read under a name.
"""
grid_configurations(specs) = sort(unique(String[s.config for s in specs]))
registered_type_rule_configurations(specs) =
    sort(unique(String[s.config for s in specs
                       if s.role in (ROLE_PRIMARY, ROLE_VARIANT, ROLE_UNIT)]))

# The two units the flag reads, split by the family of rules that can move them. Both names
# are taken from the library's own list rather than written again here, so a rename there
# cannot leave this file declaring a unit that no longer exists.
const FLAG_UNITS_TYPE = Tuple(u for u in HE.PRIMARY_UNITS if !is_chemistry_unit(u))
const FLAG_UNITS_CHEM = Tuple(u for u in HE.PRIMARY_UNITS if is_chemistry_unit(u))
const FLAG_UNITS_BOTH = Tuple(u for u in HE.PRIMARY_UNITS)

"""
    VARIANT_AXES

What an axis LETTER means, as data rather than as a naming convention: which field of a
specification the axis varies, the levels that field may take, and the units the axis can
reach.

WHY A NAME HAS TO BE BOUND TO A COMPUTATION. Everything this step publishes about a variant
is published under the variant's NAME — the breaker list beside a category, the column heads
of the written table, the levels named as excluding zero where the primary does not — and a
name is free text. A level renamed onto an axis it does not belong to, or a rule the recorded
decision deliberately keeps OUTSIDE the set running under the name of one that is inside it,
moves no count, breaks no anchor and changes no number: every check that compares a variant
with the PRIMARY still passes, and what is published is a correct number under the name of a
rule that did not produce it. So the name is checked against what the specification actually
runs, in the run, before anything is computed from it.

THE UNITS ARE DECLARED AND NOT INFERRED. A one-at-a-time grid over two orthogonal rule
families leaves cells no axis can reach, and the written table marks them `n.r.`. Derived
from the measured movement alone, that mark cannot tell an axis that was never able to touch
a unit from a variant that LOST its wiring to it — and the loss is then published as
orthogonality. Declaring the map makes those two different statements, and
`assert_axes_reach_the_units_they_declare` makes the difference fatal.
"""
const VARIANT_AXES = (
    (axis = "T", field = :tie_rule,  levels = Tuple(String(r) for r in TYPE_LEVELS),
     units = FLAG_UNITS_TYPE),      # a tie rule decides a habitat call and no chemistry flag
    (axis = "C", field = :chem_rule, levels = Tuple(String(r) for r in CHEMISTRY_LEVELS),
     units = FLAG_UNITS_CHEM),      # and a chemistry rule decides no habitat call
    (axis = "X", field = :config,    levels = TEXT_CONFIGS,
     units = FLAG_UNITS_CHEM),      # the chemistry text is read for the chemistry flags alone
    (axis = "G", field = :config,    levels = (GENERIC_CONFIG,),
     units = FLAG_UNITS_BOTH),      # which terms are usable moves both tallies
    (axis = "I", field = :config,    levels = (ISOLATION_CONFIG,),
     units = FLAG_UNITS_BOTH),      # a filter on the taxa moves both arms
    (axis = "B", field = :bg_mode,   levels = BACKGROUND_LEVELS,
     units = FLAG_UNITS_BOTH))      # and so does a different background

"""
    axis_declaration(axis) -> NamedTuple or nothing

The registered meaning of one axis letter; `nothing` for a letter that is not a registered
axis, which is a caller's error to report rather than this function's to absorb.
"""
function axis_declaration(axis::AbstractString)
    for a in VARIANT_AXES
        a.axis == axis && return a
    end
    return nothing
end

# THE SIX FIELDS THAT DECIDE WHAT A SPECIFICATION COMPUTES. The configuration of the call
# table, the two aggregation rules, the column the chemistry arm is read off, the background
# rule and the subset of taxa: change any one of them and a different computation runs.
# `chem_col` and `subset` are in the list because they are part of what a specification IS by
# this file's own definition of one — a variant quietly measured on the thin-evidence subset,
# or quietly reading the earlier pipeline's exclusive-winner chemistry column, is a different
# measurement published under a registered name, and a tuple that omitted them could not say so.
const DEFINING_FIELDS = (:config, :tie_rule, :chem_col, :chem_rule, :bg_mode, :subset)

# One field of a specification as text. `nothing` — "the level the call table already carries"
# — is the empty string, which no named rule can collide with.
field_level(s::Spec, field::Symbol) = (v = getfield(s, field); v === nothing ? "" : String(v))

"""
    defining_level(spec) -> NTuple{6,String}

The six things a specification is made of (`DEFINING_FIELDS`). Two specifications that agree
on all six run the same computation, whatever they are called — so this tuple, and not the
name, is a specification's identity.
"""
defining_level(s::Spec) = Tuple(field_level(s, f) for f in DEFINING_FIELDS)

"""
    assert_variants_vary_the_primary(specs) -> Int

Every one-at-a-time variant must differ from the registered specification in at least one of
those six. The set is declared in this file, so nothing outside it says what belongs in it;
this is the property that keeps the declaration honest. A level that agrees with the primary
on all six IS the primary, and re-running the registered specification under a variant's name
would give a category a variant that cannot break it — which inflates the flag by arithmetic
rather than by evidence, and does it while every count and every uniqueness check still
passes.
"""
function assert_variants_vary_the_primary(specs)
    p = defining_level(primary_spec())
    same = sort([s.name for s in specs if s.role == ROLE_VARIANT && defining_level(s) == p])
    isempty(same) ||
        error("the variant(s) " * join(same, ", ") * " have the registered specification's " *
              "own configuration, chemistry column, type rule, chemistry rule, background " *
              "rule and subset: they are the registered specification under another name, " *
              "and a set that counts them reports agreement it never tested")
    return count(s -> s.role == ROLE_VARIANT, specs)
end

"""
    assert_variants_differ_from_each_other(specs) -> Int

No two variants may be the same specification.

`assert_variants_vary_the_primary` compares every variant with the PRIMARY and with nothing
else, so two variants that are clones of each other pass it: each still differs from the
registered specification, each still moves estimates, the names are all present and all
unique, and every anchor is met. What such a pair does is remove a level from the set without
removing a column — the category that level broke keeps its other breakers and loses that
one, and a row the grid called fragile can be published as robust with an empty breaker list.
The count of variants, which is what a reader sees, does not move either.

The number returned is the number of DISTINCT specifications among the variants, which equals
the number of variants exactly when none of them is a clone.
"""
function assert_variants_differ_from_each_other(specs)
    seen = Dict{NTuple{length(DEFINING_FIELDS),String},Vector{String}}()
    for s in specs
        s.role == ROLE_VARIANT || continue
        push!(get!(seen, defining_level(s), String[]), s.name)
    end
    clones = sort([sort(v) for v in values(seen) if length(v) > 1])
    isempty(clones) ||
        error("these variants are the same specification under different names: " *
              join([join(c, " = ") for c in clones], "; ") * ". Each of them still differs " *
              "from the registered specification, so the comparison against the primary " *
              "passes; what a clone does is delete a level from the set while leaving its " *
              "column standing, and a category whose only other breaker was the level that " *
              "was overwritten is then published as robust")
    return length(seen)
end

"""
    assert_the_set_is_the_declared_one(enumerated, declared) -> Int

The levels this run enumerated must be exactly the levels declared OUTSIDE it, name for name,
in both directions.

WHY THE DECLARATION CANNOT LIVE IN THIS FILE. Every other check here compares one part of this
file with another part of it, and `VARIANT_AXES` reads its levels off the very tuples
`sensitivity_specs` enumerates. A level deleted from a tuple, or replaced by a rule the
recorded decision keeps outside the set, is therefore SELF-CONSISTENT: the wiring and the
vocabulary move together, the name binding still passes, the counts still agree with
themselves, and the run publishes the set it happens to hold as the set it was supposed to
hold. What that costs is a column: a category whose breakers were the levels that left is
published with a shorter breaker list, or with none, and the word beside it can change from
fragile to a colour while every number in the run is correct. Nothing inside a file can refuse
an edit to that file, so the set arrives as an input — the same input the anchor targets
arrive in, because a rule NAME is not a quantity measured on any data and this file may carry
none of those.

Both directions are refused. A level that is declared and did not run took its column with it;
a level that ran and is not declared decides a published result under a name the registered
set does not carry.

WHAT THE FAILURE MAY NOT SAY. A guard whose message names the level the declaration holds and
the run does not has written half an edit: deleting that name from the declaration makes the
failure go away, makes every count agree with itself again, and publishes the shorter grid
under the registered set's name — which is the precise outcome this guard exists to refuse.
So the difference is named, because a failure nobody can diagnose is a failure nobody fixes,
and the message then says outright that the two repairs are not equivalent and which of them
is the defect. The declared SIZE of the set is checked separately, in `read_anchor_targets`,
so that making the two lists agree by shortening one of them leaves the file inconsistent
with its own count rather than quietly consistent.
"""
function assert_the_set_is_the_declared_one(enumerated, declared)
    have = sort(unique(String[String(v) for v in enumerated]))
    want = sort(unique(String[String(v) for v in declared]))
    have == want && return length(want)
    parts = String[]
    gone  = setdiff(want, have)
    extra = setdiff(have, want)
    isempty(gone) ||
        push!(parts, "in the declared set and missing from this run: " * join(gone, ", "))
    isempty(extra) ||
        push!(parts, "run by this grid and absent from the declared set: " *
                     join(extra, ", "))
    error("this run enumerated $(length(have)) one-at-a-time level(s) against a declared set " *
          "of $(length(want)); " * join(parts, "; ") *
          ". WHICH SIDE IS WRONG IS NOT SOMETHING THIS MESSAGE CAN KNOW, and the two repairs " *
          "are not equivalent. Restoring a missing level to the run leaves the grid the " *
          "registered set describes. Striking that level from the declaration instead " *
          "publishes a SMALLER grid under the registered set's name: a category the level " *
          "broke then carries a shorter breaker list or none at all, and the word beside it " *
          "can change from fragile to a colour while every number in the run is correct and " *
          "every count in it agrees with itself. That second edit is the failure this check " *
          "exists to refuse, not the way to clear it. Reconcile both lists against the " *
          "recorded decision that fixed the set")
end

"""
    assert_variant_names_describe_what_they_run(specs) -> Int

Every variant's name must describe the computation it runs: it is `"<axis>:<level>"`, the
axis is a registered axis, the level is a registered level OF that axis, the field that axis
varies carries exactly that level, and every other defining field is the registered
specification's own — because a one-at-a-time set is one-at-a-time or it is not that set.

This is the guard a name list cannot be and a comparison with the primary cannot be. A name
list says the fourteen names are present; a comparison with the primary says each of them
differs from the registered specification. Neither says that the level called
`<axis>:<level>` is the one that ran, and every published statement about a variant is made
under its name.
"""
function assert_variant_names_describe_what_they_run(specs)
    p = primary_spec()
    bad = String[]
    for s in specs
        s.role == ROLE_VARIANT || continue
        a = axis_declaration(s.axis)
        if a === nothing
            push!(bad, "$(s.name): '$(s.axis)' is not one of the registered axes (" *
                       join([x.axis for x in VARIANT_AXES], ", ") * ")")
            continue
        end
        s.name == "$(s.axis):$(s.level)" ||
            push!(bad, "$(s.name): a variant is named '<axis>:<level>', and this one's own " *
                       "axis and level make '$(s.axis):$(s.level)'")
        s.level in a.levels ||
            push!(bad, "$(s.name): '$(s.level)' is not a registered level of axis " *
                       "$(s.axis) (" * join(a.levels, ", ") * "), so it is a rule from " *
                       "outside the set wearing a name from inside it")
        field_level(s, a.field) == s.level ||
            push!(bad, "$(s.name): axis $(s.axis) varies $(a.field), and this one runs " *
                       "'$(field_level(s, a.field))' where its name says '$(s.level)'")
        for f in DEFINING_FIELDS
            f === a.field && continue
            field_level(s, f) == field_level(p, f) ||
                push!(bad, "$(s.name): axis $(s.axis) varies $(a.field) alone, and this one " *
                           "also changes $f from the registered '$(field_level(p, f))' to " *
                           "'$(field_level(s, f))'")
        end
    end
    isempty(bad) ||
        error("the name of a variant is the only thing this step publishes about it, and " *
              "these names do not describe what they run:\n  " * join(sort(bad), "\n  "))
    return count(s -> s.role == ROLE_VARIANT, specs)
end

# What two rows have to agree on to be the same estimate. Every one is a column of the
# long-format table, so "this variant moved nothing" is a statement about the published
# numbers rather than about an intermediate nobody can see.
const COMPARED_QUANTITIES = (:a, :b, :n1, :n2, :crude, :crude_lo, :crude_hi, :std, :lo, :hi)

estimates_agree(x::Real, y::Real) = (isnan(x) && isnan(y)) || x == y
rows_agree(p, r) =
    all(c -> estimates_agree(fnum(getproperty(p, c)), fnum(getproperty(r, c))),
        COMPARED_QUANTITIES)

"""
    estimate_movement(rows, expected_variants) -> (per_cell, per_unit, per_variant)

Which variants actually moved which estimate, counted on the rows the flag reads.

`per_cell[(unit, category)]` is how many variants differ from the registered specification
there; `per_unit[(variant, unit)]` how many of that unit's categories the variant moved; and
`per_variant[variant]` how many cells it moved anywhere. A variant with no row where the
primary has one — or with one where the primary has none — counts as having moved that cell:
the two tables say different things about it, which is the only sense "moved" has here.

WHY THIS IS COUNTED AT ALL. A one-at-a-time grid over two orthogonal rule families produces
cells no axis can reach: a chemistry rule cannot move a habitat assignment and a type rule
cannot move a chemistry flag. Those cells are the primary's own numbers printed again, and a
column of them looks exactly like a variant that agreed. Counting the movement is what lets
the output say which it is, per unit and per category, instead of letting "fourteen variants
agree" stand for a number of independent checks that was never run.
"""
function estimate_movement(rows, expected_variants)
    flagged = [r for r in rows if r.feeds_flag]
    byk = Dict{Tuple{String,String,String},Any}()
    for r in flagged; byk[(r.variant, r.unit, r.category)] = r; end
    cells = sort(unique([(r.unit, r.category) for r in flagged]))
    names = String[String(v) for v in expected_variants]
    units = sort(unique([u for (u, _) in cells]))
    per_cell = Dict{Tuple{String,String},Int}()
    per_unit = Dict{Tuple{String,String},Int}((v, u) => 0 for v in names for u in units)
    per_variant = Dict{String,Int}(v => 0 for v in names)
    for (unit, cat) in cells
        p = get(byk, ("primary", unit, cat), nothing)
        n = 0
        for v in names
            r = get(byk, (v, unit, cat), nothing)
            agreed = p === nothing ? r === nothing :
                     r !== nothing && rows_agree(p, r)
            agreed && continue
            n += 1
            per_unit[(v, unit)] += 1
            per_variant[v] += 1
        end
        per_cell[(unit, cat)] = n
    end
    return (per_cell = per_cell, per_unit = per_unit, per_variant = per_variant,
            units = units)
end

"""
    assert_every_variant_moves_something(movement, expected_variants) -> Int

A variant that reproduces the registered specification on every unit it touches is refused,
by name.

An axis that cannot reach a UNIT is legitimate and expected — it is the orthogonality the
per-unit counts report — but a level that moves nothing ANYWHERE has either lost its wiring
to the rule it was supposed to vary or is the registered specification under another name.
Both produce a variant that cannot break a result, and the flag would count it as one that
agreed. This is the run-time half of `assert_variants_vary_the_primary`: that one compares
declarations, this one compares the numbers they produced.
"""
function assert_every_variant_moves_something(movement, expected_variants)
    inert = sort([String(v) for v in expected_variants
                  if get(movement.per_variant, String(v), 0) == 0])
    isempty(inert) ||
        error("the variant(s) " * join(inert, ", ") * " reproduce the registered " *
              "specification on every unit the flag reads: a level that moves no estimate " *
              "anywhere is a wiring failure or the primary under another name, and neither " *
              "is a variant that agreed. An axis that cannot reach ONE unit is expected and " *
              "is reported per unit; an axis that reaches nothing is not")
    return length(expected_variants)
end

"""
    variant_pair_differences(rows, expected_variants) -> (pairs, units)

Every unordered pair of levels, with how many of the flag's cells either of them measured and
how many of those the two disagree on — once overall and once per unit.

ONE COMPARISON, TWO READERS. The guard below needs "is this pair zero cells apart", and the
run's own disclosure needs "how far is the closest pair from being that", per unit. Computing
those separately would be two implementations of the same comparison, and the one the guard
does not use is the one that would drift into describing a different question from the one
being refused.

A cell one of them wrote a row for and the other did not counts as a difference, exactly as it
counts as movement in `estimate_movement`: the two tables say different things about it.
"""
function variant_pair_differences(rows, expected_variants)
    names = sort(unique(String[String(v) for v in expected_variants]))
    byk = Dict{Tuple{String,String,String},Any}()
    cells = Set{Tuple{String,String}}()
    for r in rows
        r.feeds_flag || continue
        byk[(r.variant, r.unit, r.category)] = r
        push!(cells, (r.unit, r.category))
    end
    cellv = sort(collect(cells))
    units = sort(unique(String[String(c[1]) for c in cellv]))
    out = NamedTuple[]
    for i in eachindex(names), j in (i + 1):lastindex(names)
        v, w = names[i], names[j]
        compared = 0; differing = 0
        cu = Dict{String,Int}(u => 0 for u in units)
        du = Dict{String,Int}(u => 0 for u in units)
        for c in cellv
            p = get(byk, (v, c[1], c[2]), nothing)
            q = get(byk, (w, c[1], c[2]), nothing)
            p === nothing && q === nothing && continue
            compared += 1; cu[String(c[1])] += 1
            if !(p !== nothing && q !== nothing && rows_agree(p, q))
                differing += 1; du[String(c[1])] += 1
            end
        end
        push!(out, (a = v, b = w, compared = compared, differing = differing,
                    compared_per_unit = cu, differing_per_unit = du))
    end
    return (pairs = out, units = units)
end

"""
    variant_pair_margin(differences, specs) -> (overall, per_unit)

How close the closest LEGITIMATE pair of levels comes to being one measurement — overall, and
on each unit separately.

WHY PER UNIT AND NOT ONLY OVERALL. The guard refuses a pair that differs on no cell the flag
reads, so a margin read against every cell the flag reads looks like far more headroom than
it is. Two levels of the chemistry family are identical on every `taxa` cell BY CONSTRUCTION
— their axis is declared unable to reach that unit — and two levels of the type family are
identical on every chemistry cell for the same reason. The cells that could have
separated them are only the ones both their axes reach, so the margin that means anything is
the one taken over those. The guard is still sound either way, because for a pair of levels on
one axis agreeing on every reachable cell implies agreeing on all of them; what changes is the
size a reader takes the headroom to be.

Pairs are restricted to those whose two axes are BOTH declared to reach the unit, which is the
same declaration `assert_axes_reach_the_units_they_declare` holds the measurement to.
"""
function variant_pair_margin(differences, specs)
    reach = Dict{String,Set{String}}()
    for s in specs
        s.role == ROLE_VARIANT || continue
        a = axis_declaration(s.axis)
        reach[s.name] = a === nothing ? Set{String}() : Set(String[String(u) for u in a.units])
    end
    best(pairs, value, denom) = isempty(pairs) ? nothing :
        (p = argmin(x -> (value(x), x.a, x.b), pairs);
         (pair = "$(p.a) vs $(p.b)", differing = value(p), compared = denom(p)))
    overall = best([p for p in differences.pairs if p.compared > 0],
                   p -> p.differing, p -> p.compared)
    per_unit = Pair{String,Any}[]
    for u in differences.units
        ps = [p for p in differences.pairs
              if u in get(reach, p.a, Set{String}()) && u in get(reach, p.b, Set{String}()) &&
                 p.compared_per_unit[u] > 0]
        m = best(ps, p -> p.differing_per_unit[u], p -> p.compared_per_unit[u])
        m === nothing || push!(per_unit, u => m)
    end
    return (overall = overall, per_unit = per_unit)
end

"""
    assert_no_two_variants_are_one_measurement(rows, expected_variants) -> Int

No two levels may produce the same numbers on every cell the flag reads.

`assert_variants_differ_from_each_other` compares DECLARATIONS, and on this set every
declaration differs by construction. What a declaration cannot see is two levels whose
declarations differ and whose NUMBERS do not — and that needs no edit to this file at all. The
five non-primary configurations of the call table are produced by one upstream script; a block
written twice, a configuration label assigned in the wrong loop iteration, or a re-export that
reused a frame makes two of them the same table, and the two levels that read them are then
one measurement published as two columns. The set has lost a level while keeping its name, the
count of variants does not move, every anchor is met, and a category both of them broke keeps
one breaker where a reader counts two.

WHAT IS REFUSED, AND WHAT IS NOT. Two levels of one axis are EXPECTED to be close, and a level
that differs from its sibling on a single taxon is a legitimate near-clone: it still separates
itself somewhere, it moves fewer cells than its sibling, and `variant_cells_moved` reports
that beside it. What is refused is a pair with no cell anywhere between them — not a small
difference but no difference at all. A cell one of them wrote a row for and the other did not
counts as a difference, for the same reason it counts as movement in `estimate_movement`: the
two tables say different things about it.

The number returned is how many distinct MEASUREMENTS the variants made, which equals the
number of variants exactly when no two of them are one measurement.
"""
function assert_no_two_variants_are_one_measurement(rows, expected_variants,
        differences = variant_pair_differences(rows, expected_variants))
    names = sort(unique(String[String(v) for v in expected_variants]))
    # Each level's representative: the earliest level it is indistinguishable from. Exact
    # agreement over the same cells is transitive, so scanning the pairs in order is enough
    # and no second pass is needed.
    group = Dict{String,String}(v => v for v in names)
    same = String[]
    for pr in differences.pairs
        # A pair with no cell in common is not two names for one measurement; it is a grid
        # with a hole in it, which `assert_grid_complete` is the guard for.
        (pr.compared > 0 && pr.differing == 0) || continue
        push!(same, "$(pr.a) = $(pr.b), on all $(pr.compared) of the cells either of them " *
                    "measured")
        group[pr.b] = group[pr.a]
    end
    isempty(same) ||
        error("these levels are two names for one measurement: " * join(same, "; ") *
              ". A level that differs from its sibling on a single taxon is a legitimate " *
              "near-clone — it separates itself somewhere, and it is reported as having " *
              "moved fewer cells; what this refuses is a pair with NO cell between them, " *
              "which is one measurement published as two columns. Their declarations all " *
              "differ, so the comparison of declarations passes; the duplicate is in what " *
              "they were fed, not in what they say they are")
    return length(unique(values(group)))
end

"""
    assert_axes_reach_the_units_they_declare(movement, specs) -> Int

The units each variant moved must be exactly the units `VARIANT_AXES` declares its axis can
reach — no fewer, and no more.

WHY EQUALITY AND NOT A DIRECTION. `assert_every_variant_moves_something` fires only for a
level that moved nothing ANYWHERE, so a variant still wired to one unit and disconnected from
the other passes it; and the written table, which derives `n.r.` from the same measured
movement, then prints the disconnection as the sentence "that rule family does not touch this
unit at all". The claim would be made by the failure itself. Declaring the reach turns a lost
wire into a contradiction between the declaration and the measurement, and a unit an axis
reached that it was declared unable to reach — a type rule that moved a chemistry flag — into
the same kind of contradiction pointed the other way.
"""
function assert_axes_reach_the_units_they_declare(movement, specs)
    bad = String[]
    for s in specs
        s.role == ROLE_VARIANT || continue
        a = axis_declaration(s.axis)
        a === nothing && continue      # the name guard has already refused this one, by name
        declared = Set(String[String(u) for u in a.units])
        measured = Set(String[u for u in movement.units
                              if get(movement.per_unit, (s.name, u), 0) > 0])
        measured == declared && continue
        for u in sort(collect(setdiff(declared, measured)))
            push!(bad, "$(s.name) moved no estimate on '$u', which axis $(s.axis) is " *
                       "declared to reach: the level has lost its wiring to that unit, and " *
                       "the table would report the loss as orthogonality")
        end
        for u in sort(collect(setdiff(measured, declared)))
            push!(bad, "$(s.name) moved an estimate on '$u', which axis $(s.axis) is " *
                       "declared unable to reach: it is varying something other than the " *
                       "rule its name names")
        end
    end
    isempty(bad) ||
        error("the grid did not reach the units its axes declare:\n  " *
              join(sort(bad), "\n  "))
    return count(s -> s.role == ROLE_VARIANT, specs)
end

"""
    declared_reachable_cells(spec, movement) -> Int

How many of the flag's cells this variant's axis is declared able to move, which is the
denominator the count of cells it DID move is read against. A variant that moves one cell of
the cells its axis can reach and reproduces the registered specification on all the rest
passes every guard here — it varies the primary, its name describes what it runs, it
is nobody's clone and it moved something — and it is what a half-applied configuration or a
partial re-export actually looks like. Nothing in the recorded rules sets a threshold on how
far a variant must move, so none is invented: the two numbers are recorded side by side, per
variant, and a one-cell level is then visible beside its five-cell sibling.
"""
function declared_reachable_cells(spec::Spec, movement)
    a = axis_declaration(spec.axis)
    a === nothing && return 0
    return count(k -> String(k[1]) in a.units, keys(movement.per_cell))
end

# ---------------------------------------------------------------------------
# Tallies, weights and the two rule seams.
# ---------------------------------------------------------------------------

"""
    parse_tally(s) -> Dict{String,Float64}

Read a `"Key:value;Key:value"` column back into the tally it was written from. The key is
everything before the LAST colon, as in the weight parser the arm builder uses, so a key
that contains a colon cannot silently lose its tail. A blank cell is an empty tally, which
is a taxon with nothing to tally rather than a taxon with a zero.
"""
function parse_tally(s::AbstractString)
    d = Dict{String,Float64}()
    isempty(strip(s)) && return d
    for part in split(s, ';')
        isempty(strip(part)) && continue
        j = findlast(==(':'), part)
        j === nothing && error("malformed tally entry '$part'")
        k = String(strip(part[1:prevind(part, j)]))
        v = tryparse(Float64, strip(part[nextind(part, j):end]))
        v === nothing && error("malformed tally value in '$part'")
        d[k] = get(d, k, 0.0) + v
    end
    return d
end

"""
    weight_string(d) -> String

A per-category weight map as the arm builder's `weight_map` reads it, at FULL precision.

Deliberately not the call table's own formatter: that one rounds for reading, and a rounded
weight is a different hit total. A grid whose whole subject is how far a hit total moves
under a rule must not move it by writing the number down. `string` of a Float64 is the
shortest text that parses back to the same value, so the round trip through the column is
exact. The ordering is the table's — heaviest first, ties by name — so a written column is
still comparable by eye with the published one.
"""
weight_string(d::AbstractDict) =
    join([string(k, ":", string(float(v)))
          for (k, v) in sort(collect(d), by = kv -> (-float(kv[2]), String(kv[1])))], ";")

"""
    combination_counts(profile, n_terms) -> Dict{String,Int}

The tally of chemistry-bearing terms per flag COMBINATION, recovered as whole numbers from
the written share map and the number of chemistry-bearing terms.

WHY NOT THE SHARES THEMSELVES. The share map is rounded for reading, and a flag's share is
the sum of the shares of every combination naming it; summing rounded shares can land a
flag on either side of a threshold that the counts decide cleanly. The counts are integers
and the rounding is far coarser than a half count, so multiplying the share back by the
denominator recovers them exactly — and `validate_combinations!` proves that on the table in
hand rather than asserting it here.
"""
function combination_counts(profile::AbstractString, n_terms::Integer)
    out = Dict{String,Int}()
    for (k, v) in parse_tally(profile)
        out[k] = round(Int, v * n_terms)
    end
    return out
end

"""
    validate_combinations!(df) -> (rows, checked)

Check the recovered combination tally against the two quantities the table wrote beside it:
it must sum to the number of chemistry-bearing terms, and the mass it gives each flag must
equal the per-flag counts. Either failure means the share map and the counts describe
different evidence, and a chemistry variant computed from the first would be reported under
the name of a rule applied to the second.
"""
function validate_combinations!(df::AbstractDataFrame)
    checked = 0
    for i in 1:nrow(df)
        n_terms = Int(coalesce(df.n_chem_terms[i], 0))
        cc = combination_counts(sv(df.chemistry_profile[i]), n_terms)
        isempty(cc) && continue
        checked += 1
        total = sum(values(cc))
        total == n_terms ||
            error("the chemistry combination tally of taxon $(df.taxid[i]) sums to $total " *
                  "and not to its $n_terms chemistry-bearing terms: the share map and the " *
                  "term count describe different evidence")
        mass = Dict{String,Int}()
        for (c, v) in cc, f in split(c, '+')
            isempty(f) && continue
            mass[String(f)] = get(mass, String(f), 0) + v
        end
        written = Dict{String,Int}(k => round(Int, v)
                                   for (k, v) in parse_tally(sv(df.chemistry_flag_counts[i])))
        mass == written ||
            error("the chemistry combinations of taxon $(df.taxid[i]) give the flags " *
                  "$(mass) where the table counted $(written): a flag's share would not be " *
                  "the share of the terms that carry it")
    end
    return (rows = nrow(df), checked = checked)
end

"""
    chemistry_weights(combos, rule; share_threshold) -> Dict{String,Float64}

The per-flag weights of one taxon under one chemistry rule.

Three of the four rules are `EnvironmentTools.resolve_chemistry`'s own. The fourth settles
the tie among the combinations at the top of the tally ALPHABETICALLY and scores every flag
that combination names — an earlier sensitivity rule this grid has to reproduce, which is
`resolve`'s alphabetical clause applied to the combination tally rather than to a type
tally. It is expressed that way, and not written out again here, so the two alphabetical
rules in this pipeline cannot come apart.
"""
function chemistry_weights(combos::AbstractDict, rule::Symbol;
                           share_threshold::Real = ET.SHARE_THRESHOLD)
    rule === COMBINATION_ALPHABETICAL ||
        return ET.resolve_chemistry(combos, rule; share_threshold = share_threshold)
    isempty(combos) && return Dict{String,Float64}()
    (label, _, _, _) = ET.resolve(combos, Dict{String,Float64}(), :alphabetical)
    return Dict{String,Float64}(String(f) => 1.0 for f in split(label, '+') if !isempty(f))
end

"""
    type_weight_column(df, rule; share_threshold) -> Vector{String}

The `type_weights` column a type rule would have written, recomputed from the table's own
`type_counts` and `type_mentions` through `EnvironmentTools.resolve`.

A taxon the rule declines to call carries an EMPTY weight string and keeps its habitat call,
so it stays one taxon of its arm and scores zero everywhere. Blanking the habitat call
instead would take it out of the arm, and a rule meant to change only the assignment would
be changing the denominator too.
"""
function type_weight_column(df::AbstractDataFrame, rule::Symbol;
                            share_threshold::Real = ET.SHARE_THRESHOLD)
    out = Vector{String}(undef, nrow(df))
    for i in 1:nrow(df)
        counts   = parse_tally(sv(df.type_counts[i]))
        mentions = parse_tally(sv(df.type_mentions[i]))
        (_, w, _, _) = ET.resolve(counts, mentions, rule; share_threshold = share_threshold)
        out[i] = weight_string(w)
    end
    return out
end

"""
    chemistry_weight_column(df, rule; share_threshold) -> Vector{String}

The per-flag weights of a chemistry rule, written into the column the arm builder scores.

WHY THE WEIGHT COLUMN AND NOT THE CHEMISTRY COLUMN. The arm builder reads its chemistry
column with `flag_weights`, which turns an already-decided LABEL into 0/1 weights and cannot
express a fraction; handed "Flag:0.5" it would read a category of that name and score it
whole. The weight column is parsed by `weight_map`, which is the parser that expresses both,
and which narrows an all-0/1 column back to the exact Bool path. The rows this produces are
therefore the chemistry rows of that specification, and the run checks that claim under the
registered rule before trusting it under any other.
"""
function chemistry_weight_column(df::AbstractDataFrame, rule::Symbol;
                                 share_threshold::Real = ET.SHARE_THRESHOLD)
    out = Vector{String}(undef, nrow(df))
    for i in 1:nrow(df)
        cc = combination_counts(sv(df.chemistry_profile[i]),
                                Int(coalesce(df.n_chem_terms[i], 0)))
        out[i] = weight_string(chemistry_weights(cc, rule; share_threshold = share_threshold))
    end
    return out
end

"""
    type_reconstruction_check(df) -> (rows, mismatches, example)

Recompute the REGISTERED type rule from the table's tallies and compare it, in the table's
own written form, with the column the table carries.

This is the licence for every type variant. The variants are recomputed from `type_counts`
and `type_mentions`; if the registered rule recomputed from those same two columns does not
reproduce the committed column, then they are not the tallies the column was built from, and
every variant built on them is a rule applied to the wrong evidence. `example` names the
first disagreement so a failure is diagnosable rather than merely fatal.
"""
function type_reconstruction_check(df::AbstractDataFrame)
    mismatches = 0; example = ""
    for i in 1:nrow(df)
        counts   = parse_tally(sv(df.type_counts[i]))
        mentions = parse_tally(sv(df.type_mentions[i]))
        (_, w, _, _) = ET.resolve(counts, mentions, :mentions_fractional)
        written = ET._weights_str(w)
        if written != sv(df.type_weights[i])
            mismatches += 1
            isempty(example) &&
                (example = "taxon $(df.taxid[i]): recomputed '$written', table '" *
                           sv(df.type_weights[i]) * "'")
        end
    end
    return (rows = nrow(df), mismatches = mismatches, example = example)
end

"""
    configuration_rows(calls, config) -> DataFrame

One configuration of the call table, refused rather than silently empty. A seam check run
over nothing reports "0 of 0 rows disagree", which is the same sentence a passing check
writes.
"""
function configuration_rows(calls::AbstractDataFrame, config::AbstractString)
    df = calls[sv.(calls.config) .== String(config), :]
    nrow(df) > 0 ||
        error("the call table has no rows for configuration '$config', which this grid " *
              "reads; it has: " * join(sort(unique(sv.(calls.config))), ", "))
    return df
end

"""
    type_reconstruction_over(calls, configs) -> (per_configuration, rows, mismatches, example)

`type_reconstruction_check` over every configuration whose committed weight column the
registered type rule wrote, not over the primary alone.

WHY THE WHOLE SET AND NOT THE ONE CONFIGURATION THE VARIANTS ARE REWRITTEN FROM. The type
variants are recomputed from the primary configuration's tallies, so checking that one is what
licenses the rewrite — and that was the whole of this check. But four of the other
configurations ARE registered levels: the text, generic-term and isolation variants are their
committed weight columns read under a level's name, and nothing recomputed those from their own
tallies. A configuration whose weight column was re-exported from a sibling — the shape a
regeneration that covered one arm and not the other leaves behind — keeps its own tallies,
its own counts and its own row count, so every other guard here stays silent: the declaration
is untouched, the name still describes the rule, the level still varies the primary and still
moves something, and it reaches exactly the units its axis declares. The tallies are the one
place the substitution shows.

`per_configuration` carries the row count each check was made over, because a check whose scope
is not recorded is a check whose scope can shrink without the report changing.
"""
function type_reconstruction_over(calls::AbstractDataFrame, configs)
    per = Pair{String,Int}[]
    rows = 0; mismatches = 0; example = ""
    for c in configs
        r = type_reconstruction_check(configuration_rows(calls, String(c)))
        push!(per, String(c) => r.rows)
        rows += r.rows; mismatches += r.mismatches
        (r.mismatches == 0 || !isempty(example)) ||
            (example = "configuration '$c', " * r.example)
    end
    return (per_configuration = per, rows = rows, mismatches = mismatches, example = example)
end

"""
    combinations_over(calls, configs) -> (per_configuration, per_configuration_with_a_flag,
                                          rows, checked)

`validate_combinations!` over every configuration the grid reads, including the earlier
pipeline's — its chemistry columns are what the audit rules are recomputed from, and a share
map that describes different evidence from the counts beside it would put those rules on the
wrong tally under the right name. The configuration is named in the failure, because the taxon
id alone does not say which of six tables holds the contradiction.
"""
function combinations_over(calls::AbstractDataFrame, configs)
    per = Pair{String,Int}[]; per_flag = Pair{String,Int}[]
    rows = 0; checked = 0
    for c in configs
        df = configuration_rows(calls, String(c))
        r = try
            validate_combinations!(df)
        catch e
            e isa ErrorException || rethrow()
            error("in configuration '$c': " * e.msg)
        end
        push!(per, String(c) => r.rows); push!(per_flag, String(c) => r.checked)
        rows += r.rows; checked += r.checked
    end
    return (per_configuration = per, per_configuration_with_a_flag = per_flag,
            rows = rows, checked = checked)
end

# ---------------------------------------------------------------------------
# Running one specification.
# ---------------------------------------------------------------------------

"""
    context(; calls, scan, ca_class_taxids, movein_taxids, min_distinct) -> NamedTuple

Everything a specification needs that is read once: the call table, the same table filtered
to the thinly-annotated-tail subset, the exclusion set of every background mode the grid
uses, the pool the move-in level moves, and the disclosed calcium arm's membership.

The exclusion sets are built HERE, once, by the enrichment step's own `exclusion_taxids`,
so every background level of this grid removes exactly what a run of that mode removes.
A mode the scan table on hand cannot answer is an error at this point rather than a level
that quietly excluded nobody.
"""
function context(; calls::DataFrame, scan, ca_class_taxids, movein_taxids::Set{Int},
                 min_distinct::Integer = MIN_DISTINCT_TERMS)
    modes = unique(vcat(String[HE.PRIMARY_MODE], collect(BACKGROUND_LEVELS),
                        collect(DISCLOSED_BACKGROUND_LEVELS), String["none"]))
    excl = Dict{String,Set{Int}}()
    for m in modes
        e = HE.exclusion_taxids(scan, m)
        e === nothing &&
            error("the background rule '$m' cannot be computed from the scan table on hand; " *
                  "a level that cannot be answered must not be run as a level that excluded " *
                  "nobody")
        excl[m] = e
    end
    keep = HS.well_supported_mask(calls.n_distinct; min_distinct = min_distinct)
    return (calls = calls, calls_min = calls[keep, :], exclusions = excl,
            movein_taxids = movein_taxids, ca_class_taxids = ca_class_taxids,
            min_distinct = min_distinct)
end

"""
    weighted_calls(calls, config, weights; share_threshold) -> DataFrame

The call table of one configuration, with the column the arm builder scores replaced by the
weights of the rule being varied — or handed back untouched when the specification runs on
the level the table already carries.
"""
function weighted_calls(calls::AbstractDataFrame, config::AbstractString, weights;
                        share_threshold::Real = ET.SHARE_THRESHOLD)
    df = calls[sv.(calls.config) .== config, :]
    nrow(df) > 0 ||
        error("the call table has no rows for configuration '$config'; it has: " *
              join(sort(unique(sv.(calls.config))), ", "))
    weights === nothing && return df
    kind, rule = weights
    col = kind === :type ? type_weight_column(df, rule; share_threshold = share_threshold) :
          kind === :chemistry ?
              chemistry_weight_column(df, rule; share_threshold = share_threshold) :
          error("a weight column is built from a type rule or a chemistry rule, not from '$kind'")
    out = copy(df)
    out.type_weights = col
    return out
end

"""
    build_arms(ctx, spec; weights, share_threshold) -> NamedTuple

One call of the arm builder, for one specification and one of its two weight axes.

The move-in level is the only one that hands taxa to the FOREGROUND, and it is the only one
that can be silently inert: its taxa carry their domain in a column the census leaves empty,
so a prokaryote test read off the census column alone moves nobody and returns the unmoved
specification under this level's name. The builder refuses that case; this function refuses
the other half of it, a pool that named nobody at all, which the builder cannot see.

AN EMPTY ARM IS NAMED HERE. A specification with nothing on one side is not a specification
with a small estimate: every interval it writes is undefined, and downstream the failure
arrives as a thin-evidence subset that "left one of the arms unchanged" — a true statement
about the wrong thing, one whole grid later. A specification with no subset twin (an anchor,
a disclosed line) would not even reach that check. So the arm that is empty is named at the
point it is built, against the specification that built it.
"""
function build_arms(ctx, spec::Spec; weights = nothing,
                    share_threshold::Real = ET.SHARE_THRESHOLD)
    source = spec.subset === :min_distinct ? ctx.calls_min : ctx.calls
    df = weighted_calls(source, spec.config, weights; share_threshold = share_threshold)
    haskey(ctx.exclusions, spec.bg_mode) ||
        error("no exclusion set was built for the background rule '$(spec.bg_mode)'")
    moved = spec.bg_mode == HE.MOVEIN_MODE ? ctx.movein_taxids : Set{Int}()
    if spec.bg_mode == HE.MOVEIN_MODE
        isempty(moved) &&
            error("the move-in level '$(spec.name)' was given an empty pool: a level that " *
                  "moves nobody is the unmoved specification wearing this level's name")
    end
    arms = HE.analyse(df, spec.config, spec.chem_col, ctx.exclusions[spec.bg_mode];
                      ca_class_taxids = ctx.ca_class_taxids, moved_in = moved)
    for (n, arm) in ((arms.n1, "FOREGROUND"), (arms.n2, "BACKGROUND"))
        n > 0 ||
            error("the specification '$(spec.name)' has an empty $arm arm: it has nothing " *
                  "to compare, so every interval it writes is undefined and every category " *
                  "it reports is a category it did not measure")
    end
    return arms
end

"""
    check_movein_growth(moved_arms, base_arms) -> NamedTuple

The move-in level must ADD taxa to the foreground, and add exactly the annotated taxa it
moved. The arm builder already refuses a pool that named frame taxa and moved none of them;
this check refuses the level that moved taxa which never reached the arm, which is the same
silent failure one step later and is invisible in every row the run writes.
"""
function check_movein_growth(moved_arms, base_arms)
    moved_arms.n_moved > 0 ||
        error("the move-in level moved no taxon at all, so its rows are the rows of the " *
              "specification it was supposed to vary")
    moved_arms.n1 == base_arms.n1 + moved_arms.n_moved_usable ||
        error("the move-in level's foreground is $(moved_arms.n1) where the unmoved " *
              "foreground of $(base_arms.n1) plus its $(moved_arms.n_moved_usable) annotated " *
              "moved taxa is $(base_arms.n1 + moved_arms.n_moved_usable): the moved taxa did " *
              "not reach the arm")
    return (n1_before = base_arms.n1, n1_after = moved_arms.n1, n2_before = base_arms.n2,
            n2_after = moved_arms.n2, moved = moved_arms.n_moved,
            moved_annotated = moved_arms.n_moved_usable)
end

"""
    check_both_arms_filtered(sub_arms, whole_arms, name) -> NamedTuple

The thin-evidence subset must take units out of BOTH arms.

A subset applied to the foreground alone is not a unit check: it is a different foreground
compared against the whole background, and it moves every estimate in one direction for a
reason that has nothing to do with how well a call is supported. The failure is invisible in
every row it writes, because a smaller foreground against the same background is exactly what
a real effect looks like — so it is checked here rather than read off the output.
"""
function check_both_arms_filtered(sub_arms, whole_arms, name::AbstractString)
    (sub_arms.n1 < whole_arms.n1 && sub_arms.n2 < whole_arms.n2) ||
        error("the thin-evidence subset of '$name' left one of the arms unchanged " *
              "(n1 $(whole_arms.n1) -> $(sub_arms.n1), n2 $(whole_arms.n2) -> " *
              "$(sub_arms.n2)): a subset that filters one arm is not a unit check, it is a " *
              "different foreground")
    return (n1 = sub_arms.n1, n2 = sub_arms.n2,
            n1_whole = whole_arms.n1, n2_whole = whole_arms.n2)
end

# The rows the flag reads: the Ln arm against the background, on the two units whose
# denominator is a taxon, under the registered specification or a one-at-a-time variant of
# it. Everything else in the table is reported and decides nothing.
feeds_flag(spec::Spec, contrast::AbstractString, unit::AbstractString) =
    (spec.role == ROLE_PRIMARY || spec.role == ROLE_VARIANT) &&
    contrast == HE.LN_VS_BG && unit in HE.PRIMARY_UNITS

"""
    grid_row(spec, r, unit) -> NamedTuple

One row of the long-format table. `lo`/`hi` are the STANDARDISED bounds, which are what every
rule in this pipeline reads; the crude bounds travel beside them under their own names,
because both estimators are reported together everywhere else and dropping one here would
make this table the only place the pair comes apart.

The unit carries the specification's suffix, which is empty for a whole-arm row and names the
subset for a thin-evidence one, so the currency of `n1` and `n2` is written into the key of
every row rather than inferred from which specification wrote it.

`n_zero_stratum_cells` and `min_bg_stratum_cell` travel too: a row whose totals are healthy
can still rest on a stratum the background barely populates, and the interval alone does not
show that.
"""
function grid_row(spec::Spec, r, unit::AbstractString)
    named = unit * spec.unit_suffix
    return (variant = spec.name, axis = spec.axis, level = spec.level, role = spec.role,
            contrast = r.contrast, unit = named, category = String(r.habitat),
            a = float(r.a), b = float(r.b), n1 = r.n1, n2 = r.n2,
            crude = r.log2_rr, crude_lo = r.ci_lo, crude_hi = r.ci_hi,
            std = r.log2_rr_std, lo = r.ci_lo_std, hi = r.ci_hi_std,
            n_strata = r.n_strata, n_zero_stratum_cells = r.n_zero_stratum_cells,
            min_bg_stratum_cell = float(r.min_bg_stratum_cell),
            tested = r.tested, q_bh = r.q_bh,
            feeds_flag = feeds_flag(spec, r.contrast, named))
end

"""
    chem_route(spec) -> Symbol
    chem_weight_rule(spec) -> Union{Symbol,Nothing}

Where a specification's chemistry rows come from, and which rule the weight column carries.

`:native` means the chemistry rows are the ones the arm builder computes from the call
table's own chemistry column — that is the registered chemistry rule, and it is what the
enrichment step publishes, so the registered specification is read off it rather than off a
recomputation of it. `:rule` means the chemistry axis is being VARIED, and the rows come
from the weight column instead.

The weight column is nonetheless built under BOTH routes, because it is the only route that
reaches the chemistry axis on the genus and per-sequence units: the arm builder emits those
two units for the weight column alone. Under the native route the weight column carries the
registered rule, so its taxon rows duplicate the native ones and are dropped — and the run
checks that they really are duplicates before believing the same seam under a varied rule.
The earlier specification's exclusive winner column is not a rule of the chemistry axis at
all, so it gets no weight column and no derived chemistry units.
"""
chem_route(spec::Spec) = spec.chem_rule === nothing ? :native : :rule
chem_weight_rule(spec::Spec) =
    spec.chem_rule !== nothing ? spec.chem_rule :
    spec.chem_col === :chemistry_flags ? :any_term : nothing

"""
    spec_rows(spec, type_arms, chem_arms) -> Vector{NamedTuple}

Every row of one specification, tagged with what produced it.

A specification whose chemistry axis is being varied carries every chemistry row in the
second run, read under the flag names; the first run's own chemistry rows are then the
REGISTERED rule and are dropped, or the table would carry the registered rule under a
variant's name. A specification on the registered chemistry rule keeps the native taxon rows
and takes only the genus and per-sequence chemistry series from the second run. A
specification that produces no row at all is an error: a variant that could not be computed
must not reach the grid as a variant that happened to have nothing to say.
"""
function spec_rows(spec::Spec, type_arms, chem_arms)
    out = NamedTuple[]
    route = chem_route(spec)
    if spec.emit !== :chemistry_only
        for r in type_arms.rows
            route === :native || !is_chemistry_unit(r.unit) || continue
            push!(out, grid_row(spec, r, r.unit))
        end
    end
    if chem_arms !== nothing
        for r in chem_arms.rows
            is_chemistry_unit(r.unit) && continue      # the call table's own column, not this rule
            route === :native && r.unit == "taxa" && continue   # the native row, already kept
            push!(out, grid_row(spec, r, chemistry_unit(r.unit)))
        end
    end
    isempty(out) &&
        error("the specification '$(spec.name)' produced no row; a variant that computed " *
              "nothing is not a variant that agrees")
    return out
end

"""
    run_grid(ctx, specs; share_threshold) -> (rows, runs, diagnostics)

Every specification, computed once. Runs are cached on what actually decides them — the
configuration, the weight column, the background rule and the subset — so the
specifications that share arms share one call of the arm builder instead of drifting apart
in the last bits of a repeated computation.
"""
function run_grid(ctx, specs::AbstractVector{Spec};
                  share_threshold::Real = ET.SHARE_THRESHOLD)
    cache = Dict{Any,Any}()
    arms(spec, weights, chem_col) = get!(cache,
        (spec.config, weights, spec.bg_mode, spec.subset, chem_col)) do
            s = Spec(spec.name, spec.axis, spec.level, spec.config, spec.tie_rule, chem_col,
                     spec.chem_rule, spec.bg_mode, spec.subset, spec.role, spec.emit,
                     spec.unit_suffix, spec.note)
            build_arms(ctx, s; weights = weights, share_threshold = share_threshold)
        end
    rows = NamedTuple[]
    runs = Dict{Tuple{String,String},Any}()
    for spec in specs
        tw = spec.tie_rule === nothing ? nothing : (:type, spec.tie_rule)
        type_arms = arms(spec, tw, spec.chem_col)
        wr = chem_weight_rule(spec)
        chem_arms = wr === nothing ? nothing : arms(spec, (:chemistry, wr), spec.chem_col)
        append!(rows, spec_rows(spec, type_arms, chem_arms))
        runs[run_key(spec)] = (type_arms = type_arms, chem_arms = chem_arms)
    end
    return (rows = rows, runs = runs)
end

"""
    chemistry_route_check(native_rows, derived_rows) -> (compared, worst)

The registered chemistry rule taken twice: once through the call table's own chemistry
column, which is what the enrichment step publishes, and once through the weight column that
every chemistry VARIANT has to cross. On the registered rule the two are the same statement,
so every taxon-unit row must agree cell for cell. If they do not, the seam does not carry the
rule it is given, and every chemistry variant computed through it is a different rule
reported under the name of the one that was asked for.

A row on one side and not the other is a failure too: the two routes must measure the same
categories on the same contrasts, or the comparison is over whichever rows happened to match.
"""
function chemistry_route_check(native_rows, derived_rows)
    unit = chemistry_unit("taxa")
    native = Dict{Tuple{String,String},Any}()
    for r in native_rows
        r.unit == unit && (native[(r.contrast, r.category)] = r)
    end
    compared = 0; worst = 0.0
    for r in derived_rows
        r.unit == unit || continue
        n = get(native, (r.contrast, r.category), nothing)
        n === nothing &&
            error("the chemistry weight column produced a $(r.contrast) row for " *
                  "'$(r.category)' that the call table's own chemistry column does not have")
        compared += 1
        for c in (:a, :b, :std, :lo, :hi)
            x = fnum(getproperty(n, c)); y = fnum(getproperty(r, c))
            (isnan(x) && isnan(y)) && continue
            worst = max(worst, abs(x - y))
        end
    end
    compared == length(native) ||
        error("the chemistry weight column produced $compared rows where the call table's " *
              "own chemistry column has $(length(native))")
    return (compared = compared, worst = worst)
end

"""
    assert_grid_complete(rows, expected)

Two properties the table has to have before anything is read off it: every expected variant
is present, and no (variant, contrast, unit, category) appears twice. A missing variant read
as an agreeing one is the failure the robustness rule exists to prevent, and a duplicated
row would let one variant vote twice.
"""
function assert_grid_complete(rows, expected)
    present = Set(String[r.variant for r in rows])
    missing_names = sort([String(v) for v in expected if !(String(v) in present)])
    isempty(missing_names) ||
        error("the grid is missing " * string(length(missing_names)) * " of the variants it " *
              "set out to run: " * join(missing_names, ", "))
    seen = Set{NTuple{4,String}}()
    for r in rows
        k = (r.variant, r.contrast, r.unit, r.category)
        k in seen && error("the grid holds two rows for " * join(k, " | "))
        push!(seen, k)
    end
    return length(seen)
end

# THE CATEGORY VOCABULARY, READ OFF THE LIBRARY THAT ASSIGNS IT. Not a list written here: a
# list written here is a second declaration of the same thing, and the copy is what drifts.
# `COARSE`'s keys are every fine habitat type the library can assign — the ontology branches
# plus the fallback pillar — and `OBT_CHEMISTRY` names every condition flag. A chemistry
# COMBINATION is deliberately not in the second set: the arm builder splits combinations into
# their flags before it counts them, so a combination reaching a written row means the column
# being read is not the column that was asked for.
const REGISTERED_TYPE_CATEGORIES =
    Set(String[String(c) for c in
               vcat(String[first(x) for x in ET.OBT_TYPE], String.(collect(keys(ET.COARSE))))])
const REGISTERED_CHEMISTRY_CATEGORIES =
    Set(String[String(first(x)) for x in ET.OBT_CHEMISTRY])

"""
    assert_categories_are_registered(rows; type_categories, chemistry_categories) -> Int

Every category this grid publishes must be one the library declares, on the kind of unit it
belongs to.

A category is not a free-text label: it is the name of an ontology branch or of a condition
flag, and the set of them is fixed by the library that assigns them. A specification that
reads the wrong column produces plausible rows under names that are not categories at all —
a combination of two flags, or a label from a column that holds something else — and every
one of them travels into the flag table, the written table and the panel as a category the
grid "measured". Refusing them here is the check that the rows describe the vocabulary the
rest of the pipeline is written against.

The vocabularies are arguments so that a caller working on fixtures can state its own; the
defaults are the library's and are what the run uses.
"""
function assert_categories_are_registered(rows;
        type_categories = REGISTERED_TYPE_CATEGORIES,
        chemistry_categories = REGISTERED_CHEMISTRY_CATEGORIES)
    bad = String[]
    seen = Set{Tuple{String,String}}()
    for r in rows
        chem = is_chemistry_unit(unit_stem(String(r.unit)))
        cat = String(r.category)
        cat in (chem ? chemistry_categories : type_categories) && continue
        k = (chem ? "chemistry" : "type", cat)
        k in seen && continue
        push!(seen, k)
        push!(bad, "'$cat' on $(r.variant) | $(r.unit)")
    end
    isempty(bad) ||
        error("the grid published " * string(length(bad)) * " name(s) that the library does " *
              "not declare as a category: " * join(sort(bad), "; ") * ". A category is an " *
              "ontology branch or a condition flag, not whatever the column being read " *
              "happened to contain")
    return length(rows)
end

# THE CONTRASTS THIS GRID MEASURES, read off the library that registers them plus the one
# disclosed comparator beside them. Every specification in this grid is run through the same
# arm builder, so every specification measures all of them; a contrast that reaches no row is
# therefore a contrast this run did not compute, not one it had nothing to say about.
const DECLARED_CONTRASTS = (HE.REGISTERED_CONTRASTS..., HE.CA_VS_BG)
# The unit series the arm builder emits, before the chemistry twin and the subset suffix are
# put on them. The taxon unit is taken from the library's own list so a rename there moves
# with it; the other two are the per-sequence and symmetric-genus series this file already
# names wherever it reports them.
const DECLARED_UNIT_SERIES = (FLAG_UNITS_TYPE..., "sequences", "genera")

"""
    assert_the_published_shape_is_the_declared_one(rows, specs; contrasts, series,
                                                  type_categories, chemistry_categories) -> Int

The table about to be published, checked against what the specifications say they measure:
every declared contrast, every declared unit series, and every declared category is there.

WHY A SHAPE CHECK AND NOT A ROW COUNT. Every other guard in this file asks whether a row that
IS there is right. None of them asks whether a row that should be there is missing, and the
two failures do not look alike: a table short of one whole contrast carries no wrong number
anywhere, meets every anchor target, reproduces every reference row and exits 0. What a reader
downstream sees is a complete-looking file, because "complete" is not a property any column of
it states. A missing input is the ordinary way this happens — the arm builder leaves a
contrast out rather than emitting a panel of zeros, which is the right decision for the arm
and the wrong one for the table — and a stale input that names nobody is the same silence with
a file still on disk.

THREE DIMENSIONS, BECAUSE THE TABLE HAS THREE. A contrast is a comparison; a unit is what the
denominator counts; a category is the name the row is about. Losing any one of them shrinks
the output without moving a number, so each is declared — from the library, never from this
run's own rows — and each is named on its way out. The per-specification check is the one that
bites hardest: a contrast absent from EVERY specification is an input that did not arrive,
while a contrast absent from one is a specification that did not run.

The vocabularies are arguments so a caller working on fixtures can state its own; the defaults
are the library's and are what the run uses.
"""
function assert_the_published_shape_is_the_declared_one(rows, specs;
        contrasts = DECLARED_CONTRASTS, series = DECLARED_UNIT_SERIES,
        type_categories = REGISTERED_TYPE_CATEGORIES,
        chemistry_categories = REGISTERED_CHEMISTRY_CATEGORIES)
    declared_contrasts = String[String(c) for c in contrasts]
    published_contrasts = Set(String[String(r.contrast) for r in rows])
    absent = [c for c in declared_contrasts if !(c in published_contrasts)]
    isempty(absent) ||
        error("the grid published no row at all on " * string(length(absent)) *
              " of the " * string(length(declared_contrasts)) * " contrast(s) it measures: " *
              join(absent, ", ") * ". A contrast with no rows is not a contrast with nothing " *
              "to report — the arms of every specification are built by one call of the arm " *
              "builder, so a contrast that reached no row is one whose evidence never " *
              "arrived. The table would go out short of that whole comparison with every " *
              "number in it correct, every anchor met and nothing in the file saying which " *
              "comparison is missing")

    # The two derived unit kinds, and the subset suffix, are what the arm builder's series
    # become on the way into this table; the stems are compared so a suffixed row counts as
    # its own series rather than as a series of its own.
    declared_stems = String[]
    for u in series
        push!(declared_stems, String(u))
        is_chemistry_unit(String(u)) || push!(declared_stems, chemistry_unit(String(u)))
    end
    published_stems = Set(String[unit_stem(String(r.unit)) for r in rows])
    missing_stems = sort(unique([u for u in declared_stems if !(u in published_stems)]))
    isempty(missing_stems) ||
        error("the grid published no row on " * join(missing_stems, ", ") *
              ", which the specifications declare they measure. A unit is what the " *
              "denominator counts, so a unit that vanished takes a whole series out of the " *
              "written table and out of the panel while every series that remains still " *
              "agrees with itself")

    # PER SPECIFICATION, because that is the difference between an input that did not arrive
    # and a specification that did not run. The subset twin carries the same name as its
    # parent and is a different run, so the key is the name and the suffix — the same key the
    # run cache is built on.
    by_spec = Dict{Tuple{String,String},Set{String}}()
    for r in rows
        u = String(r.unit)
        k = (String(r.variant), endswith(u, MIN_DISTINCT_SUFFIX) ? MIN_DISTINCT_SUFFIX : "")
        push!(get!(by_spec, k, Set{String}()), String(r.contrast))
    end
    short = String[]
    for s in specs
        got = get(by_spec, run_key(s), Set{String}())
        gone = [c for c in declared_contrasts if !(c in got)]
        isempty(gone) ||
            push!(short, "$(s.name)$(s.unit_suffix): " * join(gone, ", "))
    end
    isempty(short) ||
        error("the grid measures " * string(length(declared_contrasts)) * " contrast(s) and " *
              string(length(short)) * " specification(s) did not publish all of them: " *
              join(sort(short), "; ") * ". Every specification here is one call of the same " *
              "arm builder over the same contrasts, so a specification short of one is a " *
              "specification whose rows describe less than its name in the written table " *
              "says it does")

    # AND THE CATEGORY VOCABULARY IN THE OTHER DIRECTION. `assert_categories_are_registered`
    # refuses a published name the library does not declare; this refuses a declared name the
    # grid published nowhere, which is the same defect read from the other end and the one a
    # shrinking table produces.
    published_cats = Dict("type" => Set{String}(), "chemistry" => Set{String}())
    for r in rows
        kind = is_chemistry_unit(unit_stem(String(r.unit))) ? "chemistry" : "type"
        push!(published_cats[kind], String(r.category))
    end
    unmeasured = String[]
    for (kind, declared) in (("type", type_categories), ("chemistry", chemistry_categories))
        for c in sort(collect(declared))
            c in published_cats[kind] || push!(unmeasured, "$kind: $c")
        end
    end
    isempty(unmeasured) ||
        error("the library declares " * string(length(unmeasured)) * " category(ies) that " *
              "this grid published on no row: " * join(unmeasured, "; ") * ". A category " *
              "that left the output left it silently: the written table has no cell for it, " *
              "the panel has no line for it, and nothing in either says the name was ever " *
              "expected")
    return length(by_spec)
end

"""
    check_genus_series(rows, specifications) -> (specifications, units, cells, missing_cells)

The symmetric genus series, checked for completeness across every graded specification.

WHY A CHECK AND NOT A COUNT. A downstream decision rests on this series holding across the
whole grid, so what has to be true is a property of every CELL — each graded specification
measures each of the genus categories the registered specification measures — and not of the
number of specifications that have some genus row somewhere. A count cannot tell a complete
series from one that lost a category under one specification and gained a row from a line
outside the grid; a blank cell in the written table is all either leaves behind. The cells are
counted here so the output records what was checked rather than that something was.
"""
function check_genus_series(rows, specifications)
    names = String[String(s) for s in specifications]
    present = Set{Tuple{String,String,String}}()
    for r in rows
        r.contrast == HE.LN_VS_BG || continue
        push!(present, (String(r.variant), String(r.unit), String(r.category)))
    end
    units = String[]; checked = 0; missing_cells = String[]
    for unit in ("genera", chemistry_unit("genera"))
        cats = sort(unique([String(r.category) for r in rows
                            if r.contrast == HE.LN_VS_BG && r.unit == unit &&
                               r.variant == "primary"]))
        isempty(cats) && continue
        push!(units, unit)
        for v in names, c in cats
            checked += 1
            (v, unit, c) in present || push!(missing_cells, "$v | $unit | $c")
        end
    end
    return (specifications = length(names), units = units, cells = checked,
            missing_cells = sort(missing_cells))
end

# ---------------------------------------------------------------------------
# The robustness flag.
# ---------------------------------------------------------------------------

"""
    classify_category(primary_lo, primary_hi, variants) -> (flag, breakers, meaning,
                                                            primary_excludes_zero)

The registered robustness rule, with no relaxation: a category carries a sign only if the
PRIMARY interval excludes zero and EVERY one-at-a-time variant excludes it on the same side.
It is `fragile` when the primary excludes zero and a variant does not, and `null` whenever the
primary interval includes zero. There is no majority tier and no class of nearly-robust rows:
a rule that counts votes would let the size of the grid decide the colour of a row.

`fragile` IS A PROPERTY OF A CLAIM, AND A NULL PRIMARY MAKES NO CLAIM. The recorded rule
gives the fragile shape one meaning — a result that the registered specification states and
some variant does not reproduce — and gives a primary that includes zero a different
treatment entirely: the claim is dropped and the row stays neutral. A category whose primary
includes zero while some variant excludes it is therefore `null` and not `fragile`, however
interesting those variants are. Filing it under `fragile` would put a word shaped like "state
it tentatively" on a row that has nothing to state, and would make the flag column a third
neutral tier for anything downstream that keys colour off it.

THOSE VARIANTS ARE STILL REPORTED, because they are the informative part of a null row: a
level under which the category WOULD be distinguishable from the background is worth seeing,
and `breaking_variants_are` says which of the two senses the names carry. What changes is the
word in the column a colour is taken from, not what the table knows.

`variants` is `name => (lo, hi)`, in whatever order; a variant whose interval could not be
computed, including a category that has no row at all under that variant, is a BREAKER and
never a variant that agrees. `breakers` names them, so the flag is auditable rather than
asserted, and `meaning` says what being named there means — the three branches name three
different things, and a column of names whose sense depends on the flag beside it is a column
that will be read wrong.

`primary_excludes_zero` travels beside the flag so that anything downstream keys on the fact
rather than on the word: a reader of the word alone cannot tell a neutral row that some
variant would have coloured from one the whole grid agrees is neutral, and a reader of the
boolean never has to.
"""
function classify_category(primary_lo, primary_hi, variants)
    label(p) = String(strip(string(first(p))))
    computed = String[label(p) for p in variants
                      if OR.interval_excludes_zero(last(p)...) !== :undecidable]
    status = OR.interval_excludes_zero(primary_lo, primary_hi)
    # NO PRIMARY INTERVAL, SO NOTHING TO BE ROBUST AGAINST. This is how a category that only
    # SOME level of the grid creates arrives here: the registered specification has no row for
    # it, so it has no sign to reproduce or break. Naming the variants that do measure it says
    # which levels create it; naming all fourteen as "breakers" would report the twelve that
    # never saw the category as having disagreed about it.
    status === :undecidable &&
        return (flag = "undecidable", breakers = sort(computed),
                meaning = "the registered specification has no interval for this category; " *
                          "these variants do",
                primary_excludes_zero = false)
    if status === :excludes_zero
        depleted = fnum(primary_hi) < 0
        brk = String[]
        for p in variants
            lo, hi = last(p)
            same_side = OR.interval_excludes_zero(lo, hi) === :excludes_zero &&
                        (depleted ? fnum(hi) < 0 : fnum(lo) > 0)
            same_side || push!(brk, label(p))
        end
        isempty(brk) && return (flag = depleted ? "-" : "+", breakers = String[],
                                meaning = "", primary_excludes_zero = true)
        return (flag = "fragile", breakers = sort(brk),
                meaning = "variants that do not reproduce the primary sign",
                primary_excludes_zero = true)
    end
    # THE PRIMARY INCLUDES ZERO, so the row is neutral whatever the variants say. The names
    # are kept because they are the informative part of a neutral row; the word is not
    # changed by them, because the recorded rule reserves the fragile shape for a claim the
    # registered specification actually makes.
    exc = String[label(p) for p in variants
                 if OR.interval_excludes_zero(last(p)...) === :excludes_zero]
    isempty(exc) && return (flag = "null", breakers = String[], meaning = "",
                            primary_excludes_zero = false)
    return (flag = "null", breakers = sort(exc),
            meaning = "variants that exclude zero where the primary does not",
            primary_excludes_zero = false)
end

"""
    robustness(rows, expected_variants; movement) -> Vector{NamedTuple}

The flag of every category the grid measured, on the rows the flag reads.

A category is keyed on its unit as well as its name, because the taxon unit and the chemistry
unit are two different measurements of two different things and a category belongs to exactly
one of them. The variant intervals are also handed to the recorded tier function, which is
what makes the depletion tier decidable at all: the tier is defined by this grid, and with no
grid to consult the function can only answer that it does not know.

THE TIER IS A VERDICT ON A DEPLETION, so it is emitted only where there is one to judge. Its
three levels read "state it", "state it tentatively" and "drop the claim" about a category the
registered specification puts BELOW the background; on a category the same specification puts
above it, the same three words describe a claim nobody made, and the worst of them sits in the
column next to an enriched row looking like a verdict against it. A positive primary estimate
therefore gets an empty tier rather than the tier's own name for "nothing here".

`n_variants_that_moved_the_estimate` says how many of the expected variants actually produced
a different number for this cell. It is the honest denominator of "every variant agrees": the
rest are cells the variant's axis cannot reach, and counting them as agreement would report
perturbations that were never applied.
"""
function robustness(rows, expected_variants;
                    movement = estimate_movement(rows, expected_variants))
    flagged = [r for r in rows if r.feeds_flag]
    keys_seen = sort(unique([(r.unit, r.category) for r in flagged]))
    byname = Dict{Tuple{String,String,String},Any}()
    for r in flagged
        byname[(r.variant, r.unit, r.category)] = r
    end
    names = String[String(v) for v in expected_variants]
    out = NamedTuple[]
    for (unit, cat) in keys_seen
        p = get(byname, ("primary", unit, cat), nothing)
        plo = p === nothing ? NaN : fnum(p.lo)
        phi = p === nothing ? NaN : fnum(p.hi)
        pairs = Pair{String,Tuple{Float64,Float64}}[]
        measured = 0
        for v in names
            r = get(byname, (v, unit, cat), nothing)
            if r === nothing
                # No row means no hit on either side under that variant, so there is no
                # interval — which is a variant that did not reproduce anything, not one
                # that agreed. It is carried as an uncomputed interval on purpose.
                push!(pairs, v => (NaN, NaN))
            else
                measured += 1
                push!(pairs, v => (fnum(r.lo), fnum(r.hi)))
            end
        end
        c = classify_category(plo, phi, pairs)
        t = OR.depletion_tier(plo, phi, pairs; expected_variants = names)
        pstd = p === nothing ? NaN : fnum(p.std)
        enriched = !isnan(pstd) && pstd > 0
        push!(out, (unit = unit, category = cat, robust_sign = c.flag,
                    primary_excludes_zero = c.primary_excludes_zero,
                    n_breaking_variants = length(c.breakers),
                    breaking_variants = join(c.breakers, ";"),
                    breaking_variants_are = c.meaning,
                    depletion_tier = enriched ? "" : String(t.tier),
                    depletion_tier_breakers = enriched ? "" : join(t.breakers, ";"),
                    primary_std = pstd,
                    primary_lo = plo, primary_hi = phi,
                    primary_a = p === nothing ? NaN : fnum(p.a),
                    variants_expected = length(names), variants_with_a_row = measured,
                    n_variants_that_moved_the_estimate =
                        get(movement.per_cell, (unit, cat), 0)))
    end
    return out
end

# ---------------------------------------------------------------------------
# The anchors.
# ---------------------------------------------------------------------------

const TARGET_COLUMNS = (:kind, :anchor, :unit, :category, :quantity, :value, :tol, :note,
                        :digest)
# What a target may be taken on. Every one of them is a column of the long-format table, so a
# target names a quantity this run actually wrote rather than one it could compute on request.
const TARGET_QUANTITIES = (:a, :b, :n1, :n2, :crude, :std, :lo, :hi)

# The five kinds of row a targets file carries. Two of them are measurements of this dataset
# and three are not: a `variant` row names one level of the registered set, an `expect` row
# records a count the file makes about itself, and a `source` row names the file the targets
# were transcribed FROM together with its digest. Names, counts and digests could be written in
# this script — but the whole reason the set is declared elsewhere is that a declaration living
# beside what it describes moves with it, and the targets file is the one input the run already
# refuses to start without.
const TARGET_KINDS = ("target", "exception", "variant", "expect", "source")
# What an `expect` row must count.
#
#   targets         how many `target` rows the file holds.
#   variants        how many levels the declared set holds. NOT redundant with counting the
#                   `variant` rows: comparing the enumerated set with the declared one name for
#                   name is satisfied by ANY two lists that agree, so an edit that shortens the
#                   script's own tuples and strikes the matching names from this file is
#                   self-consistent and publishes a smaller grid under the registered set's
#                   name. A size is the one statement that edit cannot also keep consistent,
#                   because it is not a list and shortening a list does not change it.
#   reference_rows  how many rows the row-by-row comparison against the earlier run is supposed
#                   to cover. Without it a reference file regenerated at a different scope
#                   reports "N of N rows agree" for whatever N it happens to carry, and the
#                   count reads better than the honest run's rather than different from it.
const EXPECTED_COUNTS = (:targets, :variants, :reference_rows)

"""
    read_anchor_targets(path) -> (targets, exceptions, variants, expected, sources)

Everything the targets file declares: the transcribed anchor targets, the rows of an earlier
run that are known NOT to reproduce with the reason each does not, the NAMES of the registered
one-at-a-time levels, the counts the file makes about itself, and the file each target was
transcribed FROM with that file's digest.

No target is typed into this script: a target is a number measured on this dataset, and the
script carries none. The targets for the deposited data ship as a separate file,
`dat/criteria/anchor_targets.csv`, passed by `--anchor-targets`. A run given no targets file
checks nothing and says so; a run given one must meet every target in it. An EXCEPTION is a
row of the earlier run that a later, deliberate change of the analysis moved — it is excluded
from the row-by-row comparison BY NAME, with its reason travelling beside it into the output,
which is the opposite of widening a tolerance until the difference disappears.

WHY THE SET, THE COUNTS AND THE SOURCE ARE HERE. None of them is a quantity measured on
anything: one is a list of rule names, three are counts of this file's own contents, and the
last is a file name with a digest. All of them nevertheless have to reach the RUN. The set,
because everything inside the script that could contradict a deleted level is enumerated from
the same tuples the level was deleted from. The counts, because the strength of the whole
wiring check is otherwise an input: a file that lost rows is a run that checked fewer targets
and reported success in the same words as a complete one. The source, because a target read off
this grid's own output is byte-for-byte what an honestly transcribed one looks like — the two
implementations agree far below any band here — so the only thing that can distinguish them is
a machine-checkable statement of WHICH bytes the numbers were taken from, and a note in prose
is not one. A file that declares any of them is refused, so a thinned, re-transcribed or
re-sourced file is a failure rather than a quieter success.

NO TARGET MAY BE TAKEN TWICE. Two rows on one (anchor, unit, category, quantity) are one check
counted twice: they meet or miss together. A file can then carry its declared number of targets
while pinning fewer CELLS than the count suggests — and a duplicated row that cannot fail can
stand in, one for one, for a row that can, leaving every reported number identical to an honest
run's.
"""
function read_anchor_targets(path::AbstractString)
    isfile(path) || error("--anchor-targets $path does not exist")
    tbl = CSV.read(path, DataFrame)
    for c in TARGET_COLUMNS
        hascol(tbl, c) ||
            error("--anchor-targets $path needs the columns " *
                  join(String.(TARGET_COLUMNS), ", ") * " and has no '$c'")
    end
    targets = NamedTuple[]; exceptions = NamedTuple[]
    variants = String[]; expected = Dict{Symbol,Int}(); sources = NamedTuple[]
    taken = Set{NTuple{4,String}}()
    for i in 1:nrow(tbl)
        kind = sv(tbl.kind[i])
        row = (anchor = sv(tbl.anchor[i]), unit = sv(tbl.unit[i]),
               category = sv(tbl.category[i]), quantity = sv(tbl.quantity[i]),
               value = fnum(tbl.value[i]), tol = fnum(tbl.tol[i]), note = sv(tbl.note[i]),
               digest = sv(tbl.digest[i]))
        if kind == "target"
            Symbol(row.quantity) in TARGET_QUANTITIES ||
                error("a target may be taken on " * join(String.(TARGET_QUANTITIES), ", ") *
                      "; '$(row.quantity)' is not one of them")
            (isnan(row.value) || isnan(row.tol)) &&
                error("the target $(row.anchor) | $(row.unit) | $(row.category) | " *
                      "$(row.quantity) has no value or no tolerance")
            cell = (row.anchor, row.unit, row.category, row.quantity)
            cell in taken &&
                error("the target $(row.anchor) | $(row.unit) | $(row.category) | " *
                      "$(row.quantity) is taken twice. Two rows on one cell are one check " *
                      "counted twice — they meet or miss together — so the file would carry " *
                      "its declared number of targets while pinning fewer cells than that " *
                      "number says, and a copy of a row that cannot fail can stand in for " *
                      "the row that can while every reported count stays identical to an " *
                      "honest run's")
            push!(taken, cell)
            push!(targets, row)
        elseif kind == "exception"
            isempty(row.note) &&
                error("the exception $(row.anchor) | $(row.unit) | $(row.category) carries " *
                      "no reason; an unexplained exception is a widened tolerance with a " *
                      "different name")
            push!(exceptions, row)
        elseif kind == "variant"
            # The level's NAME travels in the anchor column, which is the column that names a
            # specification in every other row of this file.
            isempty(row.anchor) &&
                error("a 'variant' row declares one level of the registered set by name, in " *
                      "the anchor column, and this one names none")
            row.anchor in variants &&
                error("'$(row.anchor)' is declared twice; a file that names the set twice " *
                      "cannot say how large the set is")
            push!(variants, row.anchor)
        elseif kind == "expect"
            Symbol(row.quantity) in EXPECTED_COUNTS ||
                error("an 'expect' row counts one of " * join(String.(EXPECTED_COUNTS), ", ") *
                      "; '$(row.quantity)' is not one of them")
            haskey(expected, Symbol(row.quantity)) &&
                error("'$(row.quantity)' is expected twice, so the file states its own size " *
                      "in two places")
            (isnan(row.value) || row.value < 0 || row.value != round(row.value)) &&
                error("the expected number of $(row.quantity) is a whole count and this one " *
                      "is '$(sv(tbl.value[i]))'")
            expected[Symbol(row.quantity)] = Int(row.value)
        elseif kind == "source"
            # The file the targets were transcribed FROM, by name and by digest. The name is a
            # bare file name resolved beside this file: a path would let the declaration point
            # at something the person reading the declaration cannot see.
            isempty(row.anchor) &&
                error("a 'source' row names the file the targets were transcribed from, in " *
                      "the anchor column, and this one names none")
            basename(row.anchor) == row.anchor ||
                error("the source '$(row.anchor)' is a path; a source is a bare file name " *
                      "read from beside the targets file, so that what is declared and what " *
                      "a reader of the declaration can open are the same thing")
            is_sha256(row.digest) ||
                error("the source '$(row.anchor)' carries '$(row.digest)' where a sha256 " *
                      "digest belongs; a named source with no digest is a prose note, and a " *
                      "prose note is what a re-transcription does not update")
            any(s -> s.file == row.anchor, sources) &&
                error("the source '$(row.anchor)' is declared twice, so the file states two " *
                      "digests for one set of bytes")
            isempty(row.note) &&
                error("the source '$(row.anchor)' carries no note saying what it is and how " *
                      "it was produced; a digest pins WHICH bytes and says nothing about " *
                      "whether they were computed independently of this grid")
            push!(sources, (file = row.anchor, sha256 = lowercase(row.digest),
                            note = row.note))
        else
            error("an anchor row is one of " * join(TARGET_KINDS, ", ") * ", not a '$kind'")
        end
    end
    isempty(variants) &&
        error("$path declares no 'variant' rows, so this run would enumerate its own set and " *
              "compare it against nothing. The registered levels are named here because the " *
              "script enumerates them from its own tuples: a level deleted there is deleted " *
              "from the vocabulary that would otherwise have missed it")
    isempty(sources) &&
        error("$path declares no 'source' row, so nothing in it says which bytes its targets " *
              "were transcribed from. A target read off this grid's own output is identical " *
              "to an honestly transcribed one — the independent implementation and this grid " *
              "agree far below every band in this file — so the transcription's origin is " *
              "not something a reader, or a later run, can recover from the numbers")
    for q in EXPECTED_COUNTS
        haskey(expected, q) ||
            error("$path carries no 'expect' row for '$q'. Without it the file states one " *
                  "less fact about itself than a run needs: a file that lost rows, a " *
                  "declared set that lost a level, or a reference regenerated at another " *
                  "scope would each report success in the same words as the complete case")
    end
    expected[:targets] == length(targets) ||
        error("$path declares $(expected[:targets]) target(s) and carries " *
              "$(length(targets)): the file has been thinned, extended or merged since the " *
              "count beside it was written, and how much of the wiring is checked is not " *
              "something a run may discover from the file that was handed to it")
    # THE SIZE OF THE SET, AGAINST THE SET ITSELF. `assert_the_set_is_the_declared_one`
    # compares two lists and is satisfied by any two that agree — including two that were
    # shortened together, which is the one edit that still changes a published colour. The
    # size is not a list: shortening the names does not change it, so the three statements
    # cannot be made consistent by the same stroke.
    expected[:variants] == length(variants) ||
        error("$path declares a registered set of $(expected[:variants]) level(s) and names " *
              "$(length(variants)) of them. A set compared name for name agrees with any " *
              "list it was shortened alongside; the size is the statement that does not move " *
              "when a name is struck from both, and it disagrees here")
    return (targets = targets, exceptions = exceptions, variants = variants,
            expected = expected, sources = sources)
end

"""
    check_target_sources(targets_path, sources) -> Vector{NamedTuple}

Each declared source, digested where it is present beside the targets file.

WHAT THIS CAN AND CANNOT SAY. It cannot prove a target was transcribed honestly: a number is
a number, and the implementation this grid is checked against agrees with it far below any
band in the file. What it can do is pin WHICH bytes the file claims to have been transcribed
from, so that the claim is a fact a later run re-checks rather than a sentence in a `note`
column that no program reads and that a re-transcription does not update. A source regenerated
from a broken run's output — the shape that replaces both halves of this gate at once — no
longer has the declared digest, and the run stops.

An absent source is reported and not refused. The targets this repository ships travel with
their sources beside them, so a run on them digests both; a targets file a reader writes may
name an implementation that is not distributed with it, and a run of this script must not
require that file, or the directory it lives in, present. What is refused is a source that IS
there and is not the bytes the file names.
"""
function check_target_sources(targets_path::AbstractString, sources)
    dir = dirname(abspath(targets_path))
    out = NamedTuple[]
    for s in sources
        p = joinpath(dir, s.file)
        if !isfile(p)
            push!(out, (file = s.file, status = "not beside the targets file",
                        declared = s.sha256, found = "", ok = true))
        else
            got = bytes2hex(sha256(read(p)))
            push!(out, (file = s.file, status = got == s.sha256 ? "matches" : "differs",
                        declared = s.sha256, found = got, ok = got == s.sha256))
        end
    end
    return out
end

"""
    check_anchor_targets(rows, targets) -> (checked, met, failures)

Every transcribed target, taken against the row this run wrote for it. A target whose row is
absent fails: a target that cannot be found has not been reproduced, and reporting it as
untested would let a wiring error read as a missing input.
"""
function check_anchor_targets(rows, targets)
    byk = Dict{Tuple{String,String,String},Any}()
    for r in rows
        r.contrast == HE.LN_VS_BG || continue
        byk[(r.variant, r.unit, r.category)] = r
    end
    failures = String[]
    met = 0
    for t in targets
        r = get(byk, (t.anchor, t.unit, t.category), nothing)
        label = "$(t.anchor) | $(t.unit) | $(t.category) | $(t.quantity)"
        if r === nothing
            push!(failures, "$label: this run wrote no such row")
            continue
        end
        got = fnum(getproperty(r, Symbol(t.quantity)))
        if isnan(got) || abs(got - t.value) > t.tol
            push!(failures, @sprintf("%s: got %.5f, target %.5f (tolerance %.5f)",
                                     label, got, t.value, t.tol))
        else
            met += 1
        end
    end
    return (checked = length(targets), met = met, failures = failures)
end

"""
    levels_named_in_a_breaker_list(flags) -> Vector{String}
    unanchored_published_levels(flags, targets) -> Vector{String}

Which levels this run PUBLISHES by name, and which of those no target pins.

WHY THE REQUIREMENT IS DERIVED AND NOT LISTED. A target on the registered specification says
the primary's own row is the right number; it says nothing about the fourteen intervals that
decide, together with it, the word beside a category. The levels that reach print are the ones
a breaker list names — those names are the published claim "this rule, and not that one, is
what the result does not survive" — and a wiring error confined to one such level's inputs
rewrites that list while every count in the run agrees with itself. So the levels that must be
anchored are read off what this run actually published rather than written down: a level that
starts or stops breaking a category moves the requirement with it, and a list written here
would go stale in exactly the run where it mattered.

A target on ANY cell of the level is enough. The failure this catches — a configuration
re-exported from a sibling, a scan column half populated, two configurations that exchanged
labels upstream — moves a level's arms wholesale, not one category of them.
"""
levels_named_in_a_breaker_list(flags) =
    sort(unique(String[String(strip(v)) for f in flags
                       for v in split(f.breaking_variants, ';') if !isempty(strip(v))]))

function unanchored_published_levels(flags, targets)
    anchored = Set(String[String(t.anchor) for t in targets])
    return [v for v in levels_named_in_a_breaker_list(flags) if !(v in anchored)]
end

"""
    reference_shape(r) -> NamedTuple

One grid row in the shape the arm builder emits, so the enrichment step's own row-by-row
comparator can be CALLED on it rather than restated here. The grid names the category column
`category` because a unit and a contrast travel beside it; the comparator reads the arm
builder's `habitat`. Renaming three fields is the whole adapter, and it is written once here
instead of a second comparator being written anywhere.
"""
reference_shape(r) = (contrast = r.contrast, unit = r.unit, habitat = r.category,
                      a = r.a, b = r.b, n1 = r.n1, n2 = r.n2,
                      log2_rr = r.crude, ci_lo = r.crude_lo, ci_hi = r.crude_hi,
                      log2_rr_std = r.std, ci_lo_std = r.lo, ci_hi_std = r.hi)

"""
    scoped_reference(ref, exceptions; units) -> (table, excluded)

The earlier run's rows, scoped for the row-by-row comparison the enrichment step's own
`reference_row_diff` then makes: restricted to the units whose definition did not change,
and with the named exceptions removed.

Two things are scoped out, and both are named rather than absorbed. A unit that was
REDEFINED since the earlier run is not the same measurement under the same name, so
comparing it would report a change made on purpose as a failure to reproduce. An EXCEPTION
is a single row a later, deliberate change of the analysis moved; it leaves the comparison
by name, carrying its reason into the output, which is the opposite of widening a tolerance
until the difference disappears. The removed rows come back by name so the output can carry
them instead of losing them.

BOTH KINDS COME BACK COUNTED, and that matters more for the units than for the exception. A
reader of the provenance alone sees "N of N rows agree" and concludes the earlier run is
reproduced; if the rows a redefinition removed are not reported beside that, the count is
taken over a scope nobody can see, and the redefinition — which is the largest deliberate
difference between the two runs — is the one thing the output does not mention.
"""
function scoped_reference(ref::AbstractDataFrame, exceptions; units)
    ex = Set{Tuple{String,String}}((sv(e.unit), sv(e.category)) for e in exceptions)
    keep = BitVector(undef, nrow(ref))
    excluded = String[]; out_of_scope = String[]
    for i in 1:nrow(ref)
        u = sv(ref.unit[i]); c = sv(ref.habitat[i])
        inunit = u in units
        isex = (u, c) in ex
        keep[i] = inunit && !isex
        inunit && isex && push!(excluded, "$u | $c")
        inunit || push!(out_of_scope, "$u | $c")
    end
    return (table = ref[keep, :], excluded = sort(excluded),
            out_of_scope = sort(out_of_scope),
            out_of_scope_units = sort(unique([String(split(x, " | ")[1])
                                              for x in out_of_scope])))
end

# THE SCOPE OF THE ROW-BY-ROW COMPARISON, as a declaration rather than as an argument written
# at the call site. The genus unit is deliberately not in it, and the reason travels with the
# list so that a later widening has to delete a stated reason rather than edit a tuple.
const REFERENCE_SCOPE_UNITS = ("taxa", "sequences", "taxa (chemistry)")
const REFERENCE_OUT_OF_SCOPE_REASON =
    "these rows sit on a unit the two tables do not define the same way: a reference file " *
    "whose genus unit counts the two arms in different currencies is not comparable, row by " *
    "row, with this grid, which collapses BOTH arms to genus, so the same name does not " *
    "denote the same measurement. The earlier specification's value on the symmetric unit is " *
    "reported in this grid's own genus series instead of being compared against a number it " *
    "does not correspond to"

"""
    require_census(path) -> (census, sha256)

The census table, refused when it is not there or when it names nobody.

WHY THIS INPUT IS AN ERROR AND NOT A DISCLOSURE. The call table and the background scan both
stop a run when they are missing, and so must this one: it is the input that decides how much
of the grid gets written. The disclosed Ca-class comparator arm is defined by the census, and
an arm with no membership is left OUT of the arm builder rather than emitted as a panel of
zeros — the right decision for an arm, and the wrong one for a table, because what comes out
the other end is every other contrast, complete and correct, in a file with no column saying
that a comparison is absent. A path that was never updated, or a census not yet regenerated,
then costs a published table an entire contrast at exit 0 with no warning.

A FILE THAT NAMES NOBODY IS THE SAME FAILURE WITH BYTES ON DISK. The reader hands back
`nothing` when it cannot find the file and an empty membership when the file holds no
calcium-class taxon, and the arm builder cannot tell those apart either, so both are refused
here.

AND THE DIGEST IS TAKEN BECAUSE THE FILE IS THERE. A digest written empty for an absent file
reads in the provenance as a digest that was taken and says nothing: one blank field where
the difference between two runs is a missing contrast.
"""
function require_census(path::AbstractString)
    isfile(path) ||
        error("missing $path: the census table decides the membership of the disclosed " *
              "Ca-class comparator arm, and an arm with no members is left out of the arm " *
              "builder rather than emitted as a panel of zeros. A run without it publishes " *
              "a table short of that entire contrast, with every number that remains " *
              "correct, every anchor target met and nothing in the file saying which " *
              "comparison is missing — which is worse than no table, because every " *
              "downstream reader sees a complete-looking one")
    digest = bytes2hex(sha256(read(path)))
    is_sha256(digest) ||
        error("the census table at $path did not digest to a sha256; a provenance row that " *
              "carries no digest under the name of one is a claim about which bytes were " *
              "read that cannot be re-checked")
    census = HE.read_census_classes(path)
    (census.ca_class_taxids === nothing || isempty(census.ca_class_taxids)) &&
        error("the census table at $path gives the disclosed Ca-class comparator arm no " *
              "members, so every specification in this grid would publish one contrast " *
              "fewer than it measures. A stale or half-written census is the ordinary way " *
              "this happens, and it costs a published table a whole comparison without " *
              "moving a single number in the rows that remain")
    return (census = census, sha256 = digest)
end

"""
    reference_scope_ok(compared, expected) -> Bool

Whether the row-by-row comparison covered the scope the targets file declares it covers.

"N of N rows agree" is the same sentence at any N. A reference regenerated at another scope —
more rows, fewer rows, another vintage — therefore reports a number that reads BETTER than the
honest run's rather than different from it, and the only thing that can tell the two apart is
a count written down before the reference was on hand. `nothing` is the declined case, where
no comparison was asked for at all.
"""
reference_scope_ok(compared, expected) = expected === nothing || compared == expected

# ---------------------------------------------------------------------------
# Outputs.
# ---------------------------------------------------------------------------

fmt_est(v, lo, hi) = (isnan(v) || isnan(lo) || isnan(hi)) ? "n/a" :
                     @sprintf("%+.2f [%+.2f, %+.2f]", v, lo, hi)

# A cell whose variant's axis cannot reach the unit it is in. It is NOT an agreement: the
# number there would be the registered specification's own, printed again under another
# name, and a column of those read as a column of independent checks is exactly how a grid
# of fourteen becomes a claim of fourteen perturbations that were never applied.
#
# The mark is placed from the MEASURED movement, and it says something the measurement alone
# cannot support — that the rule family does not touch the unit at all. What licenses it is
# `assert_axes_reach_the_units_they_declare`, which the run applies before this table is
# written: measured reach and declared reach are equal there or the run stops, so a variant
# that lost its wiring to one unit cannot arrive here wearing this mark.
const NOT_REACHED = "n.r."

"""
    cell_conventions() -> String

What an unnumbered cell means. Three different facts share the unnumbered state and only one
of them is a missing measurement, so each is named where the table is read rather than left
to be inferred from a cell that happens to look empty.
"""
cell_conventions() = string(
    "A BLANK cell is a category with no hit on either side under that specification, so no ",
    "row was written.", "\n",
    "`n/a` is a row that exists but whose interval is undefined: one arm scored nothing, so ",
    "there is no", "\n", "ratio to bound.", "\n",
    "`", NOT_REACHED, "` is a cell this variant's axis cannot reach. That rule family does ",
    "not touch this unit at all,", "\n",
    "so the estimate there is the registered specification's own and is not an independent ",
    "check of it.", "\n")

"""
    columns_note() -> String

What the long-format table's columns mean, and which of them are additions to the column list
the plan registered for it.

The registered list has no room for a unit or a contrast, and without those two a genus row
and a taxon row of the same specification and category are one key with two values — so the
table would silently keep whichever was written last. It also carries one interval where this
pipeline reports two estimators everywhere else. Both additions are stated here, beside the
table they change, rather than left for a reader to infer from a column they did not expect.
"""
columns_note() = """
## the columns, and what was added to the registered list

`variant, axis, level, role, contrast, unit, category, a, b, n1, n2, crude, crude_lo,
crude_hi, std, lo, hi, n_strata, n_zero_stratum_cells, min_bg_stratum_cell, tested, q_bh,
feeds_flag`.

* `lo` and `hi` are the STANDARDISED bounds — the ones every rule in this pipeline reads.
  `crude_lo` and `crude_hi` are the crude ones, reported beside them and deciding nothing.
* `unit` and `contrast` are ADDITIONS to the registered column list. Without them a genus row
  and a taxon row of the same specification and category share one key, and so do the Ln arm
  and a comparator arm; the table would carry whichever was written last and look complete.
* `unit` carries the subset suffix, so the currency of `n1` and `n2` is part of the row key.
* `feeds_flag` is true only for the rows the robustness flag reads: the registered
  specification or a one-at-a-time variant of it, the Ln arm against the background, on the
  two units whose denominator is a taxon. Every other row is reported and decides nothing.
* `role` says what a row is for: the registered specification, a one-at-a-time variant, a
  unit check, an earlier specification reproduced, an earlier sensitivity rule reproduced, or
  a background rule disclosed beside the grid and deliberately outside it.

"""

"""
    sensitivity_table(rows, variant_order) -> String

The sensitivity table: one line per category, one column per specification, each cell the
standardised estimate on its interval. Only the rows the flag reads are tabulated; the unit
checks and the comparator arms are in the long-format table beside it, where their unit and
their contrast travel with them and cannot be read as the same measurement.
"""
function sensitivity_table(rows, variant_order; movement = nothing)
    flagged = [r for r in rows if r.feeds_flag]
    isempty(flagged) && return "no row of this run feeds the robustness flag.\n"
    order = String[String(v) for v in variant_order]
    variants = String[v for v in order if v != "primary"]
    mv = movement === nothing ? estimate_movement(rows, variants) : movement
    byk = Dict{Tuple{String,String,String},Any}()
    for r in flagged; byk[(r.variant, r.unit, r.category)] = r; end
    io = IOBuffer()
    println(io, "# Habitat sensitivity grid — standardised estimate [95% interval]")
    println(io)
    println(io, "Rows: the Ln arm against the background, on the units whose denominator is ",
                "a taxon.")
    println(io, "Columns: the registered specification and the one-at-a-time variants of it.")
    println(io)
    print(io, cell_conventions())
    println(io)
    for unit in HE.PRIMARY_UNITS
        cats = sort(unique([r.category for r in flagged if r.unit == unit]))
        isempty(cats) && continue
        reached = [v for v in variants if get(mv.per_unit, (v, unit), 0) > 0]
        println(io, "## unit: ", unit)
        println(io)
        # The honest denominator of "every variant agrees", printed where the columns are
        # read: the remaining columns are the primary's own numbers under another name, and a
        # reader who counts columns instead would count perturbations that never happened.
        println(io, "Variants whose axis reaches this unit: ", length(reached), " of ",
                    length(variants), "; the other ", length(variants) - length(reached),
                    " are marked `", NOT_REACHED, "`.")
        println(io)
        println(io, "| category | ", join(order, " | "), " |")
        println(io, "|---|", repeat("---|", length(order)))
        for c in cats
            cells = String[]
            for v in order
                r = get(byk, (v, unit, c), nothing)
                push!(cells,
                      v != "primary" && get(mv.per_unit, (v, unit), 0) == 0 ? NOT_REACHED :
                      r === nothing ? "" :
                      fmt_est(fnum(r.std), fnum(r.lo), fnum(r.hi)))
            end
            println(io, "| ", c, " | ", join(cells, " | "), " |")
        end
        println(io)
    end
    return String(take!(io))
end

"""
    unit_check_table(rows, variant_order) -> String

The unit checks, across the grid and marked as what they are. The genus series is required
across every specification, because a downstream decision rests on exactly that series;
the per-sequence series and the thin-evidence subset sit beside it. None of them feeds the
flag, and the table says so in its own text rather than leaving it to be known.
"""
function unit_check_table(rows, variant_order)
    order = String[String(v) for v in variant_order]
    byk = Dict{Tuple{String,String,String},Any}()
    for r in rows
        r.contrast == HE.LN_VS_BG || continue
        byk[(r.variant, r.unit, r.category)] = r
    end
    io = IOBuffer()
    println(io, "# Unit checks — reported across the grid, and NOT counted toward any flag")
    println(io)
    println(io, "A unit check changes what is being counted, not which rule counted it, so ",
                "it cannot")
    println(io, "reproduce or break a rule variant. These series are descriptive.")
    println(io)
    print(io, cell_conventions())
    println(io)
    units = vcat([u for u in ("genera", "sequences", "genera (chemistry)",
                              "sequences (chemistry)")],
                 [u * MIN_DISTINCT_SUFFIX
                  for u in ("taxa", "taxa (chemistry)", "genera", "sequences",
                            "genera (chemistry)", "sequences (chemistry)")])
    for unit in units
        cats = sort(unique([r.category for r in rows
                            if r.unit == unit && r.contrast == HE.LN_VS_BG]))
        isempty(cats) && continue
        println(io, "## unit: ", unit)
        println(io)
        println(io, "| category | ", join(order, " | "), " |")
        println(io, "|---|", repeat("---|", length(order)))
        for c in cats
            cells = [(r = get(byk, (v, unit, c), nothing);
                      r === nothing ? "" : fmt_est(fnum(r.std), fnum(r.lo), fnum(r.hi)))
                     for v in order]
            println(io, "| ", c, " | ", join(cells, " | "), " |")
        end
        println(io)
    end
    return String(take!(io))
end

"""
    panel_rows(rows, flags) -> Vector{NamedTuple}

What the rendered panel draws: every category the registered specification measured on the
two taxon units, including the ones too thin for the main figure, with the span of the
one-at-a-time variants beside the primary interval and the two unit series beside both.

The span is the range of the variant POINT estimates. It answers "how far does the rule
choice move the estimate", which is a different question from the interval's "how far does
the sample move it", and drawing them as two different marks is what keeps the two apart.
"""
function panel_rows(rows, flags)
    flagmap = Dict((f.unit, f.category) => f for f in flags)
    byk = Dict{Tuple{String,String,String},Any}()
    for r in rows
        r.contrast == HE.LN_VS_BG || continue
        byk[(r.variant, r.unit, r.category)] = r
    end
    out = NamedTuple[]
    for unit in HE.PRIMARY_UNITS
        cats = sort(unique([r.category for r in rows
                            if r.contrast == HE.LN_VS_BG && r.unit == unit &&
                               r.variant == "primary"]))
        for c in cats
            p = byk[("primary", unit, c)]
            f = get(flagmap, (unit, c), nothing)
            ests = Float64[]
            for r in rows
                r.feeds_flag && r.unit == unit && r.category == c || continue
                v = fnum(r.std); isnan(v) || push!(ests, v)
            end
            gunit = unit == "taxa" ? "genera" : chemistry_unit("genera")
            squnit = unit == "taxa" ? "sequences" : chemistry_unit("sequences")
            gen = get(byk, ("primary", gunit, c), nothing)
            seq = get(byk, ("primary", squnit, c), nothing)
            push!(out, (unit = unit, category = c, a = fnum(p.a),
                        est = fnum(p.std), lo = fnum(p.lo), hi = fnum(p.hi),
                        span_lo = isempty(ests) ? NaN : minimum(ests),
                        span_hi = isempty(ests) ? NaN : maximum(ests),
                        genus = gen === nothing ? NaN : fnum(gen.std),
                        sequence = seq === nothing ? NaN : fnum(seq.std),
                        flag = f === nothing ? "undecidable" : f.robust_sign))
        end
    end
    return out
end

"""
    panel_note(flags, panel) -> String

The two things a reader of the panel cannot get from the panel: what the number beside each
category name is, and which flagged categories are not on it.

A category that only some variant measures has no registered-specification row, so the panel
— which draws that specification — has nowhere to put it. It is in the flag table and in the
written table all the same, and a reader holding the three beside each other would otherwise
find rows on the page that are missing from the picture with nothing saying why. The names
come from the flag table rather than from a list written here, so a panel that drops a
different category says so without being edited.
"""
function panel_note(flags, panel)
    drawn = Set{Tuple{String,String}}((p.unit, p.category) for p in panel)
    off = [f for f in flags if !((f.unit, f.category) in drawn)]
    lines = String["The number after each category name is its weighted foreground count."]
    isempty(off) ||
        push!(lines,
              "Not drawn: " *
              join([isempty(f.breaking_variants) ? f.category :
                    f.category * " (" * replace(f.breaking_variants, ";" => ", ") * ")"
                    for f in off], ", ") *
              " — measured under those variants only, with no row under the registered " *
              "specification; both the flag table and the sensitivity table carry them.")
    return join([wrap_text(l, NOTE_WRAP) for l in lines], "\n")
end

# How wide a note line may be, in characters, before it runs off a panel of this width. The
# renderer does not wrap, and a line that overruns is not clipped visibly — it simply ends
# mid-sentence at the edge of the page, which is worse than no sentence at all.
const NOTE_WRAP = 105

"""
    wrap_text(s, width) -> String

Break a sentence into lines of at most `width` characters on word boundaries, never inside a
word: a note the renderer will not wrap has to arrive already wrapped.
"""
function wrap_text(s::AbstractString, width::Integer)
    out = String[]; line = ""
    for w in split(s, ' '; keepempty = false)
        if isempty(line)
            line = String(w)
        elseif length(line) + 1 + length(w) <= width
            line = line * " " * w
        else
            push!(out, line); line = String(w)
        end
    end
    isempty(line) || push!(out, line)
    return join(out, "\n")
end

# Okabe-Ito, the house pair: colour encodes direction and nothing else, and only a row the
# registered rule flags as robust carries any. Everything else is neutral — a row that a
# variant breaks and a row that was never distinguishable from the background are both rows
# the figure may not colour, and inventing a third mark for the first would reintroduce the
# tier the rule does not have.
const COL_ENRICHED = "#D55E00"
const COL_DEPLETED = "#0072B2"
const COL_NEUTRAL  = "#7F7F7F"
const COL_INK      = "#2E2E2E"
const COL_RULE     = "#BFBEB6"

panel_colour(flag::AbstractString) = flag == "+" ? COL_ENRICHED :
                                     flag == "-" ? COL_DEPLETED : COL_NEUTRAL

# Geometry, in the house units: a full-width panel at a 180 mm text width, sized so each
# row has room for four marks that must not touch.
const WIDTH_MM  = 180
const PT_PER_MM = 72 / 25.4
const ROW_MM    = 7.0
const FS        = 7

# The two unit series, as an offset below the row's own estimate and a marker size. They sit
# at DIFFERENT offsets because the two estimates agree closely on some rows, and two marks at
# one offset then print as a single glyph — the figure would show one measurement where it
# made two. The gap between the two offsets must therefore exceed the two markers' radii
# together, at this row height, or they touch and are separable only by looking for the
# corner of one behind the other. They are constants so that the separation is a property the
# suite can check rather than two numbers inside the renderer.
const GENUS_DY    = -0.20
const GENUS_MS    = 4.0
const SEQUENCE_DY = -0.41
const SEQUENCE_MS = 3.4
# The distance between the two marks' centres, and the distance at which their edges meet,
# both in points on the page.
unit_mark_separation_pt() = abs(GENUS_DY - SEQUENCE_DY) * ROW_MM * PT_PER_MM
unit_marks_touch_pt() = (GENUS_MS + SEQUENCE_MS) / 2

# The base name of the rendered panel. Descriptive of what the panel shows: a file name in a
# public script is read by people who have no access to whatever numbering a document gives
# its figures, and a name that is only a position in that numbering tells them nothing and
# goes stale the moment the order changes.
const PANEL_BASE = "habitat_sensitivity_panel"

"""
    cairomakie() -> Module

The renderer, loaded at the moment it is needed and not before. It costs about a minute to
load and nothing but the rendered panel uses it, so a suite that includes this file for its
rules — which is every test of everything above — never pays for it.
"""
cairomakie() = Core.eval(@__MODULE__, :(import CairoMakie; CairoMakie))

"""
    draw_panel(mk, panel, path_base; note) -> Vector{String}

The sensitivity panel: every category on one axis, the chemistry block below a rule, the
primary estimate on its interval, the span of the rule variants as a bar behind it, and the
genus and per-sequence estimates as two lighter marks.

EVERY MARK IS NAMED, including the interval. A whisker that no legend entry claims is read as
whatever the reader expects a whisker to be — and here it is the one mark the figure and the
table share, so it is the one that must not be guessed at. The two unit marks sit at
different offsets rather than one, because at the point where two of these estimates agree
they would otherwise print as a single glyph and the figure would show one measurement where
it made two.
"""
function draw_panel(mk, panel, path_base::AbstractString; note::AbstractString = "")
    hab  = sort([p for p in panel if p.unit == "taxa"], by = p -> p.est)
    chem = sort([p for p in panel if p.unit != "taxa"], by = p -> p.est)
    gap = 1.0
    ychem = collect(1.0:length(chem))
    yhab  = collect((length(chem) + gap + 1):(length(chem) + gap + length(hab)))
    yrule = length(chem) + gap / 2 + 0.5
    ymax  = (isempty(yhab) ? yrule : last(yhab)) + 0.7
    draw = vcat([(p = p, y = ychem[i]) for (i, p) in enumerate(chem)],
                [(p = p, y = yhab[i])  for (i, p) in enumerate(hab)])
    isempty(draw) && error("the panel has no row to draw")

    figw = WIDTH_MM * PT_PER_MM
    # the note is drawn under the legend, so the canvas grows by the lines it holds
    figh = (length(draw) * ROW_MM + 34 + 3.2 * count(==('\n'), note) + (isempty(note) ? 0 : 3.2)) *
           PT_PER_MM
    labels = [@sprintf("%s (%s)", x.p.category,
                       isinteger(x.p.a) ? string(Int(x.p.a)) : @sprintf("%.1f", x.p.a))
              for x in draw]

    mk.activate!(type = "pdf")
    fig = mk.Figure(size = (figw, figh), fontsize = FS)
    ax = mk.Axis(fig[1, 1], xlabel = "log₂ enrichment against the background",
                 title = "Sensitivity of every category to the aggregation and background rules",
                 titlealign = :left, titlesize = FS + 1, titlefont = :bold,
                 xlabelsize = FS, xticklabelsize = FS - 0.5, yticklabelsize = FS - 0.5,
                 yticks = ([x.y for x in draw], labels),
                 ygridvisible = false, xgridvisible = true,
                 xgridcolor = (COL_RULE, 0.5), xgridwidth = 0.4,
                 rightspinevisible = false, topspinevisible = false, spinewidth = 0.6)
    mk.vlines!(ax, [0.0]; color = (COL_INK, 0.5), linewidth = 0.8)
    length(chem) == 0 ||
        mk.hlines!(ax, [yrule]; color = (COL_RULE, 0.9), linewidth = 0.5, linestyle = :dot)

    lows = Float64[]; highs = Float64[]
    for x in draw
        p = x.p; y = x.y; c = panel_colour(p.flag)
        # the span of the rule variants, drawn behind everything as a wide pale bar
        if !isnan(p.span_lo) && !isnan(p.span_hi)
            mk.lines!(ax, [p.span_lo, p.span_hi], [y + 0.22, y + 0.22];
                      color = (COL_NEUTRAL, 0.55), linewidth = 3.2)
            push!(lows, p.span_lo); push!(highs, p.span_hi)
        end
        # the two unit series: a different visual class, so neither can be read as an estimate
        # of the same kind as the primary one
        isnan(p.genus) ||
            mk.scatter!(ax, [p.genus], [y + GENUS_DY]; marker = :diamond, color = :white,
                        strokecolor = (COL_INK, 0.6), strokewidth = 0.6,
                        markersize = GENUS_MS)
        isnan(p.sequence) ||
            mk.scatter!(ax, [p.sequence], [y + SEQUENCE_DY]; marker = :rect, color = :white,
                        strokecolor = (COL_NEUTRAL, 0.8), strokewidth = 0.6,
                        markersize = SEQUENCE_MS)
        if !isnan(p.est)
            if !isnan(p.lo) && !isnan(p.hi)
                mk.errorbars!(ax, [p.est], [y], [p.est - p.lo], [p.hi - p.est];
                              direction = :x, whiskerwidth = 3, linewidth = 1.1, color = c)
                push!(lows, p.lo); push!(highs, p.hi)
            end
            mk.scatter!(ax, [p.est], [y]; color = c, markersize = 4.6)
            push!(lows, p.est); push!(highs, p.est)
        end
        isnan(p.genus)    || (push!(lows, p.genus);    push!(highs, p.genus))
        isnan(p.sequence) || (push!(lows, p.sequence); push!(highs, p.sequence))
    end
    mk.xlims!(ax, minimum(lows) - 0.15, maximum(highs) + 0.15)
    mk.ylims!(ax, 0.3, ymax)
    length(chem) == 0 || mk.text!(ax, 0.02, (yrule - 0.35) / ymax; text = "CHEMISTRY",
                                  space = :relative, align = (:left, :top),
                                  fontsize = FS - 1.5, color = (COL_NEUTRAL, 0.95),
                                  font = :bold)
    length(hab) == 0 || mk.text!(ax, 0.02, 0.985; text = "HABITAT", space = :relative,
                                 align = (:left, :top), fontsize = FS - 1.5,
                                 color = (COL_NEUTRAL, 0.95), font = :bold)

    mk.Legend(fig[2, 1],
        [mk.MarkerElement(color = COL_ENRICHED, marker = :circle, markersize = 4.6),
         mk.MarkerElement(color = COL_DEPLETED, marker = :circle, markersize = 4.6),
         mk.MarkerElement(color = COL_NEUTRAL, marker = :circle, markersize = 4.6),
         mk.LineElement(color = (COL_INK, 0.8), linewidth = 1.1),
         mk.LineElement(color = (COL_NEUTRAL, 0.55), linewidth = 3.2),
         mk.MarkerElement(color = :white, marker = :diamond, markersize = GENUS_MS,
                          strokecolor = (COL_INK, 0.6), strokewidth = 0.6),
         mk.MarkerElement(color = :white, marker = :rect, markersize = SEQUENCE_MS,
                          strokecolor = (COL_NEUTRAL, 0.8), strokewidth = 0.6)],
        ["robust, enriched", "robust, depleted", "not robust",
         "95% interval on the estimate", "span of the rule variants",
         "genus unit", "per-sequence unit"],
        orientation = :horizontal, framevisible = false, labelsize = FS - 0.5,
        nbanks = 2, patchsize = (7.0f0, 7.0f0), colgap = 8, rowgap = 2)
    mk.rowsize!(fig.layout, 2, mk.Relative(0.10))
    isempty(note) ||
        mk.Label(fig[3, 1], note; fontsize = FS - 1.5, color = (COL_INK, 0.85),
                 halign = :left, justification = :left, lineheight = 1.15,
                 tellwidth = false)

    written = String[]
    for (ext, kw) in (("pdf", (pt_per_unit = 1,)), ("png", (px_per_unit = 4,)))
        p = path_base * "." * ext
        mk.save(p, fig; kw...)
        push!(written, p)
    end
    return written
end

# ---------------------------------------------------------------------------
# The options, and the gate between a finished run and the figure directory.
# ---------------------------------------------------------------------------

const OPTION_SPEC = Dict("--calls" => :value, "--scan" => :value, "--census" => :value,
                         "--figdir" => :value, "--anchor-reference" => :value,
                         "--anchor-targets" => :value, "--no-anchors" => :flag,
                         "--no-figure" => :flag, "--help" => :flag, "-h" => :flag)

"""
    parse_run_options(args; retired = ANCHOR_LEVELS_RETIRED) -> NamedTuple

The command line, as the settings of one run. Defined here rather than inside the program
guard because the rule it carries — that the wiring check is refused by default — is the
kind of rule that is deleted by accident and noticed by nobody, and a rule inside a program
guard cannot be tested at all.

`retired` is whether the anchor levels are left out of the grid. While they are, the grid
holds no all-legacy rows and `--anchor-reference` has nothing to be compared with, so only
`--anchor-targets` is required; with them evaluated, both are.
"""
function parse_run_options(args; retired::Bool = ANCHOR_LEVELS_RETIRED)
    cli = HE.parse_options(args, OPTION_SPEC)
    opt(flag, default) = HE.option_string(cli, flag, default)
    anchor_targets = opt("--anchor-targets", "")
    anchor_reference = opt("--anchor-reference", "")
    no_anchors = haskey(cli, "--no-anchors")
    # A WIRING CHECK THAT IS OPT-IN IS A WIRING CHECK NOBODY RUNS. Given no targets this
    # script would write every table and exit 0, and the flag table would be byte-identical
    # to a fully checked run's; the only difference would be two rows of a provenance file.
    # So the default is refused, and skipping the check is a decision that has to be typed
    # and is then recorded as one.
    #
    # BOTH HALVES WHILE THE ANCHOR LEVELS ARE EVALUATED, for the same reason and not only the
    # first. The row-by-row comparison against the earlier run is the half a run can drop by
    # losing a flag in a driver script or a shell quoting change: the run then records
    # `anchor_reference_status,not supplied` and publishes everything, and the only trace of
    # the missing check is a provenance row nobody diffs. A gate one flag removes without a
    # failure is not a gate.
    #
    # ONE HALF WHILE THEY ARE RETIRED. The grid then holds no all-legacy rows, so the
    # comparison covers none whatever file is named, and the run records it as not
    # applicable. Requiring the option would demand a path the run never opens, and any path
    # would satisfy it: that is a requirement that checks nothing.
    missing_halves = String[]
    isempty(anchor_targets) &&
        push!(missing_halves, "--anchor-targets <path>, without which nothing checks that " *
                              "this grid reproduces the measurements it is supposed to")
    !retired && isempty(anchor_reference) &&
        push!(missing_halves, "--anchor-reference <path>, without which the all-legacy " *
                              "anchor's rows are compared against nothing and the run " *
                              "publishes with the comparison recorded as 'not supplied'")
    (isempty(missing_halves) || no_anchors) ||
        error("this run would skip " *
              (retired ? "its wiring check: " :
                         string(length(missing_halves)) * " of the two wiring checks: ") *
              join(missing_halves, "; ") *
              ". Pass the missing option(s), or --no-anchors to run without " *
              (retired ? "the check" : "either half") * " deliberately — which is then " *
              "recorded as a decision rather than left as a gap a reader of the output " *
              "cannot see" *
              (retired ? ". (The anchor levels are retired, so --anchor-reference is not " *
                         "required: it would have no rows to compare.)" : ""))
    return (calls = opt("--calls", DataInputs.input_path("habitat_calls")),
            scan = opt("--scan", DataInputs.input_path("background_scan_table")),
            census = opt("--census", DataInputs.input_path("census_habitat")),
            figdir = opt("--figdir", get(ENV, "LNFINDER_FIGDIR", joinpath(REPO, "figs"))),
            anchor_reference = anchor_reference,
            anchor_targets = anchor_targets,
            no_anchors = no_anchors,
            figure = !haskey(cli, "--no-figure"))
end

"""
    record_reference_not_applicable!(rec, anchor_reference) -> true

What a run records about the row-by-row comparison while the anchor levels are retired, through
`rec` (`rec!` in a run): the comparison is not applicable, it covered no rows, and a reference
named on the command line anyway is recorded by its file name and NOT opened. Its rows would
decide nothing, and a digest recorded beside a comparison that did not happen would read as a
comparison that did. Returns the comparison's verdict, which cannot fail on zero rows.
"""
function record_reference_not_applicable!(rec, anchor_reference::AbstractString)
    rec("anchor_reference_status", "not applicable: the anchor levels are retired")
    rec("anchor_reference_rows_compared", 0)
    isempty(anchor_reference) ||
        rec("anchor_reference_named_not_opened", basename(anchor_reference))
    return true
end

# ONE staging directory, under a fixed name beside the outputs. Fixed, because a fresh name
# per run accumulates one directory per run on any filesystem that will not let go of an
# empty one; and beside the outputs, because a move within one filesystem is atomic and a
# move across two is a copy that can half-finish. Whatever a previous attempt left in it is
# cleared, so nothing from an earlier run can be published by this one.
const STAGING_DIR = "_staging"

"""
    open_staging(figdir) -> String

The staging directory for this run, empty.
"""
function open_staging(figdir::AbstractString)
    mkpath(figdir)
    stage = joinpath(figdir, STAGING_DIR)
    mkpath(stage)
    for f in readdir(stage); rm(joinpath(stage, f); force = true, recursive = true); end
    return stage
end

"""
    stale_panel_files(figdir, figure) -> Vector{String}

The panel files already in the figure directory that THIS run will not replace.

`publish!` moves what was staged and neither removes nor marks what it did not stage, so a run
given `--no-figure` publishes its tables into a directory whose panel is the LAST run's. The
result is a complete-looking set under the names a reader — and any downstream glob — takes as
one result: numbers from run n beside the picture from run n−1, with nothing in either saying
so. The panel is the artefact a reader actually looks at, which is what makes this worse than
a missing file.

Empty when the run draws its own panel, because then the panel is replaced rather than left,
and empty when no panel of the current name is there to go stale.
"""
function stale_panel_files(figdir::AbstractString, figure::Bool)
    figure && return String[]
    isdir(figdir) || return String[]
    return sort([f for f in readdir(figdir) if startswith(f, PANEL_BASE * ".")])
end

"""
    publish!(stage, figdir, staged, gates) -> Vector{String}

Move the staged outputs into the figure directory, but only if every gate passed.

NOTHING IS PUBLISHED UNTIL THE GATES HAVE PASSED, and the two halves of that sentence are one
function so that they cannot come apart. Written straight into the figure directory, a run
whose anchors failed leaves a complete-looking set of tables there — same names, same shapes,
every number present — distinguishable from a good run only by one provenance row and an exit
status nobody reads a month later. On a failure the staging directory stays, named in the
error, so the failure is still inspectable; what it is not is published.

`gates` is a list of `(ok, message)`: every one of them is evaluated by the caller and all of
their messages are reported, because a run that failed two gates should say so once.
"""
function publish!(stage::AbstractString, figdir::AbstractString, staged, gates)
    failed = String[String(g.message) for g in gates if !g.ok]
    isempty(failed) ||
        error(join(failed, "\n") * "\nNOTHING WAS PUBLISHED: the tables are in " * stage *
              " so the failure is inspectable, and " * figdir *
              " still holds the last good run")
    written = String[]
    for name in staged
        mv(joinpath(stage, name), joinpath(figdir, name); force = true)
        push!(written, joinpath(figdir, name))
    end
    # Everything is published by this point and the staging directory is now empty. A
    # filesystem that will not remove it — a sync client still holding a handle, say — leaves
    # one empty directory under a name that is reused and cleared next run, so nothing
    # accumulates and nothing is lost. Refusing to return over that would report a complete,
    # gated, published run as a failure.
    try; rm(stage; recursive = true); catch; end
    return written
end

# The gates, by name, in the order they are built. A name per gate so that the set of them is
# something a caller — the run, or a test — can enumerate and count, rather than a literal
# vector written once inside the program guard where nothing can reach it.
const PUBLISH_GATES = ("anchor targets", "published breaker lists", "target sources",
                       "reference rows", "reference scope", "stale panel")

"""
    publish_gates(; target_failures, unanchored_levels, source_checks, reference_ok,
                  reference_rows_compared, reference_rows_expected, stale_panels,
                  figdir) -> Vector{NamedTuple}

Every statement that has to be true before this run's tables may replace the last good run's,
as `(name, ok, message)` for `publish!`.

WHY THIS IS A FUNCTION AND NOT A LITERAL AT THE CALL SITE. The predicates behind these gates
are each checked on their own, and every one of them is covered. What no such check reaches
is the six lines that HAND them to `publish!`: written inline in the program body they sit
where no test can execute them, so a gate whose entry is deleted — or whose `ok` is quietly
replaced by `true` — leaves every predicate still tested, every suite still green, and the
attack that gate exists to refuse publishing at exit 0. A gate nothing tests is a gate the
next edit removes by accident. Built here, the wiring is a value: it can be enumerated,
counted, and driven from a bad input to a refusal.

Every gate is evaluated, none short-circuits, and all the failing messages are reported at
once, because a run that failed two gates should say so in one breath.
"""
function publish_gates(; target_failures, unanchored_levels, source_checks, reference_ok,
                       reference_rows_compared, reference_rows_expected, stale_panels,
                       figdir::AbstractString)
    plural(xs, one, many) = length(xs) == 1 ? one : many
    gates = [
        (name = "anchor targets",
         ok = isempty(target_failures),
         message = "$(length(target_failures)) anchor target(s) were not " *
                   "reproduced; a grid that cannot reproduce the earlier measurements is " *
                   "not wired to the data it claims to describe"),
        (name = "published breaker lists",
         ok = isempty(unanchored_levels),
         message = "this run publishes " * join(unanchored_levels, ";") * " in a breaker " *
                   "list and no target pins " *
                   plural(unanchored_levels, "it", "them") * ". A breaker list " *
                   "names the rule a result does not survive, and that name is decided by " *
                   "the level's own arms: a wiring error confined to them rewrites the " *
                   "list while every count in this run agrees with itself and every target " *
                   "on the registered specification is met. Transcribe a target on " *
                   plural(unanchored_levels, "that level", "those levels") *
                   " from the independent implementation — not from this run's output, " *
                   "which is the failure the target exists to detect"),
        (name = "target sources",
         ok = all(s -> s.ok, source_checks),
         message = "a file the targets declare as their source is present beside them and " *
                   "is not the bytes they name: " *
                   join(["$(s.file): declared $(s.declared), found $(s.found)"
                         for s in source_checks if !s.ok], ";") *
                   ". The targets were transcribed from a version of that file this run " *
                   "cannot see, so nothing here says the numbers in it were computed " *
                   "independently of this grid"),
        (name = "reference rows",
         ok = reference_ok,
         message = "the all-legacy anchor did not reproduce the earlier run's rows; see " *
                   "anchor_reference_rows_differing in the provenance table"),
        (name = "reference scope",
         ok = reference_scope_ok(reference_rows_compared, reference_rows_expected),
         message = "the row-by-row comparison against the earlier run covered " *
                   "$reference_rows_compared rows where the targets file declares " *
                   "$reference_rows_expected: the reference on hand is not the one the " *
                   "declaration was written against, and \"N of N rows agree\" is the same " *
                   "sentence at any N"),
        (name = "stale panel",
         ok = isempty(stale_panels),
         message = "--no-figure was given and " * figdir * " already holds " *
                   join(stale_panels, ";") * ": publishing would leave this run's tables " *
                   "beside an earlier run's panel, under the names a reader and any " *
                   "downstream glob take as one result. Draw the panel, or remove the " *
                   "stale one on purpose")]
    # THE SET OF GATES IS ITSELF DECLARED. A gate dropped from the vector above is a gate
    # nothing else in this file mentions, so the names are listed once beside the function and
    # compared with what it built: removing an entry then takes two edits, and the second one
    # is a bare list whose only job is to be checked.
    built = String[String(g.name) for g in gates]
    built == String[String(n) for n in PUBLISH_GATES] ||
        error("the publish gate was built from " * join(built, ", ") * " where this file " *
              "declares " * join(String.(PUBLISH_GATES), ", ") * ". Nothing is published " *
              "unless every gate passed, so a gate that is not built is not a check that " *
              "relaxed — it is a check that is no longer made, and every count in the run " *
              "would still agree with itself")
    return gates
end

# ---------------------------------------------------------------------------
# Everything above is definitions: including this file runs nothing. Everything below reads
# the inputs, computes and writes, and runs only when this file is the program.
# ---------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__

function options_or_exit()
    try
        return parse_run_options(ARGS)
    catch e
        e isa ErrorException || rethrow()
        println(stderr, "habitat_sensitivity.jl: ", e.msg, "\n")
        print(stderr, USAGE)
        exit(1)
    end
end

const OPT = options_or_exit()

# The inputs, by name, before anything is read or printed. Only the ones still at their
# REGISTRY DEFAULT are checked: a path the reader typed is theirs, and the plain messages
# further down report it as given rather than blaming the layout. Without this, a data root
# that is not there surfaces as `missing <absolute path>` followed by a stack trace out of
# Julia's include machinery — a path off one machine, printed by a step that never said
# which named input it wanted.
OPT.calls == DataInputs.input_path("habitat_calls") &&
    DataInputs.require_or_exit("habitat_calls")
OPT.scan == DataInputs.input_path("background_scan_table") &&
    DataInputs.require_or_exit("background_scan_table")

# Every reported quantity is registered here and nowhere else, so the closing report, the
# provenance file and any downstream step read the same numbers.
const PROV = Pair{String,Any}[]
rec!(k::AbstractString, v) = (push!(PROV, String(k) => v); v)
joined(xs) = join(xs, ";")

@printf("\nhabitat_sensitivity.jl\n  calls  %s\n  scan   %s\n  figdir %s\n",
        OPT.calls, OPT.scan, OPT.figdir)

isfile(OPT.calls) || error("missing $(OPT.calls) — run scripts/habitat_reaggregate.jl first")
calls_bytes = read(OPT.calls)
calls = CSV.read(calls_bytes, DataFrame)
scan_read = HE.read_scan_table(OPT.scan)
scan = scan_read.table
scan === nothing &&
    error("missing $(OPT.scan): without the background carrier scan the registered " *
          "background rule cannot be applied, and a grid whose primary is not the " *
          "registered specification measures nothing about it")

rec!("run_at", string(Dates.now()))
rec!("julia_version", string(VERSION))
rec!("calls_table", basename(OPT.calls))
rec!("calls_sha256", bytes2hex(sha256(calls_bytes)))
rec!("calls_rows", nrow(calls))
rec!("calls_configurations", joined(sort(unique(sv.(calls.config)))))
rec!("scan_table", basename(OPT.scan))
rec!("scan_sha256", scan_read.sha256)
rec!("scan_rows", nrow(scan))
# THE THIRD INPUT, REQUIRED AND DIGESTED LIKE THE OTHER TWO. The call table and the scan
# table both stop a run when they are missing, and so does this one: it is the input that
# decides how much of the grid gets WRITTEN, so a path read whether it is there or not lets a
# dropped path publish a table short of a whole contrast at exit 0. The digest goes with it,
# because a digest written empty for an absent file reads as a digest that was taken.
census_read = require_census(OPT.census)
census = census_read.census
rec!("census_table", basename(OPT.census))
rec!("census_sha256", census_read.sha256)

# ---------------------------------------------------------------- the two rule seams
# THE SPECIFICATIONS ARE BUILT HERE, before anything reads the call table, because they are
# what says which configurations this grid reads: the seam checks below run over those and not
# over one of them. `all_specs()` reads no input and decides nothing about the data.
specs = all_specs()
type_rule_configs = registered_type_rule_configurations(specs)
seam_configs = grid_configurations(specs)
rec!("grid_configurations", joined(seam_configs))

# OVER EVERY CONFIGURATION WHOSE COMMITTED COLUMN THE REGISTERED RULE WROTE, not over the
# primary alone. Four of the other configurations ARE registered levels — the text, generic-term
# and isolation variants are their own committed weight columns read under a level's name — and
# a column re-exported from a sibling keeps that configuration's tallies, its counts and its row
# count, so every other guard in this file stays silent. The tallies are where the substitution
# shows. The earlier pipeline's configuration is deliberately not in the list: its column
# carries ITS type rule, which is the thing the anchor reproduces.
tcheck = type_reconstruction_over(calls, type_rule_configs)
tcheck.mismatches == 0 ||
    error("recomputing the registered type rule from the call table's own tallies " *
          "disagrees with the column the table carries on $(tcheck.mismatches) of " *
          "$(tcheck.rows) taxa (first: $(tcheck.example)); every type variant is built " *
          "from those tallies, and a configuration whose weight column no longer follows " *
          "from them is a level published under the name of a rule that did not produce it")
rec!("type_rule_reconstruction_configurations", joined(type_rule_configs))
rec!("type_rule_reconstruction_rows_by_configuration",
     joined(["$k=$v" for (k, v) in tcheck.per_configuration]))
rec!("type_rule_reconstruction_rows", tcheck.rows)
rec!("type_rule_reconstruction_mismatches", tcheck.mismatches)
rec!("type_rule_reconstruction_configuration_not_checked", HE.ANCHOR_CONFIG)
rec!("type_rule_reconstruction_note",
     "the earlier pipeline's configuration carries its OWN type rule in the committed weight " *
     "column, which is what the all-legacy anchor reproduces, so recomputing the registered " *
     "rule over it would disagree by design and no type variant is ever built from it")
@printf("  type rule reconstruction: %d of %d taxa in %d configurations reproduce it\n",
        tcheck.rows - tcheck.mismatches, tcheck.rows, length(type_rule_configs))

ccheck = combinations_over(calls, seam_configs)
rec!("chemistry_combination_configurations", joined(seam_configs))
rec!("chemistry_combination_rows_by_configuration",
     joined(["$k=$v" for (k, v) in ccheck.per_configuration]))
rec!("chemistry_combination_rows_with_a_flag_by_configuration",
     joined(["$k=$v" for (k, v) in ccheck.per_configuration_with_a_flag]))
rec!("chemistry_combination_rows", ccheck.rows)
rec!("chemistry_combination_rows_with_a_flag", ccheck.checked)
@printf("  chemistry combinations: %d of %d rows across %d configurations carry a flag\n",
        ccheck.checked, ccheck.rows, length(seam_configs))

# ---------------------------------------------------------------- the arms' shared inputs
pool = HE.movein_pool(; scan = scan)
pool.taxids === nothing &&
    error("the move-in level has no pool: the scan table cannot answer the rule this " *
          "repository can apply alone")
# The recorded decision for this run: the move-in level runs on the scan's own lanthanide
# class call, which is a DIFFERENT measurement from the census's call on its own alignment,
# and the pool the recorded rule names is not used. The run records which rule it ran
# under, and that the other one was not taken, as data rather than as a note somewhere else.
rec!("movein_pool_rule", pool.rule)
rec!("movein_pool_rule_source", String(pool.rule_source))
rec!("movein_pool_rule_role",
     pool.rule == HE.REGISTERED_MOVEIN_POOL ? "registered" : "descriptive")
rec!("movein_pool_rule_note", HE.movein_pool_note(pool.rule))
rec!("movein_pool_rules_available", joined(HE.MOVEIN_POOL_RULES))
rec!("movein_pool_rule_not_used", HE.REGISTERED_MOVEIN_POOL)
rec!("movein_pool_rule_not_used_note",
     "the recorded rule reads the census's own pre-deduplication hit " *
     "table, which is not an input of this repository; it was NOT supplied to this run and " *
     "the level ran on the rule the committed inputs can answer alone")
rec!("movein_pool_source", pool.source)
rec!("movein_pool_taxa", length(pool.taxids))

ctx = context(; calls = calls, scan = scan, ca_class_taxids = census.ca_class_taxids,
              movein_taxids = pool.taxids, min_distinct = MIN_DISTINCT_TERMS)
rec!("min_distinct_terms", MIN_DISTINCT_TERMS)
# THE ONE ASYMMETRY BETWEEN THE PRIMARY AND ITS VARIANTS, recorded where the numbers are
# rather than in prose beside them: a reader of these three tables alone cannot otherwise know
# that part of a variant's distance from the primary is the width of a written number.
rec!("weight_column_precision",
     "the registered specification and every level that varies neither aggregation rule run " *
     "on the call table's COMMITTED type_weights column, which is rounded for reading; a " *
     "level that varies the type rule or the chemistry rule runs on a column recomputed from " *
     "the same table's tallies at full precision, so a small part of such a variant's " *
     "distance from the primary is the rounding of the published column and not the rule")
rec!("calls_rows_in_the_thin_evidence_subset", nrow(ctx.calls_min))
for m in sort(collect(keys(ctx.exclusions)))
    rec!("background_removed_$m", length(ctx.exclusions[m]))
end

# ---------------------------------------------------------------- the declared set
# READ BEFORE THE GRID, because what this file declares has to be able to stop a run before it
# computes anything. It carries two statements the run cannot make about itself: which levels
# the registered set holds, and how many targets the file is supposed to hold.
anchor_file = isempty(OPT.anchor_targets) ? nothing : read_anchor_targets(OPT.anchor_targets)

# ---------------------------------------------------------------- the grid
S = variant_names(specs)
# THE SET IS DECLARED IN THIS FILE, so nothing outside it says what belongs in it. This is the
# check that keeps the declaration honest: a level that shares the registered specification's
# configuration, type rule, chemistry rule and background rule is that specification under
# another name, and it would give every category a variant that cannot break it while the
# size of the set, the uniqueness of its names and every anchor still passed.
rec!("variants_checked_against_the_primary", assert_variants_vary_the_primary(specs))
# AND A NAME IS NOT A LABEL. Everything published about a level is published under its name,
# so the name is checked against the computation: the axis letter is a registered axis, the
# level is a registered level of THAT axis, the field the axis varies carries it, and nothing
# else moves. A level renamed onto another axis, or a rule the recorded decision keeps
# outside the set running under a name from inside it, passes every count, every anchor and
# the comparison against the primary — and publishes a correct number under the wrong rule's
# name.
rec!("variants_whose_name_describes_what_they_run",
     assert_variant_names_describe_what_they_run(specs))
# AND NO TWO OF THEM ARE THE SAME SPECIFICATION. The comparison above is with the primary and
# with nothing else, so a pair of clones passes it while deleting a level from the set and
# leaving its column standing.
rec!("distinct_variant_specifications", assert_variants_differ_from_each_other(specs))
# AND THE SET ITSELF IS DECLARED OUTSIDE THIS FILE. Every check above compares one part of this
# script with another part of it, and the axis table reads its levels off the same tuples the
# specifications are built from: delete a level and the wiring and the vocabulary move
# together, so nothing is left to disagree. The declared set is the one statement about the
# grid that an edit to the grid cannot also make.
if anchor_file === nothing
    rec!("registered_set_declared", false)
    @warn("--no-anchors was given, so nothing outside this script says which levels the " *
          "registered set holds; a level deleted from a tuple in it would publish as the set")
else
    rec!("registered_set_declared", true)
    rec!("registered_set_declared_in", basename(OPT.anchor_targets))
    rec!("variants_matching_the_declared_set",
         assert_the_set_is_the_declared_one(S, anchor_file.variants))
end
rec!("variants_in_the_registered_set", joined(S))
rec!("n_variants_in_the_registered_set", length(S))
rec!("specifications_run", length(specs))
rec!("axes", joined(sort(unique([s.axis for s in specs]))))

@printf("  running %d specifications (%d one-at-a-time variants)\n", length(specs), length(S))
t_grid = @elapsed G = run_grid(ctx, specs)
rows = G.rows
n_unique = assert_grid_complete(rows, S)
# EVERY CATEGORY IS ONE THE LIBRARY DECLARES. A specification reading a column that holds
# something other than what it was asked for writes plausible rows under names that are not
# categories at all, and every one of them travels into the flag table under a registered
# variant's name.
rec!("grid_rows_with_a_registered_category", assert_categories_are_registered(rows))
# AND NOTHING THE SPECIFICATIONS MEASURE IS SILENTLY ABSENT. Every guard above asks whether a
# row that is here is right; this is the one that asks whether a row that should be here is
# missing. The two do not look alike from the outside: a table short of a whole contrast, unit
# or category carries no wrong number, meets every target, and goes out looking complete.
rec!("published_series_checked_against_the_declaration",
     assert_the_published_shape_is_the_declared_one(rows, specs))
rec!("declared_contrasts", joined(DECLARED_CONTRASTS))
rec!("declared_unit_series", joined(DECLARED_UNIT_SERIES))
@printf("  grid: %d rows in %.1f s\n", length(rows), t_grid)
rec!("grid_rows", length(rows))
rec!("grid_unique_keys", n_unique)

# THE CHEMISTRY ROUTE, CHECKED AGAINST THE ONE IT REPLACES, on the arms the primary already
# built: the same weight column, read as if the chemistry axis were being varied.
let check_spec = Spec("check:any_term", "check", "any_term", PRIMARY_CONFIG, nothing,
                      :chemistry_flags, :any_term, HE.PRIMARY_MODE, :all, ROLE_DISCLOSED,
                      :chemistry_only, "",
                      "the registered chemistry rule through the weight column"),
    base = G.runs[("primary", "")]
    derived = spec_rows(check_spec, base.type_arms, base.chem_arms)
    c = chemistry_route_check([r for r in rows if r.variant == "primary"], derived)
    c.worst <= CHEMISTRY_ROUTE_TOL ||
        error("the chemistry weight column and the call table's own chemistry column " *
              "disagree by $(c.worst) on the registered rule; the seam every chemistry " *
              "variant crosses does not carry the rule it is supposed to carry")
    rec!("chemistry_route_rows_compared", c.compared)
    rec!("chemistry_route_worst_difference", c.worst)
    @printf("  chemistry route: %d rows agree with the call table's own column (worst %.2e)\n",
            c.compared, c.worst)
end

# THE MOVE-IN LEVEL GREW THE FOREGROUND. The arm builder refuses a pool that named frame
# taxa and moved none; this is the other half, measured against the specification the level
# varies.
# It is checked on the whole arms and again on the thin-evidence subset, because the subset is
# a second run of the same level and a level that is inert on one of them is inert.
for suffix in ("", MIN_DISTINCT_SUFFIX)
    mv = G.runs[("B:$(HE.MOVEIN_MODE)", suffix)].type_arms
    base = G.runs[("primary", suffix)].type_arms
    g = check_movein_growth(mv, base)
    tag = isempty(suffix) ? "" : "_thin_evidence"
    rec!("movein_n1_before$tag", g.n1_before); rec!("movein_n1_after$tag", g.n1_after)
    rec!("movein_n2_before$tag", g.n2_before); rec!("movein_n2_after$tag", g.n2_after)
    rec!("movein_moved_taxa$tag", g.moved)
    rec!("movein_moved_taxa_annotated$tag", g.moved_annotated)
    @printf("  move-in level%s: n1 %d -> %d, n2 %d -> %d (%d moved, %d annotated)\n",
            isempty(suffix) ? "" : suffix, g.n1_before, g.n1_after, g.n2_before, g.n2_after,
            g.moved, g.moved_annotated)
end

# The thin-evidence subset has to bite on BOTH arms, or it is a filter on the foreground
# wearing the name of a unit check. Checked for every specification it is run for, because
# the arms it filters are that specification's arms and not the registered one's.
let base_specs = vcat([primary_spec()], sensitivity_specs()), t = nothing
    for s in base_specs
        # the registered specification's own result is kept as it passes, rather than the
        # guard being run a second time on the same arms to get a value to record
        r = check_both_arms_filtered(G.runs[(s.name, MIN_DISTINCT_SUFFIX)].type_arms,
                                     G.runs[(s.name, "")].type_arms, s.name)
        s.name == "primary" && (t = r)
    end
    rec!("thin_evidence_specifications", length(base_specs))
    rec!("thin_evidence_n1", t.n1); rec!("thin_evidence_n2", t.n2)
    rec!("thin_evidence_n1_base", t.n1_whole); rec!("thin_evidence_n2_base", t.n2_whole)
end

# WHICH VARIANTS ACTUALLY MOVED SOMETHING. A one-at-a-time grid over two orthogonal rule
# families produces cells no axis can reach, and those cells hold the primary's own numbers;
# counting the movement is what separates them from agreement. A variant that moves nothing on
# any unit is refused by name — that is not orthogonality, it is a level wired to nothing, or
# the registered specification wearing a variant's name.
movement = estimate_movement(rows, S)
rec!("variants_that_moved_an_estimate_somewhere",
     assert_every_variant_moves_something(movement, S))
# AND NO TWO OF THEM ARE ONE MEASUREMENT. The clone check above compares declarations, and a
# duplicated configuration in the CALL TABLE leaves every declaration distinct while making two
# levels the same measurement — no edit to this file, fourteen columns over thirteen
# measurements, and a category both of them broke keeping one breaker where a reader counts
# two. This is the same question asked of the numbers.
pair_differences = variant_pair_differences(rows, S)
rec!("distinct_variant_measurements",
     assert_no_two_variants_are_one_measurement(rows, S, pair_differences))
# AND HOW CLOSE THE CLOSEST LEGITIMATE PAIR CAME TO IT — per unit, not only overall. The guard
# fires at zero cells apart, so a margin read against every cell the flag reads looks like far
# more headroom than it is. Two levels of one rule family agree
# on every cell of the other family's unit BY CONSTRUCTION, because their axis is declared
# unable to reach it. The cells that could have separated them are the ones both axes reach,
# and that denominator is the one a reader of the margin needs.
let m = variant_pair_margin(pair_differences, specs)
    if m.overall !== nothing
        rec!("closest_variant_pair", m.overall.pair)
        rec!("closest_variant_pair_cells_differing", m.overall.differing)
        rec!("closest_variant_pair_cells_compared", m.overall.compared)
    end
    for (u, x) in m.per_unit
        rec!("closest_variant_pair_on_$u", x.pair)
        rec!("closest_variant_pair_on_$(u)_cells_differing", x.differing)
        rec!("closest_variant_pair_on_$(u)_cells_compared", x.compared)
    end
    rec!("closest_variant_pair_note",
         "the overall row is taken over every pair of levels and every cell the flag reads; " *
         "each per-unit row is taken over the pairs whose two axes are BOTH declared able to " *
         "reach that unit, which is the only comparison in which agreeing is informative. A " *
         "pair of levels on one axis that agrees on every reachable cell agrees on all of " *
         "them, which is why the guard is sound at the smaller denominator too")
end
# AND THE UNITS EACH ONE REACHED ARE THE UNITS ITS AXIS DECLARES. Without this the "does not
# reach" lines below are read off the measurement they are supposed to describe, so a level
# that lost its wiring to one unit is published as orthogonality — in the provenance, and as
# `n.r.` in the written table, which says in words that the rule family does not touch that
# unit at all.
rec!("axis_unit_declaration",
     joined([a.axis * "=" * join(a.units, ",") for a in VARIANT_AXES]))
rec!("variants_that_reached_the_units_their_axis_declares",
     assert_axes_reach_the_units_they_declare(movement, specs))
for u in movement.units
    reached = [v for v in S if get(movement.per_unit, (v, u), 0) > 0]
    rec!("variants_whose_axis_reaches_$u", length(reached))
    rec!("variants_whose_axis_does_not_reach_$u",
         joined([v for v in S if !(v in reached)]))
end
# HOW FAR EACH LEVEL MOVED, BESIDE HOW FAR IT COULD. The guard above is a threshold at zero,
# and a variant that moves one cell of those its axis reaches passes it while
# reproducing the registered specification everywhere else — which is what a half-applied
# configuration looks like, not what a rule variant looks like. The frozen rules set no
# threshold on how far a level must move, so none is invented here: the pair of numbers is
# recorded per variant, and a one-cell level is then visible beside its siblings.
let by_name = Dict(s.name => s for s in specs if s.role == ROLE_VARIANT)
    rec!("variant_cells_moved",
         joined(["$v=$(get(movement.per_variant, v, 0))" for v in S]))
    rec!("variant_cells_reachable",
         joined(["$v=$(declared_reachable_cells(by_name[v], movement))" for v in S]))
    rec!("variant_cells_note",
         "variant_cells_moved is how many flag cells that level produced a different number " *
         "for; variant_cells_reachable is how many cells the flag has on the units its axis " *
         "is declared able to reach. The two are equal for a level that moved every cell it " *
         "could, and a level far below its own denominator while its siblings are at theirs " *
         "is a level to look at rather than a result")
end

# ---------------------------------------------------------------- the flag
flags = robustness(rows, S; movement = movement)
rec!("categories_flagged", length(flags))
for f in flags
    rec!("flag_$(f.unit)_$(f.category)", f.robust_sign)
    isempty(f.breaking_variants) ||
        rec!("breakers_$(f.unit)_$(f.category)", f.breaking_variants)
end
robust_plus  = [f.category for f in flags if f.robust_sign == "+"]
robust_minus = [f.category for f in flags if f.robust_sign == "-"]
rec!("robust_enriched", joined(sort(robust_plus)))
rec!("robust_depleted", joined(sort(robust_minus)))
rec!("null_categories", joined(sort([f.category for f in flags if f.robust_sign == "null"])))
rec!("fragile_categories",
     joined(sort([f.category for f in flags if f.robust_sign == "fragile"])))
rec!("unit_checks_feed_the_flag", false)
rec!("unit_checks", joined(("genera", "sequences", "n_distinct>=$MIN_DISTINCT_TERMS")))
# THE GENUS SERIES, CHECKED RATHER THAN COUNTED. A downstream decision rests on this series
# holding across the whole grid, so what is recorded is which specifications and which cells
# were checked — a count of specifications that have a genus row somewhere would read the same
# whether the series were complete or had lost a category under one level and gained a row
# from a line outside the grid.
let g = check_genus_series(rows, vcat(["primary"], S))
    rec!("unit_check_genus_series_specifications", g.specifications)
    rec!("unit_check_genus_series_units", joined(g.units))
    rec!("unit_check_genus_series_cells_checked", g.cells)
    rec!("unit_check_genus_series_cells_missing", length(g.missing_cells))
    # NOT AN ERROR, AND NAMED EITHER WAY. A rule variant that stops calling a category for
    # any taxon legitimately has nothing to measure there, so a missing cell is a result of
    # the grid and not a failure of it — but it is also the thing that would quietly turn
    # "every specification agrees" into "every specification that had a row agrees", so it
    # leaves this run as data with the cells named rather than as a blank in a table.
    rec!("unit_check_genus_series_cells_missing_are", joined(g.missing_cells))
    isempty(g.missing_cells) ||
        @warn("the genus series is missing $(length(g.missing_cells)) cell(s); they are " *
              "named in the provenance, and any claim over the whole series has to be read " *
              "against them")
    @printf("  genus series: %d cells across %d specifications on %s, %d missing\n",
            g.cells, g.specifications, joined(g.units), length(g.missing_cells))
end

# ---------------------------------------------------------------- the anchors
anchor_rows = [r for r in rows if r.variant == "anchor:all_legacy"]
rec!("anchor_levels_retired", ANCHOR_LEVELS_RETIRED)
targets = NamedTuple[]; exceptions = NamedTuple[]
source_checks = NamedTuple[]
expected_reference_rows = nothing
if anchor_file === nothing
    # Reached only when --no-anchors was typed: the option parser refuses the silent case.
    @warn("--no-anchors was given, so nothing checked that this grid is wired to the data " *
          "it claims to describe")
    rec!("anchor_targets_supplied", false)
    rec!("anchor_checks_declined", true)
else
    targets = anchor_file.targets; exceptions = anchor_file.exceptions
    rec!("anchor_targets_supplied", true)
    rec!("anchor_checks_declined", false)
    rec!("anchor_targets_file", basename(OPT.anchor_targets))
    # HOW MANY TARGETS THE FILE SAID IT HELD. `anchor_targets_checked` below is a count of what
    # the file contained, and a file that lost rows is a run that checked less and reported
    # success in the same words. `read_anchor_targets` refuses the mismatch; this records the
    # number the refusal was made against, so a reader of the provenance alone can see how
    # strong the wiring check was rather than only that it passed.
    rec!("anchor_targets_expected", anchor_file.expected[:targets])
    rec!("anchor_targets_on_the_registered_specification",
         count(t -> t.anchor == primary_spec().name, anchor_file.targets))
    # AND HOW MANY OF THEM SIT ON A LEVEL. The registered specification's own rows are one of
    # the fifteen SPECIFICATIONS whose intervals decide a word; the other fourteen are the
    # levels, and a target
    # count that does not separate the two says nothing about them. The names with no target at
    # all are recorded too, so the reach of the wiring check is a list rather than an
    # impression.
    let anchored = sort(unique(String[t.anchor for t in anchor_file.targets if t.anchor in S]))
        rec!("anchor_targets_on_a_registered_level",
             count(t -> t.anchor in S, anchor_file.targets))
        rec!("registered_levels_with_a_target", length(anchored))
        rec!("registered_levels_without_a_target", joined([v for v in S if !(v in anchored)]))
    end
    # WHERE THE NUMBERS CAME FROM, AS BYTES. A target transcribed from this grid's own output
    # is identical to an honestly transcribed one, so the file names the implementation it was
    # transcribed from and carries that file's digest; a source present beside the targets and
    # not matching its declared digest stops the run below.
    source_checks = check_target_sources(OPT.anchor_targets, anchor_file.sources)
    rec!("anchor_target_sources", joined([s.file for s in source_checks]))
    for s in source_checks
        rec!("anchor_target_source_status_$(s.file)", s.status)
        rec!("anchor_target_source_declared_sha256_$(s.file)", s.declared)
    end
    for s in anchor_file.sources
        rec!("anchor_target_source_note_$(s.file)", s.note)
    end
    expected_reference_rows = anchor_file.expected[:reference_rows]
    # The targets are the one check that proves this grid is wired to its data, and they arrive
    # as a file named on the command line rather than as constants of this script. Without
    # their digest a later step can confirm neither that the check ran nor that it ran against
    # the same bytes — while the call table and the scan table, which decide less, both carry
    # one.
    rec!("anchor_targets_sha256", bytes2hex(sha256(read(OPT.anchor_targets))))
end
target_result = check_anchor_targets(rows, targets)
rec!("anchor_targets_checked", target_result.checked)
rec!("anchor_targets_met", target_result.met)
rec!("anchor_targets_failed", length(target_result.failures))
isempty(target_result.failures) ||
    rec!("anchor_target_failures", joined(target_result.failures))
for e in exceptions
    rec!("anchor_named_exception_$(e.unit)_$(e.category)", e.note)
end

# THE LEVELS THIS RUN PUBLISHES BY NAME MUST BE PINNED BY SOMETHING. A breaker list is the
# published claim that a result does not survive THAT rule and does survive the others, and it
# is decided by the level's own arms — which the registered specification's targets never
# touch. The requirement is read off what this run actually published rather than listed here,
# so a level that starts or stops breaking a category carries the requirement with it; a list
# written into this file would go stale in exactly the run where it mattered.
published_levels = anchor_file === nothing ? String[] : levels_named_in_a_breaker_list(flags)
unanchored_levels = anchor_file === nothing ? String[] :
                    unanchored_published_levels(flags, targets)
rec!("levels_named_in_a_published_breaker_list", joined(published_levels))
rec!("levels_named_in_a_published_breaker_list_without_a_target", joined(unanchored_levels))

ref_ok = true
reference_rows_compared = 0
# What the targets file declares the comparison covers, carried only once a comparison was
# actually made: a run that never reached the reference has no scope to disagree with, and its
# refusal is the gate above rather than this one.
reference_rows_declared = nothing
if ANCHOR_LEVELS_RETIRED && isempty(anchor_rows)
    # NOT APPLICABLE, AND SAID SO RATHER THAN PASSED QUIETLY. The row-by-row comparison is a
    # statement about the all-legacy anchor, and that level is not evaluated against the
    # census this run reads (see ANCHOR_LEVELS_RETIRED). With no anchor rows to compare, every
    # row of an earlier run's reference would count as `absent` and the gate could never pass
    # — not because this grid disagrees with the earlier one, but because it was never asked
    # the question.
    #
    # Pointing --anchor-reference at a file holding no rows would pass for the same arithmetic
    # reason while recording nothing about WHY: a reader of the provenance would see a
    # comparison that covered zero rows and no statement that covering zero is the intended
    # outcome. So the status is written out in words, the count stays 0, and the targets
    # file's `reference_rows,0` is the declaration that agrees with it; in this branch that
    # declaration is read and not compared.
    #
    # A reference named anyway is recorded by its file name and NOT opened; see
    # `record_reference_not_applicable!`.
    ref_ok = record_reference_not_applicable!(rec!, OPT.anchor_reference)
elseif !isempty(OPT.anchor_reference)
    if !isfile(OPT.anchor_reference)
        rec!("anchor_reference_status", "missing")
        ref_ok = false
    else
        ref = HE.read_reference_rows(OPT.anchor_reference)
        if ref === nothing
            rec!("anchor_reference_status", "unreadable")
            ref_ok = false
        else
            scoped = scoped_reference(ref, exceptions; units = REFERENCE_SCOPE_UNITS)
            d = HE.reference_row_diff(reference_shape.(anchor_rows), scoped.table)
            rec!("anchor_reference_file", basename(OPT.anchor_reference))
            rec!("anchor_reference_sha256",
                 bytes2hex(sha256(read(OPT.anchor_reference))))
            rec!("anchor_reference_rows_compared", d.compared)
            # HOW MANY ROWS THE COMPARISON WAS SUPPOSED TO COVER. "N of N rows agree" is the
            # same sentence at any N, so a reference file regenerated at another scope — more
            # rows, fewer rows, a different vintage — reports a number that reads better than
            # the honest run's rather than different from it. The count is declared in the
            # targets file and compared here.
            reference_rows_compared = d.compared
            if expected_reference_rows !== nothing
                rec!("anchor_reference_rows_expected", expected_reference_rows)
                reference_rows_declared = expected_reference_rows
            end
            rec!("anchor_reference_rows_agreeing", d.agreeing)
            rec!("anchor_reference_rows_differing", joined(d.differing))
            rec!("anchor_reference_rows_absent", joined(d.absent))
            rec!("anchor_reference_worst_difference", d.worst)
            rec!("anchor_reference_named_exceptions", joined(scoped.excluded))
            # THE ROWS THAT LEFT THE COMPARISON WITH THE UNIT THEY WERE MEASURED ON. Without
            # these two rows the provenance says "every compared row agrees" over a scope
            # nobody can see, and the largest deliberate difference between the two runs is
            # the one thing the machine-readable output does not mention.
            rec!("anchor_reference_rows_out_of_scope", length(scoped.out_of_scope))
            rec!("anchor_reference_out_of_scope_units", joined(scoped.out_of_scope_units))
            rec!("anchor_reference_out_of_scope_reason", REFERENCE_OUT_OF_SCOPE_REASON)
            rec!("anchor_reference_status",
                 isempty(d.differing) && isempty(d.absent) ? "reproduces" : "differs")
            ref_ok = isempty(d.differing) && isempty(d.absent)
            @printf("  anchor vs the earlier run: %d of %d rows agree%s\n",
                    d.agreeing, d.compared,
                    isempty(scoped.excluded) ? "" :
                    " (named exceptions: " * joined(scoped.excluded) * ")")
        end
    end
else
    rec!("anchor_reference_status", "not supplied")
end

# ---------------------------------------------------------------- write
# Nothing is written into the figure directory here: everything is built in the staging
# directory beside it and moved in by `publish!`, and only if the gates below have passed.
stage = open_staging(OPT.figdir)
staged = String[]
stage!(name) = (push!(staged, name); joinpath(stage, name))

out_rows = stage!("habitat_sensitivity.csv")
CSV.write(out_rows, DataFrame(rows))

order = vcat(["primary"], S)
out_md = stage!("habitat_sensitivity_table.md")
open(out_md, "w") do io
    write(io, sensitivity_table(rows, order; movement = movement))
    write(io, "\n")
    write(io, unit_check_table(rows, order))
    write(io, "\n")
    write(io, columns_note())
end

out_flags = stage!("habitat_robustness.csv")
CSV.write(out_flags, DataFrame(flags))

panel = panel_rows(rows, flags)
note = panel_note(flags, panel)
if OPT.figure
    mk = cairomakie()
    written = Base.invokelatest(draw_panel, mk, panel,
                                joinpath(stage, PANEL_BASE); note = note)
    for p in written; push!(staged, basename(p)); end
    rec!("panel_written", joined(basename.(written)))
else
    rec!("panel_written", "")
end
# A RUN THAT DOES NOT DRAW THE PANEL MUST NOT PUBLISH BESIDE AN OLD ONE. `publish!` moves what
# was staged and leaves everything else where it is, so `--no-figure` into a directory that
# already holds a panel leaves this run's tables beside the previous run's picture, under the
# names a reader takes as one result. Measured here rather than asserted, so the number is in
# the provenance of the run that was refused.
stale_panels = stale_panel_files(OPT.figdir, OPT.figure)
rec!("stale_panel_files", joined(stale_panels))
rec!("panel_categories", length(panel))
rec!("panel_categories_not_drawn",
     joined(sort([f.category for f in flags
                  if !((f.unit, f.category) in
                       Set((p.unit, p.category) for p in panel))])))
rec!("output_columns",
     joined(String.(propertynames(DataFrame(rows[1:1])))))
rec!("output_column_note",
     "lo and hi are the STANDARDISED bounds, which every rule reads; crude_lo and crude_hi " *
     "are the crude ones. unit and contrast are additions to the planned column list, " *
     "without which a genus row and a taxon row of the same variant and category cannot be " *
     "told apart")
rec!("robustness_column_note",
     "robust_sign carries a sign only where the primary interval AND every one-at-a-time " *
     "variant exclude zero on that side; primary_excludes_zero is the fact the word rests " *
     "on, so a later step keys colour off a boolean rather than off a word. " *
     "n_variants_that_moved_the_estimate is how many of the expected variants produced a " *
     "different number for that cell: the rest are cells the variant's axis cannot reach, " *
     "and they are not independent checks. depletion_tier is a verdict on a depletion and " *
     "is left empty where the primary estimate is positive")

# ---------------------------------------------------------------- the gates
# BUILT BY A FUNCTION, NOT WRITTEN HERE, and built BEFORE the provenance is written so the run
# records which gates it made. The predicates behind them are each tested on their own; the
# wiring that hands them to `publish!` is a function rather than a literal inside this program
# guard, which no test can execute — there, a gate entry deleted, or its `ok` quietly replaced
# by `true`, leaves the whole suite green while the attack that gate exists to refuse
# publishes at exit 0.
gates = publish_gates(; target_failures = target_result.failures,
                      unanchored_levels = unanchored_levels,
                      source_checks = source_checks,
                      reference_ok = ref_ok,
                      reference_rows_compared = reference_rows_compared,
                      reference_rows_expected = reference_rows_declared,
                      stale_panels = stale_panels, figdir = OPT.figdir)
rec!("publish_gates", joined([g.name for g in gates]))
rec!("publish_gates_failed", joined([g.name for g in gates if !g.ok]))

out_prov = stage!("habitat_sensitivity_provenance.csv")
CSV.write(out_prov, DataFrame(key = first.(PROV), value = string.(last.(PROV))))

# ---------------------------------------------------------------- report
println("\nROBUSTNESS")
for f in sort(flags, by = f -> (f.unit, f.category))
    @printf("  %-18s %-14s %-9s %s\n", f.unit, f.category, f.robust_sign,
            isempty(f.breaking_variants) ? "" : "breakers: " * f.breaking_variants)
end
println("\nANCHORS")
@printf("  targets: %d checked, %d met, %d failed\n",
        target_result.checked, target_result.met, length(target_result.failures))
for m in target_result.failures; println("    MISSED ", m); end
@printf("  levels named in a published breaker list: %d, of which %d carry no target\n",
        length(published_levels), length(unanchored_levels))
for s in source_checks
    @printf("  source %s: %s\n", s.file, s.status)
end

# ---------------------------------------------------------------- publish
for p in publish!(stage, OPT.figdir, staged, gates)
    println("wrote ", p)
end

println("\ndone.")

end # program guard
