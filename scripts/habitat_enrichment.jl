#=
    habitat_enrichment.jl — where do lanthanide-dependent PQQ-ADH carriers live,
    RELATIVE TO EXPECTATION?

    Raw habitat counts measure research effort, not ecology. The habitat layer is
    text-mined from published work, so a single well-studied model organism can dominate a
    category on its own. Every number here is therefore a comparison against a denominator,
    and every denominator is computed by the run.

    NO DATA LITERALS IN THIS FILE. Counts, shares, coverages, thresholds, releases and
    p-values are computed or read at run time, written to habitat_arm_summary.csv, and the
    closing report reads them back out of that table. A number typed into a comment goes
    stale the moment an input changes and then quietly contradicts the output, so this file
    contains none. The summary table is the single place a reported number may come from.
    The one class of number this file accepts from outside is a bound the analysis plan
    fixed BEFORE the run: those arrive on the command line, are never defaulted, and are
    disclosed beside the quantity they bound (see THE RECORDED CRITERIA below).

    NO ESTIMATOR OF ITS OWN. crude(), standardised(), the Benjamini-Hochberg correction,
    the genus collapse and the row builder live in src/ModuleHabitatStats.jl, where they
    are unit-tested. This script only selects arms, builds weight vectors and reports. The
    one resampling loop it does run (the sharing ratio below) calls the module's
    standardised() on resampled indices rather than re-deriving anything.

    NO OUTCOME RULE OF ITS OWN EITHER. Every code this script emits — the robustness tier,
    the sharing grade, the host-plant codes, the category partition, the background licence
    — is decided by src/ModuleOutcomeRules.jl, where the predicates are unit-tested. A
    one-line predicate sitting here, next to the code that prints its answer, is reachable by
    no test — and an unreachable predicate that cannot fail, or was never implemented at all,
    looks exactly like one that works. That holds for the pair of ternaries that turn an
    interval into a host-plant code too: one tested mapper answers "what does this interval
    decide?" for every rule of that shape, and its third value is not this script's to
    rename, so an estimate that was never computed reports as undecidable instead of being
    folded into "includes zero".

    NO RECORD PARSER OF ITS OWN. The background scan persists three records beside its
    table, and their format belongs to the script that writes them: this one includes
    scripts/background_pqq_scan.jl behind its program-file guard (which keeps that script's
    driver from running) and calls its readers. A consumer that re-implements a producer's
    format eventually disagrees with it, and the disagreement shows up as a missing field,
    which is exactly the input every guard below has to treat as unknown.

    NO AGGREGATION OF ITS OWN EITHER. The habitat calls come pre-computed, one row per
    (configuration, taxon), from scripts/habitat_reaggregate.jl; the background carrier scan
    comes from scripts/background_pqq_scan.jl. Switching --config switches the whole tie /
    isolation / generic / chemistry-text rule set in one place, which is what makes the
    sensitivity grid downstream cheap and this script short.

    THE FRAME. Most census taxa are themselves reference proteomes, so the background is
    not an external comparator — it is the SAMPLING FRAME MINUS THE CENSUS: reference-
    proteome prokaryotes that do not carry the enzyme family. Both arms are drawn from one
    frame, so its cultivation and sequencing bias largely cancels. Reference proteomes are
    deduplicated to about one per species group, so a heavily sequenced species counts once
    — that is precisely the sampling-frequency correction wanted.

    AND MINUS THE CARRIERS IT MISSED. A background built by set subtraction still contains
    taxa that carry the enzyme but were never in the census, and those taxa carry the
    habitat signal being tested, so leaving them in biases every contrast toward the null.
    They are removed here, not argued away, on the targeted hmmsearch of the frame
    (--bg-exclusion hmmsearch, the default and the registered primary — "registered" in the
    exact sense docs/TOOLING.md §5.1 defines, and nothing stronger). That is the
    broad exclusion: it removes every frame taxon with ANY hit of the family, whatever metal
    class the hit reads, and it is the producer's own tested exclusion predicate that
    removes them, not a second copy of it. The narrower questions — which background taxa
    carry a LANTHANIDE-type enzyme specifically, or a hit whose active site reads as either
    metal — are different questions, and the modes that answer them (lntype, lnca) are
    reported as sensitivity lines, never as the headline. Which mode ran is recorded in the
    output, together with whether it was typed or inherited and whether it is the registered
    primary.

    AND ONE LEVEL THAT ADDS INSTEAD OF REMOVING. Every mode above only ever takes taxa OUT
    of the background. The recorded sensitivity grid also names a level that removes the
    carriers AND gives the lanthanide-type ones to the foreground, and that is the level
    which can break a result rather than merely move it: it asks whether the contrast is
    driven by where the carriers live or by the census's own membership, since a frame taxon
    the census never registered is a carrier by every criterion except having been counted.
    It is --bg-exclusion hmmsearch_movein, and it is a sensitivity level, never the headline.
    Which taxa it moves is not this script's opinion either: the pool arrives as a rule that
    is named in the output, with its size counted at run time, and the rule the record names
    reads the CENSUS'S OWN class for a sequence, out of the census's pre-deduplication hit
    table, which is not a repository input and so has to be supplied (--movein-pool). The
    rule this repository can answer alone — the scan's own lanthanide-class call — is a
    DIFFERENT measurement, so a run that falls back to it says so in a key and in a warning
    rather than presenting its pool size as the recorded one.

    A MOVED TAXON IS IN EXACTLY ONE ARM. It joins the foreground and leaves the frame in the
    same expression, so the background it was drawn from no longer holds it and n1 cannot
    count it twice; both properties are asserted where the arms are assembled, because a
    variant that quietly counted a taxon on both sides would report a contrast between an arm
    and part of itself. The moved taxa have NO census sequences — that is what made them
    background — so the per-sequence unit gains nothing from them and says so, and the
    comparator arms, whose membership rules are census-class rules, are left exactly as they
    were.

    WHAT CAN STILL GO WRONG WITH THAT EXCLUSION, AND WHICH HALF OF IT IS DEFINITIONAL. Two
    criteria are measured SEPARATELY, and over the arm that is ACTUALLY ANALYSED: HMM hits
    on one side, whole-frame family membership at the length cut-off — a property of the
    taxon, computed without the sequence prefilter — on the other. Under the any-hit mode
    the first of the two IS THE EXCLUSION PREDICATE, so its leftover count is zero by
    definition and can never fail, whatever the data say; it is emitted as a statement of
    that identity, not as evidence, and a key records it as definitional. The count carries
    real information only under the narrower modes and under "none". Five inputs decide the
    licence, and these four CAN fail:

      * whole-frame family membership left in the analysed arm — a full-length family
        member the model scores below the cutoff leaves its taxon unflagged while the
        whole-frame roll-up still counts it;
      * the prefilter's persisted recall verdict, RE-DERIVED from the record's own counts
        rather than read off its verdict field, and taken on the exact lower confidence
        bound of those counts rather than on their ratio;
      * a taxon-by-taxon comparison of the scanned set against the analysed frame — a
        carrier whose row is missing from the table reads as a non-carrier;
      * whether the scan's own persisted records vouch for the table that was read at all.

    A RATE AT THE BAR IS NOT A RECALL AT THE BAR. The recall check samples frame taxa, scans
    their whole proteomes, and asks what share of the carriers found there the prefilter
    also found. A sample that recovered every carrier it flagged has a rate of one at any
    size, so the rate clears the bar on three flagged taxa as readily as on three thousand.
    What the statement "the background holds no carrier the prefilter missed" needs is what
    the sample can vouch for: the exact one-sided 95 % Clopper-Pearson lower bound of the
    recovered share, held to the same bar. The scan records that bound beside the rate and
    gates its own publication on the rate, because a bound below the bar does not make its
    exclusion list wrong; the stronger statement is decided here, on the bound, and a clean
    sample too small for its bound to reach the bar leaves the licence withheld. Both the
    rate and the bound the counts give are emitted beside the record's copies of them.

    THE FIFTH INPUT EXISTS BECAUSE THREE OF THE OTHERS COME OUT OF ONE FILE. A scan table
    that covers the frame but reports no carrier anywhere satisfies every count a consumer
    can compute from it: nothing is left because nothing was ever found, and the taxon set
    still matches. Only evidence from outside the table separates "measured clean" from
    "measured nothing", and that evidence is the producer's own acceptance table, tied to
    the table by the digest of its bytes. So the module is handed the digests; the
    acceptance checks and verdicts together with the names of the gated checks the scan
    writes on every run; two counts taken from the parsed table, to hold against the
    acceptance table's; and every field of the recall record together with the same fields
    as the other records carry them. It refuses the licence when a record describes another
    table, when a check missed, when a gated check is absent or doubled, when the acceptance
    table does not count the table read, when the recall record disagrees with its own
    counts or with the other records, or when the recall verdict does not hold. Here, beside
    it, every acceptance row's verdict is re-derived from that row's own numbers by the
    producer's rule, and the recall row is held against the recall record, because a
    record whose verdict words say PASS beside numbers that fail is not evidence of
    anything. A record that is absent, unreadable or of an older shape is unknown, and
    unknown withholds the licence instead of erroring or passing.

    ONE READ OF EVERY FILE THE LICENCE RESTS ON. The digest is taken from the very bytes
    that are parsed, in one read, and each record is read from disk once. A table hashed
    from a second read, taken seconds after the first, vouches for whatever a scan run
    promoted in between, and the run then analyses one table under the records of another
    while every digest agrees. With one read, a promotion that lands between the table and
    its records shows up as a digest that does not match, which withholds the licence.

    AND A TABLE NOBODY PROMOTED IS NOT A RESULT. The producer publishes its table and its
    records together or not at all: a failing run leaves its own files as .pending and
    renames the last published set to .stale. Those names are therefore load-bearing, and
    they are checked here, because a .stale set is internally consistent — it was a passing
    run once — and would otherwise license the strong code after a later run had refused
    the frame.

    ONE ARM BUILDER, CALLED TWICE. analyse() turns a configuration, a chemistry column and
    an exclusion set into every arm, its weight vectors and the result rows, and nothing
    else in this file builds an arm. The run calls it for its own configuration, and again
    for the ANCHOR: the older configuration, winner chemistry, no background exclusion,
    which is what reproduces the earlier reference run. The anchor is what makes "this
    category excludes zero for the first time" a measurable statement rather than an
    impression — the question can only be asked against the earlier result a category would
    be new relative to. Calling one function twice is what stops the anchor from drifting
    away from the path the primary took; a second copy of the arm logic would eventually
    disagree with the first. For the same reason ONE function reads a category set off a row
    panel, and it is the one used for the primary, for the anchor and for the reference
    table the anchor is checked against.

    THE ANCHOR IS ONLY VERIFIED WHEN IT IS COMPARED. Its rows are written to their own CSV
    (--anchor-out) so the partition can be re-decided at row level without a re-run, and
    --anchor-reference compares them, row by row and set against set, with an earlier run's
    rows. Without that option the anchor is computed but unverified, and the output says so
    in a key rather than in a comment.

    THE CONFOUND, AND WHY STANDARDISATION IS NOT OPTIONAL. The two arms differ sharply in
    phylum composition (the run prints both shares). Habitat annotation tracks phylum, so
    an unadjusted ratio will report "enrichment" that is nothing but the protein family's
    phylogenetic range. Two estimates are reported for every habitat:

      CRUDE         a/n1 ÷ b/n2 — unadjusted.
      STANDARDISED  the background's habitat share REWEIGHTED to the foreground phylum
                    distribution: p0_std = Σ_k w_k (c_k/n0_k), w_k = n1_k/n1.
                    "Phylum-matched" by reweighting rather than by discarding background,
                    so it costs no data.

    The GAP between the two is the phylogenetic contribution, made legible instead of
    argued about. Trust only what survives standardisation.

    AND THE STANDARDISED BOUNDS ARE THE ONES THAT DECIDE. Where the recorded rules name an
    interval at all, they name the standardised one, so every code below is decided on the
    standardised bounds, as recorded. The crude-reading counterpart of each code is emitted
    beside it as a DESCRIPTIVE line only: it says how much of a code rests on the phylum
    adjustment, and it licenses nothing. The two readings do not always agree, which is the
    disclosure.

    AND A HEALTHY ROW CAN STILL REST ON AN EMPTY STRATUM. The standardised interval is a sum
    over phylum strata, and a stratum whose background cell is empty contributes nothing to
    the variance — so a row whose arm totals are comfortable can carry an interval decided
    by one thin cell. Each rule row therefore also carries the stratum disclosures the
    estimator returns, the run counts how many rows have an empty stratum cell, and every
    category in the partition carries the smallest background cell behind it, so a category
    resting on one thin cell cannot enter a list silently.

    UNIT = UNIQUE TAXON. Counting sequences counts how often an organism was sequenced, not
    how many organisms carry the enzyme, and a per-sequence signal can be nothing but
    sampling structure (a within-subfamily permutation test is the control for that). The
    per-sequence and per-genus series are computed as secondary so the size of that bias
    stays visible. The genus series collapses BOTH arms
    to genus, so n1 and n2 are in the same currency; putting foreground genus-by-habitat
    PAIRS against background TAXA is not a like-for-like comparison and cannot be read as a
    conservative bound either.

    THE COMPARATOR ARMS. The registered comparator for the sharing ratio is the non-Ln
    carrier arm, whose membership rule is the ABSENCE of a lanthanide-type census sequence
    — which is not the same thing as the PRESENCE of a calcium-type one. Redefining it would
    need a recorded decision, so it is not redefined: its class composition is disclosed,
    and a second arm restricted to taxa with a calcium-class census sequence and no
    lanthanide-class one is reported BESIDE it, with its own ratio, its own resampled
    interval from the same routine and seed, and its own code. Those keys are disclosures;
    they never replace the registered comparator's, and the row counters that feed the
    re-check list come in a second form counted over the registered contrasts alone.

    PROKARYOTES ONLY. The census's eukaryotic candidates are excluded here; whether any of
    them is a genuine carrier is decided outside this pipeline, by an intron/transcript
    check, and their count is in the summary table.
    Isolation sources record where a strain was obtained, not where it lives, and are
    commoner in the background, so leaving them in manufactures a fake depletion. Under the
    primary configuration they are dropped TERM by term, which recovers taxa that the older
    taxon-level drop discarded whole; the legacy configuration drops the taxon.

    THE RECORDED CRITERIA. The plan for this step fixed an upper bound on the size of the
    analysed background and on the background's share of one phylum before the step ran.
    Neither is a quantity this script may invent, so neither has a default: they arrive by
    option, and when they do the run emits the bound, the measured value, the comparison and
    whether it was met, as DISCLOSED rows that gate no code and warn when they are not met.
    Omitting them omits the rows; nothing silently passes because a bound was not supplied.

    THE SELECTED CATEGORIES ARE NOT CONSTANTS OF THIS SOFTWARE, AND THIS SCRIPT SHIPS NONE.
    Three of the recorded outcome rules speak for one habitat category each, and the
    standardisation shares are quoted for one phylum. WHICH categories those are was chosen
    from this data set's own estimates, with the estimates in hand, AFTER those estimates
    existed. The classification scheme they are drawn from (src/ModuleEnvironmentTools.jl) is
    an external ontology's mapping and was fixed before any estimate existed; picking four
    names out of it was not, and this file must not present the second as if it were the
    first. The phylum is the one entry that is DETERMINED rather than chosen — it is the
    heaviest standardisation stratum of the arms, which this script computes and records
    beside the shares, and the figure step withholds the weight identities unless the run's
    own heaviest stratum is still that phylum — so it is supplied for a different reason: the
    share keys are named after it and a reader has to be able to derive that name. A default
    here would publish
    that choice as though this software had made it, and a script-supplied default cannot be
    told apart from a selection somebody checked. So they arrive in ONE criteria file
    (--categories), with no default anywhere, and a run that does not supply it REFUSES by
    name and says which option supplies it.

    AND THE EMITTED KEYS FOLLOW THE SELECTION. The summary's rule-keyed families are named
    after the categories they are about, so a run given other categories emits keys for THOSE
    categories rather than filing other numbers under the shipped names. The prefix is derived
    from the supplied name by one documented rule (`key_slug`); the criteria file's own digest
    and every selected name are emitted as data beside them, so the table says which selection
    produced it and two selections can be told apart before any number is compared.

    A FILE A TEST CAN INCLUDE. Every function sits above the program-file guard near the end
    of this file, and the driver sits below it; the command line is parsed under the same
    guard, before any package loads, so --help and a mistyped option still cost nothing.
    Including the file therefore defines the functions and runs nothing, which is how
    test/test_habitat_enrichment.jl pins the wiring of the licence guards — the single read
    of the table, the records, the file state, the assembly of the evidence, the leftovers,
    the frame comparison, the undecidable set and the anchor-status ladder — on tiny
    fixtures. As top-level code that wiring is reachable by no test, and a one-line change to it
    licenses the strong code on a refused table with every suite still green.

    Run:  julia --project=. scripts/habitat_enrichment.jl
          julia --project=. scripts/habitat_enrichment.jl --config legacy \
                --bg-exclusion none --out <path> --summary <path>
          (--help lists every option and exits before anything is loaded or read)
=#

const REPO = normpath(joinpath(@__DIR__, ".."))

# ---------------------------------------------------------------------------
# Command line. Parsed before any package is loaded and before any file is touched, so
# --help and a mistyped option cost nothing. Both would otherwise cost a run: an ignored
# unknown option means `--bg-exlusion lntype` analyses the DEFAULT exclusion and writes its
# numbers under the name of the run that was asked for, and a `--help` read after the work
# starts runs the whole analysis and overwrites the output tables on the way. The shape of
# the parser is the producer's (scripts/background_pqq_scan.jl), for the same reason it has
# one.
# ---------------------------------------------------------------------------

# The background exclusions this script can apply. The first is the registered primary; the
# rest are descriptive, and the run records which of them was in force. The move-in level is
# LAST because it is the only one that also changes the FOREGROUND: every mode before it
# answers "which frame taxa leave the background", which is the question the sensitivity
# tally below asks of the modes that are not in force, and that tally cannot describe a level
# which hands taxa to the other arm.
const EXCLUSION_MODES = ("hmmsearch", "lnca", "lntype", "ipr017512", "none",
                         "hmmsearch_movein")
const PRIMARY_MODE    = first(EXCLUSION_MODES)
const MOVEIN_MODE     = last(EXCLUSION_MODES)

# The active-site class that names a lanthanide-type enzyme. Both the scan table and the
# census hit table spell it, and the move-in pool, the lntype exclusion and the inventory all
# turn on it, so it is written once: three spellings of one class name is three places for a
# rename to be applied twice.
const LANTHANIDE_CLASS = "Lanthanide"
# How the move-in level decides WHICH taxa are lanthanide-type. The first is the rule the
# record names — the class the CENSUS itself gave a sequence, out of its own
# pre-deduplication hit table — and it arrives by option because that table is not an input
# of this repository. The second is the only rule the committed inputs can answer alone, and
# it is a different measurement: this run's model call on the family members it fetched. A
# run states which rule it used, and whether that rule is the recorded one.
const MOVEIN_POOL_RULES      = ("census_lanthanide_hits", "scan_lanthanide_hits")
const REGISTERED_MOVEIN_POOL = first(MOVEIN_POOL_RULES)

const USAGE = """
Usage: julia --project=. scripts/habitat_enrichment.jl [options]

Options:
  --categories <path>          REQUIRED, and there is no default. The category-selection
                               criteria file: which category each recorded outcome rule
                               speaks for, which categories the grid counts are exported
                               for, and which phylum the standardisation shares are quoted
                               for. These are selections made from this data set's own
                               estimates, not constants of this software, so this step
                               neither invents them nor falls back to any. See THE SELECTED
                               CATEGORIES in the header
  --config <name>              configuration column of the taxon-call table (default primary)
  --bg-exclusion <mode>        $(join(EXCLUSION_MODES, " | "))
                               (default $PRIMARY_MODE, the registered primary; every other
                               mode is descriptive, and $MOVEIN_MODE also moves the
                               pool below into the foreground)
  --movein-pool <path>         the lanthanide-type pool the $MOVEIN_MODE level
                               moves: a census hit table with classification and TaxID
                               columns. Read only by that mode; without it the level falls
                               back to the scan's own lanthanide-class call, which is a
                               different measurement and is recorded as such (no default)
  --chem-col <flags|winner>    chemistry source; defaults by configuration
  --boot <n>                   resampling replicates for the sharing intervals
  --seed <n>                   seed of that resampling
  --calls <path>               taxon-call table
  --scan <path>                background carrier scan table
  --census <path>              census table, for the class composition of the arms
  --recall <path>              the scan's persisted recall record
  --threshold-record <path>    the scan's persisted threshold record
  --acceptance <path>          the scan's persisted acceptance table
  --anchor-reference <path>    a rows CSV from the earlier reference run; when it is given,
                               the anchor must reproduce its category set or no category
                               partition is emitted at all
  --anchor-out <path>          where the anchor's own rows go (default: --out with the
                               suffix _anchor)
  --out <path>                 the result rows
  --summary <path>             the key/value summary every reported number comes from
  --plan-bg-size-below <n>     upper bound on the analysed background size, from the plan
                               for this step; disclosed, never gating (no default)
  --plan-bg-phylum-pct-max <p>
                               the same for the background's share of the phylum the
                               criteria file names, in per cent (no default)
  -h, --help                   print this text and exit before anything is loaded or read

Exit status: 0 on a completed run; 1 on a command-line error or a missing required input.
"""

const OPTION_SPEC = Dict("--categories" => :value,
                         "--config" => :value, "--bg-exclusion" => :value,
                         "--movein-pool" => :value,
                         "--chem-col" => :value, "--boot" => :value, "--seed" => :value,
                         "--calls" => :value, "--scan" => :value, "--census" => :value,
                         "--recall" => :value, "--threshold-record" => :value,
                         "--acceptance" => :value, "--anchor-reference" => :value,
                         "--anchor-out" => :value, "--out" => :value, "--summary" => :value,
                         "--plan-bg-size-below" => :value,
                         "--plan-bg-phylum-pct-max" => :value,
                         "--help" => :flag, "-h" => :flag)

"""
    parse_options(args, spec) -> Dict{String,String}

Split a command line into option => value, a flag mapping to "". `spec` names every option
this script defines as `:value` or `:flag`. Anything else is an error — an unknown or
misspelt option, a stray word, an option given twice, a value-taking option with no value —
because an ignored token is a silently different run whose output looks exactly like the
answer to the command that was typed.
"""
function parse_options(args, spec::AbstractDict)
    out = Dict{String,String}()
    i = 1
    while i <= length(args)
        a = String(args[i])
        kind = get(spec, a, nothing)
        kind === nothing &&
            error("unknown option '$a' (defined: " * join(sort(collect(keys(spec))), ", ") * ")")
        haskey(out, a) && error("option '$a' given twice")
        if kind === :flag
            out[a] = ""
            i += 1
        else
            # A value that is itself a defined option means the value was forgotten:
            # `--scan --boot 0` must not read a scan table called "--boot".
            (i == length(args) || haskey(spec, String(args[i + 1]))) &&
                error("option '$a' needs a value")
            out[a] = String(args[i + 1])
            i += 2
        end
    end
    return out
end

"""
    cli_or_exit(f) -> the value of `f()`

Run a command-line step; an `ErrorException` from it ends the program with the message and
the usage on stderr and exit status 1. The point is that a bad command line never falls
through to a default: it stops here, before anything is loaded, read or written.
"""
function cli_or_exit(f)
    try
        return f()
    catch e
        e isa ErrorException || rethrow()
        println(stderr, "habitat_enrichment.jl: ", e.msg, "\n")
        print(stderr, USAGE)
        exit(1)
    end
end

"""
    option_string(cli, flag, default) -> String
    option_number(cli, flag, T, default; what, min_value, max_value) -> T

The value of a parsed option, or `default` when it was not given. A value that is empty,
does not parse as `T`, is not finite, or lies outside the stated range is an error naming
the option — never a fall-back to the default, which would run a different analysis under
the name of the one that was asked for.
"""
function option_string(cli::AbstractDict, flag::AbstractString, default)
    haskey(cli, flag) || return default
    v = strip(cli[flag])
    isempty(v) && error("$flag needs a value")
    return String(v)
end

function option_number(cli::AbstractDict, flag::AbstractString, ::Type{T}, default;
                       what::AbstractString = "a number", min_value = nothing,
                       max_value = nothing) where {T<:Real}
    haskey(cli, flag) || return default
    v = tryparse(T, strip(cli[flag]))
    (v === nothing || (v isa AbstractFloat && !isfinite(v))) &&
        error("$flag needs $what, got '$(cli[flag])'")
    min_value === nothing || v >= min_value ||
        error("$flag needs $what of at least $min_value, got '$(cli[flag])'")
    max_value === nothing || v <= max_value ||
        error("$flag needs $what of at most $max_value, got '$(cli[flag])'")
    return v
end

"""
    options_from(cli) -> NamedTuple

Every option the driver reads, typed and checked, in one place and in the order the usage
lists them, so the first bad value on a command line is the one reported. Nothing here
defaults a bad value: see `option_number`.

Which background exclusion is in force decides numbers that leave this repository, and the
default is the registered one — but a default that is never stated cannot be audited,
so the result records whether the choice was typed or inherited, and whether the mode it
names is the registered primary at all, and the driver prints both before it works.
"""
function options_from(cli::AbstractDict)
    # FIRST, because it is the one option with no default at all, and a run that is going to
    # be refused should be refused before anything else about the command line is judged.
    categories = option_string(cli, "--categories", "")
    isempty(categories) &&
        error("no --categories file was given, and this step has no default for what is in " *
              "it. The category each recorded outcome rule speaks for, the categories the " *
              "grid counts are exported for and the phylum the standardisation shares are " *
              "quoted for are SELECTIONS made from this data set's own estimates, not " *
              "constants of this software: a default would publish a choice as though the " *
              "software had made it, and could not be told apart from a selection somebody " *
              "checked. Pass --categories <path>; docs/DATA.md names the file and the " *
              "archived data set carries it")
    config = option_string(cli, "--config", "primary")
    bg_exclusion = option_string(cli, "--bg-exclusion", PRIMARY_MODE)
    bg_exclusion in EXCLUSION_MODES ||
        error("--bg-exclusion must be one of " * join(EXCLUSION_MODES, ", ") *
              ", got '$bg_exclusion'")
    # The pool only exists for the level that moves it. A pool handed to any other mode would
    # be read by nothing and change nothing, so the run would be a different one from the
    # command that was typed — the failure this parser exists to stop.
    movein_pool = option_string(cli, "--movein-pool", "")
    isempty(movein_pool) || bg_exclusion == MOVEIN_MODE ||
        error("--movein-pool is read only by --bg-exclusion $MOVEIN_MODE, and " *
              "--bg-exclusion $bg_exclusion was asked for")
    nboot = option_number(cli, "--boot", Int, 2000;
                          what = "a whole number of replicates", min_value = 0)
    seed = option_number(cli, "--seed", Int, 20260913; what = "a whole number")
    legacy = config == "legacy"
    # The `legacy` configuration carries the exclusive winner-chemistry column; every other
    # configuration carries the non-exclusive flag union the chemistry rows describe.
    c = option_string(cli, "--chem-col", legacy ? "winner" : "flags")
    chem_col = c == "winner" ? :chemistry : c == "flags" ? :chemistry_flags :
        error("--chem-col must be flags or winner, got '$c'")
    calls = option_string(cli, "--calls",
                          DataInputs.input_path("habitat_calls"))
    scan = option_string(cli, "--scan",
                         DataInputs.input_path("background_scan_table"))
    census = option_string(cli, "--census",
                           DataInputs.input_path("census_habitat"))
    recall = option_string(cli, "--recall",
                           DataInputs.input_path("scan_recall"))
    threshold_record = option_string(cli, "--threshold-record",
                                     DataInputs.input_path("scan_threshold"))
    acceptance = option_string(cli, "--acceptance",
                               DataInputs.input_path("scan_acceptance"))
    anchor_reference = option_string(cli, "--anchor-reference", "")
    out = option_string(cli, "--out",
                        DataInputs.input_path("enrichment_table"))
    summary = option_string(cli, "--summary",
                            DataInputs.input_path("arm_summary"))
    # The anchor's rows go next to the run's rows unless asked otherwise, so a run that
    # redirects --out cannot leave its anchor in the repository by accident.
    anchor_out = option_string(cli, "--anchor-out",
                               (b = splitext(out); b[1] * "_anchor" * (isempty(b[2]) ? ".csv" : b[2])))
    # The plan's own bounds. No default: a bound this script invented would be this script's
    # opinion, and a bound it defaulted to would be indistinguishable from one that was checked.
    plan_bg_size_below = option_number(cli, "--plan-bg-size-below", Int, nothing;
                                       what = "a whole number of taxa", min_value = 1)
    plan_bg_phylum_pct_max =
        option_number(cli, "--plan-bg-phylum-pct-max", Float64, nothing;
                      what = "a percentage", min_value = 0.0, max_value = 100.0)
    return (categories = categories, config = config, bg_exclusion = bg_exclusion,
            bg_exclusion_source = haskey(cli, "--bg-exclusion") ? :explicit : :default,
            bg_exclusion_role = bg_exclusion == PRIMARY_MODE ? :primary : :sensitivity,
            movein_pool = movein_pool,
            nboot = nboot, seed = seed, legacy = legacy, chem_col = chem_col,
            calls = calls, scan = scan, census = census, recall = recall,
            threshold_record = threshold_record, acceptance = acceptance,
            anchor_reference = anchor_reference, out = out, summary = summary,
            anchor_out = anchor_out, plan_bg_size_below = plan_bg_size_below,
            plan_bg_phylum_pct_max = plan_bg_phylum_pct_max)
end

# Only when this file is the program. test/test_habitat_enrichment.jl includes it, and must
# get the functions without the command line being read at all.
const CLI = abspath(PROGRAM_FILE) == (@__FILE__) ?
    cli_or_exit(() -> parse_options(ARGS, OPTION_SPEC)) : Dict{String,String}()
if haskey(CLI, "--help") || haskey(CLI, "-h")
    print(USAGE)
    exit(0)
end
# Where every input of this step lives, by NAME. The layout is declared once, in
# src/ModuleDataInputs.jl, and documented once, in docs/DATA.md; the option defaults
# below are those names. Loaded after --help, which therefore still answers with
# nothing but this file.
include(joinpath(REPO, "src", "ModuleDataInputs.jl"))
# Typed option values, checked at the same point, so a bad value also stops the program
# before anything is loaded or read. `nothing` when this file is included.
const OPT = abspath(PROGRAM_FILE) == (@__FILE__) ? cli_or_exit(() -> options_from(CLI)) :
    nothing

using Pkg; Pkg.activate(REPO)
using DataFrames, CSV, Printf, Statistics, Random, Dates, SHA

include(joinpath(REPO, "src", "ModuleEnvironmentTools.jl")); using .EnvironmentTools
include(joinpath(REPO, "src", "ModuleHabitatStats.jl"));     using .HabitatStats
include(joinpath(REPO, "src", "ModuleOutcomeRules.jl"));     using .OutcomeRules

# The producer of the scan table and of the three records beside it. Its program-file guard
# keeps its driver — the downloads, HMMER and every write — from running on include, so what
# this brings in is the record readers, the published/pending/stale file-name convention, the
# rule an acceptance row's verdict follows, and the exclusion predicate itself, which the
# any-hit mode calls rather than restating. See NO RECORD PARSER OF ITS OWN in the header.
module ScanRecords end
Base.include(ScanRecords, joinpath(REPO, "scripts", "background_pqq_scan.jl"))

# Provenance records that are not options: the omnicrobe cache manifest and the search
# space's UniProtKB release record, which the scan writes beside its downloads.
const MANIFEST = DataInputs.input_path("relation_manifest")
const UPRELF   = DataInputs.input_path("uniprot_release")

# The background frame's own provenance. The taxon pool is a SEEDED subsample of the fetched
# reference-proteome taxa, and the seed is the one value that makes the frame reproducible. It
# is read from the generator's own source, where the constant is located and cited by line -
# never typed here. The pool and the frame source are the two tables that generator writes,
# named in the registry like every other input.
const POOL_FILE  = DataInputs.input_path("pool_taxa")
const FRAME_FILE = DataInputs.input_path("pool_background")
const POOL_GEN   = joinpath(REPO, "scripts", "prokaryote_background.jl")

# The anchor's specification. It is not a free choice: these are the levels that reproduce
# the earlier reference run, and an anchor built any other way answers a different question.
const ANCHOR_CONFIG   = "legacy"
const ANCHOR_CHEM_COL = :chemistry
# ---------------------------------------------------------------- the category selection
# THE FOUR SELECTIONS THIS SCRIPT DOES NOT MAKE. Each of these roles is filled by a category
# name — by one or more for the grid counts, and by a phylum name for the last. Every one of
# them was chosen with this data set's estimates already in hand, so none of them is a
# constant here: the criteria file fills them, this script declares only the ROLES, and it
# refuses a file that does not fill exactly these.
#
# `repeatable` is the one role that may name more than one category. `keyed` says whether the
# role's name reaches the summary as a key prefix, which is what makes the emitted families
# follow the selection rather than the other way round.
const SELECTION_ROLES = (
    (role = "depletion_tier", repeatable = false, keyed = true,
     what = "the category whose depletion the recorded robustness-tier rule decides"),
    (role = "genus_clause", repeatable = false, keyed = true,
     what = "the category whose genus clause and category-drop the recorded rules decide"),
    (role = "sharing", repeatable = false, keyed = true,
     what = "the category whose comparator-to-estimate sharing ratio the rules decide"),
    (role = "grid_count", repeatable = true, keyed = false,
     what = "the categories whose intervals the figure step counts across every " *
            "registered specification"),
    (role = "standardisation_phylum", repeatable = false, keyed = true,
     what = "the phylum whose foreground and background shares are quoted beside the " *
            "standardisation"),
)
# The roles whose categories the recorded rules speak for BY NAME, so that they need no
# standing of their own in the category partition. The governed set is exactly the categories
# the supplied rules are about — which is the only definition of it a reader can check against
# the file in front of them, and it is not a set this script may hold an opinion about.
const GOVERNING_ROLES   = ("depletion_tier", "genus_clause", "sharing")
const SELECTION_COLUMNS = (:kind, :role, :value, :note)
const SELECTION_KINDS   = ("select", "expect")
const SELECTION_EXPECT  = "selections"

"""
    key_slug(name) -> String

The summary-key prefix a category or phylum name is filed under: lower case, with every
character that is not a letter or a digit removed. One rule, written once, because the prefix
has to be derivable by whatever READS the table as well as by what writes it. A name with no
letter or digit in it, or one that would start a key with a digit, is an error rather than a
key nobody can parse.
"""
function key_slug(name::AbstractString)
    out = String(filter(c -> ('a' <= c <= 'z') || ('0' <= c <= '9'), lowercase(String(name))))
    isempty(out) &&
        error("'$name' has no letter or digit in it, so it cannot name a key family")
    '0' <= out[1] <= '9' &&
        error("'$name' begins with a digit, so its key family would not be a readable key")
    return out
end

"""
    read_category_selection(path) -> NamedTuple

The category selection one run is held to: one value per role of `SELECTION_ROLES`, the key
prefix each keyed role's name derives to, the governed set the rule roles make, and the digest
of the bytes all of that came from.

WHY THE FILE ALSO STATES HOW MANY SELECTIONS IT HOLDS, AND WHICH ROLES. Neither is measured on
anything, and both nevertheless have to reach the run — for the reason the replication gates'
targets file states its own size. A file that lost a row is a run that filed one category's
numbers under nobody's name, and "which categories this analysis is about" is not something a
run may discover from the file it was handed. The declared roles are compared with the file's
in BOTH directions, so a file written for another version of this script — one that names a
role nobody fills, or leaves a role nobody declares unfilled — is a refusal and not a quieter
success.

WHY EVERY ROW MUST CARRY A NOTE. A category name with no record of where it came from cannot
be told apart from one transcribed out of this analysis's own output, and which of those it is
happens to be the single most important fact about a selected set. The note is not parsed and
gates nothing; its absence is a refusal.
"""
function read_category_selection(path::AbstractString)
    isfile(path) ||
        error("--categories $path does not exist. It names the category each recorded " *
              "outcome rule speaks for, the categories the grid counts are exported for " *
              "and the phylum the standardisation shares are quoted for; this step ships " *
              "no default for any of them")
    cell(x) = x === missing || x === nothing ? "" : String(strip(string(x)))
    bytes = read(path)
    tbl = CSV.read(bytes, DataFrame; types = String, pool = false)
    for c in SELECTION_COLUMNS
        c in propertynames(tbl) ||
            error("--categories $path needs the columns " *
                  join(String.(SELECTION_COLUMNS), ", ") * " and has no '$c'")
    end
    declared = Dict{String,Any}(r.role => r for r in SELECTION_ROLES)
    got   = Dict{String,Vector{String}}()
    notes = Dict{String,Vector{String}}()
    expect = nothing
    nsel = 0
    for i in 1:nrow(tbl)
        kind = cell(tbl.kind[i]); role = cell(tbl.role[i])
        value = cell(tbl.value[i]); note = cell(tbl.note[i])
        kind in SELECTION_KINDS ||
            error("'$kind' is not a row kind of $path, which holds 'select' rows and one " *
                  "'expect' row")
        if kind == "select"
            nsel += 1
            haskey(declared, role) ||
                error("'$role' is not a selection this run makes, which fills " *
                      join(sort(collect(keys(declared))), "; ") * ". A file written for " *
                      "another version of this script names a role nobody reads")
            isempty(value) &&
                error("the '$role' row of $path names nothing, and this step has no " *
                      "default to fall back to")
            isempty(note) &&
                error("the '$role' row of $path carries no note, so nothing in the file " *
                      "says where that selection came from")
            (haskey(got, role) && !declared[role].repeatable) &&
                error("'$role' is filled twice in $path, so one rule's category is named " *
                      "in two places and only one of them decides the run")
            haskey(got, role) || (got[role] = String[]; notes[role] = String[])
            value in got[role] && error("'$role' names '$value' twice in $path")
            push!(got[role], value); push!(notes[role], note)
        else
            role == SELECTION_EXPECT ||
                error("an 'expect' row counts '$SELECTION_EXPECT', not '$role'")
            expect === nothing || error("$path states its own size twice")
            v = tryparse(Int, value)
            (v === nothing || v < 0) &&
                error("the expected number of selections in $path is '$value'")
            expect = v
        end
    end
    unfilled = String[r.role for r in SELECTION_ROLES if !haskey(got, r.role)]
    isempty(unfilled) ||
        error("$path fills no " * join(unfilled, ", ") * " role, and this step ships no " *
              "default for any selection: " *
              join([r.role * " is " * r.what for r in SELECTION_ROLES
                    if r.role in unfilled], "; "))
    expect === nothing &&
        error("$path does not state how many selections it holds. An " *
              "'expect,$SELECTION_EXPECT,<n>' row is required: without it a file that lost " *
              "a row is a run that selected fewer categories and said so nowhere")
    expect == nsel || error("$path says it holds $expect selections and holds $nsel")
    sole(role) = only(got[role])
    keyed = String[r.role for r in SELECTION_ROLES if r.keyed]
    slugs = Dict{String,String}(r => key_slug(sole(r)) for r in keyed)
    # Two names that reduce to one prefix would file two selections' numbers under one family
    # of keys, and the second would overwrite the first leaving no trace in the table.
    owner = Dict{String,String}()
    for r in keyed
        s = slugs[r]
        haskey(owner, s) &&
            error("the '" * owner[s] * "' and '" * r * "' selections ('" * sole(owner[s]) *
                  "', '" * sole(r) * "') both reduce to the key prefix '" * s *
                  "', so one family of keys would carry both")
        owner[s] = r
    end
    return (path = String(path),
            digest = bytes2hex(sha256(bytes)),
            depletion_tier = sole("depletion_tier"),
            genus_clause = sole("genus_clause"),
            sharing = sole("sharing"),
            grid_count = Tuple(got["grid_count"]),
            standardisation_phylum = sole("standardisation_phylum"),
            governed = Tuple(sole(r) for r in GOVERNING_ROLES),
            key = NamedTuple{Tuple(Symbol(r) for r in keyed)}(Tuple(slugs[r] for r in keyed)),
            notes = NamedTuple{Tuple(Symbol(r.role) for r in SELECTION_ROLES)}(
                        Tuple(join(notes[r.role], " | ") for r in SELECTION_ROLES)))
end

# Read here, under the same guard the command line is read under, so that a criteria file that
# is missing, malformed, or written for another version of this script stops the run before
# any data file is opened. `nothing` when this file is included, exactly like OPT.
const SELECTION = OPT === nothing ? nothing :
    cli_or_exit(() -> read_category_selection(OPT.categories))
# The columns scripts/background_pqq_scan.jl appends for the whole-frame family membership.
# An older scan table has none of them, which is a fact about the input worth emitting.
const SCAN_APPENDED_COLS = (:n_ipr017512_entries, :n_ipr017512_ge500, :in_ipr017512_ge500,
                            :ipr017512_outside_prefilter, :hit_ipr017512, :n_hits_ln_or_ca)
# A disclosure list long enough to act on and short enough to read.
const MAX_LISTED_TAXA = 25
# The columns an anchor row is compared with the reference run's row on, and the tolerance.
# The comparison is of two computations of the same quantity, not of two measurements, so the
# tolerance is float-summation noise and nothing more.
const REF_COMPARE_COLS = (:a, :b, :n1, :n2, :log2_rr, :ci_lo, :ci_hi,
                          :log2_rr_std, :ci_lo_std, :ci_hi_std)
const REF_ROW_TOL = 1e-9

sv(x) = ismissing(x) ? "" : String(strip(string(x)))
hascol(df::AbstractDataFrame, name::Symbol) = name in propertynames(df)
# A cell of a table that may or may not hold a number. Anything that is not one reads as
# NaN, which every rule in the module treats as "not computed" rather than as a value.
fnum(x) = x === missing || x === nothing ? NaN :
          x isa Real ? float(x) : something(tryparse(Float64, sv(x)), NaN)

# ---------------------------------------------------------------- summary table
# Every reported quantity is registered here and nowhere else, so the closing report and any
# downstream text step read the same numbers the CSV carries.
const SUMMARY = Pair{String,Any}[]
rec!(k::AbstractString, v) = (push!(SUMMARY, String(k) => v); v)
function look(k::AbstractString)
    i = findfirst(p -> first(p) == k, SUMMARY)
    i === nothing ? missing : last(SUMMARY[i])
end
tally(labels) = (d = Dict{String,Int}();
                 for l in labels; d[l] = get(d, l, 0) + 1; end; d)
function tally_str(d::AbstractDict; limit::Int = 0)
    ks = sort(collect(keys(d)), by = k -> (-d[k], k))
    limit > 0 && (ks = first(ks, limit))
    join((@sprintf("%s:%d", k, d[k]) for k in ks), ";")
end
joined(xs) = join(xs, ";")
listed(xs) = isempty(xs) ? "" :
             length(xs) <= MAX_LISTED_TAXA ? joined(xs) :
             joined(first(xs, MAX_LISTED_TAXA)) * ";+$(length(xs) - MAX_LISTED_TAXA) more"
# The scan's key/value records are read with the producer's own reader, so a change of that
# format cannot leave this script parsing the old one.
read_kv(path) = ScanRecords.read_keyvals(path)

# ---------------------------------------------------------------- the background pool
# The background pool: a SEEDED subsample of the fetched reference-proteome taxa, so the seed is
# part of this step's provenance and not an implementation detail of a sibling script. It is
# read from the generator's own constant, located in its source and cited by line. Nothing is
# typed here, and a value that cannot be read is recorded as absent with the state saying why,
# because a default would be this script inventing the seed of a subsample it did not take.
"""
    pool_constant(src, name) -> Union{Nothing,NamedTuple}

The integer of a `const <name> = <n>` line in the pool generator, with the 1-based line it was
read from. `nothing` unless exactly one line defines it: a constant that cannot be located
unambiguously is not provenance.
"""
function pool_constant(src::AbstractString, name::AbstractString)
    isfile(src) || return nothing
    pat = Regex("^\\s*const\\s+" * name * "\\s*=\\s*(\\d+)\\b")
    hits = NamedTuple[]
    for (i, line) in enumerate(eachline(src))
        m = match(pat, line)
        m === nothing || push!(hits, (value = parse(Int, m[1]), line = i))
    end
    return length(hits) == 1 ? only(hits) : nothing
end
# One read of every file, digest and row count from the same bytes, for the reason spelled out
# at read_scan_table: a second read can answer differently, and then the digest describes bytes
# the count did not come from.
function pool_file!(rec, prefix::AbstractString, path::AbstractString)
    isfile(path) || return nothing
    bytes = read(path)
    tbl = CSV.read(bytes, DataFrame)
    rec("$(prefix)_file", basename(path))
    rec("$(prefix)_sha256", bytes2hex(sha256(bytes)))
    rec("$(prefix)_taxa", hascol(tbl, :taxid) ? length(unique(tbl.taxid)) : missing)
    return nrow(tbl)
end
"""
    background_pool_provenance!(rec; generator, pool, frame) -> nothing

Record the background pool's provenance through `rec` (`rec!` in a run): where the seed came
from and what it is, the subsample size and the pool size above which it applies, the
generator's digest, and the file name, digest and taxon count of the pool and of the frame
source. `background_pool_provenance_state` says how the seed was found: `generator_source`
(one `const SEED` line in the generator), `generator_absent` (no generator file) or
`seed_not_found` (a generator with no single such line). The defaults are this repository's
generator and the registry's two tables; the keywords exist so each state can be tested.
"""
function background_pool_provenance!(rec; generator::AbstractString = POOL_GEN,
                                     pool::AbstractString = POOL_FILE,
                                     frame::AbstractString = FRAME_FILE)
    src = pool_constant(generator, "SEED")
    sub = pool_constant(generator, "SUBSAMPLE")
    cap = pool_constant(generator, "MAX_TAXA")
    seed, origin, state =
        src !== nothing     ? (src.value,
                               "scripts/prokaryote_background.jl:$(src.line) (const SEED)",
                               "generator_source") :
        !isfile(generator)  ? (missing, "", "generator_absent") :
                              (missing, "", "seed_not_found")
    rec("background_pool_provenance_state", state)
    rec("background_pool_seed", seed)
    rec("background_pool_seed_origin", origin)
    rec("background_pool_subsample_size", sub === nothing ? missing : sub.value)
    rec("background_pool_subsample_above", cap === nothing ? missing : cap.value)
    if isfile(generator)
        rec("background_pool_generator", "scripts/prokaryote_background.jl")
        rec("background_pool_generator_sha256", bytes2hex(sha256(read(generator))))
    end
    npool = pool_file!(rec, "background_pool", pool)
    pool_file!(rec, "background_frame_source", frame)
    # A seed only means something if the subsample was actually taken. The pool on disk having
    # exactly SUBSAMPLE rows is the check, and it is recorded rather than asserted.
    rec("background_pool_subsample_applied",
        npool === nothing || sub === nothing ? missing : npool == sub.value)
    return nothing
end

# ---------------------------------------------------------------- weights and units
"""
    weight_map(values) -> Dict{String,AbstractVector}

Turn a `"Category:w;Category:w"` column into one per-unit weight vector per category.
A vector whose entries are all 0 or 1 is narrowed to a `BitVector`, so an undivided call
sums to an integer and takes the exact Bool path through the estimators; a fractional
call (a tie split across the tied categories) stays `Float64`.
"""
function weight_map(values::AbstractVector)
    n = length(values)
    acc = Dict{String,Vector{Float64}}()
    for (i, v) in enumerate(values)
        s = sv(v); isempty(s) && continue
        for part in split(s, ';')
            isempty(strip(part)) && continue
            j = findlast(==(':'), part)
            j === nothing && error("malformed weight entry '$part'")
            cat = String(strip(part[1:prevind(part, j)]))
            w = parse(Float64, strip(part[nextind(part, j):end]))
            haskey(acc, cat) || (acc[cat] = zeros(n))
            acc[cat][i] += w
        end
    end
    out = Dict{String,AbstractVector}()
    for (k, v) in acc
        out[k] = all(x -> x == 0.0 || x == 1.0, v) ? BitVector(x == 1.0 for x in v) : v
    end
    return out
end

expand(v::AbstractVector, idx) = v[idx]

# The census arm carries a genus from the taxonomic lineage; the reference-proteome pool
# carries only an organism name, so the genus unit would otherwise exist on one side only —
# which is the very asymmetry the symmetric collapse was introduced to remove. The name is
# therefore parsed when no lineage genus is present, and the run reports how many labels
# came from each route so the weaker provenance of the background genera stays visible.
const NOT_A_GENUS = ("uncultured", "unclassified", "unidentified", "candidate",
                     "endosymbiont", "symbiont")
"""
    genus_of(genus, organism) -> String

The lineage genus when there is one, otherwise the genus read off the organism name:
parenthetical synonyms and bracketed provisional placements are stripped, leading
non-taxonomic words are skipped, `Candidatus` keeps the following word, and anything that
does not begin with a capital (an unnamed lineage, an "endosymbiont of ...") returns blank
and is dropped from the genus unit rather than invented.
"""
function genus_of(genus, organism)
    g = sv(genus); isempty(g) || return g
    s = sv(organism); isempty(s) && return ""
    i = findfirst(==('('), s)
    i === nothing || (s = String(strip(s[1:prevind(s, i)])))
    s = replace(s, '[' => "", ']' => "")
    toks = split(s); isempty(toks) && return ""
    k = 1
    while k <= length(toks) && lowercase(String(toks[k])) in NOT_A_GENUS; k += 1; end
    k > length(toks) && return ""
    g = String(toks[k])
    lowercase(g) == "candidatus" && k < length(toks) &&
        (g = "Candidatus " * String(toks[k + 1]))
    (isempty(g) || !isuppercase(first(g))) && return ""
    return g
end

# ---------------------------------------------------------------- contrasts
# Both arms are prokaryotic (see PROKARYOTES ONLY in the header). The two names are one
# constant because they are read off two different columns: the census lineage's
# `kingdom` for taxa the census registered, and `domain` for the frame taxa the move-in
# level promotes, whose `kingdom` cell is empty.
const PROKARYOTE_NAMES = ("Bacteria", "Archaea")
const LN_VS_BG  = "Ln vs background"
const ALL_VS_BG = "all PQQ-ADH carriers vs background"
const NON_VS_BG = "non-Ln PQQ-ADH carriers vs background"
# The disclosed second comparator: carriers whose census sequences read calcium-class and
# never lanthanide-class. It stands beside NON_VS_BG and never in its place.
const CA_VS_BG  = "Ca-class carriers vs background"
# The contrasts the recorded rules are registered on. A row counter a rule's re-check reads
# must be taken over these alone, or a disclosed arm's thin row changes a registered number.
const REGISTERED_CONTRASTS = (LN_VS_BG, ALL_VS_BG, NON_VS_BG)
# The rows a category's standing is read off: the Ln arm against the background, on the two
# units whose denominator is a taxon.
const PRIMARY_UNITS = ("taxa", "taxa (chemistry)")

"""
    analyse(calls, config, chem_col, excluded; ca_class_taxids, moved_in) -> NamedTuple

Build the arms, the weight vectors and the result rows for ONE specification: a
configuration of the taxon-call table, a chemistry column, the set of frame taxids the
background exclusion removes, the census taxids whose classified sequences are
calcium-class only (`nothing` when the census table could not be read, which leaves the
disclosed calcium arm out rather than silently equating it to the non-Ln arm), and
`moved_in`, the frame taxids the move-in level hands to the FOREGROUND. `moved_in` is empty
for every other specification, which is why every other specification is unchanged by it.

This is the only arm builder in the file, and it is called twice — once for the run's own
specification and once for the anchor. That is the point: an anchor computed by a second
copy of this logic would drift away from the primary path, and the anchor's whole job is to
be the same measurement on earlier settings.

Everything it returns is derived from those arguments, so nothing about the run's
configuration can leak into the anchor. The masks come back too, because the coverage,
drop-bucket and evidence sections need the candidate sets, not just the arms.
"""
function analyse(calls::DataFrame, config::AbstractString, chem_col::Symbol,
                 excluded::Set{Int}; ca_class_taxids::Union{Set{Int},Nothing},
                 moved_in::Set{Int} = Set{Int}())
    df = calls[sv.(calls.config) .== config, :]
    nrow(df) > 0 ||
        error("no rows for config $config; have: " * joined(sort(unique(sv.(calls.config)))))
    hascol(df, chem_col) || error("the call table has no column $chem_col")

    king    = sv.(df.kingdom)
    prok    = in.(king, Ref(PROKARYOTE_NAMES))
    ht      = sv.(df.habitat_type)
    # A usable type call: a habitat that is not an isolation source. Under the primary rule
    # isolation terms never reach the tally, so this only bites on the legacy configuration.
    usable  = (ht .!= "") .& .!is_isolation_source.(ht)
    n_ln    = coalesce.(df.n_ln_seqs, 0)
    n_cen   = coalesce.(df.n_census_seqs, 0)
    in_cen  = coalesce.(df.in_census, false)
    in_pool = coalesce.(df.in_pool, false)

    # THE FRAME AS THE SCAN SAW IT, before any taxon is reassigned. Everything the scan was
    # asked about is here, which is what the scan's own coverage is compared against; the
    # move-in level reassigns part of this set rather than changing which taxa were scanned.
    frame_scanned = in_pool .& .!in_cen

    # THE MOVE-IN LEVEL, AND THE COLUMN THAT MAKES IT LOOK INERT. Frame taxa carry their
    # domain in `domain`; `kingdom` is populated from the census lineage and is EMPTY for
    # every taxon the census never registered. A prokaryote test read off `kingdom` alone
    # therefore drops every moved taxon at the mask, leaves n1 where it was, and reports
    # "this level moves nothing" — a variant that cannot break a result and so tests nothing.
    # The moved taxa are tested on either column, and the guard below refuses a run in which
    # the pool named frame taxa and none of them survived that test.
    dom        = hascol(df, :domain) ? sv.(df.domain) : fill("", nrow(df))
    named_in_pool = frame_scanned .& in.(df.taxid, Ref(moved_in))

    # PHYLUM IS THE STRATUM, so a taxon without one cannot enter the standardised contrast:
    # it belongs to no stratum, and bucketing it under "Other phyla" would state a lineage
    # where there is none. NCBI's catch-all placeholders - `bacterium` (1869227), `uncultured
    # bacterium` (77133), `uncultured archaeon` (115547) - are such taxa: they carry no phylum,
    # and their habitat relations are an aggregate over every unrelated deposit filed under
    # that label, not one organism's ecology. A taxonomy that cannot resolve those taxids
    # drops them by accident and one that resolves every taxid admits them, so the exclusion
    # is a rule here rather than a side effect of how complete the taxonomy is. It is applied
    # at the assembly point rather than in a caller, so every arm this builder returns is
    # subject to it, and the guard below stops the run, with a count per arm, if a phylum-less
    # taxon gets past.
    placed     = sv.(df.phylum) .!= ""

    moved      = named_in_pool .& (prok .| in.(dom, Ref(PROKARYOTE_NAMES))) .& placed

    ln_census  = in_cen .& (n_ln .> 0) .& prok .& placed
    ln_cand    = ln_census .| moved
    all_cand   = in_cen .& (n_cen .> 0) .& prok .& placed
    nonln_cand = all_cand .& (n_ln .== 0)
    ca_cand    = ca_class_taxids === nothing ? falses(nrow(df)) :
                 nonln_cand .& in.(df.taxid, Ref(ca_class_taxids))
    # A moved taxon leaves the frame in the same expression that gives it to the foreground,
    # so the background it was drawn from cannot still hold it whatever the exclusion mode
    # does. The comparator arms are deliberately NOT rebuilt: their membership rules are
    # rules about a taxon's CENSUS sequences (which lanthanide or calcium class it carries),
    # and a taxon with no census sequence has no standing in them. They still lose the moved
    # taxa from their shared background, which is the level's doing and is measured as such.
    frame      = frame_scanned .& .!moved
    bg_keep    = frame .& .!in.(df.taxid, Ref(excluded)) .& placed

    # THE TWO INVARIANTS OF A MOVED TAXON, HELD HERE RATHER THAN ONLY IN A TEST, because this
    # is where the arms are assembled and a later copy of the logic would be the one that
    # drifts. A taxon on both sides would put an arm against a background containing itself,
    # and a taxon entering the foreground twice would inflate n1 while every arm total still
    # looked ordinary. Both are stated as errors, not as assertions that a future build flag
    # could strip out of the one place they are checked.
    !any(ln_cand .& frame) ||
        error("the move-in level left $(count(ln_cand .& frame)) taxa in the foreground AND " *
              "in the background frame; a contrast against a background that still holds " *
              "the moved taxa compares an arm with part of itself")
    allunique(df.taxid[ln_cand]) ||
        error("the foreground holds a taxon more than once, so n1 counts it twice; the " *
              "move-in level adds taxa by mask and must never add a row")
    # The exclusion above is only honest if it actually holds on every assembled arm: a
    # phylum-less taxon reaching a contrast would be standardised against strata it is not in.
    all(placed[ln_cand]) && all(placed[all_cand]) && all(placed[bg_keep]) ||
        error("an assembled arm still holds a taxon with no phylum, so the standardised " *
              "contrast would weight it into a stratum it does not belong to: " *
              "$(count(.!placed[ln_cand])) foreground, $(count(.!placed[all_cand])) " *
              "all-carrier, $(count(.!placed[bg_keep])) background")
    count(ln_cand) == count(ln_census) + count(moved) ||
        error("the foreground is not the census carriers plus the moved taxa " *
              "($(count(ln_cand)) against $(count(ln_census)) + $(count(moved))): the two " *
              "sets overlap, and n1 counts the overlap twice")
    !(count(named_in_pool) > 0 && count(moved) == 0) ||
        error("the move-in pool names $(count(named_in_pool)) frame taxa and moved none of " *
              "them: the frame carries its domain in `domain` and leaves `kingdom` empty, " *
              "so a prokaryote test taken off `kingdom` alone silently returns the unmoved " *
              "specification under the name of the level that was asked for")

    fg  = df[ln_cand    .& usable, :]
    bgm = df[bg_keep    .& usable, :]
    ac  = df[all_cand   .& usable, :]
    nl  = df[nonln_cand .& usable, :]
    ca  = df[ca_cand    .& usable, :]

    fgw = weight_map(fg.type_weights)
    bgw = weight_map(bgm.type_weights)
    acw = weight_map(ac.type_weights)
    nlw = weight_map(nl.type_weights)
    caw = weight_map(ca.type_weights)

    fg_phy = sv.(fg.phylum); bg_phy = sv.(bgm.phylum)
    ac_phy = sv.(ac.phylum); nl_phy = sv.(nl.phylum); ca_phy = sv.(ca.phylum)

    # per-sequence unit: repeat each taxon once per lanthanide-type census sequence it carries
    #
    # A MOVED TAXON HAS NONE, AND IS NOT GIVEN ONE. Not being in the census is what put it in
    # the background to begin with, so it carries no census sequence and contributes NO
    # per-sequence unit: under the move-in level the sequences row is the census foreground
    # against the post-move background, and the run emits the count of moved taxa this leaves
    # out rather than leaving it to be inferred from an n1 that grew and an n_sequences that
    # did not. Giving each moved taxon a flat one, or the hit count its pool names, would put
    # a family hit found in a frame proteome and a census sequence into one denominator, and
    # a denominator made of two different measurements measures neither.
    seq_idx = reduce(vcat, [fill(i, coalesce(fg.n_ln_seqs[i], 0)) for i in 1:nrow(fg)];
                     init = Int[])
    sq_phy  = fg_phy[seq_idx]
    sqw     = Dict{String,AbstractVector}(k => expand(v, seq_idx) for (k, v) in fgw)

    # symmetric genus unit: BOTH arms collapsed, so n1 and n2 are both genus counts
    fg_gen = genus_of.(fg.genus, fg.organism)
    bg_gen = genus_of.(bgm.genus, bgm.organism)
    gfg = collapse_to_genus(fg_gen, fg_phy, fgw)
    gbg = collapse_to_genus(bg_gen, bg_phy, bgw)

    fgc = flag_weights(fg[!, chem_col];  splitcombo = true)
    bgc = flag_weights(bgm[!, chem_col]; splitcombo = true)
    acc = flag_weights(ac[!, chem_col];  splitcombo = true)
    nlc = flag_weights(nl[!, chem_col];  splitcombo = true)
    cac = flag_weights(ca[!, chem_col];  splitcombo = true)

    rows = NamedTuple[]
    add!(args...; kw...) = append!(rows, contrast_rows(args...; kw...))
    add!(fg_phy, fgw, bg_phy, bgw;        contrast = LN_VS_BG,  unit = "taxa")
    add!(sq_phy, sqw, bg_phy, bgw;        contrast = LN_VS_BG,  unit = "sequences")
    add!(gfg.phylum, gfg.weights, gbg.phylum, gbg.weights;
                                          contrast = LN_VS_BG,  unit = "genera")
    add!(ac_phy, acw, bg_phy, bgw;        contrast = ALL_VS_BG, unit = "taxa")
    add!(nl_phy, nlw, bg_phy, bgw;        contrast = NON_VS_BG, unit = "taxa")
    add!(fg_phy, fgc, bg_phy, bgc;        contrast = LN_VS_BG,  unit = "taxa (chemistry)")
    add!(ac_phy, acc, bg_phy, bgc;        contrast = ALL_VS_BG, unit = "taxa (chemistry)")
    add!(nl_phy, nlc, bg_phy, bgc;        contrast = NON_VS_BG, unit = "taxa (chemistry)")
    # The disclosed arm goes LAST, so adding it cannot move a row of the registered panels.
    # An empty arm is left out rather than emitted as a panel of zeros: a contrast with no
    # units on one side is not a measurement of anything.
    if nrow(ca) > 0
        add!(ca_phy, caw, bg_phy, bgw;    contrast = CA_VS_BG,  unit = "taxa")
        add!(ca_phy, cac, bg_phy, bgc;    contrast = CA_VS_BG,  unit = "taxa (chemistry)")
    end

    # n1 IS THE CENSUS ARM PLUS THE MOVED ANNOTATED TAXA, AND NOTHING ELSE. Stated as a check
    # rather than as a comment: a level whose arms are right but whose n1 does not grow by the
    # taxa it moved has silently dropped them somewhere between the mask and the arm, and the
    # rows it writes would be the unmoved specification wearing this level's name.
    nrow(fg) == count(ln_census .& usable) + count(moved .& usable) ||
        error("n1 is $(nrow(fg)), which is not the $(count(ln_census .& usable)) annotated " *
              "census carriers plus the $(count(moved .& usable)) annotated taxa this level " *
              "moved in")

    return (config = String(config), chem_col = chem_col, excluded = excluded,
            df = df, usable = usable, prok = prok,
            moved = moved, named_in_pool = named_in_pool, frame_scanned = frame_scanned,
            n_moved = count(moved), n_moved_usable = count(moved .& usable),
            n_moved_census_sequences = sum(n_ln[moved]; init = 0),
            ln_census = ln_census,
            ln_cand = ln_cand, all_cand = all_cand, nonln_cand = nonln_cand,
            ca_cand = ca_cand,
            frame = frame, bg_keep = bg_keep, in_cen = in_cen, in_pool = in_pool,
            fg = fg, bgm = bgm, ac = ac, nl = nl, ca = ca,
            fgw = fgw, bgw = bgw, acw = acw, nlw = nlw, caw = caw,
            fg_phy = fg_phy, bg_phy = bg_phy, ac_phy = ac_phy, nl_phy = nl_phy,
            ca_phy = ca_phy,
            seq_idx = seq_idx, fg_gen = fg_gen, bg_gen = bg_gen, gfg = gfg, gbg = gbg,
            rows = rows, n1 = nrow(fg), n2 = nrow(bgm))
end

pick(rows, c, u, h) = (i = findfirst(r -> r.contrast == c && r.unit == u && r.habitat == h,
                                    rows);
                       i === nothing ? nothing : rows[i])
primary_rows(rows) = [r for r in rows
                      if r.contrast == LN_VS_BG && r.unit in PRIMARY_UNITS]
registered_rows(rows) = [r for r in rows if r.contrast in REGISTERED_CONTRASTS]

# The two readings a row carries, as a bound pair. A row that does not exist has no bounds
# at all, and NaN is how the module is told that: it answers :undecidable for a pair that is
# not an interval over the reals, which is the same answer an uncomputed estimate deserves.
row_bounds(r, reading::Symbol) =
    r === nothing ? (NaN, NaN) :
    reading === :std ? (r.ci_lo_std, r.ci_hi_std) :
    reading === :crude ? (r.ci_lo, r.ci_hi) :
    error("unknown interval reading $reading")
zero_status(r, reading::Symbol) = interval_excludes_zero(row_bounds(r, reading)...)
zero_std(r)   = zero_status(r, :std)
zero_crude(r) = zero_status(r, :crude)

"""
    interval_sets(rows, reading) -> (excludes, undecidable, labels, n_undecidable_rows)

The categories of one row panel, split by what their interval decided under one reading:
those whose interval excludes zero, those whose interval could not be computed, every
category the panel holds a row for at all, and how many rows were undecidable.

ONE FUNCTION, THREE PANELS. The run's own rows, the anchor's rows and the rows of the
earlier reference run are all read by this function. Rebuilding the category set inline for
the reference table would let a change to the unit selection or to the zero test in one copy
flip the anchor's verdict between "reproduces the reference" and "differs from it" without
any test failing.

The undecidable set is what keeps an ABSENCE from being reported as a LOSS: a category
missing from the excluding-zero set is missing either because its interval includes zero or
because it has none, and only the first is a measurement.
"""
function interval_sets(rows, reading::Symbol)
    prim = primary_rows(rows)
    status = [zero_status(r, reading) for r in prim]
    labels_with(s) = sort(unique([prim[i].habitat for i in eachindex(prim) if status[i] === s]))
    return (excludes = labels_with(:excludes_zero), undecidable = labels_with(:undecidable),
            labels = sort(unique([r.habitat for r in prim])),
            n_undecidable_rows = count(s -> s === :undecidable, status))
end

"""
    undecidable_categories(primary_sets, anchor_labels) -> Vector{String}

The categories the partition must hand to the rule as undecidable in the primary run: those
whose primary interval could not be computed, and every category the anchor holds a row for
while the primary panel holds none. `primary_sets` is an `interval_sets` result.

A category missing from the primary set is a LOSS only if its primary interval was computed
and includes zero. No row means no interval, which is not a measurement of anything, so a
category the primary panel simply lacks belongs here too; without it the partition reports an
absence as a result.
"""
undecidable_categories(primary_sets, anchor_labels) =
    sort(union(primary_sets.undecidable, setdiff(anchor_labels, primary_sets.labels)))

# ---------------------------------------------------------------- background exclusion
"""
    exclusion_taxids(scan, mode) -> Union{Set{Int},Nothing}

The frame taxids a mode removes from the background, or `nothing` when the scan table on
hand cannot answer the mode — a missing table, or an older one without the column the mode
reads. `nothing` is not an empty set: an unanswerable mode must not look like a mode that
excluded nobody.

The any-hit mode is the producer's `exclude_carriers`, called on the table's own taxa,
rather than a restatement of its predicate: two copies of the registered exclusion could
drift apart, and the one that analysis ran would not be the one the scan's suite tests.
"""
function exclusion_taxids(scan, mode::AbstractString)
    mode == "none" && return Set{Int}()
    scan === nothing && return nothing
    # The move-in level REMOVES exactly what the registered primary removes — it is that
    # exclusion plus a reassignment, so it shares the predicate rather than restating it, and
    # the taxa it moves leave the background through the arm builder, not through this set.
    (mode == "hmmsearch" || mode == MOVEIN_MODE) &&
        return Set(Int.(ScanRecords.exclude_carriers(DataFrame(taxid = scan.taxid),
                                                     scan).removed.taxid))
    mode == "lntype" &&
        return Set(Int.(scan.taxid[occursin.(LANTHANIDE_CLASS, sv.(scan.hit_classes))]))
    if mode == "lnca"
        hascol(scan, :n_hits_ln_or_ca) || return nothing
        return Set(Int.(scan.taxid[coalesce.(scan.n_hits_ln_or_ca, 0) .> 0]))
    end
    if mode == "ipr017512"
        hascol(scan, :in_ipr017512_ge500) || return nothing
        return Set(Int.(scan.taxid[coalesce.(scan.in_ipr017512_ge500, false)]))
    end
    error("no exclusion rule for mode '$mode'")
end

# What each mode actually does, in the same words the emitted key carries. No mode claims a
# relationship to another mode's set: whether one mode's set contains another's is a
# property of a scan table, measured per run in the sensitivity keys, not something a note
# can assert. An unknown mode is an error and not a fallback to the most permissive note,
# which would describe a run that did not happen.
function exclusion_note(mode::AbstractString)
    mode == "hmmsearch" && return "hmmsearch: every frame taxon with at least one family " *
        "hit at the calibrated threshold"
    mode == "lnca" && return "lnca: frame taxa with at least one hit whose active site " *
        "reads Lanthanide or Calcium, at the calibrated threshold"
    mode == "lntype" && return "lntype: frame taxa with a lanthanide-classed hit only"
    mode == "ipr017512" && return "ipr017512: frame taxa with a family-signature entry at " *
        "the length cut-off, from the whole-frame query"
    mode == "none" && return "none — every frame taxon kept"
    mode == MOVEIN_MODE && return "$MOVEIN_MODE: the hmmsearch exclusion, and the " *
        "lanthanide-type pool named by the move-in pool rule joins the foreground instead " *
        "of merely leaving the background"
    error("no note for exclusion mode '$mode'")
end

# The modes the run tallies for comparison while another one is in force. Every one of them
# answers "how large would the background be under this rule", which is a count over the
# frame of the run that is executing. The move-in level is not one of them: it moves taxa to
# the other arm, so its n2 is not a property of this run's frame and its n1 is not this run's
# n1 either. Reporting it in that tally would put a number under a key whose every sibling
# means something else; it is a run of its own, and the keys it emits say so.
const SENSITIVITY_MODES = Tuple(m for m in EXCLUSION_MODES if m != MOVEIN_MODE)

"""
    scan_inventory!(scan)

Emit what the scan table holds before any of it is used: its size, how many taxa carry a
hit of each kind, the family-membership columns, and which of the appended columns an
older table lacks.
"""
function scan_inventory!(scan::AbstractDataFrame)
    nh_scan = coalesce.(scan.n_hits, 0)
    hc_scan = sv.(scan.hit_classes)
    absent  = [String(c) for c in SCAN_APPENDED_COLS if !hascol(scan, c)]
    optcount(col, f) = hascol(scan, col) ? count(f, scan[!, col]) : missing
    rec!("scan_frame_taxa", nrow(scan))
    rec!("scan_taxa_with_any_hit", count(nh_scan .> 0))
    rec!("scan_taxa_with_lanthanide_hit", count(occursin.("Lanthanide", hc_scan)))
    # in_ipr017512 is family membership of the TAXON, from a whole-frame query;
    # hit_ipr017512 is the narrower reading, an annotation on one of the taxon's hits. Both
    # are emitted under their own names, because one key that can carry either reading is
    # worse than two keys that each carry one.
    rec!("scan_taxa_in_ipr017512", optcount(:in_ipr017512, x -> coalesce(x, false)))
    rec!("scan_taxa_in_ipr017512_ge500", optcount(:in_ipr017512_ge500, x -> coalesce(x, false)))
    rec!("scan_taxa_with_hit_annotated_ipr017512",
         optcount(:hit_ipr017512, x -> coalesce(x, false)))
    rec!("scan_taxa_with_ln_or_ca_hit", optcount(:n_hits_ln_or_ca, x -> coalesce(x, 0) > 0))
    rec!("scan_taxa_ipr017512_outside_prefilter",
         optcount(:ipr017512_outside_prefilter, x -> coalesce(x, false)))
    # A zero-hit row is not the same statement as "searched and nothing found": a taxon with
    # no candidate entry was never searched at all, and that is most of what a leftover
    # count read off the hit table alone cannot see.
    rec!("scan_taxa_never_searched", optcount(:n_candidates, x -> coalesce(x, 0) == 0))
    rec!("scan_columns_missing", joined(absent))
    return nothing
end

# ---------------------------------------------------------------- census classes
taxid_of(x) = ismissing(x) ? nothing :
              x isa Real ? Int(round(x)) : tryparse(Int, sv(x))

"""
    read_census_classes(path) -> (classes, euk_seqs, euk_taxa, ca_class_taxids)

The classes of each census taxon's classified sequences, the eukaryotic disclosure, and the
disclosed calcium arm's membership: taxa with a calcium-class census sequence and no
lanthanide-class one.

Without the census table the membership is `nothing`, which leaves the arm out — an arm
whose defining evidence is missing must not fall back to the registered comparator's
membership, or a disclosed line would silently duplicate the primary one — and the two
eukaryotic counts are `missing`.
"""
function read_census_classes(path::AbstractString)
    classes = Dict{Int,Set{String}}()
    isfile(path) || return (classes = classes, euk_seqs = missing, euk_taxa = missing,
                            ca_class_taxids = nothing)
    cen = CSV.read(path, DataFrame; select = [:classification, :kingdom, :TaxID])
    euk = (sv.(cen.classification) .== "Lanthanide") .& (sv.(cen.kingdom) .== "Eukaryota")
    euk_seqs = count(euk)
    euk_taxa = length(unique(collect(skipmissing(cen.TaxID[euk]))))
    for i in 1:nrow(cen)
        t = taxid_of(cen.TaxID[i]); t === nothing && continue
        push!(get!(classes, t, Set{String}()), sv(cen.classification[i]))
    end
    ca = Set{Int}(t for (t, cs) in classes if "Calcium" in cs && !("Lanthanide" in cs))
    return (classes = classes, euk_seqs = euk_seqs, euk_taxa = euk_taxa, ca_class_taxids = ca)
end

# The share of one phylum in a list of phylum labels. WHICH phylum is a selection, so it is
# an argument: a name written in here would compute one phylum's share under a key named after
# another the moment the selection changed.
pseudo(p, phylum::AbstractString) = 100 * count(==(phylum), p) / max(length(p), 1)

# ---------------------------------------------------------------- the move-in pool
"""
    read_movein_pool(path; read_bytes = read) -> NamedTuple

The taxa a census hit table gives at least one lanthanide-class sequence, with the digest of
the bytes they were read from and the counts behind them.

WHY A HIT TABLE AND NOT A LIST OF TAXIDS. A list would carry a number and no definition, and
the one thing this pool must not be is a count that arrived from somewhere. Reading the
table the classes live in means the rule is applied HERE, by this code, to whatever table was
supplied: the size is then a measurement of that file rather than a figure to be trusted, and
a file that is not a hit table is refused instead of being read as an empty pool.

The taxon column is spelled `TaxID` in the census tables and `taxid` in this repository's
own, so both are accepted; a row whose taxon cannot be read as a whole number is counted and
skipped rather than dropped in silence.
"""
function read_movein_pool(path::AbstractString; read_bytes = read)
    isfile(path) || error("--movein-pool $path does not exist")
    bytes = read_bytes(path)
    digest = bytes2hex(sha256(bytes))
    tbl = CSV.read(bytes, DataFrame)
    have = [String(c) for c in propertynames(tbl)]
    tcol = hascol(tbl, :TaxID) ? :TaxID : hascol(tbl, :taxid) ? :taxid : nothing
    (tcol === nothing || !hascol(tbl, :classification)) &&
        error("--movein-pool $path is not a hit table: it needs a 'classification' column " *
              "and a 'TaxID' (or 'taxid') column, and has " * joined(have))
    ln = sv.(tbl.classification) .== LANTHANIDE_CLASS
    ids = Set{Int}(); unparsed = 0
    for t in tbl[ln, tcol]
        i = taxid_of(t)
        i === nothing ? (unparsed += 1) : push!(ids, i)
    end
    return (taxids = ids, sha256 = digest, rows = nrow(tbl), class_rows = count(ln),
            rows_without_a_taxon = unparsed, columns = have)
end

"""
    movein_pool(; scan, path) -> NamedTuple

The pool the move-in level moves, under the rule the run is actually able to apply: the
census hit table at `path` when one was supplied, and otherwise the scan table's own
lanthanide-class flag — which is the `lntype` set, taken from the same function that mode
uses rather than from a second copy of its predicate.

`taxids` is `nothing` when the rule cannot be answered at all (no scan table on hand), which
is not an empty pool: an empty pool is a level that moved nobody, and a level that moved
nobody is exactly the silent failure this variant exists to avoid.
"""
function movein_pool(; scan, path::AbstractString = "", read_bytes = read)
    if !isempty(path)
        p = read_movein_pool(path; read_bytes = read_bytes)
        return (taxids = p.taxids, rule = MOVEIN_POOL_RULES[1], rule_source = :explicit,
                source = basename(path), sha256 = p.sha256, rows = p.rows,
                class_rows = p.class_rows, rows_without_a_taxon = p.rows_without_a_taxon,
                columns = p.columns)
    end
    return (taxids = exclusion_taxids(scan, "lntype"), rule = MOVEIN_POOL_RULES[2],
            rule_source = :default, source = "the background carrier scan table",
            sha256 = nothing, rows = scan === nothing ? missing : nrow(scan),
            class_rows = missing, rows_without_a_taxon = missing, columns = String[])
end

# What each pool rule means, in the words the emitted key carries. As with the exclusion
# notes, no rule describes its relationship to the other: whether the two pools overlap is a
# property of a pair of tables, measured per run, not something a note may assert.
function movein_pool_note(rule::AbstractString)
    rule == MOVEIN_POOL_RULES[1] && return "$(MOVEIN_POOL_RULES[1]): frame taxa carrying at " *
        "least one sequence the census itself classified lanthanide-type, read from the " *
        "census hit table supplied with --movein-pool; this is the rule the recorded level " *
        "names, and the table is not a repository input"
    rule == MOVEIN_POOL_RULES[2] && return "$(MOVEIN_POOL_RULES[2]): frame taxa whose scan " *
        "row carries a lanthanide-classed hit — this run's own model call on the family " *
        "members it fetched, which is a DIFFERENT measurement from the census's call on its " *
        "own alignment, and is descriptive"
    error("no note for move-in pool rule '$rule'")
end

# ---------------------------------------------------------------- coverage, drops, evidence
"""
    coverage!(prefix, cand_mask, usable, label) -> (ncand, nok, label)

The biggest threat to the contrast: if one arm is annotated more often than the other,
every ratio is conditional on being annotated, and that has to be said out loud. Emits how
many candidates of a mask carry a usable call.
"""
function coverage!(prefix, cand_mask, usable, label)
    ncand = count(cand_mask); nok = count(cand_mask .& usable)
    rec!("$(prefix)_coverage_numerator", nok)
    rec!("$(prefix)_coverage_denominator", ncand)
    rec!("$(prefix)_coverage_pct", 100 * nok / max(ncand, 1))
    return (ncand = ncand, nok = nok, label = label)
end

"""
    dropped!(prefix, df, usable, cand_mask) -> NamedTuple

Why candidate taxa without a usable call were dropped. Reporting the whole dropped set as
"isolation-source" would give one cause to taxa dropped for several, so the buckets are a
partition, in priority order, and the last bucket is what is left: they sum by construction.
"""
function dropped!(prefix, df, usable, cand_mask)
    ex = df[cand_mask .& .!usable, :]
    nrel = coalesce.(ex.n_relations, 0)
    iso  = coalesce.(ex.isolation_only, false)
    gen  = coalesce.(ex.generic_only, false)
    isot = is_isolation_source.(sv.(ex.habitat_type))
    unann = nrel .== 0
    b_unann = count(unann)
    b_iso   = count(.!unann .& iso)
    b_gen   = count(.!unann .& .!iso .& gen)
    b_isot  = count(.!unann .& .!iso .& .!gen .& isot)
    b_filt  = nrow(ex) - b_unann - b_iso - b_gen - b_isot
    rec!("$(prefix)_dropped_total", nrow(ex))
    rec!("$(prefix)_dropped_unannotated", b_unann)
    rec!("$(prefix)_dropped_isolation_only", b_iso)
    rec!("$(prefix)_dropped_generic_only", b_gen)
    rec!("$(prefix)_dropped_isolation_type_taxon_level", b_isot)
    rec!("$(prefix)_dropped_filtered_other", b_filt)
    return (unann = b_unann, iso = b_iso, gen = b_gen, isot = b_isot, filt = b_filt,
            tot = nrow(ex))
end

"""
    evidence!(prefix, sub)

The evidence behind the calls of one arm: term tallies, how ties were resolved, and where
the relations came from.
"""
function evidence!(prefix, sub)
    rec!("$(prefix)_isolation_terms_total", sum(coalesce.(sub.n_isolation_terms, 0)))
    rec!("$(prefix)_generic_terms_total", sum(coalesce.(sub.n_generic_terms, 0)))
    rec!("$(prefix)_taxa_with_isolation_term", count(coalesce.(sub.n_isolation_terms, 0) .> 0))
    rec!("$(prefix)_taxa_with_generic_term", count(coalesce.(sub.n_generic_terms, 0) .> 0))
    # resolved_by records where the tie died: :none means the top distinct count was
    # unique, :mentions that the tie broke on total statements, :fractional that it
    # survived both and the unit was split across the tied categories. n_tied is the size
    # of that surviving tie, so it is a property of the :fractional taxa only.
    rb = sv.(sub.resolved_by)
    rec!("$(prefix)_taxa_tied_on_distinct", count(in.(rb, Ref(("mentions", "fractional")))))
    rec!("$(prefix)_taxa_tie_broken_by_mentions", count(rb .== "mentions"))
    rec!("$(prefix)_taxa_still_tied_after_mentions", count(rb .== "fractional"))
    rec!("$(prefix)_max_residual_tie_size",
         isempty(sub) ? 0 : maximum(coalesce.(sub.n_tied, 0)))
    rec!("$(prefix)_taxa_truncated_term_list", count(coalesce.(sub.n_truncated, 0) .> 0))
    rec!("$(prefix)_taxa_with_lt3_distinct_terms",
         count(.!well_supported_mask(sub.n_distinct)))
    rec!("$(prefix)_relations_total", sum(coalesce.(sub.n_relations, 0)))
    ntext = sum(coalesce.(sub.n_src_text, 0)); nrecd = sum(coalesce.(sub.n_src_record, 0))
    rec!("$(prefix)_relations_text_mined", ntext)
    rec!("$(prefix)_relations_database_or_culture_collection", nrecd)
    rec!("$(prefix)_relations_text_mined_pct", 100 * ntext / max(ntext + nrecd, 1))
end

# ---------------------------------------------------------------- the recorded criteria
"""
    criterion!(prefix, bound, value, relation, origin) -> Union{Bool,Nothing}

Disclose one numeric criterion the plan for this step fixed before the run: the bound, the
value measured now, the comparison between them in words, and whether it holds. Nothing
here gates a code — the criteria belong to the step, not to a rule, and a run that violates
one should still produce its numbers so the violation can be seen — but a violated
criterion is a warning, not a silence, and the keys make it diffable.

`bound === nothing` emits nothing at all. That is deliberate: an absent bound means nobody
supplied one on this run, and a default would be this script inventing the plan's number.
"""
function criterion!(prefix, bound, value, relation::Symbol, origin::AbstractString)
    bound === nothing && return nothing
    met = relation === :below ? value < bound :
          relation === :at_most ? value <= bound :
          error("unknown criterion relation $relation")
    rec!("criterion_$(prefix)_bound", bound)
    rec!("criterion_$(prefix)_value", value)
    rec!("criterion_$(prefix)_relation",
         relation === :below ? "value < bound" : "value <= bound")
    rec!("criterion_$(prefix)_met", met)
    rec!("criterion_$(prefix)_verdict", "DISCLOSED")
    rec!("criterion_$(prefix)_origin", origin)
    met || @warn("a numeric criterion recorded for this step is not met by this run",
                 criterion = prefix, measured = value, bound = bound)
    return met
end

# ---------------------------------------------------------------- rule-keyed quantities
function emit_row!(prefix, r)
    if r === nothing
        rec!("$(prefix)_present", false); return
    end
    rec!("$(prefix)_present", true)
    rec!("$(prefix)_a", r.a); rec!("$(prefix)_b", r.b)
    rec!("$(prefix)_n1", r.n1); rec!("$(prefix)_n2", r.n2)
    rec!("$(prefix)_log2_rr", r.log2_rr)
    rec!("$(prefix)_ci_lo", r.ci_lo); rec!("$(prefix)_ci_hi", r.ci_hi)
    rec!("$(prefix)_log2_rr_std", r.log2_rr_std)
    rec!("$(prefix)_ci_lo_std", r.ci_lo_std); rec!("$(prefix)_ci_hi_std", r.ci_hi_std)
    rec!("$(prefix)_q_bh", r.q_bh)
    rec!("$(prefix)_interval_excludes_zero_std", zero_std(r))
    rec!("$(prefix)_interval_excludes_zero_crude", zero_crude(r))
    # Whether the row was tested at all, and what the standardised interval rests on. A row
    # whose arm totals are healthy can still be carried by one empty stratum cell, and a
    # rule that consults the interval should not have to take that on trust.
    rec!("$(prefix)_tested", r.tested)
    rec!("$(prefix)_n_zero_stratum_cells", r.n_zero_stratum_cells)
    rec!("$(prefix)_min_bg_stratum_cell", r.min_bg_stratum_cell)
end

# The host-plant decisions: two rules keyed on one interval each, and the mapping from
# interval to code is the module's: the genus clause is on when the genus interval includes
# zero, the category is dropped when the taxa interval includes zero, and an interval that
# was never computed is undecidable — a third value this script may not rename.
genus_clause_code(r, reading) =
    interval_code(row_bounds(r, reading)...; excludes = :off, includes = :on)
category_code(r, reading) =
    interval_code(row_bounds(r, reading)...; excludes = :keep, includes = :drop)

"""
    ratio_of(num, den; standardised_scale) -> Float64

The comparator estimate as a fraction of the reference estimate, NaN when either estimate
does not exist or the reference is zero. The grading of that fraction is not decided here —
`OutcomeRules.sharing_wording_code` does it, including for the fractions that are outside
the rule's domain.
"""
function ratio_of(num, den; standardised_scale::Bool)
    (num === nothing || den === nothing) && return NaN
    a = standardised_scale ? num.log2_rr_std : num.log2_rr
    b = standardised_scale ? den.log2_rr_std : den.log2_rr
    (isnan(a) || isnan(b) || b == 0) && return NaN
    return a / b
end

"""
    ratio_interval(...) -> (lo, hi, n_ok)

A percentile interval for the ratio of two standardised estimates, by resampling units
with replacement. The two foreground arms are disjoint and are resampled independently;
the background is resampled ONCE per replicate and used by both, which keeps the
correlation the shared comparator induces instead of pretending the two ratios are
independent. No estimator is defined here — each replicate calls the module's
`standardised` on resampled indices.

Every comparator arm's interval comes out of this one function at the run's own seed, so a
disclosed arm's interval is comparable with the registered arm's rather than being a second
resampling scheme with its own behaviour.
"""
function ratio_interval(num_phy, num_w, den_phy, den_w, b_phy, b_w, cat;
                        nboot::Int, seed::Int)
    (haskey(num_w, cat) && haskey(den_w, cat) && haskey(b_w, cat)) || return (NaN, NaN, 0)
    nboot <= 0 && return (NaN, NaN, 0)
    fn = num_w[cat]; fd = den_w[cat]; fb = b_w[cat]
    rng = Xoshiro(seed)
    nn = length(num_phy); nd = length(den_phy); nb = length(b_phy)
    out = Float64[]
    for _ in 1:nboot
        i0 = rand(rng, 1:nb, nb)
        bp = view(b_phy, i0); bh_ = view(fb, i0)
        inum = rand(rng, 1:nn, nn); iden = rand(rng, 1:nd, nd)
        sn = standardised(view(num_phy, inum), view(fn, inum), bp, bh_)
        sd = standardised(view(den_phy, iden), view(fd, iden), bp, bh_)
        (isnan(sn.log2_rr) || isnan(sd.log2_rr) || sd.log2_rr == 0) && continue
        push!(out, sn.log2_rr / sd.log2_rr)
    end
    length(out) < 100 && return (NaN, NaN, length(out))
    return (quantile(out, 0.025), quantile(out, 0.975), length(out))
end

# ---------------------------------------------------------------- the scan's own evidence
"""
    scan_table_state(path) -> (state, reasons, pending_beside_it, displaced_copy_on_disk)

What the file `--scan` names IS, as a matter of file state rather than of content:
`:published` (a table a scan run promoted to its final name), `:displaced` (the `.stale`
copy a failing run moved aside), `:unpromoted` (a `.pending` file no run has promoted), or
`:absent`.

WHY A NAME IS EVIDENCE. The producer publishes its table and its three records together or
not at all, and on a refused acceptance check it renames the last published set to `.stale`.
A `.stale` set is therefore internally consistent — it was a passing run once — so every
content check is satisfied by it, and the digests match one another. What it is not is a
table any current run vouches for: a later run looked at the same frame and refused it. A
`.pending` file beside a published one says the same thing from the other side, that the most
recent attempt did not publish, so the published table may have been superseded by a run
whose outcome nobody recorded.

The reasons are joined with the module's, so the emitted reason list says which of the two
kinds of doubt applies. The two sibling flags are looked at once, here, and reported as seen.
"""
function scan_table_state(path::AbstractString)
    reasons = String[]
    pending = isfile(path * ScanRecords.PENDING_EXT)
    displaced = isfile(path * ScanRecords.STALE_EXT)
    state = !isfile(path) ? :absent :
            endswith(path, ScanRecords.STALE_EXT) ? :displaced :
            endswith(path, ScanRecords.PENDING_EXT) ? :unpromoted : :published
    state === :absent &&
        push!(reasons, "scan table: $(basename(path)) does not exist, so nothing was read")
    state === :displaced &&
        push!(reasons, "scan table: $(basename(path)) is the copy a failing scan run " *
                       "displaced; no run vouches for it")
    state === :unpromoted &&
        push!(reasons, "scan table: $(basename(path)) was never promoted, so the run that " *
                       "wrote it did not publish it")
    state === :published && pending &&
        push!(reasons, "scan table: an unpromoted $(basename(path))" *
                       "$(ScanRecords.PENDING_EXT) sits beside it, so a later scan run " *
                       "has not published and this table may be superseded")
    return (state = state, reasons = reasons, pending_beside_it = pending,
            displaced_copy_on_disk = displaced)
end

"""
    read_scan_table(path; read_bytes = read) -> (table, sha256)

The scan table and the SHA-256 digest of the very bytes it was parsed from, from ONE read of
the file; `(nothing, nothing)` when there is no file.

WHY ONE READ. The digest is the one value that lets a record be refused for describing
another table, so it may not come from anything that could itself be foreign — and a second
read of the same path is exactly that. A scan run that promotes a new table between a parse
and a later hash hands this script the old table's exclusions under the new table's digest,
and every record of the new run then agrees with a table that was never analysed.

`read_bytes` is how the file's bytes are obtained, and exists so a test can show that the
digest and the parse come from the same call: a reader that answers differently on a second
call exposes any implementation that reads twice.
"""
function read_scan_table(path::AbstractString; read_bytes = read)
    isfile(path) || return (table = nothing, sha256 = nothing)
    bytes = read_bytes(path)
    digest = bytes2hex(sha256(bytes))
    return (table = CSV.read(bytes, DataFrame), sha256 = digest)
end

"""
    read_scan_records(threshold_path, recall_path, acceptance_path, release_path;
                      read_bytes = read) -> NamedTuple

The scan's three persisted records and its search space's release record, each read from
disk ONCE: `threshold`/`recall` (typed, by the producer's readers), `threshold_kv`/
`recall_kv` (the same records as written, for the provenance keys the summary echoes
verbatim), `acceptance` with `acceptance_error`, `release_kv` (the UniProtKB release record
the scan wrote beside its search-space download), and whether each file existed at that read.

WHY ONCE. The threshold and recall records are read two ways, and two reads of one path are
two moments: a scan run promoting between them would leave typed fields that were checked
beside provenance text that was not. Each key/value file is therefore read once, through
`read_bytes`, into a private snapshot, and every reader parses the snapshot. The snapshot is
a file only because the producer's readers take a path. The acceptance table has one reader,
so its one read is that reader's. `read_bytes` exists for the same reason as in
`read_scan_table`: a test can answer differently on a second read.

WHY THE RELEASE RECORD IS HERE. The recall record names the release its search space came
from, and the licence holds that name against the release record itself; the summary reports
the same record. Reading it with the others means the release the licence was checked against
and the release the table reports are one read of one file.

`acceptance` is `nothing`, with the reader's message in `acceptance_error`, when the file is
not an acceptance table. An older record of a different shape is a record of something
else, and that refusal makes the evidence unknown, which withholds the licence, rather than
ending a run that can still report its numbers.
"""
function read_scan_records(threshold_path::AbstractString, recall_path::AbstractString,
                           acceptance_path::AbstractString, release_path::AbstractString;
                           read_bytes = read)
    mktempdir() do snapshot_dir
        function snapshot(path, slot)
            copy = joinpath(mkpath(joinpath(snapshot_dir, slot)), basename(path))
            present = isfile(path)
            present && write(copy, read_bytes(path))
            return (present = present, path = copy)
        end
        t = snapshot(threshold_path, "threshold")
        r = snapshot(recall_path, "recall")
        u = snapshot(release_path, "release")
        acceptance_present = isfile(acceptance_path)
        acc, acc_error = try
            (ScanRecords.read_acceptance(acceptance_path), "")
        catch err
            (nothing, sprint(showerror, err))
        end
        return (threshold = ScanRecords.read_threshold(t.path), threshold_kv = read_kv(t.path),
                threshold_present = t.present,
                recall = ScanRecords.read_recall(r.path), recall_kv = read_kv(r.path),
                recall_present = r.present,
                acceptance = acc, acceptance_error = acc_error,
                acceptance_present = acceptance_present,
                release_kv = read_kv(u.path), release_present = u.present)
    end
end

# The acceptance rows the scan gates on in EVERY run, by the check text its own rows carry.
# That wording belongs to scripts/background_pqq_scan.jl: the family row and the positive
# controls are spelled here from that script's own constants, and the other four are spelled
# out and must match its rows word for word — test/test_habitat_enrichment.jl holds each
# against the producer's source, so a reworded row fails a test rather than every licence.
# They are named here, from the producer, and never read off the acceptance table itself: a
# list taken from the table that was read cannot notice a row that was cut from it.
#
# Two of the producer's rows are gated only when there is something to gate — an earlier
# sample for the recall sample to extend, a download whose release was recorded — and are
# disclosures otherwise. They are not required, because a first scan run would then license
# nothing on a table nobody faults; they are still refused on a MISS, and their verdicts are
# still re-derived from their own numbers.
const FAMILY_CHECK  = "$(ScanRecords.IPR_FAMILY) frame taxa flagged " *
                      "(>= $(ScanRecords.MIN_FAMILY_LEN) aa, whole-frame membership)"
const RECALL_CHECK  = "prefilter recall"
const FLAGGED_CHECK = "flagged background taxa"
const REQUIRED_GATED_CHECKS = vcat(
    [FAMILY_CHECK, "calibration recovery", RECALL_CHECK, FLAGGED_CHECK,
     "hit sequences classed at the active-site column"],
    ["positive control $name ($t)" for (t, name) in ScanRecords.POSITIVE_CONTROLS])
# The recall row writes its rate in the shortest form that reads back to the same Float64,
# from the same two counts, so the only difference it may show is float-parsing noise.
const RECALL_ROW_TOL = 1e-12

"""
    acceptance_cell(acceptance, check, column) -> value or `nothing`

One numeric cell of the acceptance row whose check text is `check`. `nothing` when there is
no table, when that row is absent or appears more than once, or when the cell is empty: a
value picked from one of two rows is a guess, and the module reports the absence or the
duplicate on its own.
"""
function acceptance_cell(acc, check::AbstractString, column::Symbol)
    acc === nothing && return nothing
    at = findall(==(check), acc.check)
    length(at) == 1 || return nothing
    v = acc[only(at), column]
    return v === missing ? nothing : v
end

"""
    scan_table_counts(scan) -> (n_flagged, n_family)

The two counts the acceptance table has to agree with, taken from the table that was parsed:
the taxa with at least one hit, and the taxa with a family-signature member at the length
cut-off. `nothing` for a count whose column is absent or holds anything but whole counts or
flags, so an unreadable column is unknown rather than zero.
"""
function scan_table_counts(scan)
    scan === nothing && return (n_flagged = nothing, n_family = nothing)
    column_count(col, valid, f) =
        hascol(scan, col) && all(valid, scan[!, col]) ? count(f, scan[!, col]) : nothing
    return (n_flagged = column_count(:n_hits, x -> x isa Integer && !(x isa Bool), >(0)),
            n_family = column_count(:in_ipr017512_ge500, x -> x isa Bool, identity))
end

# A count for a reason string: written as a whole number only where double-precision
# arithmetic holds it exactly, so a hand-edited huge value cannot throw while the reason is
# being written.
count_text(x::Real) = (isfinite(x) && isinteger(x) && abs(x) <= maxintfloat(Float64)) ?
                      string(Int(x)) : string(x)

"""
    acceptance_row_reasons(acceptance, recall) -> Vector{String}

Two checks on the acceptance table's own numbers, made beside the module's:

  * every row's verdict is re-derived from that row's value, denominator and bar by the
    producer's own `floor_verdict`, and must equal the verdict the row prints — a row that
    says PASS beside a value below its bar, or DISCLOSED beside a bar, is refused;
  * the recall row counts the same flagged taxa as the recall record, and gives the rate the
    record's two counts give.

The module ties the recall record's bar to the recall row's and refuses a missing or doubled
row, but it is handed verdict words, not the numbers behind them, so a record whose words
disagree with its numbers would pass it. The rule used is the producer's, called, not
restated. With no acceptance table there is nothing to re-derive, and the module already
refuses that.
"""
function acceptance_row_reasons(acc, recall)
    reasons = String[]
    acc === nothing && return reasons
    for r in eachrow(acc)
        derived = ScanRecords.verdict_of(ScanRecords.floor_verdict(r.value, r.denominator, r.bar))
        r.verdict == derived || push!(reasons,
            "acceptance_rows: \"$(r.check)\" says $(r.verdict), but its own value, " *
            "denominator and bar give $derived")
    end
    nf, nr = recall.n_full_flagged, recall.n_recovered
    if nf !== nothing && nr !== nothing && 0 <= nr <= nf
        rate = nf == 0 ? NaN : nr / nf
        v = acceptance_cell(acc, RECALL_CHECK, :value)
        d = acceptance_cell(acc, RECALL_CHECK, :denominator)
        (d === nothing || d == nf) || push!(reasons,
            "acceptance_rows: the \"$RECALL_CHECK\" row counts $(count_text(d)) flagged " *
            "taxa, but the recall record counts $nf")
        (v === nothing || (isnan(rate) ? isnan(v) : (!isnan(v) && abs(v - rate) <= RECALL_ROW_TOL))) ||
            push!(reasons, "acceptance_rows: the \"$RECALL_CHECK\" row gives a rate of $v, " *
                           "but the recall record's counts give $rate")
    end
    return reasons
end

"""
    frame_comparison(frame_taxids, scan) -> (matches, missing_from_scan, extra_in_scan)

A taxon-by-taxon comparison of the analysed frame with the taxa the scan table holds a row
for. `matches` is `nothing` when there is no table. A carrier whose row is missing from the
table reads as a non-carrier to every count taken from it, which is why the two sets are
compared rather than their sizes.
"""
function frame_comparison(frame_taxids, scan)
    scan === nothing && return (matches = nothing, missing_from_scan = Int[], extra_in_scan = Int[])
    frame = Set(Int.(frame_taxids))
    scanned = Set(Int.(scan.taxid))
    missing_from_scan = sort(collect(setdiff(frame, scanned)))
    extra_in_scan = sort(collect(setdiff(scanned, frame)))
    return (matches = isempty(missing_from_scan) && isempty(extra_in_scan),
            missing_from_scan = missing_from_scan, extra_in_scan = extra_in_scan)
end

"""
    leftover_counts(scan, arm_taxids) -> (hmm_hit, ipr017512_ge500, union)

Carriers left in the analysed background arm, by two independent criteria: the HMM hit
table, and whole-frame family membership at the length cut-off, which is computed without
the sequence prefilter and so can disagree with the hit table. `arm_taxids` are the taxa
ACTUALLY ANALYSED: a count taken over the scan table with the predicate that built the
exclusion is zero whatever the data say, so it is evidence of nothing. `nothing` for a
count the table cannot answer: no table at all, or no family column.
"""
function leftover_counts(scan, arm_taxids)
    scan === nothing && return (hmm_hit = nothing, ipr017512_ge500 = nothing, union = nothing)
    hit = Dict{Int,Bool}(Int(t) => coalesce(n, 0) > 0 for (t, n) in zip(scan.taxid, scan.n_hits))
    arm = [Int(t) for t in arm_taxids]
    left_hmm = count(t -> get(hit, t, false), arm)
    hascol(scan, :in_ipr017512_ge500) ||
        return (hmm_hit = left_hmm, ipr017512_ge500 = nothing, union = nothing)
    fam = Dict{Int,Bool}(Int(t) => coalesce(f, false)
                         for (t, f) in zip(scan.taxid, scan.in_ipr017512_ge500))
    return (hmm_hit = left_hmm, ipr017512_ge500 = count(t -> get(fam, t, false), arm),
            union = count(t -> get(hit, t, false) || get(fam, t, false), arm))
end

"""
    scan_evidence_and_licence(arms; table_read, table_state, records, required_checks)
        -> NamedTuple

Everything the background licence is decided from, assembled in one place and handed to the
module: the frame comparison and the leftovers over `arms` (an `analyse` result: its `df`,
`frame` mask and analysed background `bgm`), the scan table together with the digest of the
bytes it was parsed from (`table_read`, the whole result of `read_scan_table`, so the two can
never come from different reads), its file state (`scan_table_state`), and the records and the
search space's release record (`read_scan_records`).

Returns `frame`, `leftovers`, `table_counts`, `evidence` (the module's
`scan_evidence_status`), `row_reasons` (`acceptance_row_reasons`), `ok` and `reasons` (the
scan evidence as the licence sees it), `removal` (`background_removal_complete`) and `code`.

`arms` must be an `analyse` result, or a value of the same shape: the frame comparison reads
its `frame_scanned` mask — the frame as the scan saw it — and not the analysed `frame`, which
the move-in level makes smaller.

WHY THIS IS A FUNCTION. As top-level script code this wiring is reachable by no test, and a
one-line change to it — the evidence verdict replaced by `true`, the file-state gate
dropped, the record's own pass field forwarded in place of the re-derived verdict — licenses
the strong code on a refused or displaced table while every suite stays green. Every input
the licence reads is passed here and nowhere else:

  * the recall verdict is the module's re-derivation on the lower bound, never the record's
    pass field, which is taken on the point estimate;
  * the scan evidence is `ok` only when the module's evidence holds AND the file state, the
    re-derived acceptance rows and the acceptance read raised nothing;
  * the bar the recall record judges its counts by is held against the acceptance table's
    recall row, its cutoff against the threshold record, and its release against the search
    space's release record.
"""
function scan_evidence_and_licence(arms; table_read, table_state, records,
                                   required_checks = REQUIRED_GATED_CHECKS,
                                   pre_registered_bar = ScanRecords.MIN_RECALL)
    table = table_read.table
    # The frame the scan was ASKED ABOUT, which is the frame before the move-in level hands
    # any of it to the foreground: those taxa were scanned, and counting them as "scan taxa
    # outside the analysed frame" would report that level's own definition as a gap in the
    # scan's coverage and withhold the licence on every run of it. The leftovers stay over
    # the analysed background arm, which is what they are a statement about.
    frame = frame_comparison(arms.df.taxid[arms.frame_scanned], table)
    left = leftover_counts(table, arms.bgm.taxid)
    counts = scan_table_counts(table)
    thr = records.threshold; rec = records.recall; acc = records.acceptance
    evidence = scan_evidence_status(;
        scan_sha256 = table_read.sha256,
        record_sha256s = ["threshold record" => thr.scan_sha256,
                          "recall record" => rec.scan_sha256,
                          "acceptance record" => acc === nothing ? nothing : acc.scan_sha256],
        acceptance_checks = acc === nothing ? nothing : acc.check,
        acceptance_verdicts = acc === nothing ? nothing : acc.verdict,
        required_checks = required_checks,
        flagged_row_value = acceptance_cell(acc, FLAGGED_CHECK, :value),
        table_n_flagged = counts.n_flagged,
        family_row_denominator = acceptance_cell(acc, FAMILY_CHECK, :denominator),
        table_n_family = counts.n_family,
        recall_n_sampled = rec.n_sampled,
        recall_n_full_flagged = rec.n_full_flagged,
        recall_n_recovered = rec.n_recovered,
        recall_rate_field = rec.rate,
        recall_lower_bound_field = rec.lower_bound_95,
        recall_bar = rec.bar,
        recall_pass_field = rec.pass,
        recall_threshold_bits = rec.threshold_bits,
        recall_uniprot_release = rec.uniprot_release,
        # The bar itself is evidence. Every other input is checked against a second copy, so a
        # bar that only ever comes out of the records would be the one number a two-cell edit
        # could move — and it is the number the stronger background statement turns on. The
        # scan's own constant is the second copy.
        bar_in_other_records =
            ["acceptance record, \"$RECALL_CHECK\" row" => acceptance_cell(acc, RECALL_CHECK, :bar),
             "the scan's recorded recall bar" => pre_registered_bar],
        threshold_bits_in_other_records = ["threshold record" => thr.threshold_bits],
        uniprot_release_in_other_records =
            ["search-space release record" => get(records.release_kv, "X-UniProt-Release", nothing)])
    row_reasons = acceptance_row_reasons(acc, rec)
    read_error = isempty(records.acceptance_error) ? String[] :
                 ["acceptance record: " * records.acceptance_error]
    ok = evidence.ok && isempty(table_state.reasons) && isempty(row_reasons) && isempty(read_error)
    reasons = vcat(table_state.reasons, evidence.reasons, row_reasons, read_error)
    # Two different failures land in the same `ok`, and saying so matters: a record that does not
    # describe this table is a broken provenance chain, while a recall whose interval is too wide
    # is an honest measurement that simply does not reach the bar. Only the first means the
    # records cannot be trusted, so it is reported separately from the licence.
    records_ok = isempty(filter(r -> !startswith(r, "recall_n_recovered:") &&
                                     !startswith(r, "recall_n_full_flagged:"), reasons))
    removal = background_removal_complete(; left_hmm_hit = left.hmm_hit,
                                          left_ipr_ge500 = left.ipr017512_ge500,
                                          recall_pass = evidence.recall_pass_derived,
                                          frame_matches = frame.matches,
                                          scan_evidence_ok = ok)
    return (frame = frame, leftovers = left, table_counts = counts, evidence = evidence,
            row_reasons = row_reasons, ok = ok, records_ok = records_ok, reasons = reasons,
            removal = removal,
            code = background_licence_code(removal))
end

# ---------------------------------------------------------------- the anchor
"""
    reference_row_diff(anchor_rows, ref) -> NamedTuple

Compare the anchor's rows, row by row, with the rows of the earlier reference run: how many
of the reference's rows the anchor has, how many agree cell for cell, which disagree, and
the largest cell difference. Only the columns both tables carry are compared.

WHY ROW LEVEL AS WELL AS SET LEVEL. The partition is decided on category SETS, and that is
what `anchor_status` gates on, because the recorded question is which categories excluded
zero. But a set can agree while the numbers behind it move, and an anchor whose rows have
drifted is no longer the earlier measurement even when its verdicts happen to match. The
row comparison is emitted as a disclosure for exactly that case: it decides nothing and
hides nothing.
"""
function reference_row_diff(anchor_rows, ref::AbstractDataFrame)
    have = [c for c in REF_COMPARE_COLS if hascol(ref, c)]
    rowkey(contrast, unit, habitat) = string(contrast, " | ", unit, " | ", habitat)
    byk = Dict{String,Any}()
    for r in anchor_rows
        byk[rowkey(r.contrast, r.unit, r.habitat)] = r
    end
    compared = 0; agreeing = 0
    differing = String[]; absent = String[]; worst = 0.0
    for i in 1:nrow(ref)
        k = rowkey(sv(ref.contrast[i]), sv(ref.unit[i]), sv(ref.habitat[i]))
        a = get(byk, k, nothing)
        if a === nothing
            push!(absent, k)
            continue
        end
        compared += 1
        d = 0.0; comparable = true
        for c in have
            x = fnum(ref[i, c]); y = fnum(getproperty(a, c))
            if isnan(x) && isnan(y)
                continue                       # both uncomputed: the same statement
            elseif isnan(x) || isnan(y)
                comparable = false             # one has an estimate and the other has none
            else
                d = max(d, abs(x - y))
            end
        end
        worst = max(worst, d)
        if comparable && d <= REF_ROW_TOL
            agreeing += 1
        else
            push!(differing, k)
        end
    end
    return (columns = String.(have), compared = compared, agreeing = agreeing,
            differing = sort(differing), absent = sort(absent), worst = worst)
end

"""
    read_reference_rows(path) -> Union{DataFrame,Nothing}

The earlier run's rows table, or `nothing` when the file cannot be read as one: a parse
failure, or a table without the columns a category set is read from. `nothing` becomes the
anchor status `:undecidable_reference_unreadable`, which withholds the partition. A reference
that cannot be read has verified nothing, and treating it as an empty reference would make
every anchor category look like a disagreement — or, compared the other way, like agreement
with a run that had no categories at all.
"""
function read_reference_rows(path::AbstractString)
    ref = try
        CSV.read(path, DataFrame)
    catch err
        @warn("the anchor reference could not be read as a rows table", path,
              error = sprint(showerror, err))
        return nothing
    end
    needed = (:contrast, :unit, :habitat, :ci_lo_std, :ci_hi_std)
    all(c -> hascol(ref, c), needed) && return ref
    @warn("the anchor reference is not a rows table", path,
          missing_columns = [String(c) for c in needed if !hascol(ref, c)])
    return nothing
end

"""
    anchor_status_of(; frames_match, n_undecidable_rows, reference_supplied, reference_exists,
                       reference_readable, anchor_categories, reference_categories) -> Symbol

The anchor's status, decided in a fixed order in which each rung withholds the partition
for a reason that makes every later rung meaningless:

  1. `:undecidable_frame_mismatch`       the anchor was built on another frame, so it is not
                                         the same measurement on earlier settings;
  2. `:undecidable_anchor_intervals`     a row of the anchor panel has no interval, so its
                                         category set is not the earlier result;
  3. `:computed_unverified`              no reference was supplied: the anchor stands, but
                                         nothing checked it;
  4. `:undecidable_reference_missing`    a reference was named and does not exist;
  5. `:undecidable_reference_unreadable` it exists and is not a rows table;
  6. `:reproduces_reference` or `:differs_from_reference`, on the category sets.

Only a missing reference at rung 3 lets an unverified anchor through; see
`anchor_partition_available`.
"""
function anchor_status_of(; frames_match::Bool, n_undecidable_rows::Integer,
                          reference_supplied::Bool, reference_exists::Bool,
                          reference_readable::Bool, anchor_categories, reference_categories)
    return !frames_match ? :undecidable_frame_mismatch :
           n_undecidable_rows > 0 ? :undecidable_anchor_intervals :
           !reference_supplied ? :computed_unverified :
           !reference_exists ? :undecidable_reference_missing :
           !reference_readable ? :undecidable_reference_unreadable :
           anchor_categories == reference_categories ? :reproduces_reference :
                                                       :differs_from_reference
end

# Whether an anchor status lets the partition use the anchor's category set at all.
anchor_partition_available(status::Symbol) = status in (:computed_unverified, :reproduces_reference)

# ---------------------------------------------------------------- report
cell(x) = isinteger(x) ? @sprintf("%d", Int(x)) : @sprintf("%.1f", x)

function show_panel(rows, contrast_name, unit)
    sel = [r for r in rows if r.contrast == contrast_name && r.unit == unit]
    isempty(sel) && return
    sort!(sel, by = r -> -(isnan(r.log2_rr_std) ? r.log2_rr : r.log2_rr_std))
    @printf("\n  %s (%s)   n_fg=%d  n_bg=%d\n", contrast_name, unit, sel[1].n1, sel[1].n2)
    @printf("    %-18s %6s %6s %7s %7s %9s %10s %9s %6s %s\n",
            "habitat", "fg", "bg", "%fg", "%bg", "crude", "STD", "q(BH)", "min c_k", "")
    for r in sel
        flag = !r.tested ? " (n<min)" :
               (!isnan(r.q_bh) && r.q_bh < 0.05 ?
                (isnan(r.log2_rr_std) ? "" : (r.log2_rr_std > 0 ? " ↑" : " ↓")) : "")
        @printf("    %-18s %6s %6s %6.1f%% %6.1f%% %9.2f %10s %9.3g %6s%s\n",
                first(r.habitat, 18), cell(r.a), cell(r.b), r.pct_fg, r.pct_bg,
                r.log2_rr,
                isnan(r.log2_rr_std) ? "n/a" : @sprintf("%.2f", r.log2_rr_std), r.q_bh,
                isnan(r.min_bg_stratum_cell) ? "n/a" : cell(r.min_bg_stratum_cell), flag)
    end
end

# A hit total is a sum of weights, so the column is Float64 as soon as ANY unit is
# fractional (the genus collapse always is). Writing an integral total back as an integer
# keeps whole-unit rows readable and keeps the legacy anchor textually identical to the
# run it has to reproduce; a reader still parses the column as floating point.
compact(x::Real) = isinteger(x) ? Int(round(x)) : float(x)

fmt(k) = (v = look(k); v === missing ? "n/a" : v isa AbstractFloat ? @sprintf("%.3g", v) :
          string(v))
blank(k) = (s = fmt(k); isempty(s) ? "(none)" : s)

# ---------------------------------------------------------------------------
# Everything above is definitions: including this file runs nothing. Everything below reads
# the inputs, computes and writes, and runs only when this file is the program.
# ---------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__

# The command line, bound once to the names the driver reads.
const CONFIG                 = OPT.config
const BG_EXCL                = OPT.bg_exclusion
const BG_EXCL_SOURCE         = OPT.bg_exclusion_source
const BG_EXCL_ROLE           = OPT.bg_exclusion_role
const MOVEIN_POOL_PATH       = OPT.movein_pool
const NBOOT                  = OPT.nboot
const SEED                   = OPT.seed
const LEGACY                 = OPT.legacy
const CHEM_COL               = OPT.chem_col
const CALLS                  = OPT.calls
const SCAN                   = OPT.scan
const CENSUS                 = OPT.census
const RECALLF                = OPT.recall
const THRESHF                = OPT.threshold_record
const ACCEPTF                = OPT.acceptance
const ANCHORREF              = OPT.anchor_reference
const OUT                    = OPT.out
const SUMOUT                 = OPT.summary
const ANCHOR_OUT             = OPT.anchor_out
const PLAN_BG_SIZE_BELOW     = OPT.plan_bg_size_below
const PLAN_BG_PHYLUM_PCT_MAX = OPT.plan_bg_phylum_pct_max

# The supplied selection, bound to the names the driver reads. Not one of these is a default:
# every one came out of the criteria file, which was read above before the data were opened.
const SEL       = SELECTION
const TIER_CAT  = SEL.depletion_tier
const GENUS_CAT = SEL.genus_clause
const SHARE_CAT = SEL.sharing
const PHYLUM    = SEL.standardisation_phylum
# Key prefixes, derived from those names. No family of keys below is written as a literal.
const K_TIER    = SEL.key.depletion_tier
const K_GENUS   = SEL.key.genus_clause
const K_SHARE   = SEL.key.sharing
const K_PHYLUM  = SEL.key.standardisation_phylum
# Categories the recorded rules already speak for by name: exactly the three the supplied
# rules are about.
const GOVERNED  = SEL.governed

# ---------------------------------------------------------------- the inputs, by name
# Before the banner, because a refusal about a file somebody has to go and fetch should be
# the first and only thing this run says. Only the inputs still at their REGISTRY DEFAULT
# are named here: a path the reader typed is theirs, and the plain messages further down
# report it as given rather than pretending the layout is at fault.
#
# The scan RECORDS are deliberately not in this list. Their absence is not a refusal in this
# step by design — it withholds the background licence and says so in the output — and
# turning that into a refusal here would delete the distinction the step is built around.
haskey(CLI, "--calls") ||
    DataInputs.require_or_exit("habitat_calls")
(haskey(CLI, "--scan") || BG_EXCL == "none") ||
    DataInputs.require_or_exit("background_scan_table")

# ---------------------------------------------------------------- opening banner
@printf("\nhabitat_enrichment.jl\n  config %s   chemistry column %s\n",
        CONFIG, String(CHEM_COL))
@printf("  background exclusion: %s (%s, %s)\n", BG_EXCL, BG_EXCL_SOURCE, BG_EXCL_ROLE)
BG_EXCL == MOVEIN_MODE &&
    @printf("  move-in pool: %s\n",
            isempty(MOVEIN_POOL_PATH) ? "none supplied — the scan table's own class call" :
            MOVEIN_POOL_PATH)
@printf("  calls %s\n  scan  %s\n", CALLS, SCAN)
@printf("  scan records: %s, %s, %s\n",
        basename(THRESHF), basename(RECALLF), basename(ACCEPTF))
@printf("  resampling: %d replicates, seed %d\n", NBOOT, SEED)
PLAN_BG_SIZE_BELOW === nothing ||
    @printf("  plan bound, analysed background size: < %d\n", PLAN_BG_SIZE_BELOW)
PLAN_BG_PHYLUM_PCT_MAX === nothing ||
    @printf("  plan bound, background %s share: <= %.4g%%\n",
            PHYLUM, PLAN_BG_PHYLUM_PCT_MAX)
@printf("  categories %s\n", SEL.path)
@printf("    rules speak for %s;  grid counts %s;  standardisation phylum %s\n",
        joined(GOVERNED), joined(SEL.grid_count), PHYLUM)

# ---------------------------------------------------------------- the call table
isfile(CALLS) || error("missing $CALLS — run scripts/habitat_reaggregate.jl first")
calls = CSV.read(CALLS, DataFrame)
const CONFIGS_AVAILABLE = sort(unique(sv.(calls.config)))

# ---------------------------------------------------------------- the scan and its records
# The table's file state, its bytes and its digest, then the three records and the search
# space's release record, each read once and together, before anything is computed from any
# of them. See ONE READ OF EVERY FILE THE LICENCE RESTS ON in the header.
scan_state = scan_table_state(SCAN)
scan_read = read_scan_table(SCAN)
scan = scan_read.table
scan_sha = scan_read.sha256
records = read_scan_records(THRESHF, RECALLF, ACCEPTF, UPRELF)

if scan === nothing && BG_EXCL != "none"
    # A displaced copy is named in the message and NOT read: the producer moved it aside
    # because a later run refused the frame, and falling back to it here would hand this
    # run exactly the table no run vouches for.
    stale_copy = SCAN * ScanRecords.STALE_EXT
    hint = isfile(stale_copy) ?
        "A displaced copy is on disk: " * basename(stale_copy) * ". The last scan run " *
        "refused its own acceptance table and moved the published one aside. It is not " *
        "read here, because no run vouches for it: re-run the scan, or pass the copy " *
        "explicitly with --scan and read every number as unlicensed." :
        "Run scripts/background_pqq_scan.jl, or re-run this script with --bg-exclusion " *
        "none and treat every number as a lower bound on the contrast."
    bar = repeat('!', 88)
    error("\n" * bar * "\nMISSING BACKGROUND SCAN: " * SCAN * "\n" *
          "The background frame has NOT been cleared of PQQ-ADH carriers. Every contrast " *
          "would be biased toward the null by carriers sitting in the comparator arm, and " *
          "the stronger background code would not be licensed.\n" * hint * "\n" * bar)
end

scan === nothing || scan_inventory!(scan)

bg_excluded = exclusion_taxids(scan, BG_EXCL)
bg_excluded === nothing && error("""
    --bg-exclusion $BG_EXCL cannot be computed from $SCAN: the table has no column for it.
    Re-run scripts/background_pqq_scan.jl, or choose a mode the table can answer.""")

rec!("config", CONFIG)
rec!("bg_exclusion_mode", BG_EXCL)
rec!("bg_exclusion_mode_source", BG_EXCL_SOURCE)
rec!("bg_exclusion_mode_role", BG_EXCL_ROLE)
rec!("bg_exclusion_rule", exclusion_note(BG_EXCL))
rec!("chemistry_column", String(CHEM_COL))

# ---------------------------------------------------------------- the move-in pool
# Only the level that moves taxa reads a pool, and only that level emits these keys: a run of
# any other mode moves nobody, so a movein_ key in its table would describe something the run
# did not do. The pool is read BEFORE the arms are built, so a pool that cannot be answered
# stops the run rather than turning it into the mode it falls back on.
moved_in = Set{Int}()
if BG_EXCL == MOVEIN_MODE
    pool = movein_pool(; scan = scan, path = MOVEIN_POOL_PATH)
    pool.taxids === nothing && error("""
        --bg-exclusion $MOVEIN_MODE has no pool: $SCAN cannot answer the
        $(MOVEIN_POOL_RULES[2]) rule, and no census hit table was supplied.
        Pass --movein-pool <census hit table>, or re-run the background scan.""")
    moved_in = pool.taxids
    registered = pool.rule == REGISTERED_MOVEIN_POOL
    rec!("movein_pool_rule", pool.rule)
    rec!("movein_pool_rule_source", pool.rule_source)
    rec!("movein_pool_rule_role", registered ? :registered : :descriptive)
    rec!("movein_pool_rule_note", movein_pool_note(pool.rule))
    rec!("movein_pool_rules_available", joined(MOVEIN_POOL_RULES))
    rec!("movein_pool_source", pool.source)
    rec!("movein_pool_sha256", pool.sha256 === nothing ? missing : pool.sha256)
    rec!("movein_pool_rows_read", pool.rows)
    rec!("movein_pool_class_rows", pool.class_rows)
    rec!("movein_pool_rows_without_a_taxon", pool.rows_without_a_taxon)
    rec!("movein_pool_columns", joined(pool.columns))
    rec!("movein_pool_taxa", length(moved_in))
    # The OTHER rule's size and its overlap with the pool that ran, whenever this run can
    # compute them. Two rules that name a lanthanide-type carrier differently will disagree
    # about which taxa those are, and how far apart they sit is the only comparison a single
    # run can make; under the scan rule the two are the same set, which is what the overlap
    # then says. An unanswerable count is missing, never zero.
    let alt = exclusion_taxids(scan, "lntype")
        rec!("movein_pool_taxa_under_the_scan_rule", alt === nothing ? missing : length(alt))
        rec!("movein_pool_taxa_in_both_rules",
             alt === nothing ? missing : length(intersect(moved_in, alt)))
    end
    registered ||
        @warn("the move-in level is running on the pool rule this repository can answer " *
              "alone, not the rule the recorded level names; its size is a different " *
              "measurement and the output records it as descriptive",
              rule = pool.rule, pool_taxa = length(moved_in))
end

# ---------------------------------------------------------------- census classes
# The census table is opened once, here, because the arms need it: the disclosed calcium
# comparator is defined by the classes of a taxon's census sequences. The two disclosure
# keys it also answers — how many eukaryotic candidates were left out — are emitted further
# down, with the other provenance, so this read does not move them.
census = read_census_classes(CENSUS)
census_classes = census.classes
ca_class_taxids = census.ca_class_taxids

# ---------------------------------------------------------------- the run's own arms
A = analyse(calls, CONFIG, CHEM_COL, bg_excluded; ca_class_taxids = ca_class_taxids,
            moved_in = moved_in)
df = A.df; usable = A.usable
fg = A.fg; bgm = A.bgm; ac = A.ac; nl = A.nl; ca = A.ca
fg_phy = A.fg_phy; bg_phy = A.bg_phy; ac_phy = A.ac_phy; nl_phy = A.nl_phy
ca_phy = A.ca_phy
rows = A.rows
n1 = A.n1; n2 = A.n2
gfg = A.gfg; gbg = A.gbg

# ---------------------------------------------------------------- arm inventory
rec!("n1_taxa", n1)
rec!("n1_candidate_taxa", count(A.ln_cand))
rec!("n2_taxa", n2)
rec!("frame_size", count(A.frame))
rec!("frame_taxa_with_usable_term", count(A.frame .& usable))
rec!("census_overlap_in_frame", count(A.in_pool .& A.in_cen))
rec!("bg_removed_by_scan", count(A.frame .& in.(df.taxid, Ref(bg_excluded))))
rec!("bg_removed_by_scan_with_usable_term",
     count(A.frame .& usable .& in.(df.taxid, Ref(bg_excluded))))
rec!("bg_removed_by_scan_by_phylum",
     tally_str(tally(sv.(df.phylum[A.frame .& in.(df.taxid, Ref(bg_excluded))]))))
rec!("n_allcarrier_taxa", nrow(ac))
rec!("n_nonln_carrier_taxa", nrow(nl))
rec!("n_sequences_fg", length(A.seq_idx))
rec!("n_genera_fg", length(gfg.genus)); rec!("n_genera_fg_dropped_no_genus", gfg.n_dropped)
rec!("n_genera_bg", length(gbg.genus)); rec!("n_genera_bg_dropped_no_genus", gbg.n_dropped)
rec!("fg_genus_labels_from_lineage", count(sv.(fg.genus) .!= ""))
rec!("fg_genus_labels_from_organism_name",
     count((sv.(fg.genus) .== "") .& (A.fg_gen .!= "")))
rec!("bg_genus_labels_from_lineage", count(sv.(bgm.genus) .!= ""))
rec!("bg_genus_labels_from_organism_name",
     count((sv.(bgm.genus) .== "") .& (A.bg_gen .!= "")))
# The two routes are not the same route, so the agreement between them where BOTH exist is
# reported: it is the only offline measure of how much the background genus labels can be
# trusted, and the genus unit's caveat has to be quantified rather than asserted.
let have = sv.(fg.genus) .!= "", alt = genus_of.(missing, fg.organism)
    both = have .& (alt .!= "")
    rec!("genus_route_agreement_pct",
         100 * count(both .& (alt .== sv.(fg.genus))) / max(count(both), 1))
    rec!("genus_route_agreement_denominator", count(both))
end
rec!("n_archaea_fg", count(sv.(fg.kingdom) .== "Archaea"))
rec!("n_bacteria_fg", count(sv.(fg.kingdom) .== "Bacteria"))

# ---------------------------------------------------------------- what the move moved
# The level's own arithmetic, so n1 and n2 can be taken apart from this table alone: the arm
# it started from, what it added, and what the addition cost the background. Every one of
# these is counted over the arms that were actually built.
if BG_EXCL == MOVEIN_MODE
    rec!("movein_pool_taxa_in_frame", count(A.named_in_pool))
    rec!("movein_pool_taxa_outside_frame", length(moved_in) - count(A.named_in_pool))
    rec!("movein_moved_taxa", A.n_moved)
    rec!("movein_pool_taxa_in_frame_not_prokaryotic", count(A.named_in_pool) - A.n_moved)
    rec!("movein_moved_taxa_annotated", A.n_moved_usable)
    rec!("movein_n1_before_the_move", n1 - A.n_moved_usable)
    rec!("movein_frame_size_before_the_move", count(A.frame_scanned))
    # What the move takes out of the background is not what the exclusion takes out: a moved
    # taxon the scan never flagged would otherwise have stayed in the analysed arm, and that
    # difference is the only part of n2 this level owns.
    rec!("movein_background_taxa_the_move_removed",
         count(A.moved .& usable .& .!in.(df.taxid, Ref(bg_excluded))))
    rec!("movein_moved_taxa_also_excluded_by_the_scan",
         count(A.moved .& in.(df.taxid, Ref(bg_excluded))))
    # The census-derived columns, and what this level does with taxa that have none of them.
    # The per-sequence unit counts census sequences and the moved taxa carry none, so they
    # enter the taxa and genus units and not the sequence unit; their domain is read off
    # `domain`, because `kingdom` is a census lineage column and is empty for them, which is
    # why n_archaea_fg and n_bacteria_fg above do not add up to n1 in a run of this level.
    rec!("movein_moved_census_sequences", A.n_moved_census_sequences)
    rec!("movein_sequence_units_from_moved_taxa",
         sum(coalesce.(df.n_ln_seqs[A.moved .& usable], 0); init = 0))
    rec!("movein_moved_taxa_with_a_kingdom_label", count(sv.(df.kingdom[A.moved]) .!= ""))
    rec!("movein_moved_taxa_by_domain",
         tally_str(tally(sv.(hascol(df, :domain) ? df.domain[A.moved] :
                             df.kingdom[A.moved]))))
    rec!("movein_moved_taxa_annotated_by_phylum",
         tally_str(tally(sv.(df.phylum[A.moved .& usable])), limit = MAX_LISTED_TAXA))
    rec!("movein_comparator_arms_note",
         "the all-carrier, non-Ln and Ca-class arms are NOT rebuilt by this level: their " *
         "membership rules are rules about a taxon's census sequences, and a moved taxon " *
         "has none. They do lose the moved taxa from the background they share with the Ln " *
         "arm, which is this level's doing and is counted above")
    rec!("movein_anchor_note",
         "the anchor moves nobody: it is the earlier specification, so a run of this level " *
         "reads its category partition against an unmoved anchor")
end

rec!("pct_$(K_PHYLUM)_fg", pseudo(fg_phy, PHYLUM))
rec!("pct_$(K_PHYLUM)_bg", pseudo(bg_phy, PHYLUM))
rec!("pct_$(K_PHYLUM)_bg_before_scan",
     pseudo(sv.(df.phylum[A.frame .& usable]), PHYLUM))
rec!("pct_$(K_PHYLUM)_allcarriers", pseudo(ac_phy, PHYLUM))
rec!("pct_$(K_PHYLUM)_nonln", pseudo(nl_phy, PHYLUM))
# The stratum count is a property of the UNIT, not of the category, so any drawn taxa-unit row
# carries it; the sharing category's row is the one this run is guaranteed to have asked for,
# which is why it is read there rather than at a category name of its own.
rec!("phylum_strata_used_fg_vs_bg",
     something(pick(rows, LN_VS_BG, "taxa", SHARE_CAT), (n_strata = missing,)).n_strata)

# How much of the standardised estimate rides on how little background. Removing carriers
# takes most of one phylum out of the background while leaving the foreground untouched,
# so the stratum carrying the largest foreground weight can end up matched against a small
# background cell. The estimator's variance already pays for that, but the size of the bet
# should be visible rather than buried in an interval width. The collapse rule is the
# estimator's own (its two thresholds are read from the module, not restated here).
let cf_ = tally(fg_phy), cb_ = tally(bg_phy)
    keep(p) = get(cf_, p, 0) >= HabitatStats.MIN_FG_K && get(cb_, p, 0) >= HabitatStats.MIN_BG_K
    lab(p) = keep(p) ? p : "Other phyla"
    wf = tally(lab.(fg_phy)); wb = tally(lab.(bg_phy))
    top = argmax(k -> wf[k], collect(keys(wf)))
    rec!("standardisation_heaviest_stratum", top)
    rec!("standardisation_heaviest_stratum_weight", wf[top] / max(n1, 1))
    rec!("standardisation_heaviest_stratum_bg_taxa", get(wb, top, 0))
    rec!("standardisation_strata_collapsed_to_other",
         joined(sort([p for p in keys(cf_) if !keep(p) && !isempty(p)])))
end

# ---------------------------------------------------------------- coverage
cf = coverage!("fg", A.ln_cand, usable, "Ln")
cb = coverage!("bg_frame", A.frame, usable, "background frame")
coverage!("allcarrier", A.all_cand, usable, "all carriers")
coverage!("nonln", A.nonln_cand, usable, "non-Ln carriers")
coverage!("caonly", A.ca_cand, usable, "Ca-class carriers")
# crude() is the module's Fisher exact test on the 2x2 annotated/not table; using it here
# keeps this script free of any test of its own.
covp = crude(cf.nok, cf.ncand, cb.nok, cb.ncand).fisher_p
rec!("coverage_fisher_p", covp)
rec!("coverage_differs_at_1pct", covp < 0.01)
# the same comparison after the carrier exclusion, since that is the arm actually analysed
rec!("bg_analysed_coverage_numerator", n2)
rec!("bg_analysed_coverage_denominator", count(A.bg_keep))
rec!("bg_analysed_coverage_pct", 100 * n2 / max(count(A.bg_keep), 1))
rec!("coverage_fisher_p_after_scan",
     crude(cf.nok, cf.ncand, n2, count(A.bg_keep)).fisher_p)

# ---------------------------------------------------------------- why taxa were dropped
dfg = dropped!("fg", df, usable, A.ln_cand)
dbg = dropped!("bg_frame", df, usable, A.frame)
dropped!("allcarrier", df, usable, A.all_cand)
dropped!("nonln", df, usable, A.nonln_cand)
dropped!("caonly", df, usable, A.ca_cand)

# ---------------------------------------------------------------- evidence behind a call
evidence!("fg", fg)
evidence!("bg", bgm)
evidence!("allcarrier", ac)
evidence!("nonln", nl)
evidence!("caonly", ca)

# The eukaryotic disclosure from the census read above.
rec!("eukaryotic_ln_sequences_excluded", census.euk_seqs)
rec!("eukaryotic_ln_taxa_excluded", census.euk_taxa)

# ---------------------------------------------------------------- provenance
if isfile(MANIFEST)
    mf = CSV.read(MANIFEST, DataFrame)
    rels = sort(unique(sv.(mf.omnicrobe_release)))
    obts = sort(unique(sv.(mf.obt_release)))
    rec!("omnicrobe_release", first(rels))
    rec!("omnicrobe_release_values", joined(rels))
    rec!("obt_release", first(obts))
    rec!("obt_release_values", joined(obts))
    # A single access date would have to be one of many: the cache was seeded in bulk and
    # topped up later, so the latest date describes a small minority of the bodies behind
    # these numbers. The range and the tally say what a single date cannot.
    fetched = sort(unique(sv.(mf.fetched_at)))
    rec!("omnicrobe_fetched_earliest", first(fetched))
    rec!("omnicrobe_fetched_latest", last(fetched))
    rec!("omnicrobe_fetched_date_tally", tally_str(tally(sv.(mf.fetched_at))))
    rec!("omnicrobe_cached_taxa", nrow(mf))
end
# The release and threshold records as written, from the same single read the licence is
# decided on.
let up = records.release_kv, th = records.threshold_kv
    rec!("uniprot_release", get(up, "X-UniProt-Release", "unknown"))
    rec!("uniprot_release_date", get(up, "X-UniProt-Release-Date", "unknown"))
    rec!("hmm_model", get(th, "hmm", "unknown"))
    rec!("hmm_model_sha256", get(th, "hmm_sha256", "unknown"))
    rec!("hmm_threshold_bits", get(th, "threshold_bits", "unknown"))
    rec!("hmm_threshold_rule", get(th, "rule", "unknown"))
    rec!("hmm_threshold_source", get(th, "threshold_source", "unknown"))
    rec!("hmm_active_site_column", get(th, "active_site_column", "unknown"))
    rec!("hmm_calibration_recovery", get(th, "recovery", "unknown"))
end
# The background pool the frame was drawn from: its seed, read from the generator, and the two
# tables it wrote, by digest. See background_pool_provenance!.
background_pool_provenance!(rec!)
rec!("script", "scripts/habitat_enrichment.jl")
rec!("run_at", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))

# ---------------------------------------------------------------- sensitivity on the arm
# Every mode that is NOT in force is reported here and only here, so the choice of exclusion
# can be audited without a non-primary number ever reaching a headline key. The narrow modes
# answer "which background taxa carry a lanthanide-type enzyme", or one whose active site
# reads either metal, which are different questions from "which background taxa carry this
# enzyme family at all".
let unavailable = String[]
    for mode in SENSITIVITY_MODES
        mode == BG_EXCL && continue
        set = exclusion_taxids(scan, mode)
        if set === nothing
            push!(unavailable, mode)
            rec!("sensitivity_n2_$(mode)", missing)
            rec!("sensitivity_bg_removed_$(mode)", missing)
            rec!("sensitivity_pct_$(K_PHYLUM)_bg_$(mode)", missing)
            continue
        end
        keep = A.frame .& .!in.(df.taxid, Ref(set))
        rec!("sensitivity_n2_$(mode)", count(keep .& usable))
        rec!("sensitivity_bg_removed_$(mode)", count(A.frame .& in.(df.taxid, Ref(set))))
        rec!("sensitivity_pct_$(K_PHYLUM)_bg_$(mode)",
             pseudo(sv.(df.phylum[keep .& usable]), PHYLUM))
    end
    rec!("sensitivity_modes_unavailable", joined(unavailable))
end

# ---------------------------------------------------------------- the recorded criteria
criterion!("bg_size", PLAN_BG_SIZE_BELOW, n2, :below,
           "an upper bound on the size of the analysed background, written down in the " *
           "plan for this step before it ran and supplied on the command line; disclosed " *
           "beside the value, it gates nothing here")
criterion!("bg_$(K_PHYLUM)_pct", PLAN_BG_PHYLUM_PCT_MAX, look("pct_$(K_PHYLUM)_bg"),
           :at_most,
           "an upper bound on the background's $(PHYLUM) share in per cent, from the " *
           "same plan and supplied the same way; disclosed beside the value, it gates " *
           "nothing here")

# ---------------------------------------------------------------- the supplied selection
# Emitted as DATA, so the table says which selection produced it. The digest is of the criteria
# file's own bytes: two runs whose numbers differ can be told apart by what they were asked
# about before a single estimate is compared.
rec!("selection_file", basename(SEL.path))
rec!("selection_sha256", SEL.digest)
rec!("selection_roles", joined(r.role for r in SELECTION_ROLES))
for r in SELECTION_ROLES
    v = getfield(SEL, Symbol(r.role))
    rec!("selection_$(r.role)", v isa Tuple ? joined(v) : v)
end
rec!("selection_key_prefixes",
     joined("$(r.role)=$(getfield(SEL.key, Symbol(r.role)))"
            for r in SELECTION_ROLES if r.keyed))
rec!("categories_governing_roles", joined(GOVERNING_ROLES))

# ---------------------------------------------------------------- rule-keyed quantities
# The recorded outcome rules turn these into wording elsewhere. They are emitted so each
# rule is decidable from this table alone, as NEUTRAL CODES rather than phrases: a code is
# diffable across runs and carries no prose into this repository. Nothing below adjusts a
# threshold to suit the number that came out, and every code is accompanied by its inputs.
rec!("rules_interval_reading",
     "standardised bounds decide every code below, as the recorded rules specify for the " *
     "intervals they name; each crude-reading counterpart is emitted beside its code as a " *
     "descriptive line only and decides nothing")

# -- the depletion-tier rule's category. The primary row is the chemistry contrast on the
#    taxa unit. The two strong tiers are DEFINED by the one-at-a-time sensitivity grid, which
#    is a separate step, so this run can only decide the tier that needs no grid.
rth = pick(rows, LN_VS_BG, "taxa (chemistry)", TIER_CAT)
emit_row!("$(K_TIER)_primary", rth)
emit_row!("$(K_TIER)_allcarriers", pick(rows, ALL_VS_BG, "taxa (chemistry)", TIER_CAT))
emit_row!("$(K_TIER)_nonln", pick(rows, NON_VS_BG, "taxa (chemistry)", TIER_CAT))
# This script computes no sensitivity grid, so it has no variants to offer. Handing the rule
# what it has — nothing — is what makes it withhold the two strong tiers instead of reading
# them off the primary row, and the number of variants goes into the table so the withholding
# is visible rather than looking like a missing result. The reason is keyed on the tier the
# rule returned, so it cannot describe a tier the rule did not decide.
let variants = nothing,
    (tier, breakers) = rth === nothing ? (:undecidable, String[]) :
                       depletion_tier(rth.ci_lo_std, rth.ci_hi_std, variants),
    crude_tier = rth === nothing ? :undecidable :
                 first(depletion_tier(rth.ci_lo, rth.ci_hi, variants))
    rec!("$(K_TIER)_tier_code", tier)
    rec!("$(K_TIER)_tier_code_crude", crude_tier)
    rec!("$(K_TIER)_tier_breakers", joined(breakers))
    rec!("$(K_TIER)_sensitivity_variants_supplied", 0)
    rec!("$(K_TIER)_tier_reason",
         rth === nothing ? "no $(TIER_CAT) row in this run" :
         tier === :undecidable && zero_std(rth) === :undecidable ?
             "the primary standardised interval was not computed" :
         tier === :undecidable ?
             "the primary interval excludes zero on the depleted side; the stronger tiers " *
             "are defined by the one-at-a-time sensitivity grid, which this run does not " *
             "compute" :
         tier === :t3 ?
             "decided by the primary interval alone, which needs no sensitivity grid" :
             "decided by the sensitivity grid supplied with the primary interval")
end

# -- the genus-clause rule's category. See genus_clause_code and category_code above the
#    guard.
rhp_t = pick(rows, LN_VS_BG, "taxa", GENUS_CAT)
rhp_g = pick(rows, LN_VS_BG, "genera", GENUS_CAT)
emit_row!("$(K_GENUS)_taxa", rhp_t)
emit_row!("$(K_GENUS)_genera", rhp_g)
rec!("$(K_GENUS)_genus_clause_code", genus_clause_code(rhp_g, :std))
rec!("$(K_GENUS)_genus_clause_code_crude", genus_clause_code(rhp_g, :crude))
rec!("$(K_GENUS)_category_code", category_code(rhp_t, :std))
rec!("$(K_GENUS)_category_code_crude", category_code(rhp_t, :crude))

# -- the sharing rule's category. The registered comparator is the non-Ln carrier arm; the
#    all-carrier ratio has one arm inside the other, and is labelled so; the calcium-class
#    arm is a disclosure beside the registered one.
rte_ln  = pick(rows, LN_VS_BG,  "taxa", SHARE_CAT)
rte_non = pick(rows, NON_VS_BG, "taxa", SHARE_CAT)
rte_all = pick(rows, ALL_VS_BG, "taxa", SHARE_CAT)
rte_ca  = pick(rows, CA_VS_BG,  "taxa", SHARE_CAT)
emit_row!("$(K_SHARE)_ln", rte_ln)
emit_row!("$(K_SHARE)_nonln", rte_non)
emit_row!("$(K_SHARE)_allcarriers", rte_all)
emit_row!("$(K_SHARE)_caonly", rte_ca)

const CI_METHOD = "percentile interval over $(NBOOT) unit resamples, seed $(SEED); " *
                  "shared background resample"

r_non = ratio_of(rte_non, rte_ln; standardised_scale = true)
r_all = ratio_of(rte_all, rte_ln; standardised_scale = true)
lo, hi, nok = ratio_interval(nl_phy, A.nlw, fg_phy, A.fgw, bg_phy, A.bgw, SHARE_CAT;
                             nboot = NBOOT, seed = SEED)
cmp_std = zero_std(rte_non)
share_std = sharing_wording_code(r_non, lo, hi; comparator_excludes_zero = cmp_std)
rec!("sharing_ratio_nonln_over_ln", r_non)
rec!("sharing_ratio_ci_lo", lo)
rec!("sharing_ratio_ci_hi", hi)
rec!("sharing_ratio_ci_replicates", nok)
rec!("sharing_ratio_ci_method", CI_METHOD)
rec!("sharing_comparator_excludes_zero", cmp_std)
rec!("sharing_code", share_std.code)
rec!("sharing_rule_used", share_std.rule_used)
# The crude-reading counterpart. There is no resampled interval on the crude scale, so the
# point form decides it; the two readings need not agree, and that is the disclosure.
let r_crude = ratio_of(rte_non, rte_ln; standardised_scale = false),
    share_crude = sharing_wording_code(r_crude, NaN, NaN;
                                       comparator_excludes_zero = zero_crude(rte_non))
    rec!("sharing_ratio_nonln_over_ln_crude", r_crude)
    rec!("sharing_comparator_excludes_zero_crude", zero_crude(rte_non))
    rec!("sharing_code_crude", share_crude.code)
    rec!("sharing_rule_used_crude", share_crude.rule_used)
end
rec!("sharing_ratio_allcarriers_over_ln", r_all)
rec!("sharing_ratio_allcarriers_note",
     "one arm inside the other — the all-carrier arm contains the Ln arm, so this fraction " *
     "compares an arm with a superset of itself")
rec!("sharing_ratio_allcarriers_ci",
     "not computed — the two arms overlap, so resampling them independently would " *
     "understate their correlation and return an interval that is too narrow")

# The disclosed calcium-class comparator: the same ratio, the same resampling routine and
# seed, the same grading function, on an arm whose membership is the PRESENCE of a
# calcium-class census sequence rather than the absence of a lanthanide-class one. It is
# emitted beside the registered comparator's keys and never instead of them.
let r_ca = ratio_of(rte_ca, rte_ln; standardised_scale = true),
    (calo, cahi, cank) = ratio_interval(ca_phy, A.caw, fg_phy, A.fgw, bg_phy, A.bgw,
                                        SHARE_CAT; nboot = NBOOT, seed = SEED),
    cmp_ca = zero_std(rte_ca)
    share_ca = sharing_wording_code(r_ca, calo, cahi; comparator_excludes_zero = cmp_ca)
    rec!("sharing_ratio_caonly_over_ln", r_ca)
    rec!("sharing_ratio_caonly_ci_lo", calo)
    rec!("sharing_ratio_caonly_ci_hi", cahi)
    rec!("sharing_ratio_caonly_ci_replicates", cank)
    rec!("sharing_ratio_caonly_ci_method", CI_METHOD)
    rec!("sharing_comparator_caonly_excludes_zero", cmp_ca)
    rec!("sharing_code_caonly", share_ca.code)
    rec!("sharing_rule_used_caonly", share_ca.rule_used)
    let r_ca_crude = ratio_of(rte_ca, rte_ln; standardised_scale = false),
        share_ca_crude = sharing_wording_code(r_ca_crude, NaN, NaN;
                                    comparator_excludes_zero = zero_crude(rte_ca))
        rec!("sharing_ratio_caonly_over_ln_crude", r_ca_crude)
        rec!("sharing_comparator_caonly_excludes_zero_crude", zero_crude(rte_ca))
        rec!("sharing_code_caonly_crude", share_ca_crude.code)
        rec!("sharing_rule_used_caonly_crude", share_ca_crude.rule_used)
    end
    rec!("sharing_caonly_role",
         "disclosed beside the registered comparator, never in its place: its arm is " *
         "defined by the presence of a calcium-class census sequence, which is a different " *
         "membership rule from the registered comparator's absence of a lanthanide-class one")
end

# What the comparator arms are made of. The registered arm's membership rule is the ABSENCE
# of a lanthanide-type census sequence, which is not the same thing as the presence of a
# calcium-type one, and a code above depends on the arm being what its name says. Redefining
# that arm would need a recorded decision, so the composition is disclosed instead of
# changed, and the calcium-class arm is reported separately above.
let ids = [taxid_of(t) for t in nl.taxid]
    classes = [i === nothing ? Set{String}() : get(census_classes, i, Set{String}())
               for i in ids]
    rec!("comparator_arm_taxa", nrow(nl))
    rec!("comparator_arm_taxa_with_calcium_class", count(c -> "Calcium" in c, classes))
    rec!("comparator_arm_taxa_gap_or_unknown_class_only",
         count(c -> !isempty(c) && issubset(c, Set(["Gap", "Unknown"])), classes))
    rec!("comparator_arm_taxa_without_classified_sequence", count(isempty, classes))
    rec!("comparator_arm_class_source", isfile(CENSUS) ? basename(CENSUS) : "unavailable")
    rec!("caonly_arm_taxa", nrow(ca))
    rec!("caonly_arm_candidate_taxa", count(A.ca_cand))
    rec!("caonly_arm_census_taxa_available",
         ca_class_taxids === nothing ? missing : length(ca_class_taxids))
    rec!("pct_$(K_PHYLUM)_caonly", pseudo(ca_phy, PHYLUM))
end

# -- the background code: is the comparator arm really free of the carriers it is supposed
#    to be free of? Five inputs, each measured over the arm that is ACTUALLY ANALYSED or
#    read from the persisted records, all assembled by scan_evidence_and_licence and all
#    emitted so the code can be re-decided from the table.
L = scan_evidence_and_licence(A; table_read = scan_read, table_state = scan_state,
                              records = records)
let fc = L.frame, left = L.leftovers, ev = L.evidence,
    thr_rec = records.threshold, rec_rec = records.recall, acc_tbl = records.acceptance,
    rk = records.recall_kv
    rec!("background_scan_frame_matches_analysis_frame",
         fc.matches === nothing ? missing : fc.matches)
    rec!("background_analysis_frame_taxa_missing_from_scan", length(fc.missing_from_scan))
    rec!("background_analysis_frame_taxa_missing_from_scan_list", listed(fc.missing_from_scan))
    rec!("background_scan_taxa_outside_analysis_frame", length(fc.extra_in_scan))
    rec!("background_scan_taxa_outside_analysis_frame_list", listed(fc.extra_in_scan))

    rec!("background_carriers_left_hmm_hit",
         left.hmm_hit === nothing ? missing : left.hmm_hit)
    # Under the any-hit mode the exclusion predicate and this count read the same column, so
    # the count is an identity and not a check. Saying so in the table is the difference
    # between a guard that cannot fail and a guard that is known not to be one. The move-in
    # level removes by the same predicate, so the identity holds there too — and saying so
    # matters more in a level whose arms moved, where a zero could be read as evidence that
    # the move cleaned the background.
    rec!("background_carriers_left_hmm_hit_definitional",
         BG_EXCL == PRIMARY_MODE || BG_EXCL == MOVEIN_MODE)
    rec!("background_carriers_left_ipr017512_ge500",
         left.ipr017512_ge500 === nothing ? missing : left.ipr017512_ge500)
    rec!("background_carriers_left_union", left.union === nothing ? missing : left.union)

    # The persisted records, and what they vouch for. A recall that is only printed cannot
    # be checked by anything downstream; one that is read but believed is no better, so the
    # verdict the licence uses is the module's RE-DERIVATION from the record's own counts, on
    # their lower confidence bound, and the record's own fields are emitted as written beside
    # the rate and bound those counts give, as claims to be compared.
    rec!("background_prefilter_recall_record_present", !isempty(rk))
    rec!("background_prefilter_recall_pass",
         rec_rec.pass === nothing ? missing : rec_rec.pass)
    rec!("background_prefilter_recall_rate", get(rk, "rate", missing))
    rec!("background_prefilter_recall_rate_derived",
         ev.recall_rate === nothing ? missing : ev.recall_rate)
    rec!("background_prefilter_recall_bar", get(rk, "bar", missing))
    rec!("background_prefilter_recall_sampled", get(rk, "n_sampled", missing))
    rec!("background_prefilter_recall_full_flagged", get(rk, "n_full_flagged", missing))
    rec!("background_prefilter_recall_recovered", get(rk, "n_recovered", missing))
    rec!("background_prefilter_recall_never_searched_in_sample",
         get(rk, "n_never_searched_in_sample", missing))
    rec!("background_prefilter_recall_missed_taxids", get(rk, "missed_taxids", missing))
    rec!("background_prefilter_recall_lower_bound", get(rk, "lower_bound_95", missing))
    rec!("background_prefilter_recall_lower_bound_derived",
         ev.recall_lower_bound === nothing ? missing : ev.recall_lower_bound)
    rec!("background_recall_pass_derived",
         ev.recall_pass_derived === nothing ? missing : ev.recall_pass_derived)
    rec!("background_recall_pass_rule",
         "met only when the full-proteome search flagged at least one sampled taxon and the " *
         "exact one-sided 95 % Clopper-Pearson lower bound of recovered over flagged, taken " *
         "from the record's counts, reaches the record's bar; the record's own pass field is " *
         "the scan's gate on the rate and is checked against the rate only")

    # The scan evidence, input by input, so a reader can re-decide it: the digest of the
    # bytes parsed here, the digest each record claims to describe, the acceptance verdicts,
    # the counts the acceptance table was held to, and the file state of the table itself.
    rec!("background_scan_table_state", scan_state.state)
    rec!("background_scan_table_pending_beside_it", scan_state.pending_beside_it)
    rec!("background_scan_table_displaced_copy_on_disk", scan_state.displaced_copy_on_disk)
    rec!("background_scan_sha256", scan_sha === nothing ? missing : scan_sha)
    rec!("background_scan_record_threshold_present", records.threshold_present)
    rec!("background_scan_record_threshold_sha256",
         thr_rec.scan_sha256 === nothing ? missing : thr_rec.scan_sha256)
    rec!("background_scan_record_threshold_run_id",
         thr_rec.run_id === nothing ? missing : thr_rec.run_id)
    rec!("background_scan_record_recall_present", records.recall_present)
    rec!("background_scan_record_recall_sha256",
         rec_rec.scan_sha256 === nothing ? missing : rec_rec.scan_sha256)
    rec!("background_scan_record_recall_run_id",
         rec_rec.run_id === nothing ? missing : rec_rec.run_id)
    rec!("background_scan_record_acceptance_present", records.acceptance_present)
    rec!("background_scan_record_acceptance_read_error", records.acceptance_error)
    rec!("background_scan_acceptance_rows", acc_tbl === nothing ? missing : nrow(acc_tbl))
    for (label, word) in (("pass", "PASS"), ("miss", "MISS"), ("disclosed", "DISCLOSED"))
        rec!("background_scan_acceptance_$(label)",
             acc_tbl === nothing ? missing : count(==(word), acc_tbl.verdict))
    end
    rec!("background_scan_acceptance_other_verdicts",
         acc_tbl === nothing ? missing :
         count(v -> !(v in ("PASS", "MISS", "DISCLOSED")), acc_tbl.verdict))
    rec!("background_scan_acceptance_required_gated_checks", length(REQUIRED_GATED_CHECKS))
    rec!("background_scan_acceptance_flagged_row_value",
         something(acceptance_cell(acc_tbl, FLAGGED_CHECK, :value), missing))
    rec!("background_scan_table_taxa_flagged",
         something(L.table_counts.n_flagged, missing))
    rec!("background_scan_acceptance_family_row_denominator",
         something(acceptance_cell(acc_tbl, FAMILY_CHECK, :denominator), missing))
    rec!("background_scan_table_taxa_in_family_ge500",
         something(L.table_counts.n_family, missing))
    rec!("background_scan_record_acceptance_sha256",
         acc_tbl === nothing ? missing : joined(sort(unique(acc_tbl.scan_sha256))))
    rec!("background_scan_record_acceptance_run_id",
         acc_tbl === nothing ? missing : joined(sort(unique(acc_tbl.run_id))))
    # One identifier per run is the producer's contract, and three records that agree on the
    # digest but not on the run_id describe one table written by two invocations. The module
    # is not asked about it, because tying a record to the TABLE is what licenses the code;
    # the disagreement is disclosed so it cannot pass unnoticed either.
    rec!("background_scan_record_run_ids_agree",
         let ids = filter(!isempty,
                          vcat(something(thr_rec.run_id, ""), something(rec_rec.run_id, ""),
                               acc_tbl === nothing ? String[] : unique(acc_tbl.run_id)))
             isempty(ids) ? missing : length(unique(ids)) == 1
         end)
    rec!("background_scan_evidence_ok", L.ok)
    rec!("background_scan_records_vouch_for_table", L.records_ok)
    rec!("background_scan_evidence_reasons", joined(L.reasons))

    rec!("background_removal_complete", L.removal.complete)
    rec!("background_removal_incomplete_reasons", joined(L.removal.reasons))
    rec!("background_licence_code", L.code)
    if fc.matches === false
        @warn("the scanned taxon set is not the analysed frame: the leftover counts " *
              "can only see the taxa the scan holds a row for",
              frame_taxa_absent_from_scan = length(fc.missing_from_scan),
              scan_taxa_outside_the_frame = length(fc.extra_in_scan))
    end
    if !L.records_ok
        @warn("the scan's persisted records do not vouch for the table that was read",
              reasons = L.reasons)
    elseif !L.ok
        @info("the scan's records vouch for the table, but the recall verdict derived from " *
              "them does not reach the bar, so the stronger background statement is withheld",
              reasons = L.reasons)
    end
end

# -- categories whose primary interval excludes zero, and which of them are new. "New" is
#    measurable only against the earlier result a category would be new relative to, so the
#    ANCHOR is computed in this same run — same builder, earlier settings, no exclusion —
#    and the partition is taken against it. An expression comparing the run against a
#    hard-coded list of names nominates a category even on the anchor itself, where by
#    construction nothing can be new.
primary_std   = interval_sets(rows, :std)
primary_crude = interval_sets(rows, :crude)
sig_std_primary   = primary_std.excludes
sig_crude_primary = primary_crude.excludes

anchor_present = ANCHOR_CONFIG in CONFIGS_AVAILABLE
anchor_status = :undecidable_config_absent
anchor_sig = nothing              # what the partition may use: nothing = no anchor at all
anchor_sig_seen = String[]        # what the anchor computed, emitted either way
anchor_sig_crude = String[]
anchor_labels = String[]
anchor_rows = NamedTuple[]
ref_sig = String[]
ref_rows_read = missing
if anchor_present
    B = analyse(calls, ANCHOR_CONFIG, ANCHOR_CHEM_COL, Set{Int}();
                ca_class_taxids = ca_class_taxids)
    anchor_rows = B.rows
    anchor_std = interval_sets(B.rows, :std)
    anchor_crude = interval_sets(B.rows, :crude)
    anchor_sig_seen  = anchor_std.excludes
    anchor_sig_crude = anchor_crude.excludes
    anchor_labels = anchor_std.labels
    # An anchor is only an anchor if it is the same measurement on the same frame. If the
    # frame moved, or if a row of the anchor panel has no interval, its category set is not
    # the earlier result and must not be treated as one.
    # The frame the two runs drew from, not the arm each of them kept: the anchor applies no
    # exclusion and moves nobody, so comparing it with the run's ANALYSED frame would call
    # every move-in run a frame mismatch and withhold its partition for doing exactly what
    # the level says it does. What must still agree is the pool both specifications sampled.
    frames_match = Set(Int.(B.df.taxid[B.frame_scanned])) ==
                   Set(Int.(df.taxid[A.frame_scanned]))
    ref_supplied = !isempty(ANCHORREF)
    ref_exists = ref_supplied && isfile(ANCHORREF)
    ref_readable = true
    if ref_exists
        ref = read_reference_rows(ANCHORREF)
        if ref === nothing
            ref_readable = false
        else
            # The reference's rows in the shape a row panel has, so ONE function reads its
            # category set — the same one that read the run's and the anchor's.
            ref_rows = [(contrast = sv(ref.contrast[i]), unit = sv(ref.unit[i]),
                         habitat = sv(ref.habitat[i]),
                         ci_lo_std = ref.ci_lo_std[i], ci_hi_std = ref.ci_hi_std[i],
                         ci_lo = hascol(ref, :ci_lo) ? ref.ci_lo[i] : NaN,
                         ci_hi = hascol(ref, :ci_hi) ? ref.ci_hi[i] : NaN)
                        for i in 1:nrow(ref)]
            ref_sig = interval_sets(ref_rows, :std).excludes
            ref_rows_read = nrow(ref)
            rd = reference_row_diff(B.rows, ref)
            rec!("anchor_reference_columns_compared", joined(rd.columns))
            rec!("anchor_reference_rows_matched_in_anchor", rd.compared)
            rec!("anchor_reference_rows_agreeing", rd.agreeing)
            rec!("anchor_reference_rows_disagreeing", length(rd.differing))
            rec!("anchor_reference_rows_disagreeing_list", listed(rd.differing))
            rec!("anchor_reference_rows_absent_from_anchor", length(rd.absent))
            rec!("anchor_reference_rows_absent_from_anchor_list", listed(rd.absent))
            rec!("anchor_reference_max_abs_cell_difference", rd.worst)
        end
    end
    anchor_status = anchor_status_of(; frames_match = frames_match,
                                     n_undecidable_rows = anchor_std.n_undecidable_rows,
                                     reference_supplied = ref_supplied,
                                     reference_exists = ref_exists,
                                     reference_readable = ref_readable,
                                     anchor_categories = anchor_sig_seen,
                                     reference_categories = ref_sig)
    anchor_partition_available(anchor_status) && (anchor_sig = anchor_sig_seen)
    rec!("anchor_n1_taxa", B.n1)
    rec!("anchor_n2_taxa", B.n2)
    rec!("anchor_rows_total", length(B.rows))
    rec!("anchor_rows_interval_undecidable", anchor_std.n_undecidable_rows)
    rec!("anchor_frame_matches_run_frame", frames_match)
    rec!("anchor_pct_$(K_PHYLUM)_bg", pseudo(B.bg_phy, PHYLUM))
end
# The categories handed to the rule as undecidable in the primary run: see
# undecidable_categories above the guard for why a category with no row belongs there.
undecidable_std = undecidable_categories(primary_std, anchor_labels)
undecidable_crude = undecidable_categories(primary_crude, anchor_labels)
cats = category_sets(sig_std_primary, anchor_sig, GOVERNED;
                     primary_undecidable = undecidable_std)
cats_crude = category_sets(sig_crude_primary,
                           anchor_sig === nothing ? nothing : anchor_sig_crude, GOVERNED;
                           primary_undecidable = undecidable_crude)
rec!("anchor_config", ANCHOR_CONFIG)
rec!("anchor_chemistry_column", String(ANCHOR_CHEM_COL))
rec!("anchor_bg_exclusion_mode", "none")
rec!("anchor_status", anchor_status)
rec!("anchor_verified", anchor_status === :reproduces_reference)
rec!("anchor_reference_supplied", !isempty(ANCHORREF))
rec!("anchor_reference_rows", ref_rows_read)
rec!("anchor_reference_categories_excluding_zero", joined(ref_sig))
rec!("anchor_rows_out", basename(ANCHOR_OUT))
rec!("anchor_rows_written", !isempty(anchor_rows))
rec!("categories_governed", joined(GOVERNED))
rec!("categories_interval_excludes_zero_primary", joined(sig_std_primary))
rec!("categories_interval_excludes_zero_anchor", joined(anchor_sig_seen))
rec!("categories_newly_excluding_zero", joined(cats.newly))
rec!("categories_preexisting_excluding_zero", joined(cats.preexisting))
rec!("categories_no_longer_excluding_zero", joined(cats.no_longer))
rec!("categories_undecidable_in_primary", joined(cats.undecidable_in_primary))
rec!("categories_legend_only", joined(cats.legend_only))
rec!("categories_anchor_available", cats.anchor_available)
rec!("categories_interval_excludes_zero_primary_crude", joined(sig_crude_primary))
rec!("categories_interval_excludes_zero_anchor_crude", joined(anchor_sig_crude))
rec!("categories_newly_excluding_zero_crude", joined(cats_crude.newly))
rec!("categories_undecidable_in_primary_crude", joined(cats_crude.undecidable_in_primary))
# A category resting on a thin or untested row should not enter a list quietly, so the rows
# behind the partition carry their own disclosures: which rows were not tested, which have
# no interval, which rest on an empty background stratum cell, and how small the smallest
# cell behind each category in the set is.
rec!("categories_primary_rows", length(primary_rows(rows)))
rec!("categories_primary_rows_untested",
     joined(sort([r.habitat for r in primary_rows(rows) if !r.tested])))
rec!("categories_primary_rows_interval_undecidable", joined(primary_std.undecidable))
rec!("categories_primary_undecidable_supplied", joined(undecidable_std))
rec!("categories_primary_rows_with_a_zero_stratum_cell",
     joined(sort(unique([r.habitat for r in primary_rows(rows)
                         if r.n_zero_stratum_cells > 0]))))
rec!("categories_excluding_zero_min_bg_stratum_cell",
     joined([@sprintf("%s:%.4g", r.habitat, r.min_bg_stratum_cell)
             for r in primary_rows(rows) if r.habitat in sig_std_primary]))

# -- re-check list: the row-level properties a downstream text step reads. The plain form
#    describes every row the rows CSV carries, the disclosed arm's rows included; the
#    _registered form counts the registered contrasts alone, which is the form a re-check
#    of a registered rule reads, so a disclosed arm's thin row cannot move it.
reg_rows = registered_rows(rows)
rec!("rows_total", length(rows))
rec!("rows_continuity_corrected", count(r -> r.continuity_corrected, rows))
rec!("rows_continuity_corrected_registered", count(r -> r.continuity_corrected, reg_rows))
rec!("rows_rounded_for_fisher", count(r -> r.rounded_for_fisher, rows))
# Two different quantities, named so they cannot be confused: an empty WHOLE-ARM hit total,
# which is what the continuity correction reacts to, and an empty BACKGROUND STRATUM cell
# inside a standardised estimate, which the arm totals cannot reveal.
rec!("rows_with_a_zero_arm_total", count(r -> r.a == 0 || r.b == 0, rows))
rec!("rows_with_a_zero_arm_total_registered", count(r -> r.a == 0 || r.b == 0, reg_rows))
rec!("rows_with_a_zero_stratum_cell", count(r -> r.n_zero_stratum_cells > 0, rows))
rec!("rows_with_a_zero_stratum_cell_registered",
     count(r -> r.n_zero_stratum_cells > 0, reg_rows))
rec!("rows_min_bg_stratum_cell",
     (m = [r.min_bg_stratum_cell for r in rows if !isnan(r.min_bg_stratum_cell)];
      isempty(m) ? NaN : minimum(m)))
rec!("rows_tested", count(r -> r.tested, rows))
rec!("rows_registered_contrasts", length(reg_rows))
rec!("rows_disclosed_contrasts", count(r -> r.contrast == CA_VS_BG, rows))

# ---------------------------------------------------------------- report
@printf("\nCONFIG %s   background exclusion: %s (%s, %s)\n",
        CONFIG, exclusion_note(BG_EXCL), BG_EXCL_SOURCE, BG_EXCL_ROLE)
@printf("Ln (prokaryotic, habitat-annotated): %d sequences → %d taxa, %d genera\n",
        length(A.seq_idx), n1, length(gfg.genus))
@printf("Background frame: %d taxa; %d dropped as census members; %d removed by the scan; %d analysed\n",
        count(A.frame), look("census_overlap_in_frame"), look("bg_removed_by_scan"), n2)
@printf("All carriers: %d taxa;  non-Ln carriers: %d taxa;  Ca-class carriers: %d taxa\n",
        nrow(ac), nrow(nl), nrow(ca))
println("\nWhy candidate taxa were dropped (a partition, not a single bucket):")
@printf("  Ln          total %4d = unannotated %4d + isolation-only %3d + generic-only %3d + isolation-type %3d + other filtered %3d\n",
        dfg.tot, dfg.unann, dfg.iso, dfg.gen, dfg.isot, dfg.filt)
@printf("  background  total %4d = unannotated %4d + isolation-only %3d + generic-only %3d + isolation-type %3d + other filtered %3d\n",
        dbg.tot, dbg.unann, dbg.iso, dbg.gen, dbg.isot, dbg.filt)
@printf("\nPhylum composition: %s %.1f%% in the Ln arm vs %.1f%% in the background (%.1f%% before the scan)\n",
        PHYLUM, look("pct_$(K_PHYLUM)_fg"), look("pct_$(K_PHYLUM)_bg"),
        look("pct_$(K_PHYLUM)_bg_before_scan"))
@printf("\nHabitat-annotation coverage: Ln %.1f%% (%d/%d) vs background frame %.1f%% (%d/%d) (Fisher p = %.3g)%s\n",
        look("fg_coverage_pct"), cf.nok, cf.ncand,
        look("bg_frame_coverage_pct"), cb.nok, cb.ncand, covp,
        covp < 0.01 ? "\n  ⚠ coverage DIFFERS — annotated organisms are the well-studied ones;\n" *
                      "    the contrast is conditional on being annotated. Report it with every ratio." : "")

println("\n" * "="^92)
println("HABITAT TYPE — Ln carriers vs non-carrier prokaryote reference proteomes")
println("="^92)
show_panel(rows, LN_VS_BG, "taxa")
show_panel(rows, LN_VS_BG, "sequences")
show_panel(rows, LN_VS_BG, "genera")
show_panel(rows, ALL_VS_BG, "taxa")
show_panel(rows, NON_VS_BG, "taxa")
show_panel(rows, CA_VS_BG, "taxa")

println("\n" * "="^92)
println("HABITAT CHEMISTRY — non-exclusive flags, denominator = taxa")
println("="^92)
show_panel(rows, LN_VS_BG, "taxa (chemistry)")
show_panel(rows, ALL_VS_BG, "taxa (chemistry)")
show_panel(rows, NON_VS_BG, "taxa (chemistry)")
show_panel(rows, CA_VS_BG, "taxa (chemistry)")

outdf = DataFrame(rows)
outdf.a = compact.(outdf.a); outdf.b = compact.(outdf.b)
CSV.write(OUT, outdf)
# The anchor's own rows, so the partition can be re-decided at row level without a re-run.
# When there is no anchor the file is written EMPTY rather than left alone: a rows table from
# an earlier run, sitting beside a summary that says the anchor was unavailable, is exactly
# the kind of mismatched pair the scan records now guard against.
anchordf = isempty(anchor_rows) ? similar(outdf, 0) : DataFrame(anchor_rows)
if nrow(anchordf) > 0
    anchordf.a = compact.(anchordf.a); anchordf.b = compact.(anchordf.b)
end
CSV.write(ANCHOR_OUT, anchordf)
CSV.write(SUMOUT, DataFrame(key = first.(SUMMARY), value = last.(SUMMARY)))
@printf("\n✅ wrote %s (%d rows)\n", OUT, length(rows))
@printf("✅ wrote %s (%d anchor rows)\n", ANCHOR_OUT, nrow(anchordf))
@printf("✅ wrote %s (%d keys)\n", SUMOUT, length(SUMMARY))

println("""

READING IT. `crude` is unadjusted; `STD` reweights the background to the foreground phylum
mix. Every number below is read back from $(basename(SUMOUT)); none is typed into this file.
Every decision below is a CODE, not a phrase: the step that writes prose maps codes to words.

  Ln arm $(fmt("n1_taxa")) taxa ($(fmt("n_archaea_fg")) archaeal, $(fmt("n_bacteria_fg")) bacterial) against $(fmt("n2_taxa")) background taxa.
  Selection $(fmt("selection_file")) [$(SEL.digest[1:12])]: the rules speak for $(joined(GOVERNED));
    grid counts $(joined(SEL.grid_count)); standardisation phylum $(PHYLUM). None of those is a default.
  $(PHYLUM): $(fmt("pct_$(K_PHYLUM)_fg"))% of the Ln arm vs $(fmt("pct_$(K_PHYLUM)_bg"))% of the background, so where crude
  and STD disagree the crude number is phylogeny, not habitat. Report STD.
  Background carrier removal: code $(fmt("background_licence_code")), complete = $(fmt("background_removal_complete")).
    left in the analysed arm: HMM hit $(fmt("background_carriers_left_hmm_hit")) (definitional here: $(fmt("background_carriers_left_hmm_hit_definitional"))),
      family signature $(fmt("background_carriers_left_ipr017512_ge500")), either $(fmt("background_carriers_left_union"));
    scanned set == analysed frame: $(fmt("background_scan_frame_matches_analysis_frame")) (frame taxa absent from the scan $(fmt("background_analysis_frame_taxa_missing_from_scan")));
    prefilter recall re-derived on the lower bound: $(fmt("background_recall_pass_derived")) (bound $(fmt("background_prefilter_recall_lower_bound_derived")), rate $(fmt("background_prefilter_recall_rate_derived")); the record says $(fmt("background_prefilter_recall_pass")), rate $(fmt("background_prefilter_recall_rate")), bar $(fmt("background_prefilter_recall_bar")));
    scan records vouch for the table read: $(fmt("background_scan_evidence_ok")) (table state $(fmt("background_scan_table_state")),
      acceptance $(fmt("background_scan_acceptance_pass")) PASS / $(fmt("background_scan_acceptance_miss")) MISS / $(fmt("background_scan_acceptance_disclosed")) DISCLOSED);
    evidence reasons: $(blank("background_scan_evidence_reasons"));
    licence reasons: $(blank("background_removal_incomplete_reasons")).
  $(SHARE_CAT) sharing, registered comparator / Ln = $(fmt("sharing_ratio_nonln_over_ln")) [$(fmt("sharing_ratio_ci_lo")), $(fmt("sharing_ratio_ci_hi"))]
    → code $(fmt("sharing_code")) by $(fmt("sharing_rule_used")); crude reading $(fmt("sharing_code_crude")) by $(fmt("sharing_rule_used_crude")).
    comparator arm: $(fmt("comparator_arm_taxa")) taxa, $(fmt("comparator_arm_taxa_with_calcium_class")) with a calcium-class census sequence,
    $(fmt("comparator_arm_taxa_gap_or_unknown_class_only")) with none but gap/unknown-class sequences.
  DISCLOSED, beside it: Ca-class arm $(fmt("caonly_arm_taxa")) taxa, ratio $(fmt("sharing_ratio_caonly_over_ln")) [$(fmt("sharing_ratio_caonly_ci_lo")), $(fmt("sharing_ratio_caonly_ci_hi"))]
    → code $(fmt("sharing_code_caonly")) by $(fmt("sharing_rule_used_caonly")); it does not replace the registered comparator.
  $(TIER_CAT) tier code: $(fmt("$(K_TIER)_tier_code")) — $(fmt("$(K_TIER)_tier_reason")).
  $(GENUS_CAT): category $(fmt("$(K_GENUS)_category_code")), genus clause $(fmt("$(K_GENUS)_genus_clause_code")).
  Categories excluding zero — primary: $(blank("categories_interval_excludes_zero_primary"))
    anchor ($(fmt("anchor_config")), $(fmt("anchor_status")), verified $(fmt("anchor_verified"))): $(blank("categories_interval_excludes_zero_anchor"))
    newly $(blank("categories_newly_excluding_zero")) | pre-existing $(blank("categories_preexisting_excluding_zero")) | no longer $(blank("categories_no_longer_excluding_zero")) | undecidable now $(blank("categories_undecidable_in_primary")) | legend-only $(blank("categories_legend_only"))
    smallest background stratum cell behind each: $(blank("categories_excluding_zero_min_bg_stratum_cell"))
  Rows with an empty background stratum cell: $(fmt("rows_with_a_zero_stratum_cell")) of $(fmt("rows_total")) ($(fmt("rows_with_a_zero_stratum_cell_registered")) of $(fmt("rows_registered_contrasts")) registered); smallest cell $(fmt("rows_min_bg_stratum_cell")).
  Coverage differs between the arms: $(fmt("coverage_differs_at_1pct")) — every ratio is
  conditional on being annotated.""")
BG_EXCL == MOVEIN_MODE &&
    println("""
  The move-in level — a sensitivity level, never the headline:
    pool rule $(fmt("movein_pool_rule")) ($(fmt("movein_pool_rule_role")), $(fmt("movein_pool_rule_source"))), from $(fmt("movein_pool_source"))
    pool $(fmt("movein_pool_taxa")) taxa, $(fmt("movein_pool_taxa_in_frame")) of them in the frame; moved $(fmt("movein_moved_taxa")), of which $(fmt("movein_moved_taxa_annotated")) carry a usable call
    n1 $(fmt("movein_n1_before_the_move")) → $(fmt("n1_taxa")); the move took $(fmt("movein_background_taxa_the_move_removed")) taxa out of the analysed background that the exclusion had left in
    the moved taxa carry $(fmt("movein_moved_census_sequences")) census sequences, so they add $(fmt("movein_sequence_units_from_moved_taxa")) units to the per-sequence series""")
PLAN_BG_SIZE_BELOW === nothing && PLAN_BG_PHYLUM_PCT_MAX === nothing ||
    println("""
  Recorded criteria for this step, disclosed and gating nothing:
    analysed background size $(fmt("criterion_bg_size_value")) against bound $(fmt("criterion_bg_size_bound")): met = $(fmt("criterion_bg_size_met"))
    background $(PHYLUM) % $(fmt("criterion_bg_$(K_PHYLUM)_pct_value")) against bound $(fmt("criterion_bg_$(K_PHYLUM)_pct_bound")): met = $(fmt("criterion_bg_$(K_PHYLUM)_pct_met"))""")

end # driver
