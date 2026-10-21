#=
    habitat_reaggregate.jl — one habitat call per taxon per rule variant, computed offline
    from a replayable relation cache.

    THE FAILURE THIS PREVENTS. A habitat number taken from a live re-fetch is not
    reproducible. The call for a taxon then depends on what OmniCrobe answers that
    afternoon, on whether the fetch succeeds at all (a failed fetch enters the analysis as
    the ecological statement "this organism has no known habitat"), and — where two
    categories tie — on Julia's Dict iteration order. Two runs of "the same" pipeline are
    not comparable that way, so a change in the output cannot be attributed to the rule
    change rather than to the refetch.

    The remedy is to freeze the INPUTS, not the outputs: the raw response bodies live on disk
    under habitat/omnicrobe_cache/relations in the data root, the ontology climb lives in
    habitat/obt_path_cache.csv beside it, and every rule variant is computed from those same
    bytes in one pass. A run that holds those bytes reproduces the table exactly with the
    network unplugged.

    WHERE THAT CACHE COMES FROM. A run writes it; the data set does not carry it. The response
    bodies arrive from the service and the ontology climb from the ontology endpoint, and both
    stay on disk for every later run. What the data set carries is the derived call table this
    step writes, so a reader holding only the data set enters the habitat route at that table.
    To run this step itself, seed the cache from directories of already-fetched <taxid>.json
    bodies, given as bare arguments or in $LNFINDER_RELATION_SEED. A run with neither cache nor
    seed fetches, and the service is re-mined and unversioned, so it answers with a different
    data set rather than with these bytes: such a run measures the habitat layer again rather
    than replaying it.

    WHAT IT WRITES. habitat/habitat_taxon_calls.csv, under the data root: one row per (config, taxid)
    carrying the winner, the raw distinct-term and mention counts per type, the 1/k weights
    a tie was split into, the non-exclusive chemistry flags and the isolation/generic
    bookkeeping. Everything a downstream sensitivity rule needs is in that row, so no
    variant of the rule costs a refetch.

    CONFIGS. `primary` is the rule this release computes with. `legacy` reproduces the rule
    the committed CSVs were computed under and exists only so those CSVs can be re-derived
    and diffed field by field — never for new analysis. The other four move exactly one axis
    off `primary`.

    REPLICATION GATES. The script ends by re-deriving the committed calls under `legacy` and
    comparing them taxon by taxon, and by checking coverage and isolation-term recovery
    under `primary`. Every gate is an assertion, not a print: this table feeds every
    downstream estimate, so it must not be allowed to drift quietly. A gate that fails
    prints the measured value next to the expected one and then aborts the script. The
    four comparisons against the committed calls are made only when the targets file
    declares, by its sha256, the census they were written for and the census this run
    reads has that digest. Otherwise - `--no-gates`, a targets file that declares no
    census, or a census with another digest - each of the four prints its measurement as
    NOT APPLICABLE with the reason, and the closing line says so: those calls belong to
    one census, and neither a default nor a path says which census a file holds.

    NOTHING IS PUBLISHED UNTIL THE GATES HAVE PASSED. The gates can only be computed from a
    table that has been built, and building that table straight over the committed one and
    gating it afterwards would mean a run whose gates FAILED had already replaced the input
    of every estimate downstream of it, with a failure message about a number while the
    damage was to a file. The outputs are therefore written
    into a `_staging` directory beside their destinations and moved into place by `publish!`,
    and only if every gate passed. A failed run leaves what it built in `_staging`, so the
    failure is inspectable, and leaves the destination holding the last good run.

    INCLUDING THIS FILE RUNS NOTHING. Everything above the two `PROGRAM_FILE` guards is
    definitions, so the self-test can call the option reader, the targets reader and the
    publish gate directly. The command line is read in the FIRST guard, before the packages
    are loaded, because a bad command line should cost nothing and must never fall through
    to a default.

    NO MEASURED NUMBER IN THIS FILE. A gate's expected value is one of two things: a
    requirement of the METHOD — zero, where zero means "nothing may be missing", or a count
    this run computes for itself — or a quantity MEASURED on this dataset. This file carries
    only the first kind. The second kind lives in the targets file `--gate-targets` names,
    one row per gate, each with the note that says where its number came from. That split is
    not a formality: a measured number typed in here is a result of this dataset published in
    code, and it also goes stale the moment an input changes and then quietly contradicts the
    table beside it. Every other script in this pipeline computes or reads its numbers at run
    time, and so does this one.

    A run given no targets file is REFUSED. The gates are the only thing standing between a
    silently drifted call table and every estimate downstream of it, so a gate set that is
    skipped by default is a gate set nobody runs: skipping it has to be typed as
    `--no-gates`, and is then said out loud in the report instead of looking like a pass.

    SEEDING. The relation cache can be primed from directories of already-fetched
    <taxid>.json response bodies, given as bare arguments or in LNFINDER_RELATION_SEED (a
    semicolon-separated path list). Whatever is still missing is fetched over the network,
    one taxon at a time with a short pause after each fetch that succeeds; a transient
    failure (status 429, 500, 502, 503 or 504, or a dropped connection) is asked again up
    to three times with a growing wait before the taxon is recorded as a failed fetch. No
    address and no account is sent: the relation requests name this client in their
    User-Agent header and say nothing about who runs it.

    Run:  julia --project=. scripts/habitat_reaggregate.jl \
              --gate-targets <path> [seed-dir ...]
=#

const REPO = normpath(joinpath(@__DIR__, ".."))

# Every path below is a NAME in the one registry, resolved under the data root. No
# location is spelled out here: src/ModuleDataInputs.jl declares the layout, docs/DATA.md
# documents it, and a rename has to pass through those two and nothing else.
include(joinpath(REPO, "src", "ModuleDataInputs.jl"))

const CENSUS     = DataInputs.input_path("census_habitat")
const POOL       = DataInputs.input_path("pool_taxa")
const BACKGROUND = DataInputs.input_path("pool_background")
const CACHE_DIR  = DataInputs.input_path("relation_cache")
const REL_DIR    = joinpath(CACHE_DIR, "relations")
const MANIFEST   = DataInputs.input_path("relation_manifest")
const OBT_CSV    = DataInputs.input_path("obt_path_cache")
const OUT        = DataInputs.input_path("habitat_calls")

# Stamped into the manifest next to every cached body. The call depends on both releases,
# and a cache with no release stamp cannot be told apart from one taken years later.
const OMNICROBE_RELEASE = "OmniCrobe 1.7"
const OBT_RELEASE       = "OBT 2022-04-28"

# The cap on distinct terms per taxon, held at the library default so the `legacy` config
# reproduces which terms that cap drops.
const MAX_DISTINCT = 50

# ---- the command line, read before anything is loaded, fetched or written ---------
# A bad command line costs nothing here and never falls through to a default: a run under
# settings nobody asked for produces a table that looks exactly like the answer to the
# command that was typed. What follows is the reader; the program guard below it is where a
# run actually reads ARGS.

const USAGE = """
Usage: julia --project=. scripts/habitat_reaggregate.jl [options] [seed-dir ...]

A bare argument is a seed directory of already-fetched <taxid>.json response bodies; so is
every entry of LNFINDER_RELATION_SEED (a semicolon-separated path list).

Options:
  --gate-targets <path>  the expected value, tolerance and comparison of every replication
                         gate whose expected value is a MEASUREMENT of this dataset, as
                         kind,gate,expected,tol,note with kind in
                         equal|at_least|expect|census. A `census` row carries, in its
                         expected column, the sha256 of the census the four comparisons
                         against the committed calls were written for; they are made only
                         when the census this run reads has that digest
  --no-gates             run without those gates at all, and say so in the report. This has
                         to be typed: the gates are the only thing between a silently
                         drifted call table and every estimate downstream of it, so a gate
                         set skipped by default is a gate set nobody runs. With no targets
                         file nothing declares a census, so the four comparisons against
                         the committed calls are not made either
  --out <path>           write the call table here instead of over the committed one. A
                         verification re-run needs this, because the committed table is the
                         input of every downstream step: a re-run that overwrites it has
                         destroyed the thing it was re-running to check
  --census <path>        read the census from here instead of the deposited one (the
                         census_habitat input under the data root). The census fixes the
                         taxid universe and every per-taxon label, so a run against a
                         re-annotated census is a different measurement: the path and the
                         sha256 are recorded in the report, and the four gates that compare
                         this run against the committed calls are made only if the targets
                         file declares that digest
  -h, --help             print this text and exit

Nothing is published unless every gate passed. The outputs are built in a `_staging`
directory beside each destination and moved into place at the very end, so a run that fails a
gate leaves what it built there and leaves the destination holding the last good run.

Exit status: 0 on a completed run with every gate met; 1 on a command-line error, a missing
or malformed input, or a failed replication gate.
"""

const OPTION_SPEC = Dict("--gate-targets" => :value, "--no-gates" => :flag,
                         "--out" => :value, "--census" => :value,
                         "--help" => :flag, "-h" => :flag)

"""
    parse_options(args, spec) -> (options, positional)

Split a command line into option => value (a flag mapping to "") and the bare words left
over, which are this script's seed directories. A token beginning with `-` is an option and
must be one this script defines; anything else is a seed directory. An unknown or misspelt
option, an option given twice, or a value-taking option with no value is an error, because
an ignored token is a silently different run.
"""
function parse_options(args, spec::AbstractDict)
    out = Dict{String,String}(); positional = String[]
    i = 1
    while i <= length(args)
        a = String(args[i])
        if !startswith(a, "-")
            push!(positional, a); i += 1; continue
        end
        kind = get(spec, a, nothing)
        kind === nothing &&
            error("unknown option '$a' (defined: " * join(sort(collect(keys(spec))), ", ") * ")")
        haskey(out, a) && error("option '$a' given twice")
        if kind === :flag
            out[a] = ""
            i += 1
        else
            # A value that is itself a defined option means the value was forgotten:
            # `--gate-targets --no-gates` must not read a targets file called "--no-gates".
            (i == length(args) || haskey(spec, String(args[i + 1]))) &&
                error("option '$a' needs a value")
            out[a] = String(args[i + 1])
            i += 2
        end
    end
    return (options = out, positional = positional)
end

"""
    cli_or_exit(f) -> the value of `f()`

Run a command-line step; an `ErrorException` from it ends the program with the message and
the usage on stderr and exit status 1. The point is that a bad command line stops here,
before anything is loaded, fetched or written.
"""
function cli_or_exit(f)
    try
        return f()
    catch e
        e isa ErrorException || rethrow()
        println(stderr, "habitat_reaggregate.jl: ", e.msg, "\n")
        print(stderr, USAGE)
        exit(1)
    end
end

"""
    option_string(cli, flag, default) -> String

The value of a parsed option, or `default` when it was not given. An empty value is an error
naming the option rather than a fall-back to the default, which would run under settings
other than the ones that were asked for.
"""
function option_string(cli::AbstractDict, flag::AbstractString, default)
    haskey(cli, flag) || return default
    v = strip(cli[flag])
    isempty(v) && error("$flag needs a value")
    return String(v)
end

"""
    gate_options(options, positional, default_out) -> NamedTuple

The settings a run needs, from an already-parsed command line: which targets file the gates
come from, whether they were declined, where the call table goes, which census is read (the
deposited `census_habitat` input unless `--census` names another) and whether that was typed,
and the seed directories. Separate from the parser so the suite can ask what a given command
line means without running anything.

`census_given` decides only how a missing census is reported: through the registry, which
names the step that writes it, or by the path that was typed. It does not decide which gates
are made; that is settled by the census's digest (`reproduction_scope`).
"""
function gate_options(options::AbstractDict, positional, default_out::AbstractString)
    targets  = option_string(options, "--gate-targets", "")
    no_gates = haskey(options, "--no-gates")
    # A GATE SET THAT IS OPT-IN IS A GATE SET NOBODY RUNS. Given no targets file this script
    # would write the whole call table and exit 0, and the table would be byte-for-byte what
    # a fully gated run writes; the only difference would be a few lines of a report nobody
    # diffs. So the default is refused, and declining the gates is a decision that has to be
    # typed and is then recorded as one.
    isempty(targets) && !no_gates &&
        error("no --gate-targets were given, so the gates whose expected value is a " *
              "measurement of this dataset would not be checked at all, and this table is " *
              "the input of every estimate downstream of it. Pass --gate-targets <path>, " *
              "or --no-gates to run without them deliberately")
    # Both together is a contradiction, and resolving it either way would run something other
    # than what was typed.
    !isempty(targets) && no_gates &&
        error("--gate-targets and --no-gates were both given; one names the gates to check " *
              "and the other declines to check any")
    return (gate_targets = targets, no_gates = no_gates,
            out = option_string(options, "--out", default_out),
            census = option_string(options, "--census", CENSUS),
            census_given = haskey(options, "--census"),
            seed_dirs = collect(String, positional))
end

# THE FIRST PROGRAM GUARD: the command line, and nothing else. A run's command line is read
# here, before the packages below are loaded, because a bad one costs nothing at this point
# and a run under settings nobody asked for produces a table that looks exactly like the
# answer to the command that was typed. Including this file reads no command line at all.
const IS_PROGRAM = abspath(PROGRAM_FILE) == @__FILE__

if IS_PROGRAM
    const CLI = cli_or_exit(() -> parse_options(ARGS, OPTION_SPEC))
    if haskey(CLI.options, "--help") || haskey(CLI.options, "-h")
        print(USAGE)
        exit(0)
    end
    const OPT = cli_or_exit(() -> gate_options(CLI.options, CLI.positional, OUT))
end

using Pkg; Pkg.activate(REPO)
using DataFrames, CSV, Printf, Dates, ProgressMeter, SHA

include(joinpath(REPO, "src", "ModuleEnvironmentTools.jl")); using .EnvironmentTools

# One axis moved at a time off `primary`, so a downstream sensitivity row can attribute a
# change to the axis that caused it. `legacy` moves all four at once on purpose: it is the
# rule the committed CSVs carry, and so the only config they can be compared against.
const CONFIGS = (
    (name = "primary",      kw = (isolation = :term,  drop_generic = true,
                                  chem_text = :forms_union, tie_rule = :mentions_fractional)),
    (name = "legacy",       kw = (isolation = :taxon, drop_generic = false,
                                  chem_text = :first_form,  tie_rule = :legacy_dict)),
    (name = "generic_none", kw = (isolation = :term,  drop_generic = false,
                                  chem_text = :forms_union, tie_rule = :mentions_fractional)),
    (name = "iso_taxon",    kw = (isolation = :taxon, drop_generic = true,
                                  chem_text = :forms_union, tie_rule = :mentions_fractional)),
    (name = "text_first",   kw = (isolation = :term,  drop_generic = true,
                                  chem_text = :first_form,  tie_rule = :mentions_fractional)),
    (name = "text_object",  kw = (isolation = :term,  drop_generic = true,
                                  chem_text = :object_name, tie_rule = :mentions_fractional)),
)

# `string`, not `String`: the CSV reader types a column by what it holds, so a re-read of the
# manifest hands back a `Date` where the first run wrote text, and a String-only conversion
# turns that into a crash on the second run only.
sv(x) = ismissing(x) ? "" : strip(string(x))
fnum(x) = x === missing || x === nothing ? NaN :
          x isa Real ? float(x) : something(tryparse(Float64, sv(x)), NaN)
iso_date(t::Real) = Dates.format(Dates.unix2datetime(t), "yyyy-mm-dd")

# ---- the gates whose expected value is a measurement of this dataset --------------
# Read here, before a single body is seeded, fetched or aggregated, so a malformed or
# thinned targets file stops the run in a second rather than after the whole pass.

const GATE_COLUMNS = (:kind, :gate, :expected, :tol, :note)
# Two kinds of row compare a measurement with an expected value and a third says how many of
# those the file holds. `equal` is |measured - expected| <= tol. `at_least` is
# measured >= expected - tol, which is how an ALLOWANCE for a known, named set of rows that
# do not reproduce is written down — as a stated number of rows, rather than folded invisibly
# into the expected value or hidden by widening a two-sided band. A fourth, `census`, names by
# sha256 the census the REPRODUCTION_GATES were written for (see `reproduction_scope`).
const GATE_KINDS = ("equal", "at_least", "expect", "census")
# What an `expect` row may count. Only the gate rows: the point of it is that the file cannot
# lose rows without saying so.
const EXPECTED_COUNTS = ("gates",)

# The gates whose expected value is a MEASUREMENT of this dataset rather than a requirement
# of the method. A measurement of the data is never a constant of this code, so each of these
# takes its expected value, tolerance and comparison from the targets file, by name. The
# NAMES are here because a rule's name is not a quantity measured on anything, and because a
# file naming a gate this run does not have — or a gate this run has and the file does not —
# is then a refusal in both directions rather than a check that quietly did not happen.
#
# `from_run` marks a gate measured against a number THIS RUN computes: the size of an arm it
# has just built. Such a gate takes no expected value from the file, because a second
# statement of a number the run computes is a number that can disagree with it silently. What
# the file supplies for it is the tolerance — the allowance for rows a later, deliberate
# change of the analysis is known to have moved.
#
# `kind` is the SHAPE of the comparison, and it is declared here as well as in the file,
# which must agree with it. The shape is not a quantity measured on anything, and leaving it
# to the file alone would make it an input: one word changed from `equal` to `at_least` turns
# a gate that catches drift in both directions into one that catches it in a single
# direction, without touching a number and without failing anything.
const MEASURED_GATES = (
    (gate = "Ln chemistry reproduced",                from_run = true,  kind = "at_least"),
    (gate = "Ln taxa with a usable type term",        from_run = false, kind = "equal"),
    (gate = "frame taxa with a usable type term",     from_run = false, kind = "equal"),
    (gate = "Ln, at least one non-isolation term",    from_run = false, kind = "equal"),
    (gate = "frame, at least one non-isolation term", from_run = false, kind = "equal"),
)

# THE FOUR GATES THAT COMPARE THIS RUN AGAINST THE COMMITTED CALLS: the `legacy` config's
# habitat type and chemistry against the calls the census and the background carry in their
# own columns, taxon by taxon. Three take their expected value from the arm this run builds;
# the fourth is also a measured gate above, because its allowance is a number measured on the
# data. They are named here because whether they are made at all is decided per run
# (`reproduction_scope`), and a gate that is declined still has to be reached: the accounting
# guard at the end holds every name here, made or declined, so a deleted call site cannot hide
# behind a decline.
const REPRODUCTION_GATES = ("Ln habitat_type reproduced", "background habitat_type reproduced",
                            "Ln chemistry reproduced", "background chemistry reproduced")
# Every gate the accounting guard holds a run to: the measured gates and the reproduction
# gates, each name once.
const ACCOUNTED_GATES = unique(vcat(String[g.gate for g in MEASURED_GATES],
                                    collect(String, REPRODUCTION_GATES)))

"""
    read_gate_targets(path, declared) -> (gates, census)

`gates`: every gate the targets file declares, by name, with its comparison, its expected
value where the file carries one, its tolerance, and the note recording where the number came
from. `census`: the census the file says the REPRODUCTION_GATES were written for, as
`(sha256, note)`, or `nothing` when it declares none — which is not a refusal, because a file
can carry the measured gates of a census the committed calls do not belong to; it is what
`reproduction_scope` declines those four gates on.

A census is declared by content, not by name: one `census` row, no gate named, the census's
sha256 in the expected column, no tolerance (a digest matches or it does not), and a note
saying which census that is.

WHY THE FILE ALSO STATES HOW MANY GATES IT HOLDS, AND WHICH. Neither is a quantity measured
on anything, and both nevertheless have to reach the run. The count, because a file that lost
rows is a run that checked fewer gates and printed "All replication gates passed" in the same
words as a complete one — how much of this table is checked is not something a run may
discover from the file it was handed. The set, because a file written for an older version of
this script would otherwise declare gates nobody evaluates and stay silent about the ones
nobody declared. Both are compared with `declared` in both directions, so a thinned,
half-merged or stale file is a failure rather than a quieter success.
"""
function read_gate_targets(path::AbstractString, declared)
    isfile(path) || error("--gate-targets $path does not exist")
    # The expected column is read as text: it holds counts and, on a census row, a digest, and
    # a digest that happens to be all decimal digits must not be read back as a number. A file
    # with no such column is refused by name below rather than by the CSV reader.
    tbl = CSV.read(path, DataFrame; types = Dict(:expected => String), validate = false)
    for c in GATE_COLUMNS
        c in propertynames(tbl) ||
            error("--gate-targets $path needs the columns " *
                  join(String.(GATE_COLUMNS), ", ") * " and has no '$c'")
    end
    want = Dict{String,NamedTuple}(g.gate => (from_run = g.from_run, kind = g.kind)
                                   for g in declared)
    gates = Dict{String,NamedTuple}(); declared_n = nothing; census = nothing
    for i in 1:nrow(tbl)
        kind = sv(tbl.kind[i]); name = sv(tbl.gate[i]); note = sv(tbl.note[i])
        if kind == "equal" || kind == "at_least"
            isempty(name) &&
                error("a '$kind' row names the gate its numbers belong to, and row $i of " *
                      "$path names none")
            haskey(gates, name) &&
                error("'$name' is declared twice in $path, so one gate's expected value is " *
                      "stated in two places and only one of them decides the run")
            haskey(want, name) ||
                error("'$name' is not a gate of this run, which measures " *
                      join(sort(collect(keys(want))), "; ") * ". A file written for another " *
                      "version of this script declares numbers nobody checks")
            kind == want[name].kind ||
                error("the gate '$name' is declared here as a '$(want[name].kind)' " *
                      "comparison and in $path as a '$kind' one. The shape of the " *
                      "comparison is not a number measured on anything, so it is stated in " *
                      "both places and has to agree: one word changed in the file alone " *
                      "would narrow what the gate can catch without moving a number")
            # A number with no record of where it came from cannot be told apart from one
            # transcribed out of this script's own output, which would make the gate a
            # statement that the run reproduces itself.
            isempty(note) &&
                error("the gate '$name' carries no note, so nothing in $path says where its " *
                      "expected value came from")
            tol = fnum(tbl.tol[i])
            (isnan(tol) || tol < 0) &&
                error("the gate '$name' has no tolerance. An absent tolerance is not a " *
                      "tolerance of zero; it is a gate whose width nobody decided")
            raw = sv(tbl.expected[i])
            if want[name].from_run
                isempty(raw) ||
                    error("the gate '$name' is measured against a number THIS RUN computes, " *
                          "and $path also carries an expected value for it ('$raw'): two " *
                          "statements of one number can disagree, and this one would never " *
                          "be checked against the run's")
                gates[name] = (kind = kind, expected = nothing, tol = tol, note = note)
            else
                v = fnum(tbl.expected[i])
                (isnan(v) || v < 0 || v != round(v)) &&
                    error("the gate '$name' is a count of taxa, and its expected value in " *
                          "$path is '$raw'")
                gates[name] = (kind = kind, expected = Int(v), tol = tol, note = note)
            end
        elseif kind == "expect"
            name in EXPECTED_COUNTS ||
                error("an 'expect' row counts one of " * join(EXPECTED_COUNTS, ", ") *
                      ", not '$name'")
            declared_n === nothing ||
                error("'$name' is expected twice, so $path states its own size in two places")
            v = fnum(tbl.expected[i])
            (isnan(v) || v < 0 || v != round(v)) &&
                error("the expected number of $name is a whole count and this one is " *
                      "'$(sv(tbl.expected[i]))'")
            declared_n = Int(v)
        elseif kind == "census"
            census === nothing ||
                error("$path declares its census twice, so which census the reproduction " *
                      "gates were written for is stated in two places")
            isempty(name) ||
                error("a 'census' row names no gate: it states, by the sha256 in its expected " *
                      "column, which census the reproduction gates were written for, and row " *
                      "$i of $path names '$name'")
            raw = sv(tbl.expected[i])
            digest = lowercase(raw)
            (length(digest) == 64 && all(c -> c in "0123456789abcdef", digest)) ||
                error("the census row of $path carries '$raw' where the sha256 of a census " *
                      "belongs (64 hexadecimal digits)")
            isempty(sv(tbl.tol[i])) ||
                error("the census row of $path carries a tolerance. A census is matched by " *
                      "its digest exactly, and a width on that match is one nobody can decide")
            # Which census a digest belongs to cannot be read off the digest, so a declaration
            # with no note is a number nobody can check against anything.
            isempty(note) &&
                error("the census row of $path carries no note, so nothing in it says which " *
                      "census that digest belongs to or where it was taken from")
            census = (sha256 = digest, note = note)
        else
            error("a row of $path is one of " * join(GATE_KINDS, ", ") * ", not a '$kind'")
        end
    end
    declared_n === nothing &&
        error("$path carries no 'expect' row for 'gates', so a file that lost rows would be " *
              "a run that checked fewer gates and reported success in the same words as a " *
              "complete one")
    declared_n == length(gates) ||
        error("$path declares $declared_n gate(s) and carries $(length(gates)): the file has " *
              "been thinned, extended or merged since the count beside it was written, and " *
              "how much of this table is checked is not something a run may discover from " *
              "the file it was handed")
    absent = sort([g for g in keys(want) if !haskey(gates, g)])
    isempty(absent) ||
        error("$path declares no expected value for " * join(absent, "; ") * ", so " *
              (length(absent) == 1 ? "that gate" : "those gates") * " would be measured " *
              "against nothing. Every gate whose expected value is a measurement of this " *
              "dataset takes it from this file, because the code carries none")
    return (gates = gates, census = census)
end

"""
    census_digest(bytes) -> String

The sha256 of a census, lowercase hex, taken from the very bytes the run parses: a digest
from a second read can describe a file that changed in between, and then the declaration it is
compared with vouches for a census this run did not read.
"""
census_digest(bytes::AbstractVector{UInt8}) = bytes2hex(sha256(bytes))

"""
    reproduction_scope(targets, census_sha256) -> (evaluated, reason, detail)

Whether this run makes the four REPRODUCTION_GATES, and if not, why.

REPRODUCING THE COMMITTED CALLS IS A QUESTION ABOUT ONE CENSUS: the one those gates were
written for. Asked of another it is not a weaker test but a different one, which an honest
run can fail and whose pass would read as evidence that two censuses agree. Neither a default
nor a path says which census a file holds, so the answer is taken from content: the gates are
made only when a targets file was read (`targets`, `nothing` under --no-gates), it declares a
census, and `census_sha256` — the digest of the census this run reads — is that census's.
`reason` is the short form printed beside each declined gate; `detail` says the same in full,
once.
"""
function reproduction_scope(targets, census_sha256::AbstractString)
    targets === nothing &&
        return (evaluated = false, reason = "--no-gates",
                detail = "--no-gates was given, so no targets file says which census these " *
                         "comparisons were written for")
    targets.census === nothing &&
        return (evaluated = false, reason = "no census declared",
                detail = "the targets file declares no census, so nothing says which census " *
                         "these comparisons were written for")
    got = lowercase(census_sha256)
    targets.census.sha256 == got ||
        return (evaluated = false, reason = "another census",
                detail = "the targets file declares them for the census with sha256 " *
                         targets.census.sha256 * ", and this run reads one with sha256 " * got)
    return (evaluated = true, reason = "",
            detail = "the targets file declares them for the census with sha256 " * got *
                     ", which is the census this run reads")
end

"""
    decline_gate!(reached, label, measured, reason) -> true

Print a gate this run does not make as NOT APPLICABLE, with its measurement and the reason,
and register it in `reached`.

A GATE DECLINED IS STILL A GATE REACHED. The declared-against-evaluated cross-check exists so
that a deleted call site cannot hide behind --no-gates; a gate declined for a stated reason
must therefore still register, or that check fires on the decline instead of on the deletion
it was written to catch. The measurement is printed either way, so declining reports a number
and withholds only the verdict.
"""
function decline_gate!(reached::Vector{String}, label::AbstractString, measured::Integer,
                       reason::AbstractString)
    reach_gate!(reached, label)
    @printf("  %-46s %6d   NOT APPLICABLE (%s)\n", label, measured, reason)
    return true
end

"""
    reach_gate!(reached, label) -> reached

Register `label` as reached, once. A gate evaluated twice under one name has two measurements
and only one of them is reported.
"""
function reach_gate!(reached::Vector{String}, label::AbstractString)
    label in reached &&
        error("the gate '$label' was evaluated twice under one name, so one of its two " *
              "measurements is reported and the other is not")
    return push!(reached, String(label))
end

# ---- nothing is published until the gates have passed -----------------------------

# ONE staging directory per destination directory, under a fixed name beside the file it
# will become. Beside, because a move within one filesystem is atomic while a move across
# two is a copy that can half-finish; under a fixed name, because a fresh name per run
# accumulates one directory per run on any filesystem that will not let go of an empty one.
const STAGING_DIR = "_staging"

"""
    open_staging(destinations) -> Vector{NamedTuple}

One `(staged, dest)` pair per destination, with every staging directory created and emptied
of whatever an earlier attempt left in it, so nothing from a previous run can be published
by this one.

The staged path is never the destination: that is the whole property, and it is what lets a
run build its outputs in full and still decide afterwards whether they may be seen.
"""
function open_staging(destinations)
    dests = String[abspath(String(d)) for d in destinations]
    length(unique(dests)) == length(dests) ||
        error("two outputs of this run would be published to one path: " *
              join(sort(dests), ", "))
    for dir in unique(String[joinpath(dirname(d), STAGING_DIR) for d in dests])
        mkpath(dir)
        for f in readdir(dir)
            rm(joinpath(dir, f); force = true, recursive = true)
        end
    end
    return NamedTuple[(staged = joinpath(dirname(d), STAGING_DIR, basename(d)), dest = d)
                      for d in dests]
end

"""
    publish!(plan, failures) -> Vector{String}

Move every staged file onto its destination, but ONLY if `failures` is empty.

NOTHING IS PUBLISHED UNTIL THE GATES HAVE PASSED, and the two halves of that sentence are
one function so that they cannot come apart. Written straight to its destination, a run
whose gates failed leaves a complete-looking call table there — same name, same shape, every
column present — distinguishable from a good run only by an exit status nobody reads a month
later; and this table is not an output anybody reads, it is the INPUT of every estimate
downstream of it. On a failure the staged files stay where they are, named in the error, so
the failure is inspectable; what they are not is published.

`failures` is the list of gate failure messages, all of them, because a run that failed two
gates should say so once.
"""
function publish!(plan, failures)
    if !isempty(failures)
        staged = join(String[p.staged for p in plan], "\n  ")
        dests  = join(String[p.dest for p in plan], "\n  ")
        error(join(failures, "\n") *
              "\nNOTHING WAS PUBLISHED. What this run built is here, so the failure is " *
              "inspectable:\n  " * staged *
              "\nand these still hold the last good run:\n  " * dests)
    end
    written = String[]
    for p in plan
        isfile(p.staged) ||
            error("the gates passed and $(p.staged) is not there to publish, so this run " *
                  "would leave $(p.dest) holding another run's numbers while reporting success")
        mv(p.staged, p.dest; force = true)
        push!(written, p.dest)
    end
    # Everything is published by this point and the staging directories are now empty. A
    # filesystem that will not remove one — a sync client still holding a handle, say —
    # leaves an empty directory under a name that is reused and cleared next run, so nothing
    # accumulates and nothing is lost. Refusing to return over that would report a complete,
    # gated, published run as a failure.
    for dir in unique(String[dirname(p.staged) for p in plan])
        try; rm(dir; recursive = true); catch; end
    end
    return written
end

"""
    accounting_failures(declared, reached) -> Vector{String}

EVERY DECLARED GATE WAS ACTUALLY EVALUATED, AND NO OTHER. The targets file and this script
name the same gates and `read_gate_targets` compares those two lists — but a name in both is
not a check. A call site deleted or renamed in an edit leaves the gate declared, no longer
evaluated, and mentioned by nothing, and the run would end with the same sentence as a
complete one. This is the statement that the run REACHED each of them, and it is made
whether or not the gates were checked, so `--no-gates` cannot hide a deleted call site
either. A run is held to ACCOUNTED_GATES: the measured gates and the four reproduction gates,
each of which registers whether it is made or declined.
"""
function accounting_failures(declared, reached)
    d = sort(String[String(x) for x in declared])
    r = sort(String[String(x) for x in reached])
    out = String[]
    absent = setdiff(d, r)
    isempty(absent) ||
        push!(out, "gate(s) this script declares and this run never evaluated: " *
                   join(absent, "; "))
    extra = setdiff(r, d)
    isempty(extra) ||
        push!(out, "evaluated under a name nothing declares: " * join(extra, "; "))
    return out
end

# ---------------------------------------------------------------------------
# THE SECOND PROGRAM GUARD. Everything above is definitions: including this file runs
# nothing. Everything below reads the inputs, builds the table, gates it and publishes it,
# and runs only when this file is the program.
# ---------------------------------------------------------------------------

if IS_PROGRAM

# `nothing` when --no-gates was typed; the option parser refuses the silent case.
const TARGETS = OPT.no_gates ? nothing :
                cli_or_exit(() -> read_gate_targets(OPT.gate_targets, MEASURED_GATES))
const GATES = TARGETS === nothing ? nothing : TARGETS.gates
if GATES === nothing
    @warn("--no-gates was given, so the gates whose expected value is a measurement of " *
          "this dataset are NOT checked; nothing in this run compares the call table " *
          "against any number measured before it")
else
    @printf("Gate targets: %d gate(s) read from %s\n", length(GATES),
            basename(OPT.gate_targets))
    TARGETS.census === nothing ||
        @printf("  census the reproduction gates were written for: sha256 %s\n",
                TARGETS.census.sha256)
end

# Both files this run REPLACES are staged, and neither is touched again until the gates at
# the bottom have passed. The relation bodies and the OBT path cache are not here: those are
# accumulating caches whose entries a later run reuses, and a run that adds to them and then
# fails a gate has cost nothing and destroyed nothing.
# Both inputs are named before anything is staged: a missing one is a plain statement of
# what is absent and which step writes it, not a CSV parse error and a stack trace, and no
# _staging directory is created for a run that cannot start. A census named by --census is
# not a registry input, so it is checked by the path that was typed, and the deposited one is
# then not required: this run does not read it.
OPT.census_given ? DataInputs.require_or_exit("pool_taxa", "pool_background") :
                   DataInputs.require_or_exit("census_habitat", "pool_taxa", "pool_background")
if OPT.census_given && !isfile(OPT.census)
    println(stderr, "habitat_reaggregate.jl: --census names $(OPT.census), and there is no " *
                    "file there")
    exit(1)
end

# WHICH CENSUS THIS RUN READS, BY CONTENT. Read once, here, before anything is staged or
# fetched: its digest decides whether the four comparisons against the committed calls are
# made (`reproduction_scope`), so a run that will decline them says so now rather than after
# the whole pass, and the census is parsed below from these same bytes.
const CENSUS_BYTES  = read(OPT.census)
const CENSUS_SHA256 = census_digest(CENSUS_BYTES)
const SCOPE         = reproduction_scope(TARGETS, CENSUS_SHA256)
if SCOPE.evaluated
    println("The four gates that compare this run against the committed calls are made: " *
            SCOPE.detail)
else
    @warn("The four gates that compare this run against the COMMITTED calls are not " *
          "applicable and are not checked: " * SCOPE.detail * ". Nothing here says this run " *
          "reproduces the committed calls")
end

# The seed directories, read here rather than below, because whether this run brings its own
# response bodies decides whether the next check applies.
seed_dirs = String[]
append!(seed_dirs, OPT.seed_dirs)
haskey(ENV, "LNFINDER_RELATION_SEED") &&
    append!(seed_dirs, filter(!isempty, strip.(split(ENV["LNFINDER_RELATION_SEED"], ";"))))

# THE TWO CACHES ARE REQUIRED TOO, AND THIS IS THE CHECK THAT KEEPS THE RUN OFFLINE. A cache
# MISS in this step is a LIVE FETCH, and a missing ontology cache sends the ontology climb
# back to the network as well. Without this check a cache put in place one level wrong would
# pass and then issue one request per taxon — on a run this repository documents as
# offline, with nothing in the output saying the bodies came from the service rather than
# from disk. `relation_bodies` names the BODIES, not the directory holding them, because
# `omnicrobe_cache/` existing with nothing underneath it is exactly the shape of that
# mistake. A run that was given seed directories is left alone: it is bootstrapping the
# cache on purpose, which is the one case where not having it yet is the point.
if isempty(seed_dirs)
    for n in ("relation_bodies", "relation_manifest", "obt_path_cache")
        try
            DataInputs.require_input(n)
        catch e
            e isa DataInputs.MissingInput || rethrow()
            showerror(stderr, e)
            println(stderr)
            println(stderr,
                    "  Prime the cache from directories of already-fetched <taxid>.json\n" *
                    "  bodies, given as bare arguments or in \$LNFINDER_RELATION_SEED: every\n" *
                    "  taxon whose body is there is read from disk.\n" *
                    "  Whatever the cache does not hold, this step fetches - one request per\n" *
                    "  uncached taxon, and one more for every ontology term the OBT cache is\n" *
                    "  missing. A run that fetches measures the habitat layer as the service\n" *
                    "  answers today rather than replaying the call table the data set\n" *
                    "  carries, and is not the offline run this step otherwise reports itself\n" *
                    "  as.")
            exit(1)
        end
    end
end

const PLAN  = open_staging((OPT.out, MANIFEST))
const STAGE = (calls = PLAN[1].staged, manifest = PLAN[2].staged)
@printf("Staging this run's outputs in %s beside each destination\n", STAGING_DIR)

# ---- 1. seed the relation cache by local copy, no network -------------------------
mkpath(REL_DIR)

# taxid => (source, fetched_at). Filled by the seeding and the live fetch below; an earlier
# run's manifest supplies the rest, so a re-run does not restamp old bodies with today.
prov = Dict{Int,Tuple{String,String}}()
if isfile(MANIFEST)
    for r in eachrow(CSV.read(MANIFEST, DataFrame))
        prov[Int(r.taxid)] = (sv(r.source), sv(r.fetched_at))
    end
end

n_seeded = 0; n_seen = 0
for d in seed_dirs
    if !isdir(d)
        @warn "seed directory not found; skipped" dir = d
        continue
    end
    for f in readdir(d)
        endswith(f, ".json") || continue
        t = tryparse(Int, f[1:end-5]); t === nothing && continue
        global n_seen += 1
        src = joinpath(d, f); dst = joinpath(REL_DIR, f)
        # The fetch date travels with the body: `cp` does not carry mtime across, and a body
        # dated "whenever it happened to be copied" is provenance that states nothing.
        stamp = iso_date(mtime(src))
        if isfile(dst)
            get!(prov, t, ("seed", stamp))
            continue
        end
        cp(src, dst)
        global n_seeded += 1
        prov[t] = ("seed", stamp)
    end
end
isempty(seed_dirs) ? println("No seed directory given; the cache is used as it stands.") :
    @printf("Seeded relation cache: %d of %d offered bodies copied (rest already present)\n",
            n_seeded, n_seen)

# ---- 2. taxid universe, and a network fetch for whatever the seed did not cover ----
census = CSV.read(CENSUS_BYTES, DataFrame)
pool   = CSV.read(POOL, DataFrame)
prok   = in.(sv.(census.kingdom), Ref(("Bacteria", "Archaea")))
universe = sort!(collect(union(Set(Int.(skipmissing(census.TaxID[prok]))),
                              Set(Int.(pool.taxid)))))
@printf("Census file: %s\n", OPT.census)
@printf("Census sha256: %s\n", CENSUS_SHA256)
@printf("Census: %d rows, %d with a TaxID (%d unique, %d prokaryotic)\n", nrow(census),
        count(!ismissing, census.TaxID), length(unique(skipmissing(census.TaxID))),
        length(unique(skipmissing(census.TaxID[prok]))))
@printf("Taxid universe: %d (census prokaryotic, union the %d-taxon pool)\n",
        length(universe), nrow(pool))
n_todo = count(t -> !isfile(joinpath(REL_DIR, "$t.json")), universe)
@printf("Not yet cached: %d taxa - these are fetched over the network\n", n_todo)

relations = Dict{Int,Vector{EnvironmentTools.Relation}}()
fetch_failed = Int[]
n_live = 0
@showprogress desc = "Reading relations " for t in universe
    was_cached = isfile(joinpath(REL_DIR, "$t.json"))
    # THIS is the run the retry exists for: one transient 502 among thousands of GETs
    # otherwise lands its taxon in fetch_failed, which a tolerance-zero gate reads as a
    # changed dataset -- after the whole fetch has been paid for. The pause keeps a long
    # serial run polite to a public academic service.
    r = EnvironmentTools.fetch_relations(t; cache_dir = REL_DIR, retries = 3,
                                         retry_delay = 2.0, backoff = 2.0, pause = 0.2)
    # `nothing` is a FAILED fetch. It is recorded as such and kept out of the aggregation:
    # a blank call from a lost request must never be read as "this taxon has no habitat".
    if r === nothing
        push!(fetch_failed, t)
        continue
    end
    if !was_cached
        global n_live += 1
        prov[t] = ("live", Dates.format(Dates.today(), "yyyy-mm-dd"))
    end
    relations[t] = r
end
@printf("Relations held: %d taxa, %d relations; fetched live this run: %d; fetch FAILED: %d\n",
        length(relations), sum(length, values(relations); init = 0), n_live, length(fetch_failed))
isempty(fetch_failed) || @warn "taxa whose relation fetch failed are left out of the table, " *
                              "not recorded as unannotated" n = length(fetch_failed)

# ---- manifest: what is in the cache, where each body came from, and when ----------
mrows = NamedTuple[]
for t in universe
    haskey(relations, t) && isfile(joinpath(REL_DIR, "$t.json")) || continue
    s, d = get(prov, t, ("cache", iso_date(mtime(joinpath(REL_DIR, "$t.json")))))
    push!(mrows, (taxid = t, n_relations = length(relations[t]), source = s, fetched_at = d,
                  omnicrobe_release = OMNICROBE_RELEASE, obt_release = OBT_RELEASE))
end
CSV.write(STAGE.manifest, DataFrame(mrows))
@printf("Staged the manifest for %s (%d cached bodies)\n", MANIFEST, length(mrows))

# ---- 3. fill the OBT path cache, then run strictly offline ------------------------
# Only the terms the walk actually REACHES need a path, and which terms those are depends on
# the paths already known: a term that resolves to a generic node or to a culture medium is
# dropped, which lets the next term into the max_distinct budget. So the discovery pass is
# repeated until it reaches nothing new - a fixpoint, not a single sweep.
n_before = EnvironmentTools.load_obt_cache!(OBT_CSV)
@printf("Primed OBT path cache with %d entries\n", n_before)

reached = Set{String}()
record_path(oid) = (push!(reached, oid); EnvironmentTools.get_obt_path(oid))
unresolved = String[]
for attempt in 1:3
    empty!(reached)
    @showprogress desc = "Resolving OBT paths (pass $attempt) " for t in universe
        rels = get(relations, t, nothing)
        rels === nothing && continue
        for cfg in CONFIGS
            aggregate_habitats(rels, record_path; max_distinct = MAX_DISTINCT, cfg.kw...)
        end
    end
    global unresolved = sort!([o for o in reached if !haskey(EnvironmentTools.OBT_CACHE, o)])
    @printf("  pass %d: %d terms reached, %d still unresolved\n",
            attempt, length(reached), length(unresolved))
    isempty(unresolved) && break
end
@printf("Saved OBT path cache (%d entries, %d added this run)\n",
        EnvironmentTools.save_obt_cache(OBT_CSV),
        length(EnvironmentTools.OBT_CACHE) - n_before)

# From here the lookup is strict: an id that never made it into the cache must raise, not
# quietly return "" and type its term as "Other".
path_strict(oid) = EnvironmentTools.get_obt_path(oid; strict = true)

# ---- 4. aggregate every taxon under every config ----------------------------------
# Per-taxon labels. The census supplies kingdom/phylum/genus, the pool phylum/domain; where
# both know a taxon the census wins, since its lineage comes from the taxonomy database
# rather than from a lineage-string heuristic.
c_king = Dict{Int,String}(); c_phy = Dict{Int,String}(); c_gen = Dict{Int,String}()
c_spec = Dict{Int,String}(); c_nseq = Dict{Int,Int}(); c_nln = Dict{Int,Int}()
c_type = Dict{Int,String}(); c_chem = Dict{Int,String}()
for r in eachrow(census)
    ismissing(r.TaxID) && continue
    t = Int(r.TaxID)
    c_nseq[t] = get(c_nseq, t, 0) + 1
    sv(r.classification) == "Lanthanide" && (c_nln[t] = get(c_nln, t, 0) + 1)
    get!(c_type, t, sv(r.habitat_type)); get!(c_chem, t, sv(r.habitat_chemistry))
    for (d, v) in ((c_king, r.kingdom), (c_phy, r.phylum), (c_gen, r.genus), (c_spec, r.species))
        s = sv(v)
        (isempty(s) || haskey(d, t)) && continue
        d[t] = s
    end
end
p_phy = Dict{Int,String}(); p_dom = Dict{Int,String}(); p_org = Dict{Int,String}()
for r in eachrow(pool)
    t = Int(r.taxid)
    p_phy[t] = sv(r.phylum); p_dom[t] = sv(r.domain); p_org[t] = sv(r.organism)
end
phylum_of(t) = (s = get(c_phy, t, ""); isempty(s) ? get(p_phy, t, "") : s)
organism_of(t) = (s = get(p_org, t, ""); isempty(s) ? get(c_spec, t, "") : s)

# A text-mined sentence and a curated strain record are not the same evidence, and the split
# cannot be recovered from the aggregate afterwards, so it is carried per taxon here.
const RECORD_SOURCES = ("GenBank", "DSMZ", "CIRM")
is_record(s) = any(p -> startswith(s, p), RECORD_SOURCES)

blank = aggregate_habitats(EnvironmentTools.Relation[], identity)
calls = Dict{String,Dict{Int,Any}}(c.name => Dict{Int,Any}() for c in CONFIGS)
rows = NamedTuple[]
@showprogress desc = "Aggregating " for t in universe
    rels = get(relations, t, nothing)
    ok = rels !== nothing
    nrec = ok ? count(r -> is_record(r.source), rels) : 0
    for cfg in CONFIGS
        a = ok ? aggregate_habitats(rels, path_strict; max_distinct = MAX_DISTINCT, cfg.kw...) :
                 blank
        calls[cfg.name][t] = a
        push!(rows, (config = cfg.name, taxid = t,
                     in_census = haskey(c_nseq, t), in_pool = haskey(p_dom, t),
                     n_census_seqs = get(c_nseq, t, 0), n_ln_seqs = get(c_nln, t, 0),
                     kingdom = get(c_king, t, ""), domain = get(p_dom, t, ""),
                     phylum = phylum_of(t), genus = get(c_gen, t, ""), organism = organism_of(t),
                     fetch_ok = ok,
                     n_src_record = nrec, n_src_text = (ok ? length(rels) : 0) - nrec,
                     name = ismissing(a.name) ? "" : String(a.name),
                     category = ismissing(a.category) ? "" : String(a.category),
                     confidence = round(a.confidence, digits = 4),
                     habitat_type = ismissing(a.habitat_type) ? "" : String(a.habitat_type),
                     type_confidence = round(a.type_confidence, digits = 4),
                     chemistry = a.chemistry,
                     chemistry_share = round(a.chemistry_share, digits = 4),
                     chemistry_coverage = round(a.chemistry_coverage, digits = 4),
                     n_relations = a.n_relations, n_distinct = a.n_distinct,
                     n_truncated = a.n_truncated,
                     profile = a.profile, type_profile = a.type_profile,
                     chemistry_profile = a.chemistry_profile,
                     type_counts = a.type_counts, type_mentions = a.type_mentions,
                     type_weights = EnvironmentTools._weights_str(a.type_weights),
                     n_tied = a.n_tied, resolved_by = String(a.resolved_by),
                     chemistry_flags = a.chemistry_flags,
                     chemistry_flag_counts = a.chemistry_flag_counts,
                     n_chem_terms = a.n_chem_terms,
                     n_isolation_terms = a.n_isolation_terms,
                     n_generic_terms = a.n_generic_terms,
                     isolation_only = a.isolation_only, generic_only = a.generic_only))
    end
end
CSV.write(STAGE.calls, DataFrame(rows))
@printf("Staged the call table for %s (%d rows = %d configs x %d taxa)\n", OPT.out,
        length(rows), length(CONFIGS), length(universe))
# A verification re-run writes somewhere else on purpose, and the gates below still run
# against the table this run just built. Saying which file that is keeps a report from being
# read as a statement about the committed table when it is not one.
OPT.out == OUT ||
    @printf("  (--out was given: the committed table at %s is not this run's destination)\n",
            OUT)

# ---- 5. replication report --------------------------------------------------------
# The committed arms, rebuilt exactly as the enrichment script defines them, so the gates
# below compare like with like: a usable committed call is a non-empty habitat type that is
# not an isolation source.
committed_usable(ty) = !isempty(ty) && !is_isolation_source(ty)
ln_all = sort!(unique(Int.(skipmissing(
    census.TaxID[(string.(census.classification) .== "Lanthanide") .& prok]))))
ln_used = filter(t -> committed_usable(get(c_type, t, "")), ln_all)
ln_iso  = filter(t -> is_isolation_source(get(c_type, t, "")), ln_all)

bg = unique(CSV.read(BACKGROUND, DataFrame), :taxid)
census_tids = Set(Int.(skipmissing(census.TaxID)))
b_type = Dict{Int,String}(); b_chem = Dict{Int,String}()
for r in eachrow(bg)
    b_type[Int(r.taxid)] = sv(r.habitat_type); b_chem[Int(r.taxid)] = sv(r.habitat_chemistry)
end
frame   = sort!([Int(t) for t in bg.taxid if !(Int(t) in census_tids)])
bg_used = filter(t -> committed_usable(b_type[t]), frame)
bg_iso  = filter(t -> is_isolation_source(b_type[t]), frame)

got(cfg, t, f) = (v = getfield(calls[cfg][t], f); ismissing(v) ? "" : String(v))
has_type(cfg, t) = !ismissing(calls[cfg][t].habitat_type)

failures = String[]
# Which of the declared measured gates this run actually reached. A call site deleted or
# renamed in an edit would otherwise leave its gate declared in the targets file, never
# evaluated, and unmentioned — and the report would still end with "All replication gates
# passed" over a check that no longer exists.
gates_reached = String[]

# The two comparisons a gate can make, written once, so the gate whose numbers come from the
# targets file and the gate whose numbers are the run's own cannot drift apart.
gate_ok(kind, measured, expected, tol) =
    kind == "at_least" ? measured >= expected - tol : abs(measured - expected) <= tol
# A whole tolerance prints as a whole number: the file's tolerances arrive as floats, and
# "+/-5.0 taxa" reads as a precision the count does not have.
gate_tol(tol) = tol == round(tol) ? string(Int(tol)) : string(tol)
gate_qual(kind, tol) = kind == "at_least" ? (tol == 0 ? ">= " : ">=-$(gate_tol(tol)) ") :
                                            (tol == 0 ? "   " : "+/-$(gate_tol(tol)) ")

function report_gate(label, measured, kind, expected, tol)
    ok = gate_ok(kind, measured, expected, tol)
    qual = gate_qual(kind, tol)
    @printf("  %-46s %6d   expected %s%-6d %s\n",
            label, measured, qual, expected, ok ? "OK" : "MISMATCH")
    ok || push!(failures,
                @sprintf("%s: measured %d, expected %s%d", label, measured, qual, expected))
    return ok
end

"""
    gate(label, measured, expected; tol, atleast) -> Bool

A gate whose expected value is a requirement of the METHOD: zero, where zero means "nothing
may be missing", or a count this run computes for itself. Neither is a quantity measured on
this dataset, so both may be written in this file.
"""
gate(label, measured, expected; tol = 0, atleast = false) =
    report_gate(label, measured, atleast ? "at_least" : "equal", expected, tol)

"""
    measured_gate(label, measured; baseline = nothing) -> Bool

A gate whose expected value is a MEASUREMENT of this dataset. Its expected value, tolerance
and comparison come from the targets file by name, and nothing about the number is written
in this repository. `baseline` is the run's own count, for the gates declared `from_run` in
MEASURED_GATES: those take their tolerance from the file and their expected value from the
run, and exactly one of the two sources must supply it — a gate given both is a number
stated twice that can disagree, and a gate given neither is a gate measured against nothing.

Under --no-gates the gate is announced as NOT CHECKED rather than printed as met: a skipped
check that reads like a passed one is worse than no check at all.
"""
function measured_gate(label, measured; baseline = nothing)
    reach_gate!(gates_reached, label)
    if GATES === nothing
        @printf("  %-46s %6d   NOT CHECKED (--no-gates)\n", label, measured)
        return true
    end
    # read_gate_targets has already refused a file that does not declare this gate.
    g = GATES[label]
    (baseline === nothing) == (g.expected === nothing) &&
        error("the gate '$label' takes its expected value either from this run or from the " *
              "targets file, and takes it from " *
              (baseline === nothing ? "neither" : "both"))
    return report_gate(label, measured, g.kind,
                       g.expected === nothing ? baseline : g.expected, g.tol)
end

# A declined gate registers as reached; see `decline_gate!`.
declined_gate(label, measured, reason) = decline_gate!(gates_reached, label, measured, reason)

println("\n", "="^78)
println("REPLICATION REPORT")
println("="^78)
# REPRODUCING THE COMMITTED CALLS IS A QUESTION ABOUT ONE CENSUS. Asked of a run that reads a
# different census it is not a weaker test but a different test: `ln_used` and the background
# frame are different sets of taxa, and any disagreement says only "these are two datasets",
# which is already known. Asked anyway it would fail an honest re-annotation run, or - worse -
# pass one and be read as evidence that the census this run reads reproduces the one the
# committed calls belong to. So the four gates below are made only when SCOPE says the targets
# file declares this run's census by digest; otherwise each is declined explicitly and said
# out loud with its reason, and what the run replicates is left to the gates that measure THIS
# dataset: the coverage and recovery gates after them, whose expected values the targets file
# must then carry as measurements of the census the run reads. The targets file is read the
# same way either way, so it still declares "Ln chemistry reproduced" with its tolerance; a
# run that declines that gate prints its measurement and applies no tolerance to it.
# A reproduction gate the run makes registers as reached exactly as a declined one does, so
# the accounting guard holds all four whichever way SCOPE went.
function reproduction_gate(label, measured, expected; kw...)
    SCOPE.evaluated || return declined_gate(label, measured, SCOPE.reason)
    reach_gate!(gates_reached, label)
    return gate(label, measured, expected; kw...)
end

println("legacy config against the committed calls, taxon by taxon")
reproduction_gate("Ln habitat_type reproduced",
                  count(t -> got("legacy", t, :habitat_type) == c_type[t], ln_used),
                  length(ln_used))
reproduction_gate("background habitat_type reproduced",
                  count(t -> got("legacy", t, :habitat_type) == b_type[t], bg_used),
                  length(bg_used))
# The Ln arm is measured against its own size, which this run computed; what the targets file
# supplies is the ALLOWANCE the recorded criterion carries for committed chemistry calls the
# legacy config does not reproduce. The background arm has no such allowance, so its expected
# value is the arm size and nothing about it is a number measured on this dataset.
if SCOPE.evaluated
    measured_gate("Ln chemistry reproduced",
                  count(t -> got("legacy", t, :chemistry) == c_chem[t], ln_used);
                  baseline = length(ln_used))
else
    declined_gate("Ln chemistry reproduced",
                  count(t -> got("legacy", t, :chemistry) == c_chem[t], ln_used),
                  SCOPE.reason)
end
reproduction_gate("background chemistry reproduced",
                  count(t -> got("legacy", t, :chemistry) == b_chem[t], bg_used),
                  length(bg_used); atleast = true)

println("\nprimary config: taxa with at least one usable type term")
@printf("  %-46s %6d\n", "Ln prokaryotic taxa in the census", length(ln_all))
measured_gate("Ln taxa with a usable type term", count(t -> has_type("primary", t), ln_all))
@printf("  %-46s %6d\n", "background frame taxa", length(frame))
measured_gate("frame taxa with a usable type term", count(t -> has_type("primary", t), frame))

# A taxon is RECOVERED at term level when it holds at least one term that is not an isolation
# source: the taxon-level filter threw that evidence away wholesale, the term-level one keeps
# it. Whether the taxon then ends with a usable type is a second, stricter question - the
# generic-term filter can still blank it - so that count is printed beside the gate rather
# than folded into it.
println("\nisolation-source taxa recovered at term level")
recovered(t) = !calls["primary"][t].isolation_only
@printf("  %-46s %6d\n", "Ln taxa typed as an isolation source", length(ln_iso))
measured_gate("Ln, at least one non-isolation term", count(recovered, ln_iso))
@printf("  %-46s %6d\n", "  ... of those, ending with a usable type",
        count(t -> has_type("primary", t), ln_iso))
@printf("  %-46s %6d\n", "frame taxa typed as an isolation source", length(bg_iso))
measured_gate("frame, at least one non-isolation term", count(recovered, bg_iso))
@printf("  %-46s %6d\n", "  ... of those, ending with a usable type",
        count(t -> has_type("primary", t), bg_iso))

# The analysed arms are the taxa a usable type term survives for; a taxon with no habitat at
# all never reaches the contrast, so a missing phylum there cannot bias it.
#
# A MISSING PHYLUM IS REPORTED HERE, NOT GATED, AND THE RULE IT BEARS ON IS ENFORCED WHERE IT
# CAN BE: where the contrast is built. This table is the record of what OmniCrobe says and
# must keep every taxon it answered for; whether a taxon may enter a contrast is a question
# about the CONTRAST, so habitat_enrichment.jl excludes phylum-less taxa at arm assembly and
# holds that as an error on every assembled arm. A gate here would refuse to publish an
# honest call table over a downstream policy. Reported, not gated.
ln_arm = filter(t -> has_type("primary", t), ln_all)
bg_arm = filter(t -> has_type("primary", t), frame)

println("\ncache and ontology completeness")
gate("OBT paths still missing after the fill", length(unresolved), 0)
gate("relation fetches that failed", length(fetch_failed), 0)
@printf("  %-46s %6d   excluded at arm assembly\n", "phylum missing, Ln analysed arm",
        count(t -> isempty(phylum_of(t)), ln_arm))
@printf("  %-46s %6d   excluded at arm assembly\n", "phylum missing, frame analysed arm",
        count(t -> isempty(phylum_of(t)), bg_arm))
@printf("  %-46s %6d\n", "phylum missing, all Ln prokaryotic taxa",
        count(t -> isempty(phylum_of(t)), ln_all))

println("\nfor orientation, not gated")
# `nrec_all`, not `nrec`: the per-taxon count above is a loop local of the same shape, and
# one name for both is a soft-scope ambiguity the whole program block is lowered under.
ntot = sum(length, values(relations); init = 0)
nrec_all = sum(v -> count(r -> is_record(r.source), v), values(relations); init = 0)
@printf("  %-46s %6d\n", "taxa in the table", length(universe))
@printf("  %-46s %6d\n", "relations behind it", ntot)
@printf("  %-46s %6d\n", "distinct OBT terms the walk reached", length(reached))
@printf("  %-46s %5.1f%% text-mined / %.1f%% strain record\n", "relation source split",
        100(ntot - nrec_all) / max(ntot, 1), 100nrec_all / max(ntot, 1))
# Where each primary arm came from, so its size reads as an accounting identity rather than
# as a claim: kept + regained - blanked has to land on the arm size printed above.
for (lab, used, iso) in (("Ln", ln_used, ln_iso), ("frame", bg_used, bg_iso))
    keep = count(t -> has_type("primary", t), used)
    @printf("  %-6s arm: %5d of %5d committed-usable kept, %4d blanked generic, %4d from isolation\n",
            lab, keep, length(used), length(used) - keep,
            count(t -> has_type("primary", t), iso))
end
for cfg in CONFIGS
    @printf("  %-14s ln typed %5d   frame typed %5d   ln ties %5d\n", cfg.name,
            count(t -> has_type(cfg.name, t), ln_all),
            count(t -> has_type(cfg.name, t), frame),
            count(t -> calls[cfg.name][t].n_tied > 1, ln_all))
end

# The accounting guard: every gate this script declares was reached, and no other. See
# `accounting_failures` for what a name in both lists does not prove.
append!(failures, accounting_failures(ACCOUNTED_GATES, gates_reached))

println("\n", "="^78)
if isempty(failures)
    # A run that declined half the gates must not sign off in the words of one that made them
    # all. The closing sentence is the only part of this report most readers keep, so it is
    # the part that has to say which run this was.
    println(GATES === nothing ?
            "Every gate this run MADE passed. The gates whose expected value is a " *
            "measurement of this dataset were NOT made (--no-gates)." :
            SCOPE.evaluated ? "All replication gates passed." :
            "Every gate this run MADE passed.")
    SCOPE.evaluated ||
        println("The four gates that compare this run against the committed calls were NOT " *
                "made ($(SCOPE.reason)).")
else
    println("REPLICATION GATES FAILED - the numbers above are the measured ones, unchanged:")
    for f in failures
        println("  - ", f)
    end
end

# THE PUBLISH GATE. Nothing this run built has reached its destination yet; `publish!` moves
# the staged files only when `failures` is empty and raises otherwise, so the gates and the
# writing cannot come apart in a later edit. This is the last statement of the run.
for p in publish!(PLAN, failures)
    @printf("Published %s\n", p)
end

end  # the second program guard
