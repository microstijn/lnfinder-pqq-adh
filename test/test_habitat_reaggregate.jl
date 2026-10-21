#=
    Offline self-test for scripts/habitat_reaggregate.jl — the step that rebuilds one habitat
    call per taxon per rule variant, and the checks that decide whether that table is allowed
    to replace the committed one. No network and no repository data file: every input below is
    built here, at the size of the property under test.

    The script is included behind its own `PROGRAM_FILE` guards, so including it defines the
    functions and runs nothing: no command line is read, no cache is touched, nothing is
    written. One testset starts it as its own program, offline, against a data root built
    in a temp directory, because which gates a run makes is decided where no include
    reaches.

    EVERY FIXTURE IS ABSTRACT. The gate names in the fixtures are role names such as
    "gate-stated", the expected values and tolerances are small arbitrary integers chosen so
    that one branch is decided, and the tables that get staged and published are two-line
    strings. Nothing here is a quantity measured on the data this script aggregates: a suite
    that encodes a run's own numbers passes only until the next run, and it publishes what the
    run found.

    THE FAILURES THIS SUITE EXISTS TO PREVENT.

    1. A FAILED RUN THAT PUBLISHED ANYWAY. The call table this script builds is not an output
       anybody reads; it is the INPUT of every estimate downstream of it — the arms, the
       enrichment, the grid, the figure. Built straight over its destination and gated
       afterwards, a run whose replication gates FAILED has already replaced that input, and
       the failure message is about a number while the damage is to a file. The staging tests
       below require that a failing gate leaves the previous table byte-unchanged and moves
       nothing, and that what the failed run built stays where it can be inspected.

    2. A GATE SET NOBODY RAN. The gates are the only thing between a silently drifted call
       table and every estimate downstream of it. A run given no targets file must be refused,
       and declining the gates must be something a person typed rather than something a
       default did.

    3. A GATE MEASURED AGAINST NOTHING. A targets file that is missing a gate this script
       declares must be a refusal, not a check that quietly did not happen — and so must a
       file naming a gate this script does not have, a gate stated twice, a gate whose
       comparison shape disagrees with the script's, and a file whose own declared size no
       longer matches what it carries.

    4. A CALL SITE DELETED IN AN EDIT. The targets file and the script naming the same gate
       does not mean the run evaluated it. A gate whose call site was deleted or renamed
       leaves the name declared in both places, evaluated by neither, and mentioned by
       nothing, and the run would end with the same closing sentence as a complete one.

    5. A CONFIGURATION SET THAT STOPS BEING ONE-AT-A-TIME. Every config but `legacy`
       moves exactly one axis off `primary`; that is what lets a downstream sensitivity row
       attribute a change to the axis that caused it. A second axis moved by accident makes
       every row under that name a two-variable result reported as a one-variable one.

    6. A REPRODUCTION ASKED OF THE WRONG CENSUS. The four comparisons against the committed
       calls belong to one census. Made on another they fail an honest run or, worse, pass
       and read as evidence that two censuses agree. Which census a file holds is not
       something its path or a default can say, so the comparisons are made only when the
       targets file declares the census by digest and the run reads that digest; every
       other case — no targets file, no declaration, another digest — must decline them, and
       a declined gate must still count as reached.

    Run:  julia --project=. test/test_habitat_reaggregate.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Test, DataFrames, CSV

# The script goes into a module of its own, the way the grid suite includes the grid: the
# aggregate runner includes every suite into one process, and a script's constants must not
# collide with another suite's.
module HabitatReagg end
Base.include(HabitatReagg, joinpath(@__DIR__, "..", "scripts", "habitat_reaggregate.jl"))
const RG = HabitatReagg

# ---------------------------------------------------------------------------
# Fixtures. A declared gate set of two, one of each kind the script supports: one whose
# expected value comes from the file, one measured against a number the run computes for
# itself and for which the file supplies only a tolerance.
# ---------------------------------------------------------------------------

const DECLARED = ((gate = "gate-stated",   from_run = false, kind = "equal"),
                  (gate = "gate-computed", from_run = true,  kind = "at_least"))

"""
    targets_file(dir, rows; name) -> path

A gate-targets CSV built row by row. Each row is a `(kind, gate, expected, tol, note)` tuple
of strings; an empty string is an empty field, which is how the file says "no expected value
here" and how a malformed row leaves one out.
"""
function targets_file(dir, rows; name = "targets.csv")
    path = joinpath(dir, name)
    CSV.write(path, DataFrame(kind     = String[r[1] for r in rows],
                              gate     = String[r[2] for r in rows],
                              expected = String[r[3] for r in rows],
                              tol      = String[r[4] for r in rows],
                              note     = String[r[5] for r in rows]))
    return path
end

# A complete, well-formed file for DECLARED: two gates and the count that says so.
good_rows() = [("equal",    "gate-stated",   "7", "1", "a fixture note"),
               ("at_least", "gate-computed", "",  "2", "a fixture note"),
               ("expect",   "gates",         "2", "",  "")]

"""
    fails_with(f, fragment) -> Bool

`f()` raises an `ErrorException` whose message carries `fragment`. Testing the message and
not only the type is deliberate: two of the refusals below differ from each other by which
side of a comparison is wrong, and a suite that accepted any error would pass with the two
swapped.
"""
function fails_with(f, fragment::AbstractString)
    try
        f()
    catch e
        e isa ErrorException || rethrow()
        occursin(fragment, e.msg) && return true
        @info "refused, but not with the expected wording" expected = fragment got = e.msg
        return false
    end
    return false
end

@testset "habitat_reaggregate.jl" begin

# ---------------------------------------------------------------------------
@testset "the configuration set this script declares" begin
    names = String[c.name for c in RG.CONFIGS]
    @test length(unique(names)) == length(names)
    @test "primary" in names
    @test "legacy" in names

    # The four axes a config may move. Named here as well as in the script so that an axis
    # added to one config and forgotten in another is a failure rather than a config that
    # silently runs under the library default.
    axes = (:isolation, :drop_generic, :chem_text, :tie_rule)
    for c in RG.CONFIGS
        @test Set(keys(c.kw)) == Set(axes)
    end

    primary = first(c for c in RG.CONFIGS if c.name == "primary")
    for c in RG.CONFIGS
        c.name == "primary" && continue
        moved = count(a -> getfield(c.kw, a) != getfield(primary.kw, a), axes)
        # ONE AXIS AT A TIME, except `legacy`, which moves all four on purpose:
        # it is the only config the committed CSVs can be compared against.
        @test moved == (c.name == "legacy" ? length(axes) : 1)
    end

    # Two configs that differ in no axis are one measurement reported under two names.
    kws = [c.kw for c in RG.CONFIGS]
    @test length(unique(kws)) == length(kws)
end

# ---------------------------------------------------------------------------
@testset "the gate set this script declares" begin
    names = String[g.gate for g in RG.MEASURED_GATES]
    @test length(unique(names)) == length(names)
    @test all(!isempty, names)
    # The SHAPE of each comparison is declared in the script as well as in the file, so a
    # single word changed in the file alone cannot narrow what a gate catches.
    @test all(g -> g.kind in ("equal", "at_least"), RG.MEASURED_GATES)
    @test all(g -> g.from_run isa Bool, RG.MEASURED_GATES)
    @test Set(RG.GATE_COLUMNS) == Set((:kind, :gate, :expected, :tol, :note))
    @test "expect" in RG.GATE_KINDS
    @test "census" in RG.GATE_KINDS
    @test RG.EXPECTED_COUNTS == ("gates",)

    # The four comparisons against the committed calls, each named once; one of them is also a
    # measured gate, because its allowance is a number measured on the data.
    @test length(RG.REPRODUCTION_GATES) == 4
    @test allunique(RG.REPRODUCTION_GATES)
    @test intersect(Set(RG.REPRODUCTION_GATES), Set(names)) == Set(["Ln chemistry reproduced"])
    # The accounting guard holds a run to every measured and every reproduction gate, once.
    @test allunique(RG.ACCOUNTED_GATES)
    @test Set(RG.ACCOUNTED_GATES) == union(Set(names), Set(RG.REPRODUCTION_GATES))
end

# ---------------------------------------------------------------------------
@testset "a targets file is required unless the no-gates flag is typed" begin
    # The default is the dangerous one, so the default is refused.
    @test fails_with("no --gate-targets were given") do
        RG.gate_options(Dict{String,String}(), String[], "out.csv")
    end

    typed = RG.gate_options(Dict("--no-gates" => ""), String[], "out.csv")
    @test typed.no_gates
    @test isempty(typed.gate_targets)
    @test typed.out == "out.csv"

    given = RG.gate_options(Dict("--gate-targets" => "t.csv"), ["seed-a", "seed-b"], "out.csv")
    @test !given.no_gates
    @test given.gate_targets == "t.csv"
    @test given.seed_dirs == ["seed-a", "seed-b"]

    # Both at once is a contradiction; resolving it either way runs something other than what
    # was typed.
    @test fails_with("both given") do
        RG.gate_options(Dict("--gate-targets" => "t.csv", "--no-gates" => ""), String[], "o.csv")
    end

    redirected = RG.gate_options(Dict("--gate-targets" => "t.csv", "--out" => "elsewhere.csv"),
                                 String[], "out.csv")
    @test redirected.out == "elsewhere.csv"

    # THE CENSUS IS AN INPUT, NOT A SETTING. A run against a re-annotated census measures a
    # different dataset, so the path has to travel with the run rather than being edited into
    # a const where the report cannot see it. A command line that names no census reads the
    # deposited one, the `census_habitat` input under the data root.
    default_census = RG.gate_options(Dict("--no-gates" => ""), String[], "out.csv")
    @test default_census.census == RG.CENSUS
    @test RG.CENSUS == RG.DataInputs.input_path("census_habitat")
    @test !default_census.census_given

    elsewhere = RG.gate_options(Dict("--no-gates" => "", "--census" => "reannotated.csv"),
                                String[], "out.csv")
    @test elsewhere.census == "reannotated.csv"
    @test elsewhere.out == "out.csv"
    # Typed, even when the path typed is the deposited one: whether it was typed decides only
    # how a missing census is reported, never which gates are made.
    @test elsewhere.census_given
    @test RG.gate_options(Dict("--no-gates" => "", "--census" => RG.CENSUS),
                          String[], "out.csv").census_given

    @test fails_with("needs a value") do
        RG.gate_options(Dict("--no-gates" => "", "--census" => "  "), String[], "out.csv")
    end
    @test fails_with("given twice") do
        RG.parse_options(["--census", "a.csv", "--census", "b.csv"], RG.OPTION_SPEC)
    end
    @test fails_with("needs a value") do
        RG.parse_options(["--census"], RG.OPTION_SPEC)
    end

    # An option present but empty names the option rather than falling back to the default.
    @test fails_with("needs a value") do
        RG.gate_options(Dict("--gate-targets" => "   "), String[], "out.csv")
    end

    # An ignored token is a silently different run, so the parser refuses every shape of one.
    @test fails_with("unknown option") do
        RG.parse_options(["--gate-target", "t.csv"], RG.OPTION_SPEC)
    end
    @test fails_with("given twice") do
        RG.parse_options(["--out", "a.csv", "--out", "b.csv"], RG.OPTION_SPEC)
    end
    @test fails_with("needs a value") do
        RG.parse_options(["--gate-targets"], RG.OPTION_SPEC)
    end
    # A value that is itself a defined option means the value was forgotten.
    @test fails_with("needs a value") do
        RG.parse_options(["--gate-targets", "--no-gates"], RG.OPTION_SPEC)
    end
    parsed = RG.parse_options(["--gate-targets", "t.csv", "seed-dir"], RG.OPTION_SPEC)
    @test parsed.options["--gate-targets"] == "t.csv"
    @test parsed.positional == ["seed-dir"]
end

# ---------------------------------------------------------------------------
@testset "a gate with no expected value in the file is refused, not skipped" begin
    mktempdir() do dir
        read_back = RG.read_gate_targets(targets_file(dir, good_rows()), DECLARED)
        gates = read_back.gates
        # A file may declare no census: that declines the reproduction gates, it is not a
        # malformed file.
        @test read_back.census === nothing
        @test Set(keys(gates)) == Set(("gate-stated", "gate-computed"))
        @test gates["gate-stated"].expected == 7
        @test gates["gate-stated"].tol == 1
        @test gates["gate-stated"].kind == "equal"
        # The run computes this one's expected value, so the file carries none.
        @test gates["gate-computed"].expected === nothing
        @test gates["gate-computed"].tol == 2

        # THE PROPERTY: a declared gate the file does not carry is an error. Removing the row
        # and lowering the count with it is the self-consistent thinning, and it is the one
        # that has to fail.
        @test fails_with("declares no expected value for gate-stated") do
            RG.read_gate_targets(
                targets_file(dir, [("at_least", "gate-computed", "", "2", "a fixture note"),
                                   ("expect",   "gates",         "1", "",  "")];
                             name = "thinned.csv"), DECLARED)
        end

        # The same refusal against the set the SCRIPT declares, so the property is pinned on
        # the real declaration and not only on the fixture's.
        real_rows = [(g.kind, g.gate, g.from_run ? "" : "3", "1", "a fixture note")
                     for g in RG.MEASURED_GATES]
        push!(real_rows, ("expect", "gates", string(length(real_rows)), "", ""))
        real = RG.read_gate_targets(targets_file(dir, real_rows; name = "real.csv"),
                                    RG.MEASURED_GATES).gates
        @test Set(keys(real)) == Set(String[g.gate for g in RG.MEASURED_GATES])

        dropped = RG.MEASURED_GATES[end]
        short_rows = [(g.kind, g.gate, g.from_run ? "" : "3", "1", "a fixture note")
                      for g in RG.MEASURED_GATES if g.gate != dropped.gate]
        push!(short_rows, ("expect", "gates", string(length(short_rows)), "", ""))
        @test fails_with("declares no expected value for " * dropped.gate) do
            RG.read_gate_targets(targets_file(dir, short_rows; name = "real_short.csv"),
                                 RG.MEASURED_GATES)
        end
    end
end

# ---------------------------------------------------------------------------
@testset "the rest of the targets file's refusals" begin
    mktempdir() do dir
        rows = good_rows()

        @test fails_with("does not exist") do
            RG.read_gate_targets(joinpath(dir, "absent.csv"), DECLARED)
        end

        # A file missing a column the reader needs is a file for another script.
        noteless = DataFrame(kind = ["expect"], gate = ["gates"], expected = ["2"],
                             tol = [""])
        CSV.write(joinpath(dir, "nocol.csv"), noteless)
        @test fails_with("has no 'note'") do
            RG.read_gate_targets(joinpath(dir, "nocol.csv"), DECLARED)
        end

        # A gate this run does not have: a file written for another version of this script.
        @test fails_with("is not a gate of this run") do
            RG.read_gate_targets(
                targets_file(dir, vcat(rows[1:2],
                                       [("equal", "gate-retired", "1", "0", "a fixture note"),
                                        ("expect", "gates", "3", "", "")]);
                             name = "extra.csv"), DECLARED)
        end

        # One gate's expected value stated in two places: only one of them decides the run.
        @test fails_with("declared twice") do
            RG.read_gate_targets(
                targets_file(dir, [rows[1], rows[1], rows[2],
                                   ("expect", "gates", "3", "", "")];
                             name = "twice.csv"), DECLARED)
        end

        # The comparison's shape is declared in both places and has to agree.
        @test fails_with("shape of the") do
            RG.read_gate_targets(
                targets_file(dir, [("at_least", "gate-stated", "7", "1", "a fixture note"),
                                   rows[2], rows[3]]; name = "shape.csv"), DECLARED)
        end

        # A number with no record of where it came from cannot be told from one transcribed
        # out of this script's own output.
        @test fails_with("carries no note") do
            RG.read_gate_targets(
                targets_file(dir, [("equal", "gate-stated", "7", "1", ""), rows[2], rows[3]];
                             name = "nonote.csv"), DECLARED)
        end

        # An absent tolerance is not a tolerance of zero; it is a width nobody decided.
        @test fails_with("no tolerance") do
            RG.read_gate_targets(
                targets_file(dir, [("equal", "gate-stated", "7", "", "a fixture note"),
                                   rows[2], rows[3]]; name = "notol.csv"), DECLARED)
        end

        # A gate measured against a number the run computes must not also be given one.
        @test fails_with("also carries an expected value") do
            RG.read_gate_targets(
                targets_file(dir, [rows[1],
                                   ("at_least", "gate-computed", "5", "2", "a fixture note"),
                                   rows[3]]; name = "both.csv"), DECLARED)
        end

        # These gates are counts of taxa.
        @test fails_with("is a count of taxa") do
            RG.read_gate_targets(
                targets_file(dir, [("equal", "gate-stated", "7.5", "1", "a fixture note"),
                                   rows[2], rows[3]]; name = "fractional.csv"), DECLARED)
        end

        # How much of this table is checked is not something a run may discover from the file
        # it was handed, so the file states its own size and the two have to agree.
        @test fails_with("carries no 'expect' row") do
            RG.read_gate_targets(targets_file(dir, rows[1:2]; name = "nocount.csv"), DECLARED)
        end
        @test fails_with("declares 3 gate(s) and carries 2") do
            RG.read_gate_targets(
                targets_file(dir, [rows[1], rows[2], ("expect", "gates", "3", "", "")];
                             name = "miscount.csv"), DECLARED)
        end
        @test fails_with("states its own size in two places") do
            RG.read_gate_targets(
                targets_file(dir, [rows[1], rows[2], rows[3], rows[3]]; name = "twocount.csv"),
                DECLARED)
        end
        @test fails_with("an 'expect' row counts one of") do
            RG.read_gate_targets(
                targets_file(dir, [rows[1], rows[2], ("expect", "taxa", "2", "", "")];
                             name = "othercount.csv"), DECLARED)
        end
        @test fails_with("not a 'equals'") do
            RG.read_gate_targets(
                targets_file(dir, [("equals", "gate-stated", "7", "1", "a fixture note"),
                                   rows[2], rows[3]]; name = "badkind.csv"), DECLARED)
        end
    end
end

# ---------------------------------------------------------------------------
@testset "every declared gate was actually evaluated, and no other" begin
    @test isempty(RG.accounting_failures(["gate-a", "gate-b"], ["gate-b", "gate-a"]))
    @test isempty(RG.accounting_failures(String[], String[]))

    # A call site deleted or renamed in an edit: the name is declared, nothing evaluates it.
    absent = RG.accounting_failures(["gate-a", "gate-b"], ["gate-a"])
    @test length(absent) == 1
    @test occursin("never evaluated", absent[1])
    @test occursin("gate-b", absent[1])

    extra = RG.accounting_failures(["gate-a"], ["gate-a", "gate-c"])
    @test length(extra) == 1
    @test occursin("nothing declares", extra[1])
    @test occursin("gate-c", extra[1])

    both = RG.accounting_failures(["gate-a", "gate-b"], ["gate-a", "gate-c"])
    @test length(both) == 2
end

# ---------------------------------------------------------------------------
# The census the reproduction gates were written for. A fixture digest is the sha256 of a
# fixture's own bytes, never of any real census.
fixture_digest(s::AbstractString) = RG.census_digest(Vector{UInt8}(s))
census_row(digest; gate = "", tol = "", note = "a fixture census") =
    ("census", gate, digest, tol, note)

@testset "a targets file declares its census by digest, or not at all" begin
    mktempdir() do dir
        d = fixture_digest("census bytes\n")
        t = RG.read_gate_targets(targets_file(dir, vcat(good_rows(), [census_row(d)]);
                                              name = "declared.csv"), DECLARED)
        @test t.census.sha256 == d
        @test t.census.note == "a fixture census"
        # The census row is not a gate, and the count the file states of its gates leaves it
        # out.
        @test Set(keys(t.gates)) == Set(("gate-stated", "gate-computed"))
        # A digest written in capitals is the same digest.
        up = RG.read_gate_targets(targets_file(dir, vcat(good_rows(), [census_row(uppercase(d))]);
                                               name = "upper.csv"), DECLARED)
        @test up.census.sha256 == d
        # A digest made only of decimal digits is still read as the text it is, not as a number.
        digits = repeat("0123456789", 7)[1:64]
        num = RG.read_gate_targets(targets_file(dir, vcat(good_rows(), [census_row(digits)]);
                                                name = "digits.csv"), DECLARED)
        @test num.census.sha256 == digits

        # Each malformed declaration is refused on its own.
        bad(rows, name) = RG.read_gate_targets(targets_file(dir, rows; name = name), DECLARED)
        @test fails_with("declares its census twice") do
            bad(vcat(good_rows(), [census_row(d), census_row(d)]), "twice.csv")
        end
        @test fails_with("names no gate") do
            bad(vcat(good_rows(), [census_row(d; gate = "gate-stated")]), "named.csv")
        end
        @test fails_with("64 hexadecimal digits") do
            bad(vcat(good_rows(), [census_row(d[1:63])]), "short.csv")
        end
        @test fails_with("64 hexadecimal digits") do
            bad(vcat(good_rows(), [census_row("g" * d[2:end])]), "nothex.csv")
        end
        @test fails_with("64 hexadecimal digits") do
            bad(vcat(good_rows(), [census_row("")]), "empty.csv")
        end
        @test fails_with("carries a tolerance") do
            bad(vcat(good_rows(), [census_row(d; tol = "0")]), "tol.csv")
        end
        @test fails_with("carries no note") do
            bad(vcat(good_rows(), [census_row(d; note = "")]), "nonote.csv")
        end
        # A declaration does not excuse a missing gate: the census is not a gate.
        @test fails_with("declares no expected value for gate-stated") do
            bad([good_rows()[2], ("expect", "gates", "1", "", ""), census_row(d)], "thin.csv")
        end
    end
end

@testset "the reproduction gates are made only for the census the targets file declares" begin
    mktempdir() do dir
        census_path = joinpath(dir, "census.csv")
        write(census_path, "TaxID,habitat_type\n1,cat-a\n")
        run_digest = RG.census_digest(read(census_path))
        declared = RG.read_gate_targets(
            targets_file(dir, vcat(good_rows(), [census_row(run_digest)]); name = "t.csv"),
            DECLARED)

        # DECLARED AND MATCHING: the one case in which the four comparisons are made.
        made = RG.reproduction_scope(declared, run_digest)
        @test made.evaluated
        @test isempty(made.reason)
        @test occursin(run_digest, made.detail)
        # The digest is compared as a digest, whatever case it is written in.
        @test RG.reproduction_scope(declared, uppercase(run_digest)).evaluated

        # --NO-GATES: no targets file, so nothing declares a census.
        no_gates = RG.reproduction_scope(nothing, run_digest)
        @test !no_gates.evaluated
        @test no_gates.reason == "--no-gates"

        # DECLARED, AND ANOTHER CENSUS READ: one byte changed is another census.
        write(census_path, "TaxID,habitat_type\n1,cat-b\n")
        other_digest = RG.census_digest(read(census_path))
        @test other_digest != run_digest
        other = RG.reproduction_scope(declared, other_digest)
        @test !other.evaluated
        @test other.reason == "another census"
        @test occursin(run_digest, other.detail) && occursin(other_digest, other.detail)

        # A TARGETS FILE THAT DECLARES NO CENSUS declines them as well.
        silent = RG.read_gate_targets(targets_file(dir, good_rows(); name = "silent.csv"),
                                      DECLARED)
        undeclared = RG.reproduction_scope(silent, run_digest)
        @test !undeclared.evaluated
        @test undeclared.reason == "no census declared"

        # THE PATH DECIDES NOTHING. The same bytes under another name are the declared census;
        # the declared name holding other bytes is not.
        write(census_path, "TaxID,habitat_type\n1,cat-a\n")
        copy_path = joinpath(dir, "renamed_copy.csv")
        cp(census_path, copy_path)
        @test RG.reproduction_scope(declared, RG.census_digest(read(copy_path))).evaluated
        write(census_path, "TaxID,habitat_type\n1,cat-a\n2,cat-b\n")
        @test !RG.reproduction_scope(declared, RG.census_digest(read(census_path))).evaluated
    end
end

@testset "a declined gate is printed, withholds its verdict, and still counts as reached" begin
    # What a declined gate prints, captured.
    function declined_output(reached, label, measured, reason)
        mktemp() do path, io
            ok = redirect_stdout(io) do
                RG.decline_gate!(reached, label, measured, reason)
            end
            flush(io)
            return (ok = ok, text = read(path, String))
        end
    end
    reached = String[]
    d = declined_output(reached, "Ln chemistry reproduced", 3, "--no-gates")
    @test d.ok === true
    @test reached == ["Ln chemistry reproduced"]
    @test occursin("NOT APPLICABLE (--no-gates)", d.text)
    @test occursin("Ln chemistry reproduced", d.text) && occursin(" 3 ", d.text)

    # Declining one name twice is the same fault as evaluating it twice.
    @test fails_with("evaluated twice") do
        declined_output(reached, "Ln chemistry reproduced", 3, "--no-gates")
    end
    @test fails_with("evaluated twice") do
        RG.reach_gate!(reached, "Ln chemistry reproduced")
    end

    # THE PROPERTY: a run that declines all four reproduction gates and makes every measured
    # gate passes the accounting guard, so declining is not mistaken for a deleted call site.
    all_declined = String[]
    for g in RG.REPRODUCTION_GATES
        declined_output(all_declined, g, 0, "another census")
    end
    for g in RG.MEASURED_GATES
        g.gate in all_declined || RG.reach_gate!(all_declined, g.gate)
    end
    @test isempty(RG.accounting_failures(RG.ACCOUNTED_GATES, all_declined))

    # And the guard still catches a reproduction gate whose call site is gone, declined or not.
    for dropped in RG.REPRODUCTION_GATES
        partial = filter(!=(dropped), all_declined)
        f = RG.accounting_failures(RG.ACCOUNTED_GATES, partial)
        @test length(f) == 1 && occursin("never evaluated", f[1]) && occursin(dropped, f[1])
    end
end

# ---------------------------------------------------------------------------
@testset "nothing is published until the gates have passed" begin
    mktempdir() do dir
        dest  = joinpath(dir, "calls.csv")
        other = joinpath(dir, "cache", "manifest.csv")
        mkpath(dirname(other))
        write(dest,  "previous,table\n1,2\n")
        write(other, "previous,manifest\n")
        was_dest  = read(dest)
        was_other = read(other)

        plan = RG.open_staging((dest, other))
        @test length(plan) == 2
        # THE PROPERTY THIS RESTS ON: what the run builds is not where the run publishes.
        @test all(p -> p.staged != p.dest, plan)
        @test all(p -> basename(dirname(p.staged)) == RG.STAGING_DIR, plan)
        @test all(p -> basename(p.staged) == basename(p.dest), plan)
        # Beside the destination, so the move that publishes is within one filesystem.
        @test all(p -> dirname(dirname(p.staged)) == dirname(p.dest), plan)
        @test [p.dest for p in plan] == [abspath(dest), abspath(other)]

        write(plan[1].staged, "new,table\n3,4\n")
        write(plan[2].staged, "new,manifest\n")

        # A FAILING GATE PUBLISHES NOTHING, and says so.
        failed = try
            RG.publish!(plan, ["a replication gate failed"])
            nothing
        catch e
            e
        end
        @test failed isa ErrorException
        @test occursin("NOTHING WAS PUBLISHED", failed.msg)
        @test occursin("a replication gate failed", failed.msg)
        @test occursin(plan[1].staged, failed.msg)
        @test occursin(plan[1].dest, failed.msg)
        # Byte-unchanged: the previous table is still the input every downstream step reads.
        @test read(dest)  == was_dest
        @test read(other) == was_other
        # And what the failed run built is still there, so the failure is inspectable.
        @test isfile(plan[1].staged)
        @test isfile(plan[2].staged)
        # Two failed gates are reported once, both of them.
        two = try
            RG.publish!(plan, ["first failure", "second failure"])
            nothing
        catch e
            e
        end
        @test occursin("first failure", two.msg) && occursin("second failure", two.msg)
        @test read(dest) == was_dest

        # With no failures the same call publishes, and only then.
        written = RG.publish!(plan, String[])
        @test written == [abspath(dest), abspath(other)]
        @test read(dest)  == Vector{UInt8}("new,table\n3,4\n")
        @test read(other) == Vector{UInt8}("new,manifest\n")
        @test !isdir(dirname(plan[1].staged))
        @test !isdir(dirname(plan[2].staged))
    end
end

# ---------------------------------------------------------------------------
@testset "the staging directory carries nothing from an earlier run" begin
    mktempdir() do dir
        dest = joinpath(dir, "calls.csv")
        write(dest, "previous\n")
        was = read(dest)

        first_try = RG.open_staging((dest,))
        write(first_try[1].staged, "half-finished\n")
        write(joinpath(dirname(first_try[1].staged), "stray.csv"), "left behind\n")

        # A second attempt must not be able to publish what the first one left.
        second = RG.open_staging((dest,))
        @test !isfile(second[1].staged)
        @test isempty(readdir(dirname(second[1].staged)))

        # The gates passed and there is nothing to publish: that is a run that would leave the
        # destination holding another run's numbers while reporting success.
        @test fails_with("is not there to publish") do
            RG.publish!(second, String[])
        end
        @test read(dest) == was

        # Two outputs that resolve to one destination: one of them would be lost silently.
        @test fails_with("would be published to one path") do
            RG.open_staging((dest, joinpath(dir, ".", "calls.csv")))
        end
    end
end

# ---------------------------------------------------------------------------
# THE PROGRAM, RUN. Which gates a run makes is decided inside the second program guard, which
# no include reaches, so the decision is measured from outside: the script is started as its
# own program against a data root built here, offline — every relation body is on disk and
# empty, so no request is made. The committed calls in the fixture are ones the `legacy`
# config cannot reproduce from empty bodies, so a run that MAKES the two habitat_type
# comparisons fails them and publishes nothing, and a run that declines them publishes. The
# three runs start together and are read back once all have ended.
# ---------------------------------------------------------------------------
@testset "the program makes the reproduction gates only for the declared census" begin
    script = joinpath(@__DIR__, "..", "scripts", "habitat_reaggregate.jl")
    census_csv = "TaxID,kingdom,classification,habitat_type,habitat_chemistry," *
                 "phylum,genus,species\n" *
                 "101,Bacteria,Lanthanide,cat-a,,phy-1,gen-1,Org 101\n" *
                 "102,Bacteria,Calcium,,,phy-1,gen-2,Org 102\n"
    # Every fixture file goes where the registry puts it under the fixture's own data root.
    function data_root()
        root = mktempdir()
        at(name) = RG.DataInputs.input_path(name; env = Dict("LNFINDER_DATA" => root))
        put(name, text) = (mkpath(dirname(at(name))); write(at(name), text))
        put("census_habitat", census_csv)
        put("pool_taxa", "taxid,phylum,domain,organism\n" *
                         "101,phy-1,Bacteria,Org 101\n201,phy-2,Bacteria,Org 201\n")
        put("pool_background", "taxid,habitat_type,habitat_chemistry\n201,cat-a,\n101,cat-a,\n")
        put("obt_path_cache", "obtid,path\n")
        mkpath(at("relation_bodies"))
        manifest = "taxid,n_relations,source,fetched_at,omnicrobe_release,obt_release\n"
        for t in (101, 102, 201)
            write(joinpath(at("relation_bodies"), "$t.json"), "[]")
            manifest *= "$t,0,cache,2000-01-01,fixture,fixture\n"
        end
        put("relation_manifest", manifest)
        return root
    end
    targets(dir, digest) = targets_file(dir,
        vcat([(g.kind, g.gate, g.from_run ? "" : "0", "0", "a fixture note")
              for g in RG.MEASURED_GATES],
             [("expect", "gates", string(length(RG.MEASURED_GATES)), "", "")],
             digest === nothing ? Tuple{String,String,String,String,String}[] :
                                  [census_row(digest)]); name = "targets.csv")
    function start(args)
        root = data_root()
        out = joinpath(root, "out", "calls.csv")
        # Two files, because one path opened twice is two writers over the same bytes.
        logs = (joinpath(root, "stdout.log"), joinpath(root, "stderr.log"))
        all_args = vcat(args, ["--out", out])
        cmd = addenv(ignorestatus(`$(Base.julia_cmd()) --startup-file=no $script $all_args`),
                     "LNFINDER_DATA" => root, "LNFINDER_RELATION_SEED" => "")
        return (process = run(pipeline(cmd; stdout = logs[1], stderr = logs[2]); wait = false),
                root = root, out = out, logs = logs)
    end
    function finish(s)
        wait(s.process)
        return (code = s.process.exitcode, out = s.out,
                log = join((isfile(f) ? read(f, String) : "" for f in s.logs), "\n"))
    end
    digest = RG.census_digest(Vector{UInt8}(census_csv))
    other = RG.census_digest(Vector{UInt8}(census_csv *
                                           "103,Bacteria,Calcium,,,phy-1,gen-3,Org 103\n"))
    tdir = mktempdir()
    declared, other_dir = joinpath(tdir, "declared"), joinpath(tdir, "other")
    mkpath(declared); mkpath(other_dir)
    runs = (start(["--gate-targets", targets(declared, digest)]),
            start(["--no-gates"]),
            start(["--gate-targets", targets(other_dir, other)]))
    made, no_gates, elsewhere = finish.(runs)
    reproduction_lines(log) = filter(l -> any(g -> occursin(g, l), RG.REPRODUCTION_GATES),
                                     split(log, '\n'))

    # DECLARED AND MATCHING: the comparisons are made, the fixture fails them, nothing is
    # published.
    @test made.code == 1
    @test occursin("are made: the targets file declares them for the census with sha256 " *
                   digest, made.log)
    @test occursin(r"Ln habitat_type reproduced\s+0\s+expected\s+1\s+MISMATCH", made.log)
    @test !occursin("NOT APPLICABLE", made.log)
    @test occursin("NOTHING WAS PUBLISHED", made.log)
    @test !isfile(made.out)
    # Every gate it made registered: the only failures are the fixture's mismatches.
    @test !occursin("never evaluated", made.log)

    # --NO-GATES: all four declined with the reason, the run publishes, and says what it did not
    # make.
    @test no_gates.code == 0
    @test count(l -> occursin("NOT APPLICABLE (--no-gates)", l),
                reproduction_lines(no_gates.log)) == 4
    @test occursin("Census sha256: " * digest, no_gates.log)
    @test occursin("were NOT made (--no-gates)", no_gates.log)
    @test isfile(no_gates.out)

    # DECLARED FOR ANOTHER CENSUS: all four declined, and the run still publishes, which it
    # can only do if the declined gates registered with the accounting guard.
    @test elsewhere.code == 0
    @test count(l -> occursin("NOT APPLICABLE (another census)", l),
                reproduction_lines(elsewhere.log)) == 4
    @test occursin(other, elsewhere.log) && occursin(digest, elsewhere.log)
    @test occursin("Every gate this run MADE passed.", elsewhere.log)
    @test !occursin("never evaluated", elsewhere.log)
    @test isfile(elsewhere.out)
    for r in (made, no_gates, elsewhere)
        r.code in (0, 1) || println(r.log)
    end
end

end
