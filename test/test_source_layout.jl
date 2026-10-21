#=
    test_source_layout.jl — every file this repository names must exist, spelled exactly.

    WHY THIS SUITE EXISTS. Julia resolves `include` through the operating system, and
    Windows and a default macOS volume are case-insensitive: there, `include("Foo.jl")`
    happily opens `foo.jl` and nothing ever complains. On a case-sensitive filesystem —
    every Linux clone, every ubuntu-latest CI runner — the same line is a fatal
    `SystemError`, and because it fires at load time it takes the package, the whole test
    suite and every script down with it. A suite that runs only on a case-insensitive
    filesystem can never see the difference, so this one does not ask the filesystem
    whether a file exists. It reads the directory and compares names byte for byte.

    WHAT IT CHECKS.
      1. Every `include` target in src/, scripts/ and test/ resolves to a file that is
         really there, spelled the way the `include` spells it.
      2. The data files the code names by path are there, spelled the same way.
      3. No two files in the tree have names that differ only by case, which would make
         a checkout ambiguous on a case-insensitive filesystem.
      4. Every `.jl` file under src/, scripts/ and test/ parses. Nothing in this suite
         loads or runs a script, so parsing is the only reading its syntax gets.

    Comments and docstrings are stripped before the scan, so a usage example in a header
    is documentation and is not held to the same standard as code.

    Run:  julia --project=. test/test_source_layout.jl
=#

using Test

const LAYOUT_REPO = normpath(joinpath(@__DIR__, ".."))
const LAYOUT_DIRS = ("src", "scripts", "test")

"""
    strip_julia_comments(src) -> String

`src` with `#= … =#` blocks (nested), `#` line comments and `\"\"\"` docstrings blanked
out, newlines preserved so reported line numbers stay meaningful. Crude but sufficient:
its only job is to keep a usage example inside a header comment from being read as code.
"""
function strip_julia_comments(src::AbstractString)
    out = IOBuffer()
    i, n = 1, lastindex(src)
    depth = 0                 # nesting depth of #= =#
    in_doc = false            # inside a \"\"\" … \"\"\" string
    while i <= n
        rest = SubString(src, i)
        if depth > 0
            if startswith(rest, "#=")
                depth += 1; i += 2; continue
            elseif startswith(rest, "=#")
                depth -= 1; i += 2; continue
            end
            src[i] == '\n' ? write(out, '\n') : write(out, ' ')
            i = nextind(src, i); continue
        end
        if in_doc
            if startswith(rest, "\"\"\"")
                in_doc = false; i += 3; continue
            end
            src[i] == '\n' ? write(out, '\n') : write(out, ' ')
            i = nextind(src, i); continue
        end
        if startswith(rest, "#=")
            depth = 1; i += 2; continue
        elseif startswith(rest, "\"\"\"")
            in_doc = true; i += 3; continue
        elseif src[i] == '#'
            while i <= n && src[i] != '\n'
                i = nextind(src, i)
            end
            continue
        end
        write(out, src[i])
        i = nextind(src, i)
    end
    return String(take!(out))
end

"""
    exists_case_exact(root, relpath) -> Bool

Whether `relpath` exists under `root` with exactly that spelling. Every component is
checked against `readdir` rather than against `isfile`, because on Windows and on a
default macOS volume `isfile` answers yes for the wrong case — which is the entire bug
this suite exists to catch.
"""
function exists_case_exact(root::AbstractString, relpath::AbstractString)
    parts = filter(!isempty, split(replace(relpath, '\\' => '/'), '/'))
    here = root
    for (k, part) in enumerate(parts)
        part == "." && continue
        if part == ".."
            here = dirname(here); continue
        end
        isdir(here) || return false
        part in readdir(here) || return false
        here = joinpath(here, part)
        k < length(parts) && !isdir(here) && return false
    end
    return ispath(here)
end

# `include("name.jl")` and `include(joinpath(...))`, optionally with a module as the first
# argument. Built from strings so this file does not match its own patterns.
const RE_PLAIN = Regex("include\\(\\s*(?:[A-Za-z_][A-Za-z0-9_]*\\s*,\\s*)?\"([^\"\\n]+)\"\\s*\\)")
const RE_JOIN  = Regex("include\\(\\s*(?:[A-Za-z_][A-Za-z0-9_]*\\s*,\\s*)?joinpath\\(([^()\\n]*)\\)\\s*\\)")

"""
    include_targets(file) -> Vector{Tuple{String,String}}

Every path `file` includes, as `(as written, resolved absolute path)`. `@__DIR__` resolves
to the including file's directory and `REPO` to the repository root, which is how every
script in this repository spells it. A `joinpath` built from anything else is returned
with an empty resolution and reported rather than silently skipped.
"""
function include_targets(file::AbstractString)
    src = strip_julia_comments(read(file, String))
    here = dirname(file)
    found = Tuple{String,String}[]
    for m in eachmatch(RE_PLAIN, src)
        push!(found, (m.captures[1], normpath(joinpath(here, m.captures[1]))))
    end
    for m in eachmatch(RE_JOIN, src)
        raw = strip(m.captures[1])
        parts = strip.(split(raw, ','))
        resolved = String[]
        ok = true
        for (k, p) in enumerate(parts)
            if p == "@__DIR__"
                k == 1 || (ok = false)
                push!(resolved, here)
            elseif p == "REPO"
                k == 1 || (ok = false)
                push!(resolved, LAYOUT_REPO)
            elseif startswith(p, '"') && endswith(p, '"') && length(p) >= 2
                push!(resolved, p[2:end-1])
            else
                ok = false
            end
        end
        push!(found, ("joinpath(" * raw * ")",
                      ok ? normpath(joinpath(resolved...)) : ""))
    end
    return found
end

"""
    parse_error(file) -> Union{Nothing,String}

A one-line summary of the first syntax error in `file`, or `nothing` when it parses.
`Meta.parseall` reports a syntax error by putting an `:error` or an `:incomplete` node
into the expression it returns rather than by raising, so that expression is walked for
one. Only the text is read: nothing in `file` is evaluated and nothing it includes is
opened, so a file that resolves inputs or reads a data root at top level is as safe to
check as a module.
"""
function parse_error(file::AbstractString)
    stack = Any[Meta.parseall(read(file, String); filename = relpath(file, LAYOUT_REPO))]
    while !isempty(stack)
        e = pop!(stack)
        e isa Expr || continue
        if e.head === :error || e.head === :incomplete
            detail = isempty(e.args) ? "syntax error" : sprint(showerror, e.args[1])
            return first(join(split(detail), " "), 120)   # the report is one line
        end
        append!(stack, e.args)
    end
    return nothing
end

# Every directory the repository ships, for the case-collision scan. `include` targets only
# ever live in LAYOUT_DIRS, but a case collision anywhere breaks a checkout: LICENSE against
# license, README.md against Readme.md, dat/models against dat/Models.
const CASE_SCAN_DIRS = (".", "src", "scripts", "test", "docs", "dat",
                        joinpath("dat", "models"), joinpath("dat", "control_set"),
                        joinpath("dat", "criteria"),
                        joinpath(".github", "workflows"))

# Shipped paths that prose, the ignore rules or a document spell out by name. No `include`
# guards any of these, so a rename that misses one is invisible until a reader hits it.
const PROSE_PATHS = ("dat/models/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm",
                     "dat/models/PROVENANCE.md",
                     "dat/control_set/known_cofactor_mdhs.tsv",
                     "dat/control_set/control_validation_predictions.csv",
                     "dat/criteria/anchor_targets.csv",
                     "dat/criteria/impl_4.csv",
                     "dat/criteria/impl_w1.csv",
                     "docs/TOOLING.md",
                     "docs/ANALYSIS.md",
                     "docs/REPRODUCING.md",
                     "docs/DATA.md",
                     "docs/data_manifest.sha256",
                     "scripts/verify_data.jl",
                     "README.md", "LICENSE", "CITATION.cff",
                     ".gitignore", ".gitattributes",
                     "Project.toml", "Manifest.toml")

# Prose and CI files that carry commands a reader copies and runs. docs/DATA.md belongs
# here for the same reason README.md does, and more so: it is the document a reader builds
# a data root from, so every command in it is typed by somebody who has not run anything
# from this repository yet.
const PROSE_FILES = let acc = String[]
    for rel in ("README.md", joinpath("docs", "TOOLING.md"),
                joinpath("docs", "ANALYSIS.md"),
                joinpath("docs", "REPRODUCING.md"),
                joinpath("docs", "DATA.md"),
                joinpath("dat", "models", "PROVENANCE.md"),
                joinpath(".github", "workflows", "ci.yml"))
        f = joinpath(LAYOUT_REPO, rel)
        isfile(f) && push!(acc, f)
    end
    acc
end

const LAYOUT_FILES = let acc = String[]
    for d in LAYOUT_DIRS
        dir = joinpath(LAYOUT_REPO, d)
        isdir(dir) || continue
        for f in sort(readdir(dir))
            endswith(f, ".jl") && push!(acc, joinpath(dir, f))
        end
    end
    acc
end

@testset "source layout" begin

    @testset "the tree is not empty" begin
        # A scan that found nothing would pass every check below for the wrong reason.
        @test length(LAYOUT_FILES) > 20
    end

    @testset "every include resolves, case-exactly" begin
        for file in LAYOUT_FILES
            rel_file = relpath(file, LAYOUT_REPO)
            for (written, resolved) in include_targets(file)
                if isempty(resolved)
                    @test ("$(rel_file): cannot resolve $(written)", false) == ("", true)
                    continue
                end
                rel = relpath(resolved, LAYOUT_REPO)
                hit = exists_case_exact(LAYOUT_REPO, rel)
                hit || @error "include target missing or misspelled" file = rel_file target = written resolves_to = rel
                @test hit
            end
        end
    end

    @testset "the data files the code names are present, case-exactly" begin
        # Named by src/scriptLn.jl, by the classifier's control-set workflow and, for the
        # anchor targets and their two sources, by the registry in src/ModuleDataInputs.jl.
        # All five are inputs this repository ships. test/test_data_inputs.jl asks the
        # filesystem for the registered ones, and on a case-insensitive volume the filesystem
        # answers yes for the wrong case, so their spelling is checked here.
        for rel in ("dat/models/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm",
                    "dat/control_set/known_cofactor_mdhs.tsv",
                    "dat/criteria/anchor_targets.csv",
                    "dat/criteria/impl_4.csv",
                    "dat/criteria/impl_w1.csv")
            @test exists_case_exact(LAYOUT_REPO, rel)
        end
    end

    @testset "no two files differ only by case" begin
        # Such a pair cannot both be checked out on Windows or on a default macOS volume:
        # one silently overwrites the other.
        # The hazard is not confined to the directories that hold `include` targets, so
        # every directory this repository ships is scanned.
        for d in CASE_SCAN_DIRS
            dir = d == "." ? LAYOUT_REPO : joinpath(LAYOUT_REPO, d)
            isdir(dir) || continue
            names = filter(n -> n != ".git", readdir(dir))
            folded = lowercase.(names)
            length(unique(folded)) == length(folded) ||
                @error "two names differ only by case" directory = d names = names
            @test length(unique(folded)) == length(folded)
        end
    end

    @testset "the paths named in prose are the paths on disk" begin
        # A shipped file can be renamed out from under a document or an ignore rule, where
        # no `include` is watching. Every shipped path that prose or .gitignore spells out
        # is checked here against the directory listing, byte for byte.
        for rel in PROSE_PATHS
            hit = exists_case_exact(LAYOUT_REPO, rel)
            hit || @error "a path named in prose is missing or misspelled" path = rel
            @test hit
        end
        # ...and each one is still named somewhere, so a stale entry in the list above is
        # caught rather than passing forever on a file nothing refers to any more.
        prose = join([isfile(joinpath(LAYOUT_REPO, f)) ?
                      read(joinpath(LAYOUT_REPO, f), String) : ""
                      for f in ("README.md", ".gitignore",
                                joinpath("docs", "TOOLING.md"),
                                joinpath("docs", "DATA.md"),
                                joinpath("dat", "models", "PROVENANCE.md"))], "\n")
        for rel in PROSE_PATHS
            occursin(basename(rel), prose) ||
                @error "a path in PROSE_PATHS is named in no prose file any more" path = rel
            @test occursin(basename(rel), prose)
        end
    end

    @testset "no invocation depends on a version manager" begin
        # A `+channel` argument in front of `--project` is understood by the juliaup
        # shim, not by Julia itself. A reader who installed Julia any other way gets
        # `SystemError: opening file ".../+channel"` from the very first command they
        # copy out of this repository, so no invocation here may carry one.
        # The scan covers prose and CI as well as code: a reader copies a command out of
        # README.md or docs/TOOLING.md far more often than out of a source file, and a
        # workflow that needs a shim fails on a runner that has none.
        offenders = String[]
        for file in vcat(LAYOUT_FILES, PROSE_FILES)
            for (i, line) in enumerate(eachline(file))
                occursin(Regex("julia \\+[A-Za-z0-9._-]+"), line) &&
                    push!(offenders, "$(relpath(file, LAYOUT_REPO)):$(i)")
            end
        end
        isempty(offenders) || @error "version-manager-only invocation" offenders
        @test isempty(offenders)
    end

    @testset "every source file parses" begin
        # A syntax error in a script is invisible to every other check here: `using LnFinder`
        # loads only what src/LnFinder.jl includes, the scans above read these files as text,
        # and src/scriptLn.jl resolves its inputs against a data root at top level, so no test
        # loads it. A reader meets the scripts the moment they try the pipeline, and they are
        # long, so every .jl file the repository ships is parsed. Parsing executes none of it.
        for file in LAYOUT_FILES
            rel_file = relpath(file, LAYOUT_REPO)
            err = parse_error(file)
            err === nothing || @error "source file does not parse" file = rel_file detail = err
            @test (rel_file, err) == (rel_file, nothing)
        end
    end
end
