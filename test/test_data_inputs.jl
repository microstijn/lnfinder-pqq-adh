#=
    Offline self-test for src/ModuleDataInputs.jl — the one place this repository says
    where its inputs live. No network and nothing written outside a temporary directory:
    every root below is an empty `mktempdir` this suite fills itself. The only repository
    data files it reads are the shipped anchor targets and the two sources they name (6).

    THE FAILURES THIS SUITE EXISTS TO PREVENT.

    1. A CONTRACT NOBODY CAN SATISFY. The whole point of one data root is that a reader
       sets one variable, unpacks one archive and has answered every path at once. That is
       only true if EVERY declared input really does resolve under the root. The first
       testset below builds a root from the registry itself — one file or directory per
       declared entry — and then requires every one of them, so an entry whose layout
       nobody can satisfy, or that quietly resolves somewhere else, fails here.

    2. A MISSING FILE REPORTED AS A STACK TRACE. Without a registry, a missing input
       surfaces as a `SystemError` naming an absolute path on the machine that ran it, or
       as a bare CSV parse error: a path no reader can have, produced by a step that never
       said what it wanted. A refusal has to name the input, say where the layout puts it
       and say which variable moves it — and it has to be a type of its own, so that a
       caller can tell "you have not got this file" apart from "this file is corrupt".

    3. A PATH PRINTED OUT OF SOMEBODY'S HOME DIRECTORY. The refusal quotes the layout,
       never the resolved absolute path, so a message copied into an issue or a paper
       carries no drive letter and no home directory. The message is scanned for both.

    4. A VARIABLE THAT EXISTS ONLY IN A SCRIPT. Every LNFINDER_* name this repository
       understands is declared in the registry; one invented in a script and documented
       nowhere is a second contract. src/ and scripts/ are scanned for the pattern and
       every hit has to be a declared name.

    5. DOCUMENTATION THAT HAS DRIFTED OFF THE CODE. docs/DATA.md is what a reader builds
       their data root from. Every declared input and every declared layout path must
       appear there, or the document describes a layout the code does not use.

    6. A SOURCE CHECK THAT SILENTLY NEVER RUNS. The shipped anchor targets name the two
       files their cells were transcribed from, by file name and sha256, and step 4 checks
       a source only when it finds one beside the targets, reporting a missing one without
       refusing. A clone that lost a source, or a checkout that rewrote a source's line
       endings, would leave that check not run or failing on a file nobody edited. Every
       declared source must be a shipped input, registered beside the targets file, whose
       bytes are the declared digest.

    Run:  julia --project=. test/test_data_inputs.jl
=#

using Test, SHA

# Into a module of its own: the aggregate runner includes every suite into one process, and
# a pipeline script that this registry also reaches brings its own copy of it into Main.
module DataInputsHost end
Base.include(DataInputsHost, joinpath(@__DIR__, "..", "src", "ModuleDataInputs.jl"))

const DI = DataInputsHost.DataInputs
const DI_REPO = normpath(joinpath(@__DIR__, ".."))

"A drive letter or a UNC root, in any of the forms an error message could carry."
const ABSOLUTE_PATH_RE = r"([A-Za-z]:[\\/])|(\\\\[A-Za-z0-9])|(/(home|Users|mnt)/)"

"""
    populate!(root) -> Vector{String}

Create, under `root`, exactly what the registry says a complete data root holds: one
empty file or one directory per declared entry that is not shipped with the code. Returns
the names it made something for.
"""
function populate!(root::AbstractString)
    env = Dict("LNFINDER_DATA" => root)
    made = String[]
    for i in DI.INPUTS
        i.origin === :repo && continue
        p = DI.input_path(i.name; env = env)
        if i.kind === :dir
            mkpath(p)
        else
            mkpath(dirname(p))
            write(p, "")
        end
        push!(made, i.name)
    end
    return made
end

@testset "data inputs" begin

    @testset "the registry is a set, not a list with duplicates" begin
        names = [i.name for i in DI.INPUTS]
        @test length(unique(names)) == length(names)
        layouts = [DI.layout(i) for i in DI.INPUTS]
        @test length(unique(layouts)) == length(layouts)
        @test length(DI.INPUTS) > 20                 # a registry that lost its entries
        # Every entry says what it is, who reads it and what shape it has.
        for i in DI.INPUTS
            @test !isempty(i.what)
            @test !isempty(i.read_by)
            @test i.kind in (:file, :dir)
            @test i.origin in (:repo, :deposit, :regenerate, :refetch, :run,
                               :notdistributed, :unreleased)
        end
    end

    @testset "every documented input resolves under a temporary root" begin
        mktempdir() do root
            env = Dict("LNFINDER_DATA" => root)
            made = populate!(root)
            @test !isempty(made)
            for i in DI.INPUTS
                p = DI.require_input(i.name; env = env)       # throws if it does not resolve
                if i.origin === :repo
                    # The shipped inputs travel with the CODE and must not move when the
                    # data root does, or a clone would stop being able to find its own files.
                    @test startswith(normpath(p), normpath(DI_REPO))
                    @test !startswith(normpath(p), normpath(root))
                else
                    @test startswith(normpath(p), normpath(root))
                    # ...and at exactly the place the layout advertises.
                    @test normpath(p) == normpath(joinpath(root, i.segments...))
                end
            end
        end
    end

    @testset "the root moves everything at once" begin
        mktempdir() do a
            mktempdir() do b
                pa = DI.input_path("census_habitat"; env = Dict("LNFINDER_DATA" => a))
                pb = DI.input_path("census_habitat"; env = Dict("LNFINDER_DATA" => b))
                @test startswith(normpath(pa), normpath(a))
                @test startswith(normpath(pb), normpath(b))
                @test pa != pb
            end
        end
        # Unset, the root is dat/ in the clone — stated, so that a reader who unpacks the
        # data set in place needs no variable at all.
        @test normpath(DI.data_root(; env = Dict{String,String}())) ==
              normpath(joinpath(DI_REPO, "dat"))
        # An empty value is an unset value; joined literally, a variable exported as ""
        # gives a path that begins with a separator and resolves at the filesystem root.
        @test normpath(DI.data_root(; env = Dict("LNFINDER_DATA" => ""))) ==
              normpath(joinpath(DI_REPO, "dat"))
    end

    @testset "the two per-input overrides win over the root, and only for their own input" begin
        mktempdir() do root
            mktempdir() do elsewhere
                for i in DI.INPUTS
                    isempty(i.env) && continue
                    env = Dict("LNFINDER_DATA" => root, i.env => joinpath(elsewhere, i.name))
                    @test DI.input_path(i.name; env = env) == joinpath(elsewhere, i.name)
                    # every other input is unmoved
                    for j in DI.INPUTS
                        j.name == i.name && continue
                        j.origin === :repo && continue
                        @test startswith(normpath(DI.input_path(j.name; env = env)),
                                         normpath(root))
                    end
                end
            end
        end
    end

    @testset "a missing input is a named refusal, not an ENOENT" begin
        mktempdir() do root
            env = Dict("LNFINDER_DATA" => root)              # an EMPTY root
            for name in ("census_habitat", "relation_cache", "hits_fasta", "taxdump")
                err = nothing
                try
                    DI.require_input(name; env = env)
                catch e
                    err = e
                end
                @test err isa DI.MissingInput
                @test !(err isa SystemError)
                msg = sprint(showerror, err)
                i = DI.input(name)
                @test occursin(name, msg)                    # which input
                @test occursin(DI.layout(name), msg)         # where the layout puts it
                @test occursin("LNFINDER_DATA", msg)         # what moves it
                @test occursin("docs/DATA.md", msg)          # where the rest is written down
                @test occursin(i.what, msg)                  # what it is, in words
                # ...and nothing off this machine. A refusal is copied into issues and
                # papers; it must carry no drive letter and no home directory.
                @test !occursin(ABSOLUTE_PATH_RE, msg)
                @test !occursin(root, msg)
            end
        end
    end

    @testset "a refusal through an override says the override is the thing that is wrong" begin
        mktempdir() do root
            env = Dict("LNFINDER_DATA" => root,
                       "LNFINDER_AFA" => joinpath(root, "nowhere", "absent.afa"))
            err = try DI.require_input("model_alignment"; env = env) catch e; e end
            @test err isa DI.MissingInput
            msg = sprint(showerror, err)
            @test occursin("LNFINDER_AFA", msg)
            @test !occursin(ABSOLUTE_PATH_RE, msg)
        end
    end

    @testset "an input that is there is returned, and a directory is not a file" begin
        mktempdir() do root
            env = Dict("LNFINDER_DATA" => root)
            populate!(root)
            @test isfile(DI.require_input("census_habitat"; env = env))
            @test isdir(DI.require_input("relation_cache"; env = env))
            # A file where a directory belongs is a refusal, not a silent pass: the cache is
            # read with readdir, and "it exists" is not the property that matters.
            rm(DI.input_path("relation_cache"; env = env); recursive = true)
            write(DI.input_path("relation_cache"; env = env), "")
            @test_throws DI.MissingInput DI.require_input("relation_cache"; env = env)
        end
    end

    @testset "a name nobody declared is refused, with the list of names there are" begin
        err = try DI.input("no_such_input") catch e; e end
        @test err isa ErrorException
        @test occursin("no_such_input", err.msg)
        @test occursin("census_habitat", err.msg)            # the list is really there
    end

    @testset "every LNFINDER_ variable in the tree is a declared one" begin
        declared = Set{String}(collect(DI.ALL_VARS))
        seen = Dict{String,Vector{String}}()
        for dir in ("src", "scripts")
            d = joinpath(DI_REPO, dir)
            isdir(d) || continue
            for f in readdir(d)
                endswith(f, ".jl") || continue
                for m in eachmatch(r"LNFINDER_[A-Z_]+", read(joinpath(d, f), String))
                    push!(get!(seen, m.match, String[]), joinpath(dir, f))
                end
            end
        end
        @test !isempty(seen)                                  # the scan found something
        for (v, files) in seen
            v in declared || @error "an LNFINDER_ variable nobody declared" variable = v files = files
            @test v in declared
        end
        # ...and the converse: a declared variable nothing reads is a contract with no
        # implementation. The two switches are not LNFINDER_ names and are checked by the
        # testset below instead, so they are skipped here.
        switches = Set{String}(v for (v, _) in DI.SWITCH_VARS)
        for v in declared
            v in switches && continue
            haskey(seen, v) || @error "a declared variable nothing reads" variable = v
            @test haskey(seen, v)
        end
    end

    @testset "every environment variable read by name is a declared one" begin
        # THE PREFIX SCAN ABOVE CANNOT SEE A NAME OUTSIDE THE PREFIX. Two switches the
        # census driver reads — RUN_HABITAT and EXPECTED_CENTROIDS — carry no LNFINDER_
        # prefix, so a scan for that prefix alone leaves the claim "these are the whole
        # set" checked by nothing. This scan is for the ACCESS rather than for the name:
        # every `ENV["X"]`, `get(ENV, "X", …)` and `haskey(ENV, "X")` in src/ and
        # scripts/, whatever X is called.
        declared = Set{String}(collect(DI.ALL_VARS))
        pattern = r"(?:ENV|env)\s*(?:\[|,\s*)\"([A-Z][A-Z0-9_]*)\""
        seen = Dict{String,Vector{String}}()
        for dir in ("src", "scripts")
            d = joinpath(DI_REPO, dir)
            isdir(d) || continue
            for f in readdir(d)
                endswith(f, ".jl") || continue
                for m in eachmatch(pattern, read(joinpath(d, f), String))
                    push!(get!(seen, m.captures[1], String[]), joinpath(dir, f))
                end
            end
        end
        @test !isempty(seen)
        for (v, files) in seen
            v in declared ||
                @error "an environment variable nobody declared" variable = v files = files
            @test v in declared
        end
        # Both switches are really read, so the declaration is not a contract with no
        # implementation either.
        for (v, _) in DI.SWITCH_VARS
            haskey(seen, v) || @error "a declared switch nothing reads" variable = v
            @test haskey(seen, v)
        end
    end

    @testset "no second spelling of a shipped input's path" begin
        # The classifier keeps its own fallback path to the profile HMM rather than taking
        # a second copy of the registry into a namespace two suites already share. That is
        # a deliberate exception, and this is what keeps it honest: the two spellings must
        # resolve to the same file, so a rename in the registry cannot leave the classifier
        # pointing at the old name.
        src = read(joinpath(DI_REPO, "src", "validate_classifier.jl"), String)
        m = match(r"_DEFAULT_HMM\s*=\s*normpath\(joinpath\(@__DIR__,\s*([^\)]*)\)\)", src)
        @test m !== nothing
        parts = [strip(p, ['"', ' ', '\n', '\r']) for p in split(m.captures[1], ",")]
        parts = filter(p -> !isempty(p) && p != "\"", parts)
        spelled = normpath(joinpath(joinpath(DI_REPO, "src"), parts...))
        @test normpath(spelled) == normpath(DI.input_path("profile_hmm"))
    end

    @testset "the shipped anchor targets and the sources they name are the same bytes" begin
        # A 'source' row is `source,<file>,...,<sha256>`; the note between them may carry
        # commas, so the file name is the second field and the digest the last.
        targets = DI.input_path("anchor_targets")
        declared = Dict{String,String}()
        source_rows = 0
        for line in eachline(targets)
            startswith(line, "source,") && (source_rows += 1)
            m = match(r"^source,([^,]+),.*,([0-9a-f]{64})\s*$", line)
            m === nothing || (declared[m.captures[1]] = m.captures[2])
        end
        @test length(declared) >= 1
        # Every 'source' row was read, so none escapes the check below by a shape this
        # pattern does not match.
        @test length(declared) == source_rows
        shipped = Dict(basename(DI.layout(i)) => i.name
                       for i in DI.INPUTS if i.origin === :repo)
        for (file, digest) in declared
            haskey(shipped, file) ||
                @error "a source the anchor targets name is not a shipped input" source = file
            @test haskey(shipped, file)
            haskey(shipped, file) || continue
            # Step 4 looks for a source BESIDE the targets file and nowhere else, so a source
            # registered in any other directory would be shipped and still never checked.
            @test dirname(DI.input_path(shipped[file])) == dirname(targets)
            got = open(io -> bytes2hex(sha256(io)), DI.input_path(shipped[file]))
            got == digest ||
                @error("a shipped source is not the bytes the anchor targets name",
                       source = file, declared = digest, found = got)
            @test got == digest
        end
    end

    @testset "docs/DATA.md documents every declared input" begin
        doc = joinpath(DI_REPO, "docs", "DATA.md")
        @test isfile(doc)
        text = read(doc, String)
        for i in DI.INPUTS
            occursin(i.name, text) || @error "an input the data document never names" input = i.name
            @test occursin(i.name, text)
            occursin(DI.layout(i), text) ||
                @error "a layout path the data document never names" path = DI.layout(i)
            @test occursin(DI.layout(i), text)
        end
    end
end
