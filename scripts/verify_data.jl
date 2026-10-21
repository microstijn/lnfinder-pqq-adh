#=
    verify_data.jl — is the data root the one the analysis was run against?

    WHAT IT IS FOR. The analyses in this repository read files this repository does not
    carry. Those files are archived separately, and an archive can arrive truncated, half
    unpacked, line-ending-mangled or simply out of date. This checks a data root against a
    checksum manifest and says, per file, whether it is there and whether it is the same
    bytes — and exits non-zero if anything is missing or different, so it can stand in
    front of a long run in a script.

    WHY JULIA RATHER THAN A SHELL SCRIPT. Julia is the only prerequisite this repository
    declares, so a reader who can run the analysis can run this. A shell version would
    need `sha256sum` on Linux, `shasum -a 256` on macOS and PowerShell's `Get-FileHash` on
    Windows — three variants to keep in step — and a Python version would add a runtime
    nothing else here needs. This file uses only the standard library (`SHA`, `Printf`), so
    it also runs BEFORE `Pkg.instantiate`, which is exactly when a reader is setting up.

    The manifest is plain sha256sum format, so `sha256sum -c ../docs/data_manifest.sha256`
    from inside the data root is a second, independent way to run the same check.

    Run:  julia scripts/verify_data.jl [options]
          (--help lists them and exits; nothing is ever written, read-only throughout)
=#

const USAGE = """
Usage: julia scripts/verify_data.jl [options]

Checks a data root against a checksum manifest: every file the manifest lists is looked
for, hashed and compared. Nothing is written and nothing is fetched.

Options:
  --root <path>      the data root to check
                     (default: \$LNFINDER_DATA, else dat/ in this clone)
  --manifest <path>  the checksum manifest to check against
                     (default: MANIFEST.sha256 in the data root if there is one,
                      else docs/data_manifest.sha256 in this clone)
  --quiet            print only the summary and whatever failed
  -h, --help         print this text and exit

Exit status: 0 when every listed file is present and matches; 1 when a required file is
missing, any digest differs, or the two records disagree; 2 on a command-line error.

A manifest has two kinds of entry. REQUIRED entries are the archived files: a missing one
is a failure. OPTIONAL entries are files the archive does not carry because a reader
regenerates or re-fetches them from a public source; a missing one is reported and is not
a failure, a mismatching one is. docs/DATA.md says which is which and why.

TWO RECORDS, NOT ONE. An archive carries its own MANIFEST.sha256, and that copy is the
complete one, so it is what a run is checked against when it is there. On its own it
would only establish that the archive agrees with itself, so whenever some other record
is used the repository's docs/data_manifest.sha256 is ALSO consulted, for the paths it
covers, and a disagreement between the two is a failure in its own right.
"""

# ---- command line, before anything is loaded or opened ----------------------------
function parse_args(args)
    root = ""
    manifest = ""
    quiet = false
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-h", "--help")
            print(USAGE); exit(0)
        elseif a == "--quiet"
            quiet = true; i += 1
        elseif a in ("--root", "--manifest")
            i + 1 <= length(args) ||
                (println(stderr, "verify_data.jl: $(a) needs a path\n"); print(stderr, USAGE); exit(2))
            a == "--root" ? (root = args[i+1]) : (manifest = args[i+1])
            i += 2
        else
            println(stderr, "verify_data.jl: unknown option '$(a)'\n")
            print(stderr, USAGE)
            exit(2)
        end
    end
    return (root = root, manifest = manifest, quiet = quiet)
end

const OPT = abspath(PROGRAM_FILE) == (@__FILE__) ? parse_args(ARGS) :
            (root = "", manifest = "", quiet = false)

using SHA, Printf

const REPO = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "src", "ModuleDataInputs.jl"))

# ---- the manifest ------------------------------------------------------------------

"""
Section markers. A manifest's entries are required or optional, and which one they are is
carried in a comment line, so that the file is still plain sha256sum format. The markers
are EXACT, and a comment that opens with `# ---` and is not one of them is an error rather
than a line to guess at: the alternative is deciding from a word that happens to appear in
a sentence, and a manifest written by the archive could then demote a required file to
optional with nothing but prose.
"""
const SECTION_REQUIRED = "# --- deposit: every one of these must be present in a complete data root ---"
const SECTION_OPTIONAL = "# --- optional: checked only if present"

"""
    read_manifest(path) -> Vector{NamedTuple}

`(digest, path, required)` per entry. Plain sha256sum format — `<64 hex>  <path>` — with
`#` comment lines, which is what makes the same file readable by `sha256sum -c`. One of the
two section markers above switches what follows between required and optional entries; a
manifest with no marker at all is read as all-required, which is the safe reading.
"""
function read_manifest(path::AbstractString)
    entries = NamedTuple[]
    required = true
    for (n, raw) in enumerate(eachline(path))
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, "#")
            if startswith(line, SECTION_OPTIONAL)
                required = false
            elseif startswith(line, SECTION_REQUIRED)
                required = true
            elseif startswith(line, "# ---")
                error("$(basename(path)) line $(n): a section comment this tool does not " *
                      "know. A manifest says which of its files are required by one of two " *
                      "exact markers; guessing from the words in a comment is how a " *
                      "required file becomes an optional one by accident. The markers are:" *
                      "\n  " * SECTION_REQUIRED * "\n  " * SECTION_OPTIONAL * " ---…")
            end
            continue
        end
        m = match(r"^([0-9a-fA-F]{64})\s+\*?(.+)$", line)
        m === nothing &&
            error("$(basename(path)) line $(n): not a sha256sum entry and not a comment")
        push!(entries, (digest = lowercase(m.captures[1]),
                        path = String(strip(m.captures[2])),
                        required = required))
    end
    isempty(entries) && error("$(basename(path)) lists no files")
    return entries
end

"Hex sha256 of a file, streamed, so a multi-hundred-megabyte input is not read into memory."
file_digest(p::AbstractString) = open(p) do io
    bytes2hex(sha256(io))
end

human(bytes::Integer) = bytes >= 1 << 30 ? @sprintf("%.2f GB", bytes / (1 << 30)) :
                        bytes >= 1 << 20 ? @sprintf("%.1f MB", bytes / (1 << 20)) :
                        @sprintf("%d B", bytes)

# ---- the run -----------------------------------------------------------------------

function main(opt)
    root = isempty(opt.root) ? DataInputs.data_root() : abspath(opt.root)
    isdir(root) || begin
        println(stderr, "verify_data.jl: the data root does not exist: $(root)")
        println(stderr, "Pass --root <path>, or set \$", DataInputs.ROOT_VAR,
                        ". docs/DATA.md says what belongs under it.")
        return 1
    end

    repo_record = joinpath(REPO, "docs", "data_manifest.sha256")
    archive_record = joinpath(root, "MANIFEST.sha256")

    manifest = if !isempty(opt.manifest)
        abspath(opt.manifest)
    elseif isfile(archive_record)
        archive_record                        # the archive's own, and the more complete one
    else
        repo_record
    end
    isfile(manifest) || begin
        println(stderr, "verify_data.jl: no manifest at $(manifest)")
        return 1
    end

    # THE ARCHIVE'S OWN RECORD CANNOT BE THE ONLY ONE CONSULTED. It travels with the files it
    # describes, so on its own it says no more than "this archive agrees with itself": bytes
    # and record could have moved together and every line would still read ok. The
    # repository's copy is the second, independent location, which is the whole reason it is
    # shipped — so whenever the run is against some other record, this one is cross-checked
    # too, for the paths it covers, and a disagreement between the two is a failure.
    cross = (!isempty(repo_record) && manifest != repo_record &&
             isfile(repo_record)) ? repo_record : ""

    # A manifest that cannot be read is a sentence, not a stack trace: it is a file the reader
    # was handed, and the message has to say what is wrong with it.
    entries = try
        read_manifest(manifest)
    catch e
        e isa ErrorException || rethrow()
        println(stderr, "verify_data.jl: ", e.msg)
        return 1
    end
    println("data root: ", root)
    println("manifest:  ", manifest, "  (", length(entries), " files listed)")
    println(isempty(cross) ?
            "checked against this repository's own record." :
            "cross-checked against this repository's own record: " * repo_record)
    println()

    ok, bad, absent_required, absent_optional, bytes = 0, String[], String[], String[], 0
    seen = Dict{String,String}()          # layout path => the digest measured in this run
    for e in entries
        p = joinpath(root, split(e.path, '/')...)
        if !isfile(p)
            e.required ? push!(absent_required, e.path) : push!(absent_optional, e.path)
            opt.quiet || @printf("  %-9s %s\n", e.required ? "MISSING" : "absent", e.path)
            continue
        end
        d = file_digest(p)
        seen[e.path] = d
        bytes += filesize(p)
        if d == e.digest
            ok += 1
            opt.quiet || @printf("  %-9s %s\n", "ok", e.path)
        else
            push!(bad, e.path)
            @printf("  %-9s %s\n", "DIFFERS", e.path)
            @printf("            expected %s\n            found    %s\n", e.digest, d)
        end
    end

    # The second record. Files it lists that the first one also listed are compared against
    # the digest already measured, so nothing is hashed twice; files only it lists are
    # hashed here. A file it lists that is absent is NOT reported again — the first record
    # has already said so, or does not claim it — because the question this section answers
    # is only whether the two records agree about the same bytes.
    disagree, crosschecked = String[], 0
    if !isempty(cross)
        cross_entries = try
            read_manifest(cross)
        catch e
            e isa ErrorException || rethrow()
            println(stderr, "verify_data.jl: ", e.msg)
            return 1
        end
        for e in cross_entries
            p = joinpath(root, split(e.path, '/')...)
            d = get(seen, e.path, "")
            if isempty(d)
                isfile(p) || continue
                d = file_digest(p)
            end
            crosschecked += 1
            d == e.digest && continue
            push!(disagree, e.path)
            @printf("  %-9s %s\n", "DISAGREE", e.path)
            @printf("            this repository records %s\n            the data root holds    %s\n",
                    e.digest, d)
        end
    end

    # The manifest is one record; the code's own registry is another. A data root can pass
    # the digests and still be missing something a step reads, so the contract itself is
    # checked too — this is the list of archived inputs the pipeline will refuse on.
    # Inputs a reader has to obtain elsewhere are outside it; docs/DATA.md section 5
    # lists those, and each names itself when a step asks for it and it is not there.
    unresolved = String[]
    for i in DataInputs.deposit_inputs()
        env = Dict{String,String}(k => v for (k, v) in ENV)
        env[DataInputs.ROOT_VAR] = root
        p = DataInputs.input_path(i.name; env = env)
        (i.kind === :dir ? isdir(p) : isfile(p)) || push!(unresolved, i.name)
    end

    println()
    @printf("%d verified (%s), %d differing, %d required missing, %d optional absent\n",
            ok, human(bytes), length(bad), length(absent_required), length(absent_optional))
    if isempty(cross)
        println("The record used was this repository's own; there was no second one to ",
                "cross-check against.")
    else
        @printf("%d of those also checked against this repository's own record, %d disagreeing\n",
                crosschecked, length(disagree))
    end
    if !isempty(bad)
        println("\nDIFFERENT from the recorded bytes:")
        for p in bad; println("  ", p); end
    end
    if !isempty(disagree)
        println("\nThe two records DISAGREE about these files. One of them is not the record ",
                "of\nthe bytes in this data root, and which one cannot be decided from here:")
        for p in disagree; println("  ", p); end
    end
    if !isempty(absent_required)
        println("\nMISSING and required:")
        for p in absent_required; println("  ", p); end
        println("docs/DATA.md says what each of these is and where it comes from.")
    end
    if !isempty(absent_optional)
        println("\nAbsent, and not a failure — regenerate or re-fetch these only if you need",
                "\nthe stage that reads them (docs/DATA.md gives the command for each):")
        for p in absent_optional; println("  ", p); end
    end
    if !isempty(unresolved)
        println("\nDeclared inputs this data root does not satisfy:")
        for n in unresolved
            println("  ", n, "  ->  ", DataInputs.layout(n))
        end
    end

    failed = !isempty(bad) || !isempty(disagree) || !isempty(absent_required) ||
             !isempty(unresolved)
    println(failed ? "\nFAILED" : "\nOK — this data root matches the record.")
    return failed ? 1 : 0
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    exit(main(OPT))
end
