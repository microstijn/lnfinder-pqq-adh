#=
    ModuleDataInputs.jl — the one place this repository says where its inputs live.

    THE CONTRACT, IN ONE PARAGRAPH. Every file the analysis code reads that this
    repository does not carry lives under a single directory, the DATA ROOT. The root is
    the directory named by the environment variable LNFINDER_DATA, and when that is not
    set it is `dat/` inside the clone. Underneath the root the layout is fixed and is
    declared in INPUTS below, so a reader who unpacks the archived data set at the root
    has satisfied every path at once. Two inputs are not in that archive because a reader
    makes or fetches them from a public source — the alignment and the taxonomy dump —
    and those two, and only those two, also answer to a variable of their own, so that a
    large file already on the machine does not have to be copied under the root.

    WHY ONE ROOT. A pipeline of this size otherwise accumulates one absolute root per
    stage, each of them valid on one machine. Every input here is a path under
    LNFINDER_DATA, so there is one thing to set, one layout to document, and one place —
    this file — that a rename has to pass through.

    WHAT THIS FILE IS NOT. It is not a downloader and it opens no network connection. It
    resolves names to paths, and it refuses by name when the file a name points at is not
    there. `scripts/verify_data.jl` checks a data root against a checksum manifest;
    `docs/DATA.md` says what each entry is, how it was produced and how to obtain it.

    This module deliberately EXPORTS NOTHING. Several scripts include it, some of them
    into the same namespace as each other, and a name that arrived from two copies of one
    module would be ambiguous at the point of use. Everything is called qualified:
    `DataInputs.input_path("census_habitat")`.

    Run:  julia --project=. test/test_data_inputs.jl
=#

module DataInputs

# ---------------------------------------------------------------------------------
# The variables. Every LNFINDER_* name this repository understands is declared here
# and nowhere else; test/test_data_inputs.jl scans src/ and scripts/ and fails if a
# name appears there that is not in one of these three lists.
# ---------------------------------------------------------------------------------

"The one variable that moves the whole data root."
const ROOT_VAR = "LNFINDER_DATA"

"""
Per-input overrides. Kept for the two inputs that the archived data set does NOT
contain, because a reader regenerates one and re-fetches the other from a public
archive, and both are large enough that they are commonly kept outside a project tree.
Every other input resolves under the root and has no variable of its own: a second way
to name a file is a second thing to get wrong.
"""
const OVERRIDE_VARS = ("LNFINDER_AFA", "LNFINDER_TAXDUMP")

"""
Variables that name where a run WRITES, or that seed a cache. They are not inputs and
they have no place under the data root, so they are declared here and resolved by the
step that owns them.
"""
const OTHER_VARS = (
    ("LNFINDER_CENSUS_OUT", "the directory a census run writes its table and manifest into"),
    ("LNFINDER_FIGDIR",     "the figure directory of habitat steps 4 and 5"),
    ("LNFINDER_RELATION_SEED",
     "a semicolon-separated list of directories of already-fetched <taxid>.json bodies, " *
     "used to prime step 1's relation cache without the network"),
)

"""
Environment variables that are neither the root, an input override nor a write location:
two switches the census driver reads. They do NOT carry the `LNFINDER_` prefix, and they
are declared here for exactly that reason — a name outside the prefix is a name the
prefix scan cannot see, so it is listed instead. `test/test_data_inputs.jl` scans `src/`
and `scripts/` for every environment variable read by name and fails on one that appears
in none of these four lists.
"""
const SWITCH_VARS = (
    ("RUN_HABITAT",
     "false turns off the census driver's habitat layer, which is its only per-taxon " *
     "network route besides accession resolution"),
    ("EXPECTED_CENTROIDS",
     "a whole number the census driver checks its centroid count against, so a clustering " *
     "that silently changed size stops the run"),
)

"Every environment variable name this repository understands, in one set."
const ALL_VARS = Tuple(vcat([ROOT_VAR], collect(OVERRIDE_VARS),
                            [v for (v, _) in OTHER_VARS],
                            [v for (v, _) in SWITCH_VARS]))

# ---------------------------------------------------------------------------------
# The registry
# ---------------------------------------------------------------------------------

"""
    DataInput

One named input. `segments` is its path relative to the data root — or, when `origin`
is `:repo`, relative to the repository root, because the shipped inputs travel with the
code and never move when the data root does.

`origin` says where a reader gets it:

  `:repo`            shipped in this repository
  `:deposit`         in the archived data set
  `:regenerate`      rebuilt from a deposited file with a command in docs/DATA.md
  `:refetch`         downloaded from a public archive by identifier or release
  `:run`             written by an earlier step of the pipeline, not archived
  `:notdistributed`  obtainable only by re-querying a live service, which returns
                     different data

`env` is the per-input override variable, or `""` when the input has none.
"""
struct DataInput
    name::String
    segments::Tuple{Vararg{String}}
    kind::Symbol           # :file or :dir
    origin::Symbol
    env::String
    what::String
    read_by::String
end

const NO_STAGE = "no stage in this release"

"""
    INPUTS

Every input of every analysis entry point in this repository, and every file the
archived data set carries. The order is the order docs/DATA.md lists them in.
"""
const INPUTS = DataInput[

    # ---- shipped with the code -------------------------------------------------
    DataInput("profile_hmm",
              ("dat", "models", "methanol_dehydrogenase_pqq_xoxF_mxaF.hmm"), :file, :repo, "",
              "the profile HMM the cofactor call is calibrated against",
              "the census driver, habitat step 2 and control-set validation"),
    DataInput("control_set",
              ("dat", "control_set", "known_cofactor_mdhs.tsv"), :file, :repo, "",
              "the curated set of enzymes whose cofactor is known from the literature",
              "control-set validation"),
    # The anchor targets ship with the code, and so do the two files they were transcribed
    # from. Step 4 digests a declared source only when it finds it beside the targets file
    # and reports one it cannot find without refusing, so a source that did not travel with
    # the targets is a check that silently never runs.
    DataInput("anchor_targets", ("dat", "criteria", "anchor_targets.csv"), :file, :repo, "",
              "the transcribed anchor cells, the registered level names and the digests of " *
              "the two files the cells were transcribed from, which step 4 checks its grid " *
              "against",
              "habitat step 4, as --anchor-targets"),
    DataInput("anchor_source_primary", ("dat", "criteria", "impl_4.csv"), :file, :repo, "",
              "the independent implementation the registered specification's anchor cells " *
              "are transcribed from",
              "habitat step 4, which checks it against the digest anchor_targets.csv declares"),
    DataInput("anchor_source_variants", ("dat", "criteria", "impl_w1.csv"), :file, :repo, "",
              "the independent implementation the registered levels' anchor cells are " *
              "transcribed from",
              "habitat step 4, which checks it against the digest anchor_targets.csv declares"),

    # ---- the habitat route -----------------------------------------------------
    DataInput("census_habitat", ("habitat", "census_habitat.csv"), :file, :deposit, "",
              "the census table carrying the habitat layer, one row per sequence",
              "habitat steps 1, 2, 3 and 4"),
    DataInput("pool_taxa", ("habitat", "prokaryote_taxa.csv"), :file, :deposit, "",
              "the reference-proteome taxon pool the comparison arm is drawn from",
              "habitat steps 1 and 2, and step 3 for its provenance when present; " *
              "written by step 0"),
    DataInput("pool_background", ("habitat", "prokaryote_background.csv"), :file, :deposit, "",
              "the habitat annotation of that pool",
              "habitat step 1, and step 3 for its provenance when present; written by step 0"),
    DataInput("obt_path_cache", ("habitat", "obt_path_cache.csv"),
              :file, :notdistributed, "",
              "the habitat ontology climb, cached, so the ontology is walked once",
              "habitat steps 0 and 1"),
    DataInput("relation_cache", ("habitat", "omnicrobe_cache"), :dir, :notdistributed, "",
              "the directory the stored response bodies and their manifest live in; the archive " *
              "carries the manifest and not the bodies",
              "habitat step 1"),
    # The bodies are a separate entry from the directory that holds them, and deliberately.
    # A cache miss in step 1 is a live fetch, so "the cache directory exists" is not the
    # property that has to hold before an offline run: "the bodies are there" is. A cache
    # put in place one level wrong satisfies the parent and nothing else.
    DataInput("relation_bodies", ("habitat", "omnicrobe_cache", "relations"),
              :dir, :notdistributed, "",
              "the stored response bodies themselves, one <taxid>.json per taxon",
              "habitat step 1"),
    DataInput("relation_manifest", ("habitat", "omnicrobe_cache", "manifest.csv"),
              :file, :deposit, "",
              "what that cache holds, where each body came from and when",
              "habitat steps 1 and 3"),
    DataInput("habitat_calls", ("habitat", "habitat_taxon_calls.csv"), :file, :deposit, "",
              "one habitat call per taxon per rule variant",
              "habitat steps 3, 4 and 5; written by step 1"),
    DataInput("background_scan_dir", ("habitat", "background_scan"), :dir, :deposit, "",
              "the background scan's own records: its threshold, its recall and its acceptance rows",
              "habitat step 3; written by step 2"),
    DataInput("background_scan_table", ("habitat", "background_pqq_scan.csv"),
              :file, :deposit, "",
              "which taxa of the comparison arm carry the enzyme at the calibrated threshold",
              "habitat steps 3 and 4; written by step 2"),
    DataInput("scan_recall", ("habitat", "background_scan", "recall.txt"), :file, :deposit, "",
              "the measured recall of the scan's prefilter, with the counts it was measured on",
              "habitat step 3"),
    DataInput("scan_threshold", ("habitat", "background_scan", "threshold.txt"),
              :file, :deposit, "",
              "the calibrated score threshold and how it was set",
              "habitat step 3"),
    DataInput("scan_acceptance", ("habitat", "background_scan", "acceptance.tsv"),
              :file, :deposit, "",
              "the scan's acceptance rows, one per check it had to pass",
              "habitat step 3"),
    DataInput("uniprot_release", ("habitat", "background_scan", "uniprot_release.txt"),
              :file, :deposit, "",
              "the sequence-database release the comparison arm was built from",
              "habitat step 3"),
    DataInput("enrichment_table", ("habitat", "habitat_enrichment.csv"), :file, :run, "",
              "the estimates, one row per contrast",
              "habitat step 5; written by step 3"),
    DataInput("arm_summary", ("habitat", "habitat_arm_summary.csv"), :file, :run, "",
              "what each arm was made of, and every provenance key the panel echoes",
              "habitat step 5; written by step 3"),

    # ---- the criteria a run is held to, supplied on the command line ------------
    DataInput("category_selection", ("criteria", "category_selection.csv"), :file, :deposit, "",
              "which category each recorded outcome rule speaks for, which categories the " *
              "grid counts are exported for, and which phylum the standardisation shares " *
              "are quoted for, each with a note saying where that selection came from",
              "habitat steps 3 and 5, as --categories"),

    # ---- the census route ------------------------------------------------------
    DataInput("hits_fasta", ("census", "firstTry.fa"), :file, :deposit, "",
              "the unaligned sequences of the database search the census was built from",
              "the census driver"),
    DataInput("uc_clusters", ("census", "clusters.uc"), :file, :deposit, "",
              "the clustering of those sequences that defines the census centroids",
              "the census driver"),
    DataInput("uc_clusters_50", ("census", "clusters_50.uc"), :file, :deposit, "",
              "the coarser clustering used as the sub-family unit of the phylogenetic control",
              NO_STAGE),
    DataInput("census_table", ("census", "final_deduplicated_data.csv"), :file, :deposit, "",
              "the census as the driver in this repository writes it; the habitat input " *
              "carries its identifiers, calls and taxonomy row for row",
              NO_STAGE),
    DataInput("model_alignment", ("derived", "aligned_to_model.afa"), :file, :regenerate,
              "LNFINDER_AFA",
              "those sequences aligned to the profile, in model coordinates",
              "the census driver"),
    DataInput("taxdump", ("derived", "taxdump"), :dir, :refetch, "LNFINDER_TAXDUMP",
              "an unpacked NCBI taxonomy dump; the loader opens nodes, names, merged and delnodes",
              "the census driver"),

    # ---- the determinant and structure route -----------------------------------
    DataInput("pocket_labels", ("determinant", "pocket_labels_real.csv"), :file, :deposit, "",
              "the per-accession metal label set of the whole census",
              NO_STAGE),
    DataInput("structure_accessions", ("structures", "structure_accessions.txt"),
              :file, :deposit, "",
              "the identifiers of every predicted and experimental structure the analysis read, " *
              "so the structures themselves can be fetched from their own archives",
              NO_STAGE),
    DataInput("structures", ("derived", "structures"), :dir, :refetch, "",
              "the structures those identifiers name, once fetched",
              NO_STAGE),
]

const BY_NAME = Dict{String,DataInput}(i.name => i for i in INPUTS)

# A registry with two entries under one name would resolve whichever came first, and the
# duplicate would look declared. Checked when this file loads, not by a suite that may or
# may not run.
length(BY_NAME) == length(INPUTS) ||
    error("ModuleDataInputs: two entries share a name")

# ---------------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------------

"The repository this file ships in."
repo_root() = normpath(joinpath(@__DIR__, ".."))

"""
    default_root()

Where the data root is when `LNFINDER_DATA` is not set: `dat/` inside the clone. That
directory already exists, already holds the shipped inputs, and every format a run
writes under it is ignored by `.gitignore`, so a reader who simply unpacks the data set
in place cannot commit it by accident.
"""
default_root() = joinpath(repo_root(), "dat")

"""
    data_root(; env = ENV) -> String

The data root: `\$LNFINDER_DATA` when it is set to something, else `default_root()`.
Never fails, and never asks whether anything is there — path composition has to work
before a command line has been read, and a step that cannot even print its usage because
a variable is unset would be worse than a hard-coded absolute root.

**Always absolute.** A relative `\$LNFINDER_DATA` would otherwise mean a different
directory in every step that was started from a different working directory, and two
steps of one pipeline would disagree about where the data set is while both reported the
same variable. `default_root()` is absolute already, so this only ever moves a value a
reader typed relative.
"""
function data_root(; env = ENV)
    v = strip(get(env, ROOT_VAR, ""))
    return isempty(v) ? default_root() : abspath(String(v))
end

"""
    input(name) -> DataInput

The registry entry called `name`, or an error naming every entry there is. A typo in a
name would otherwise become a path nobody declared.
"""
function input(name::AbstractString)
    haskey(BY_NAME, name) && return BY_NAME[name]
    error("no input is called \"$(name)\". The declared inputs are: " *
          join(sort([i.name for i in INPUTS]), ", "))
end

"""
    layout(name) -> String

Where `name` sits, written the way docs/DATA.md writes it: relative to the data root,
with forward slashes, and never as an absolute path. This is what error messages quote,
so that a message says what the reader has to put where rather than echoing a path off
one machine.
"""
layout(name::AbstractString) = layout(input(name))
layout(i::DataInput) = join(i.segments, "/")

"""
    input_path(name; env = ENV) -> String

The absolute path of `name`, without asking whether anything is there. Resolution order:

  1. the input's own override variable, when it has one and it is set;
  2. the repository, for the inputs that ship with the code;
  3. the data root, for everything else.

Absolute in all three cases, for the reason `data_root` gives: a relative value in a
variable means a different file in every step started from a different directory.
"""
function input_path(name::AbstractString; env = ENV)
    i = input(name)
    if !isempty(i.env)
        v = strip(get(env, i.env, ""))
        isempty(v) || return abspath(String(v))
    end
    base = i.origin === :repo ? repo_root() : data_root(; env = env)
    return joinpath(base, i.segments...)
end

"""
    MissingInput

A required input is not there. Carried as a type of its own, rather than as whatever
`open` or the CSV reader would have raised, so that the failure says which named input
is missing, where that name expects the file, and which variable moves it — and says it
without quoting an absolute path off anybody's machine.
"""
struct MissingInput <: Exception
    name::String
    layout::String
    kind::Symbol
    origin::Symbol
    env::String
    what::String
    read_by::String
    via_override::Bool
end

function Base.showerror(io::IO, e::MissingInput)
    thing = e.kind === :dir ? "directory" : "file"
    print(io, "LnFinder: the required input \"", e.name, "\" is not there.")
    print(io, "\n  what it is:  ", e.what)
    print(io, "\n  read by:     ", e.read_by)
    if e.via_override
        print(io, "\n  looked for:  the ", thing, " named by \$", e.env, ", which is set")
    elseif e.origin === :repo
        print(io, "\n  looked for:  ", e.layout, ", in this repository — which should carry it, ",
                  "so this checkout is incomplete")
    else
        print(io, "\n  looked for:  ", e.layout, ", under the data root")
        print(io, "\n  data root:   \$", ROOT_VAR, " when that is set; dat/ in the clone when it is not")
        isempty(e.env) ||
            print(io, "\n  or set:      \$", e.env, ", to name this one ", thing, " somewhere else")
    end
    if e.origin === :notdistributed
        print(io, "\n  That ", thing, " holds a live service's own content, which is not",
                  "\n  redistributed here or with the archived data set: a fresh query",
                  "\n  answers with different data rather than these bytes.",
                  "\n  docs/DATA.md section 5 says what a run with it replays and what a",
                  "\n  run without it measures instead.")
    elseif e.origin !== :repo
        print(io, "\n  This repository does not carry that ", thing, ";",
                  "\n  docs/DATA.md says which stage reads it and how to obtain or rebuild it.")
    end
    return nothing
end

"""
    require_input(name; env = ENV) -> String

The absolute path of `name`, or `MissingInput`. This is the only way an analysis step
should reach a file it does not carry: it is what turns "no such file or directory,
followed by a stack trace" into a sentence naming the input, its place in the layout and
the variable that moves it.
"""
function require_input(name::AbstractString; env = ENV)
    i = input(name)
    p = input_path(name; env = env)
    ok = i.kind === :dir ? isdir(p) : isfile(p)
    ok && return p
    via = !isempty(i.env) && !isempty(strip(get(env, i.env, "")))
    throw(MissingInput(i.name, layout(i), i.kind, i.origin, i.env, i.what, i.read_by, via))
end

"""
    require_or_exit(names...; env = ENV, io = stderr)

Require every one of `names`, in order; on the first one that is not there, print the
refusal and stop the process with status 1.

This is for use inside a program guard, and it exists for one reason: a `throw` at the top
level of a script prints the sentence and then several screens of Julia's own include
machinery, and the stack trace is the part a reader sees first. A refusal about a file
somebody has to go and fetch should be a sentence. Anything that INCLUDES a step calls
`require_input` instead and gets the exception.
"""
function require_or_exit(names::AbstractString...; env = ENV, io = stderr)
    for n in names
        try
            require_input(n; env = env)
        catch e
            e isa MissingInput || rethrow()
            showerror(io, e)
            println(io)
            exit(1)
        end
    end
    return nothing
end

"Every entry the archived data set carries."
deposit_inputs() = filter(i -> i.origin === :deposit, INPUTS)

end # module
