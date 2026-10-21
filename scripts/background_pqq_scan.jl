#=
    background_pqq_scan.jl — which BACKGROUND taxa actually carry a PQQ-ADH?

    The habitat contrast compares taxa that carry a PQQ-dependent alcohol dehydrogenase
    against a background of reference-proteome prokaryotes. That background is built by
    taxonomy alone: every reference-proteome taxon that is not a census taxon. It has
    never been checked for the enzyme itself.

    THE FAILURE THIS PREVENTS. The census is a UniProtKB hmmsearch snapshot, so a taxon
    lands outside it for two very different reasons:
      (a) it has no PQQ-ADH — a true background organism, and
      (b) it has one, but the census run did not see it (a proteome added after the census
          release, an entry whose TaxID moved, a sequence below the census threshold).
    Case (b) organisms are carriers sitting in the comparator arm. They drag the background
    habitat profile towards the carrier profile, which biases every contrast towards the
    null — and the bias is invisible, because nothing in the pipeline ever asks whether a
    background taxon has the gene. This script asks, with the census HMM itself.

    METHOD — targeted scan, not a full-proteome scan. Downloading every background proteome
    is several GB for a question about one domain family, so the search space is narrowed
    first, to the UniProtKB entries InterPro places in EITHER of two entries:
      IPR011047  Quinoprotein alcohol dehydrogenase-like superfamily (homologous superfamily)
      IPR017512  PQQ-dependent dehydrogenase, methanol/ethanol family (family)
    The names are InterPro's own. IPR017512 is the whole methanol/ethanol family, Ca- and
    Ln-dependent members alike; it is not a lanthanide signature, and nothing here treats
    it as one.

    WHY THE UNION. The completeness criterion this step is held to names both entries, and
    the recall clause it carries can only be evaluated over a search space that covers
    both. The superfamily alone does not: InterPro annotates some family members without
    the superfamily entry, and under a superfamily-only prefilter a taxon whose only family
    members are of that kind never reaches hmmsearch, so it comes out with zero hits
    however permissive the threshold. The union puts those entries in front of the model;
    whether they then reach the cutoff is measured, not assumed — the run reports how many
    family entries and frame taxa sit outside the superfamily and how many of those
    entries are hits. The superfamily pages keep their original directory names, so a
    download already on disk is reused rather than re-paged; the family pages go into
    their own directories and the union is taken by accession, locally. Each download
    keeps its own release record, and the run stops if the two came from different
    UniProtKB releases: a union of two snapshots is not a snapshot.

    RECALL. The union is still a prefilter, and a prefilter is an assumption, so the
    recall stage measures it: --recall-n frame taxa get their reference proteome
    downloaded and scanned at the same threshold, and of the taxa the full scan flags, a
    share of at least MIN_RECALL must also be flagged by the prefilter. A TaxID UniProt
    lists no reference proteome under gets every UniProtKB entry of that TaxID instead;
    the release record beside the download says which, and the acceptance table counts
    them. The sample is the first N draws of one seeded permutation of the frame, so a
    larger N keeps every taxon a smaller one drew — and the run CHECKS that against the
    sample column of an earlier scan table instead of taking it on faith. A sample in
    which the full scan flags nothing is a MISS, not a pass: it measured no recall at all.

    A RATE IS NOT A BOUND, AND THE TWO ANSWER DIFFERENT QUESTIONS. The recall row this
    scan gates on is the rate: it asks whether the prefilter finds the carriers the
    full-proteome scan finds, which is what decides whether the exclusion list can be
    published. A small sample with every carrier recovered answers that question, and a
    wide interval around the answer does not make the list wrong. The stronger statement,
    that the background left after the exclusion holds no carrier the prefilter missed,
    needs more than an observed rate: a perfect score on a small number of flagged taxa is
    still compatible with a true recall below the bar. So recall.txt also records
    lower_bound_95, the exact one-sided 95 % Clopper-Pearson lower bound of the recovered
    share, and the acceptance table discloses it without a bar of its own here. The
    habitat analysis downstream decides that stronger statement on the bound, held to the
    same recall bar, rather than on the rate row's PASS, so a scan can pass its own gate
    and still leave that statement unlicensed.

    Both sides assign a sequence to an organism by the entry's own TaxID — the OX= token
    of the proteome FASTA on one side, the organism column of the search-space pages on
    the other — and not by the proteome a sequence was downloaded under. A proteome entry
    whose TaxID is not a sampled taxon is therefore counted for neither. Whether any such
    entry exists is not assumed either way: the run counts them and their hits. The one
    exception is a header with no OX= token at all, which has no TaxID of its own to go
    by; it is attributed to the taxon it was downloaded under, and the run counts those
    headers too, so the exception is a number rather than an assumption.

    WHOLE-FRAME FAMILY MEMBERSHIP. The per-taxon IPR017512 columns are read from the same
    two downloads the search space is built from, which carry xref_interpro and the
    sequence length for every entry and whose releases the run requires to agree. The
    membership is a property of the TAXON, taken from every family entry the taxon has,
    independent of the hit table. That independence is what lets the family acceptance
    check fail: a full-length family member the model scores below the cutoff leaves its
    taxon unflagged while the roll-up still counts it. The other hole this check answers
    for, a frame family entry outside the search space, cannot arise while the search
    space is the union: every family entry of a frame taxon that the download holds is a
    candidate. A carrier annotated in neither entry is invisible to the check and to the
    search alike; that is the recall sample's question, not this check's.

    TWO READINGS OF IPR017512 MEMBERSHIP, EACH UNDER ITS OWN NAME. `in_ipr017512` is
    n_ipr017512_entries > 0, a property of the TAXON. The narrower reading — true when some
    HIT row of the taxon is annotated IPR017512 — is false for every zero-hit taxon by
    construction, so a guard built on it can never see a carrier the scan missed; it is
    carried separately as hit_ipr017512. Columns appended after the first nine:
      n_ipr017512_entries          family entries of this exact TaxID, whole-frame query
      n_ipr017512_ge500            of those, length >= 500 aa
      in_ipr017512_ge500           n_ipr017512_ge500 > 0
      ipr017512_outside_prefilter  some family entry of the taxon lacks the superfamily
                                   annotation, so a superfamily-only prefilter never saw it
      hit_ipr017512                the narrow reading: some hit row is annotated IPR017512
      n_hits_ln_or_ca              hit sequences classed Lanthanide or Calcium

    THRESHOLD. The HMM carries no GA/TC/NC cutoffs, so hmmsearch has no threshold of its
    own and -T is mandatory: left to itself it falls back to an E-value default, and an
    E-value moves with the size of the database searched, so the small candidate set and
    the large proteome set would be scanned at two different stringencies. A bit score
    does not move, but a hand-picked one would be unreproducible, so the threshold is set
    EMPIRICALLY — the same hmmsearch is run over the census accessions present in the
    downloaded search space and the cutoff is read off their full-sequence bit scores.
    That puts the background scan on the census's own scale, as permissive as the census
    it is compared against up to the low tail the rule gives up — the share of census
    accessions below the cutoff is not assumed small, it is the recovery the acceptance
    table gates on. That is the property the comparison needs. `--threshold <bits>`
    replaces the calibrated cutoff with a recorded one for a cheap re-run; threshold.txt
    records which of the two was used.

    SPECIFICITY IS DISCLOSED, NOT GATED. Every acceptance criterion bounds the flagged set
    from one side: a carrier left in the background fails a check, a non-carrier removed
    from it fails none. A hit at the census cutoff is not necessarily a PQQ-ADH — the
    census itself keeps sequences whose alignment never reaches the active-site column —
    so the acceptance table reports how many flagged taxa have no hit classed Lanthanide
    or Calcium. There is no pre-registered bar for that number, and the table says so
    instead of inventing one.

    THE POSITIVE CONTROLS CHECK THE WIRING, NOT THE FINDING. Five named frame taxa must come
    out flagged. They were not drawn from the literature: each was selected, before this
    search was first run, because a UniProtKB entry of that exact TaxID is annotated to the
    methanol/ethanol family at full length. So the check asks whether this prefilter and this
    profile recover what public annotation already claims, and a PASS is not evidence about
    any of those organisms' cofactors. POSITIVE_CONTROLS names the five entries.

    ONE RUN, ONE IDENTITY. The per-taxon table and the three records that describe it —
    threshold.txt, recall.txt and acceptance.tsv — mean something only together, so they
    are published together or not at all. The table is written first, as
    background_pqq_scan.csv.pending; the sha256 of those bytes is the run's scan_sha256,
    and it goes, with one run_id stamp, into all three records, each written as
    <name>.pending too. Only when no acceptance verdict is MISS are the four renamed to
    their final names. On a MISS nothing is promoted: the failing run's files stay
    .pending for inspection, a table and records left by an earlier run are renamed to
    .stale, and the run exits 1.
    The failure this prevents: written in run order, a failing run would overwrite
    threshold.txt and recall.txt and leave the previous run's table in place, so a
    consumer would read one run's cutoff and recall beside another run's exclusion list
    with no file saying they did not belong together. Under the rule above there is no table
    left to misread, and a consumer that hashes the table it reads can refuse any record
    that names another.

    acceptance.tsv columns: check, value, denominator, bar, verdict, note, run_id,
    scan_sha256. value and bar are numbers on the same scale, and a row with a bar is
    gated: its verdict is value >= bar, so it can be re-derived from the row. denominator
    is the population the value was measured over; a gated row whose denominator is zero,
    or whose value is missing or NaN, is a MISS, because a check over nothing measured
    nothing. A row without a bar is DISCLOSED and never blocks. Every bar is a fixed
    parameter of this script or computed by this run, and the note says which.

    CACHES ARE CHECKED, NOT TRUSTED. A derived file is reused only when the parameters
    that produced it still hold: a subset FASTA must hold exactly the requested accession
    set, a hmmsearch table must have been produced at the requested -T (read back from
    the table's own option line), the reference FASTA must hold the requested pair, and
    the recall concatenation must carry a manifest naming exactly the sampled taxa. A
    stale derived file is rebuilt, loudly, and the search built on it is dropped. A
    download that disagrees with the server's own X-Total-Results cannot be repaired
    locally, so that stops the run. One kind of file is reused on its name alone: a
    per-taxon proteome download. It is written through a rename, so it is complete as
    fetched, but nothing can re-check it against a release that has since moved on — so
    every new download records the UniProtKB release it came from, an acceptance row
    fails if a recorded release differs from the search space's, and a download made
    before those records existed is counted as unrecorded rather than assumed current.

    ASYMMETRY, ON PURPOSE. Carriers found here are removed from the BACKGROUND only. The
    census arm is left untouched: this step cleans the comparator, it does not recruit new
    members into the foreground.

    Everything downloaded lands under habitat/background_scan/ in the data root and every stage
    resumes from what is already on disk, so a dropped connection costs one page, not the
    run. Requests are anonymous.

    Run:  julia --project=. scripts/background_pqq_scan.jl [options]
          (--help lists them and exits before anything is loaded, fetched or written)
=#

const REPO = normpath(joinpath(@__DIR__, ".."))

# ---------------------------------------------------------------------------
# Command line. Parsed before any package is loaded and before any file is touched, so
# --help and a mistyped option cost nothing. Both would otherwise cost a run: an ignored
# unknown option means `--thresold 50` scans at the calibrated cutoff and rewrites every
# record as if the flag had been honoured, and a `--help` read after the work starts begins
# the downloads.
# ---------------------------------------------------------------------------

# The sample size lives with the usage text because the usage text quotes it.
const RECALL_N_DEFAULT = 150

const USAGE = """
Usage: julia --project=. scripts/background_pqq_scan.jl [options]

Options:
  --threshold <bits>        scan at this hmmsearch bit score instead of the cutoff
                            calibrated on the census accessions
  --recall-n <n>            frame taxa in the prefilter recall sample (default $RECALL_N_DEFAULT);
                            the first n draws of one seeded permutation of the frame
  --recall-reference <csv>  an earlier background_pqq_scan.csv whose recall sample this
                            run's sample must extend (default: the published table, or
                            the .stale one a failing run displaced)
  -h, --help                print this text and exit

Exit status: 0 when the table and its records were promoted; 1 on an acceptance MISS
(nothing promoted), a failed stage or a command-line error.
"""

const OPTION_SPEC = Dict("--threshold" => :value, "--recall-n" => :value,
                         "--recall-reference" => :value, "--help" => :flag, "-h" => :flag)

"""
    parse_options(args, spec) -> Dict{String,String}

Split a command line into option => value, a flag mapping to "". `spec` names every
option the script defines as `:value` or `:flag`. Anything else is an error — an unknown
or misspelt option, a stray word, an option given twice, a value-taking option with no
value — because an ignored token is a silently different run: the output looks exactly
like the answer to the command that was typed.
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
            # `--threshold --recall-n 20` must not scan at a threshold called "--recall-n".
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

Run a command-line step; an `ErrorException` from it ends the program with the message
and the usage on stderr and exit status 1. The point is that a bad command line never
falls through to a default: it stops here, before anything is loaded or written.
"""
function cli_or_exit(f)
    try
        return f()
    catch e
        e isa ErrorException || rethrow()
        println(stderr, "background_pqq_scan.jl: ", e.msg, "\n")
        print(stderr, USAGE)
        exit(1)
    end
end

# Only when this file is the program. test/test_background_scan.jl includes it with its
# own ARGS, and must get the helpers without the command line being read at all.
const CLI = abspath(PROGRAM_FILE) == (@__FILE__) ?
    cli_or_exit(() -> parse_options(ARGS, OPTION_SPEC)) : Dict{String,String}()
if haskey(CLI, "--help") || haskey(CLI, "-h")
    print(USAGE)
    exit(0)
end

"""
    option_number(cli, flag, T, default; what, min_value = nothing) -> T

The typed value of a parsed option, or `default` when it was not given. A value that does
not parse as `T`, is not finite, or lies below `min_value` is an error naming the option —
never a fall-back to the default, which would run a different scan under the name of the
one that was asked for.
"""
function option_number(cli::AbstractDict, flag::AbstractString, ::Type{T}, default;
                       what::AbstractString = "a number", min_value = nothing) where {T<:Real}
    haskey(cli, flag) || return default
    v = tryparse(T, cli[flag])
    (v === nothing || (v isa AbstractFloat && !isfinite(v))) &&
        error("$flag needs $what, got '$(cli[flag])'")
    min_value === nothing || v >= min_value ||
        error("$flag needs $what of at least $min_value, got '$(cli[flag])'")
    return v
end

# Typed option values, checked at the same point, so a bad value also stops the program
# before anything is loaded, fetched or written. `nothing` when this file is included.
const OPT = abspath(PROGRAM_FILE) == (@__FILE__) ? cli_or_exit(() -> (
    threshold = option_number(CLI, "--threshold", Float64, nothing; what = "a bit score"),
    recall_n = option_number(CLI, "--recall-n", Int, RECALL_N_DEFAULT;
                             what = "a whole number of taxa", min_value = 1),
    recall_reference = get(CLI, "--recall-reference", ""))) : nothing

using Pkg; Pkg.activate(REPO)
using DataFrames, CSV, HTTP, CodecZlib, Printf, Statistics, Random, Dates, SHA

# Names in the one registry (src/ModuleDataInputs.jl). The profile HMM ships with the
# code and resolves in the repository; everything else resolves under the data root.
include(joinpath(REPO, "src", "ModuleDataInputs.jl"))

const SCAN_DIR   = DataInputs.input_path("background_scan_dir")
const HMM        = DataInputs.input_path("profile_hmm")
const CENSUS_CSV = DataInputs.input_path("census_habitat")
const POOL_CSV   = DataInputs.input_path("pool_taxa")
const OUT_CSV    = DataInputs.input_path("background_scan_table")
const RECALL_TXT = DataInputs.input_path("scan_recall")
const ACCEPT_TSV = DataInputs.input_path("scan_acceptance")
const THRESH_TXT = DataInputs.input_path("scan_threshold")

const IPR_SUPERFAMILY = "IPR011047"   # Quinoprotein alcohol dehydrogenase-like superfamily
const IPR_FAMILY      = "IPR017512"   # PQQ-dependent dehydrogenase, methanol/ethanol family
# The calibration pair is the one the census itself calibrated on (a recorded project
# decision). Record order matters: calibrate_active_site reads record 1 as the Ln-type
# reference and record 2 as the Ca-type one.
const LN_REF      = "C5B120"
const CA_REF      = "P16027"
const PAGE        = 500           # UniProt REST maximum page size
const RECALL_SEED = 20260913      # recorded so the recall sample is reproducible
const HMM_CPU     = 8
# A family member shorter than this is treated as a fragment rather than as a carrier; the
# column names n_ipr017512_ge500 / in_ipr017512_ge500 record the floor.
const MIN_FAMILY_LEN = 500
# Acceptance bars, fixed before the scan was run. Both are floors on a share; the flagged-
# count floor is not a constant at all but computed by the run (see the acceptance stage),
# and specificity has no bar, for the reason the header gives.
const MIN_RECOVERY = 0.99
const MIN_RECALL   = 0.98
# The default level of clopper_pearson_lower. The bound recall.txt records does not read
# it: that key is named lower_bound_95, so its level is written where the key is filled.
const RECALL_ALPHA = 0.05
const LN_OR_CA     = ("Lanthanide", "Calcium")
# THE POSITIVE CONTROLS, AND WHERE THEY CAME FROM. Five taxa of the frame this step searches
# — present in the reference-proteome pool, absent from the census — each of which must come
# out flagged or the run publishes nothing. They were NOT taken from the literature, and
# nothing here claims a published cofactor measurement for any of them. They were selected
# from PUBLIC UNIPROT ANNOTATION: each is a frame taxon that has a UniProtKB entry of at least
# MIN_FAMILY_LEN residues annotated to the methanol/ethanol-family PQQ-dependent dehydrogenase
# entry (InterPro IPR017512) under that exact TaxID. The entry that put each taxon on this
# list is named beside it below, so the selection is reproducible from public records alone.
#
# WHAT THAT MAKES THE CHECK. The criterion is annotation, and the list was fixed before this
# search was first run, so flagging all five asks whether this prefilter and this profile
# recover what public annotation already claims. That is a wiring check on the search, not a
# measurement of what the search discovers, and not evidence about any of these organisms'
# cofactors — all five of those annotations are inferred from homology, not measured.
# The list cannot have been read off this scan's output: it does not order or top the scan's
# own hit counts or bit scores, and no rule over the scan table selects these five.
#
# ONE ORGANISM HAS PUBLISHED COFACTOR EVIDENCE, AND THIS LIST DOES NOT INHERIT IT. A
# lanthanide resolved in the active site of a purified methanol dehydrogenase is published for
# Methylacidiphilum fumariolicum STRAIN SolV (Pol et al. 2014). That is a different protein
# and a different TaxID from the one used here: the control below is the species-level TaxID,
# and the annotated entry behind it comes from a later, separately submitted genome whose
# strain is not recorded, so it is not established to be that strain. The strain-level
# evidence is carried where it belongs, against its own accession, in
# dat/control_set/known_cofactor_mdhs.tsv. The name string below is UniProt's spelling for
# that TaxID; the validly published epithet is fumariolicum, as the control set writes it.
#
# WHAT THE LITERATURE DOES AND DOES NOT ESTABLISH, TAXON BY TAXON — recorded so that
# "positive control" is never read as a claim about the organism's metabolism. None of it
# changes a taxon's standing under the annotation criterion above; a PASS means this
# prefilter and this profile recover what public annotation claims.
#
#   591154  The only one of the five with a published cofactor determination, and it does
#           not attach to the entry selected here (the ONE ORGANISM paragraph above). Pol
#           et al. 2014 (doi:10.1111/1462-2920.12249) show growth strictly dependent on
#           lanthanides and a lanthanide in the crystallised XoxF-type dehydrogenase — for
#           strain SolV, whose protein sits under a different TaxID. 591154 is the species
#           taxon, so species-level support holds, but the entry that makes it a carrier
#           here comes from an unpublished direct genome submission that names no strain.
#   511746  Genome evidence only: an mxaF/xoxF homologue (Minf_0992) alongside the PQQ
#           biosynthesis genes, Hou et al. 2008 (doi:10.1186/1745-6150-3-26). No cofactor
#           determination is published for any dehydrogenase of this species.
#   418708  Listed with mxaF, xoxF and exaF in a published genome screen, Diep et al. 2026
#           (doi:10.1038/s41589-026-02176-3). Its species description reports growth on a
#           complex medium containing 0.5 % methanol, Weon et al. 2008
#           (doi:10.1099/ijs.0.65047-0) — methanol tolerance, not growth on methanol as the
#           sole carbon and energy source.
#   217069  Published evidence AGAINST methylotrophy, and it is about the right strain: the
#           type strain NCIMB 1785T = DSM 15561T, the genome behind the entry selected here,
#           is described as using neither methanol nor methylamine as a sole carbon and
#           energy source, Fritz et al. 2004 (doi:10.1099/ijs.0.02746-0). Two limits travel
#           with that. Failure to grow on methanol does not prove the absence of a
#           xoxF-family gene, and the 2004 medium was calcium-supplemented with no
#           lanthanides, at a date when lanthanide-dependent dehydrogenases were unknown, so
#           a lanthanide-dependent enzyme could not have been detected by that assay. The
#           taxon was described as Angulomicrobium amanitiforme and later reclassified.
#   980347  Nothing published on methanol, a PQQ-dehydrogenase or lanthanides for this
#           species; a moderately halophilic salt-lake isolate, Chen et al. 2011
#           (doi:10.1007/s10482-011-9591-0). That reading rests on the description's
#           abstract, not its full text, and the full text's physiology table is where a
#           methanol entry would appear if one existed.
const POSITIVE_CONTROLS = [(591154, "Methylacidiphilum fumarolicum"),   # A0ABM9ICZ5
                           (511746, "Methylacidiphilum infernorum"),    # A0ABX7PVQ1
                           (418708, "Methylobacterium aerolatum"),      # A0ABU0I1Z8
                           (217069, "Ancylobacter amanitiformis"),      # A0ABU0LS59
                           (980347, "Litchfieldella qijiaojingensis")]  # A0ABQ2ZAD6
# A record is written to <name>.pending and renamed to <name> only when the run publishes;
# a record an unpublished run displaced is renamed to <name>.stale.
const PENDING_EXT = ".pending"
const STALE_EXT   = ".stale"

# ---------------------------------------------------------------------------
# Pure helpers. No network, no HMMER — unit-tested in test/test_background_scan.jl.
# ---------------------------------------------------------------------------

"""
    uniprot_accession(target) -> String

Recover the plain accession from a hmmsearch target name or a FASTA header. UniProt
writes `db|ACC|ENTRY_NAME`; anything else is returned as its first whitespace token, so
a locally written `>ACC` header works too.
"""
function uniprot_accession(target::AbstractString)
    s = strip(target)
    isempty(s) && return ""
    t = String(first(split(s)))
    parts = split(t, '|')
    return length(parts) >= 2 ? String(parts[2]) : t
end

"""
    parse_tblout(path) -> DataFrame

Read a `hmmsearch --tblout` table into (target, accession, full_evalue, full_bits,
best_dom_evalue, best_dom_bits). The table LOOKS column-aligned but is not fixed-width —
a long target name pushes the numbers right — so fields are taken by whitespace
splitting, never by character offset. Comment lines start with '#', and the free-text
description from field 19 on is dropped.
"""
function parse_tblout(path::AbstractString)
    df = DataFrame(target = String[], accession = String[], full_evalue = Float64[],
                   full_bits = Float64[], best_dom_evalue = Float64[],
                   best_dom_bits = Float64[])
    isfile(path) || return df
    for line in eachline(path)
        (isempty(strip(line)) || startswith(lstrip(line), '#')) && continue
        f = split(line)
        length(f) >= 9 || continue
        tgt = String(f[1])
        push!(df, (tgt, uniprot_accession(tgt),
                   parse(Float64, f[5]), parse(Float64, f[6]),
                   parse(Float64, f[8]), parse(Float64, f[9])))
    end
    return df
end

"""
    tblout_threshold(path) -> Union{Float64,Nothing}

The `-T` a `hmmsearch --tblout` table was produced at, read from the `# Option settings:`
line HMMER writes into the table's own trailer. `nothing` when the table has no such
line or no `-T` in it. This is what lets a cached table be checked against the cutoff a
run asks for: the file name does not change with the threshold, the trailer does.
"""
function tblout_threshold(path::AbstractString)
    isfile(path) || return nothing
    for line in eachline(path)
        startswith(line, "# Option settings:") || continue
        m = match(r"(?:^|\s)-T\s+([-+0-9.eE]+)", line)
        return m === nothing ? nothing : tryparse(Float64, m.captures[1])
    end
    return nothing
end

"""
    header_taxid(header) -> Union{Int,Nothing}

Pull the organism TaxID out of a UniProt FASTA header's `OX=` token. The recall stage
reads it because its downloads are keyed by proteome while the search space is keyed by
each entry's own organism; comparing the two on the per-entry TaxID keeps both sides on
the exact-TaxID rule the frame uses. Whether a proteome ever holds an entry of another
TaxID is not assumed: `proteome_taxon_map` counts such entries and `recall_flagged` the
hits it sets aside for them.
"""
function header_taxid(header::AbstractString)
    m = match(r"\bOX=(\d+)", header)
    return m === nothing ? nothing : parse(Int, m.captures[1])
end

"""
    frame_taxids(pool_taxids, census_taxids) -> Vector{Int}

The background frame: reference-proteome taxa that are not census taxa. Exact-TaxID
matching on both sides, with no lineage climbing, because the census TaxIDs are exact
organism ids too; mixing the two rules would silently drop or keep whole genera.
"""
frame_taxids(pool_taxids, census_taxids) =
    sort(collect(setdiff(Set{Int}(pool_taxids), Set{Int}(census_taxids))))

"""
    exclude_carriers(frame, scan) -> (kept = ..., removed = ...)

The primary exclusion rule, stated once as a tested function: drop every frame taxon this
scan flagged (`n_hits > 0`). `frame` needs a `taxid` column, `scan` needs `taxid` and
`n_hits`. This script uses it for its closing report, and the downstream habitat analysis
calls it for its any-hit exclusion, so the rule is defined and tested in one place.

What a zero means has to be read carefully. A frame taxon with `n_hits == 0` — or one
absent from `scan` altogether, which is KEPT the same way — was not necessarily searched:
hmmsearch only ever saw the taxa that had at least one entry in the downloaded search
space (`n_candidates > 0`). For a taxon with no candidate the zero says "UniProtKB
annotates no member of either InterPro entry for this TaxID", not "the HMM found
nothing"; only the recall sample reads such a taxon's whole proteome. The join keeps both
kinds, because this step removes carriers it found and does not guess at the rest.
"""
function exclude_carriers(frame::AbstractDataFrame, scan::AbstractDataFrame)
    flagged = Set{Int}(scan.taxid[scan.n_hits .> 0])
    mask = [!(t in flagged) for t in frame.taxid]
    return (kept = frame[mask, :], removed = frame[.!mask, :])
end

"""
    calibration_threshold(scores; rule = :p1) -> Float64

Turn the census accessions' full-sequence bit scores into the hmmsearch `-T` cutoff.
`:min` keeps every calibration accession by construction but is hostage to the single
most truncated fragment in the set; `:p1`, the 1st percentile, tolerates a 1 % tail of
such fragments. Both are computed so the choice is visible rather than assumed.
"""
function calibration_threshold(scores::AbstractVector{<:Real}; rule::Symbol = :p1)
    isempty(scores) && error("no calibration scores — the calibration hmmsearch returned nothing")
    s = sort(Float64.(scores))
    rule === :min && return first(s)
    rule === :p1 && return quantile(s, 0.01)
    error("unknown threshold rule: $rule")
end

"""
    recovery(scores, threshold, n_total) -> Float64

Fraction of the calibration accessions scoring at or above `threshold`. `n_total` counts
every census accession present in the family, INCLUDING those hmmsearch returned no
score for at all: leaving those out of the denominator would turn a miss into a pass.
"""
recovery(scores, threshold::Real, n_total::Integer) =
    n_total == 0 ? 0.0 : count(>=(threshold), scores) / n_total

"""
    ipr017512_by_taxon(entries, frame; min_len = MIN_FAMILY_LEN) -> DataFrame

Roll the downloaded UniProtKB entries up to one row per FRAME taxon: how many family
(IPR017512) entries the exact TaxID has, how many of those reach `min_len`, and whether
any of them lacks the superfamily (IPR011047) annotation, which is what puts an entry
outside a superfamily-only prefilter.

Membership is read off the InterPro annotation of each entry, not off which download the
entry arrived in, and `entries` covers the whole frame, so every column here is a
property of the TAXON. Nothing in it depends on the hit table: that is the whole point,
because the acceptance check it feeds exists to catch a carrier the model missed, and a
column derived from hits cannot do that.

`entries` needs `taxid`, `len` and `interpro`, one row per entry, already unique by
accession — the roll-up counts rows and cannot tell a duplicate from a second entry.
Frame taxa with no family entry get a zero row, so the result covers `frame` exactly.
"""
function ipr017512_by_taxon(entries::AbstractDataFrame, frame::AbstractVector{<:Integer};
                            min_len::Integer = MIN_FAMILY_LEN,
                            family::AbstractString = IPR_FAMILY,
                            prefilter::AbstractString = IPR_SUPERFAMILY)
    n_all, n_long, outside = Dict{Int,Int}(), Dict{Int,Int}(), Dict{Int,Bool}()
    for r in eachrow(entries)
        occursin(family, r.interpro) || continue
        t = Int(r.taxid)
        n_all[t] = get(n_all, t, 0) + 1
        r.len >= min_len && (n_long[t] = get(n_long, t, 0) + 1)
        occursin(prefilter, r.interpro) || (outside[t] = true)
    end
    df = DataFrame(taxid = Int[], n_ipr017512_entries = Int[], n_ipr017512_ge500 = Int[],
                   in_ipr017512 = Bool[], in_ipr017512_ge500 = Bool[],
                   ipr017512_outside_prefilter = Bool[])
    for t in frame
        na, nl = get(n_all, Int(t), 0), get(n_long, Int(t), 0)
        push!(df, (Int(t), na, nl, na > 0, nl > 0, get(outside, Int(t), false)))
    end
    return df
end

"""
    ipr017512_unflagged(rollup, flagged_taxids) -> Vector{Int}

The frame taxa that carry a full-length family member (`in_ipr017512_ge500`) but were
not flagged by the scan. The denominator comes from the whole-frame roll-up, never from
the prefilter's own candidate rows: a denominator drawn from the search space can only
ever contain taxa the search space reached, and a check built on it passes by
construction exactly when the search space has a hole.
"""
function ipr017512_unflagged(rollup::AbstractDataFrame, flagged_taxids)
    flagged = Set{Int}(flagged_taxids)
    return sort([t for (t, c) in zip(rollup.taxid, rollup.in_ipr017512_ge500)
                 if c && !(t in flagged)])
end

"""
    specificity_disclosure(n_hits, n_hits_ln_or_ca) -> (n_flagged, n_without_ln_or_ca, share)

Among flagged taxa (`n_hits > 0`), how many have no hit sequence classed Lanthanide or
Calcium. Reported, never gated: no bar was fixed for it before the scan, and one made up
afterwards would be tuned to the answer. `share` is NaN when nothing is flagged.
"""
function specificity_disclosure(n_hits::AbstractVector{<:Integer},
                                n_hits_ln_or_ca::AbstractVector{<:Integer})
    flagged = n_hits .> 0
    nf = count(flagged)
    nw = count(flagged .& (n_hits_ln_or_ca .== 0))
    return (n_flagged = nf, n_without_ln_or_ca = nw, share = nf == 0 ? NaN : nw / nf)
end

"""
    recall_sample(frame, n, seed) -> (taxa = sorted Vector{Int}, draw_order = Vector{Int})

The recall sample: the first `n` draws of ONE seeded permutation of `frame` (all of it
when `n` exceeds the frame). `draw_order` keeps the draws in permutation order, which is
what makes samples of different sizes comparable: a sample of `n` is exactly the first
`n` draws of any larger one, so enlarging the sample adds taxa and never swaps any out.
Drawing a fresh subset per size instead would silently replace the taxa an earlier,
smaller recall was measured on.
"""
function recall_sample(frame::AbstractVector{<:Integer}, n::Integer, seed::Integer)
    n >= 1 || error("the recall sample needs at least one taxon, got n = $n")
    perm = randperm(MersenneTwister(seed), length(frame))
    drawn = Int.(frame[perm[1:min(n, length(frame))]])
    return (taxa = sort(drawn), draw_order = drawn)
end

"""
    sample_extends(earlier, draw_order) -> (n_earlier, n_matched, ok)

Does an earlier run's recall sample survive as the first draws of this run's sample?
`ok` is true when the `k = length(earlier)` taxa of `earlier` are exactly the first `k`
entries of `draw_order`. `n_matched` counts the earlier taxa among those first
`min(k, length(draw_order))` draws, so `ok == (n_matched >= n_earlier)` and an acceptance
row can carry the two numbers and re-derive the verdict from them. `ok` is `nothing` when
there is no earlier sample to compare. The failure this catches is a frame or a seed that
moved between runs: the sample then has the requested size and different taxa, and a
recall said to extend an earlier one in fact replaces it.
"""
function sample_extends(earlier, draw_order::AbstractVector{<:Integer})
    e = Set{Int}(earlier)
    k = length(e)
    k == 0 && return (n_earlier = 0, n_matched = 0, ok = nothing)
    prefix = Set{Int}(draw_order[1:min(k, length(draw_order))])
    matched = length(intersect(prefix, e))
    return (n_earlier = k, n_matched = matched, ok = matched >= k)
end

"""
    clopper_pearson_lower(k, n; alpha = RECALL_ALPHA) -> Float64

The exact one-sided Clopper-Pearson lower confidence bound for `k` successes in `n`
trials: the share `p` at which `k` or more successes have probability exactly `alpha`.
It is 0 for `k == 0`, `alpha^(1/n)` for `k == n`, and NaN for `n == 0` — an empty sample
bounds nothing, and a number invented for it would read as evidence. The binomial upper
tail rises monotonically in `p`, so the bound is found by bisection on that tail, summed
in log space so that no binomial coefficient overflows.
"""
function clopper_pearson_lower(k::Integer, n::Integer; alpha::Real = RECALL_ALPHA)
    (0 < alpha < 1) || error("alpha must lie in (0, 1), got $alpha")
    (0 <= k <= n) || error("need 0 <= k <= n, got k = $k, n = $n")
    n == 0 && return NaN
    k == 0 && return 0.0
    k == n && return Float64(alpha)^(1 / n)
    logfact = cumsum([0.0; log.(1:n)])                  # logfact[m + 1] = log(m!)
    tail(p) = sum(exp(logfact[n + 1] - logfact[i + 1] - logfact[n - i + 1] +
                      i * log(p) + (n - i) * log1p(-p)) for i in k:n)
    lo, hi = 0.0, 1.0
    for _ in 1:200
        mid = (lo + hi) / 2
        tail(mid) < alpha ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

"""
    recall_summary(; n_sampled, full_flagged, pre_flagged, n_never_searched, bar, seed)

Score the recall sample: of the sampled taxa the full-proteome scan flags, how many the
prefilter flags too. `missed` lists the ones it did not.

An EMPTY full-proteome result is a MISS. With nothing to recover the rate is undefined
(NaN), and reading "nothing to recover" as "everything recovered" would let a broken or
empty proteome download certify the prefilter. `n_never_searched` counts sampled taxa
with no candidate entry at all, the stratum the prefilter never scans, so the file says
how much of the sample actually probed the hole a prefilter can have.

`pass` is decided on the rate, because this scan's gate asks whether the prefilter finds
carriers. `lower_bound_95` is the exact one-sided 95 % Clopper-Pearson lower bound of the
same two counts, NaN when the full scan flags nothing, and it decides nothing here: it is
recorded for the consumer that needs the stronger statement the header describes. Its
level is written here and not read from RECALL_ALPHA, because the level is part of the
name the record gives it.
"""
function recall_summary(; n_sampled::Integer, full_flagged, pre_flagged,
                          n_never_searched::Integer, bar::Real, seed::Integer)
    ff, pf = Set{Int}(full_flagged), Set{Int}(pre_flagged)
    recovered = length(intersect(ff, pf))
    rate = isempty(ff) ? NaN : recovered / length(ff)
    return (n_sampled = Int(n_sampled), n_full_flagged = length(ff),
            n_recovered = recovered, rate = rate, bar = Float64(bar),
            pass = !isempty(ff) && rate >= bar,
            n_never_searched_in_sample = Int(n_never_searched),
            missed = sort(collect(setdiff(ff, pf))), seed = Int(seed),
            lower_bound_95 = clopper_pearson_lower(recovered, length(ff); alpha = 0.05))
end

"""
    recall_flagged(hit_accessions, acc2tax, sample_set) -> (flagged, n_outside)

The sampled taxa the full-proteome scan flags. Each hit counts for the TaxID its entry
names (`acc2tax`, read from the OX= tokens), and only when that TaxID is a sampled taxon —
the rule the prefilter side follows, since its hits carry each entry's own organism.
`n_outside` counts the hits that rule sets aside, so a proteome holding entries of another
TaxID shows up as a number instead of shrinking the flagged set in silence.
"""
function recall_flagged(hit_accessions, acc2tax::AbstractDict, sample_set)
    flagged = Set{Int}()
    outside = 0
    for a in hit_accessions
        t = get(acc2tax, String(a), 0)
        if t in sample_set
            push!(flagged, t)
        else
            outside += 1
        end
    end
    return (flagged = flagged, n_outside = outside)
end

"""
    recall_keyvals(summary; threshold_bits, uniprot_release) -> Vector{Pair{String,String}}

The recall.txt record up to its identity, keys in a fixed order. `missed_taxids` is
`;`-joined and empty when nothing was missed; `rate` is written as NaN when it is
undefined. The last three keys — run_id, scan_sha256, derived_at — are appended by
`publish_scan`, which is the only place that knows the hash of the table the record
describes.

`lower_bound_95` comes after the keys earlier records already had, so none of those moved.
It is written in the shortest form that reads back to the same Float64 (NaN as "NaN"),
not to the rate's six decimals: a consumer recomputes the bound from the two counts and
refuses a record that disagrees with it beyond rounding error, and a bound rounded for
display would disagree by up to half a millionth.
"""
recall_keyvals(s; threshold_bits::Real, uniprot_release::AbstractString) = [
    "n_sampled" => string(s.n_sampled),
    "n_full_flagged" => string(s.n_full_flagged),
    "n_recovered" => string(s.n_recovered),
    "rate" => @sprintf("%.6f", s.rate),
    "bar" => @sprintf("%.6f", s.bar),
    "pass" => string(s.pass),
    "n_never_searched_in_sample" => string(s.n_never_searched_in_sample),
    "missed_taxids" => join(s.missed, ";"),
    "seed" => string(s.seed),
    "threshold_bits" => @sprintf("%.4f", threshold_bits),
    "uniprot_release" => String(uniprot_release),
    "lower_bound_95" => string(Float64(s.lower_bound_95)),
]

# ---------------------------------------------------------------------------
# Records. Offline, but they read and write files; tested on temp directories.
# ---------------------------------------------------------------------------

"""
    write_keyvals(path, pairs); read_keyvals(path) -> Dict{String,String}

Tab-separated key/value records (recall.txt, threshold.txt, the release records).
Written through a `.part` file and renamed, so a reader never sees half a record and
mistakes a missing `pass` line for an absent check.
"""
function write_keyvals(path::AbstractString, pairs)
    open(path * ".part", "w") do io
        for (k, v) in pairs
            println(io, k, '\t', v)
        end
    end
    mv(path * ".part", path; force = true)
    return path
end

function read_keyvals(path::AbstractString)
    d = Dict{String,String}()
    isfile(path) || return d
    for line in eachline(path)
        p = split(line, '\t'; limit = 2)
        length(p) == 2 && (d[String(strip(p[1]))] = String(strip(p[2])))
    end
    return d
end

# Typed field readers for a key/value record. Absent, empty or unparseable is `nothing`,
# never a default: a consumer that fills a missing `pass` with `false` is safe by luck, and
# one that fills a missing count with 0 turns "not recorded" into "measured nothing".
_kv_str(d, k)   = (s = get(d, k, ""); isempty(s) ? nothing : s)
_kv_int(d, k)   = tryparse(Int, get(d, k, ""))
_kv_float(d, k) = tryparse(Float64, get(d, k, ""))
_kv_bool(d, k)  = (s = get(d, k, ""); s == "true" ? true : s == "false" ? false : nothing)
_is_sha256(s)   = s !== nothing && occursin(r"^[0-9a-f]{64}$", s)

"""
    read_recall(path) -> NamedTuple

recall.txt parsed back into typed fields, every one `nothing` when absent or unparseable.
`missed_taxids` is `nothing` when any listed id fails to parse, rather than the ids that
did: a shortened miss list is a quietly better recall. `scan_sha256` is `nothing` unless it
is a well-formed sha256, so a truncated hash cannot match anything. `lower_bound_95` reads
"NaN" back as NaN, the bound of a sample in which the full scan flagged nothing, and an
absent key as `nothing`, the state of a record written before the key existed.
"""
function read_recall(path::AbstractString)
    d = read_keyvals(path)
    ids = strip(get(d, "missed_taxids", ""))
    parsed = isempty(ids) ? Int[] : [tryparse(Int, strip(t)) for t in split(ids, ';')]
    sha = _kv_str(d, "scan_sha256")
    return (n_sampled = _kv_int(d, "n_sampled"), n_full_flagged = _kv_int(d, "n_full_flagged"),
            n_recovered = _kv_int(d, "n_recovered"), rate = _kv_float(d, "rate"),
            bar = _kv_float(d, "bar"), pass = _kv_bool(d, "pass"),
            n_never_searched_in_sample = _kv_int(d, "n_never_searched_in_sample"),
            missed_taxids = any(isnothing, parsed) ? nothing : Vector{Int}(parsed),
            seed = _kv_int(d, "seed"), threshold_bits = _kv_float(d, "threshold_bits"),
            uniprot_release = _kv_str(d, "uniprot_release"),
            lower_bound_95 = _kv_float(d, "lower_bound_95"), run_id = _kv_str(d, "run_id"),
            scan_sha256 = _is_sha256(sha) ? sha : nothing, derived_at = _kv_str(d, "derived_at"))
end

"""
    read_threshold(path) -> NamedTuple

threshold.txt parsed back into the fields a consumer decides or reports on, typed the
same way as `read_recall`.
"""
function read_threshold(path::AbstractString)
    d = read_keyvals(path)
    sha = _kv_str(d, "scan_sha256")
    return (rule = _kv_str(d, "rule"), threshold_bits = _kv_float(d, "threshold_bits"),
            threshold_source = _kv_str(d, "threshold_source"),
            threshold_min_bits = _kv_float(d, "threshold_min_bits"),
            search_floor_bits = _kv_float(d, "search_floor_bits"),
            n_calibration_accessions = _kv_int(d, "n_calibration_accessions"),
            n_scored = _kv_int(d, "n_scored"), recovery = _kv_float(d, "recovery"),
            active_site_column = _kv_int(d, "active_site_column"),
            hmm = _kv_str(d, "hmm"), hmm_sha256 = _kv_str(d, "hmm_sha256"),
            run_id = _kv_str(d, "run_id"), scan_sha256 = _is_sha256(sha) ? sha : nothing,
            derived_at = _kv_str(d, "derived_at"))
end

"""
    verdict_of(ok) -> "PASS" | "MISS" | "DISCLOSED"

`nothing` marks a quantity that is reported without a bar. It can never read as a pass,
and it never blocks the run.
"""
verdict_of(ok::Bool) = ok ? "PASS" : "MISS"
verdict_of(::Nothing) = "DISCLOSED"

"""
    floor_verdict(value, denominator, bar) -> Union{Bool,Nothing}

The verdict of one acceptance row, derived from the row's own numbers so the table can
never print PASS beside numbers that fail. `bar` missing: a disclosure, `nothing`.
Otherwise `value >= bar`, except that a missing or NaN value and a denominator of zero are
false — a check that measured nothing must not pass, and "no taxon to recover" is not
"every taxon recovered".
"""
function floor_verdict(value, denominator, bar)
    (bar === missing || bar === nothing) && return nothing
    (value === missing || value === nothing) && return false
    value isa AbstractFloat && isnan(value) && return false
    denominator === missing || denominator === nothing || denominator > 0 || return false
    return value >= bar
end

"""
    num_cell(x) -> String

A machine-readable number for acceptance.tsv: "" for a missing quantity, an integer as
itself, anything else in the shortest form that reads back to the same Float64, NaN as
"NaN". Rounding for display would let a value print as the bar while failing it.
"""
num_cell(::Missing) = ""
num_cell(::Nothing) = ""
num_cell(x::Integer) = string(x)
num_cell(x::Real) = string(Float64(x))

"""
    acceptance_table() -> DataFrame

An empty acceptance record: check, value, denominator, bar, verdict, note, run_id,
scan_sha256, all strings as written. See the header for what the columns mean; run_id and
scan_sha256 stay empty until `publish_scan` stamps them.
"""
acceptance_table() = DataFrame(check = String[], value = String[], denominator = String[],
                               bar = String[], verdict = String[], note = String[],
                               run_id = String[], scan_sha256 = String[])

"""
    accept!(table, check; value, denominator, bar, note, ok) -> table

Append one row. `ok` defaults to the verdict `floor_verdict` derives from the row's own
numbers; a caller passes it only for a check whose verdict is not a floor, and then the
numbers it passes must still re-derive it. A check or note may not contain a tab or a
newline, which would split the row, nor a double quote, which makes the writer quote the
cell, so that a reader splitting on tabs takes the quotes for part of the text.
"""
function accept!(table::AbstractDataFrame, check::AbstractString; value = missing,
                 denominator = missing, bar = missing, note::AbstractString = "",
                 ok = floor_verdict(value, denominator, bar))
    any(c -> c == '\t' || c == '\n' || c == '\r' || c == '"', check * note) &&
        error("acceptance row '$check': a tab, newline or double quote would not survive the record")
    push!(table, (String(check), num_cell(value), num_cell(denominator), num_cell(bar),
                  verdict_of(ok), String(note), "", ""))
    return table
end

"""
    acceptance_blocks(table) -> Bool

True when any acceptance row is a MISS. Nothing is promoted on it.
"""
acceptance_blocks(table::AbstractDataFrame) = any(==("MISS"), table.verdict)

"""
    recall_acceptance!(table, summary) -> table

The two recall rows of the acceptance table, from one `recall_summary`. The first is gated
on the rate against the summary's bar: it is this scan's own question, whether the
prefilter finds carriers, and a MISS on it stops publication. The second discloses the
exact one-sided 95 % lower bound with no bar, because a bound below the bar does not make
the exclusion list wrong; the stronger background statement is decided on that bound
downstream. Both rows live here, above the driver, so a suite can check that the gate
stays on the rate and that the bound never gains a bar that would block a correct list.
"""
function recall_acceptance!(table::AbstractDataFrame, s)
    accept!(table, "prefilter recall"; value = s.rate, denominator = s.n_full_flagged,
            bar = s.bar,
            note = "fixed bar; $(s.n_recovered) of $(s.n_full_flagged) taxa the full-proteome scan flags, seeded sample of $(s.n_sampled); an empty full scan measures nothing and fails")
    accept!(table, "prefilter recall, exact one-sided 95 % lower bound";
            value = s.lower_bound_95, denominator = s.n_full_flagged,
            note = "no pre-registered bar, disclosed; Clopper-Pearson on $(s.n_recovered) of $(s.n_full_flagged); this scan gates on the rate, and the stronger background statement downstream is held to the recall bar on this bound")
    return table
end

"""
    write_acceptance(path, table); read_acceptance(path) -> DataFrame

acceptance.tsv out and back. The reader returns value, denominator and bar as
`Union{Missing,Float64}` and the text columns as strings; a missing file is an EMPTY table
— no verdicts, which a consumer must read as "not accepted", never as "nothing failed".
A file without the expected columns is an error: it is a record of some other shape, and
guessing its columns would turn an unknown into a verdict.
"""
function write_acceptance(path::AbstractString, table::AbstractDataFrame)
    CSV.write(path * ".part", table; delim = '\t')
    mv(path * ".part", path; force = true)
    return path
end

const ACCEPTANCE_COLUMNS = ["check", "value", "denominator", "bar", "verdict", "note",
                            "run_id", "scan_sha256"]

function read_acceptance(path::AbstractString)
    out = DataFrame(check = String[], value = Union{Missing,Float64}[],
                    denominator = Union{Missing,Float64}[], bar = Union{Missing,Float64}[],
                    verdict = String[], note = String[], run_id = String[],
                    scan_sha256 = String[])
    isfile(path) || return out
    raw = CSV.read(path, DataFrame; delim = '\t', types = String, missingstring = nothing)
    names(raw) == ACCEPTANCE_COLUMNS ||
        error("$(basename(path)) does not have the acceptance columns: $(names(raw))")
    num(s) = isempty(s) ? missing : something(tryparse(Float64, s), missing)
    for r in eachrow(raw)
        push!(out, (r.check, num(r.value), num(r.denominator), num(r.bar), r.verdict,
                    r.note, r.run_id, r.scan_sha256))
    end
    return out
end

"""
    csv_bytes(table) -> Vector{UInt8}

The per-taxon table serialised exactly as the promoted CSV will hold it. `publish_scan`
writes these bytes and hashes these bytes, so scan_sha256 is the hash of the file a
consumer reads and not of some other rendering of the same table.
"""
function csv_bytes(table::AbstractDataFrame)
    io = IOBuffer()
    CSV.write(io, table)
    return take!(io)
end

function write_bytes(path::AbstractString, bytes::AbstractVector{UInt8})
    open(io -> write(io, bytes), path * ".part", "w")
    mv(path * ".part", path; force = true)
    return path
end

"""
    publish_scan(; out_csv, csv, threshold_path, threshold_pairs, recall_path, recall_pairs,
                   acceptance_path, acceptance, run_id) -> NamedTuple

Write one run's per-taxon table and its three records, then promote all four or none.

1. The table's bytes go to `<out_csv>.pending`; their sha256 is the run's scan_sha256.
2. threshold.txt and recall.txt get `run_id`, `scan_sha256` and `derived_at` (= run_id)
   appended, acceptance.tsv gets both identity columns filled, and each is written as
   `<name>.pending`. The hash is taken HERE and nowhere else, so no record can carry a
   hash other than that of the table beside it.
3. No verdict MISS: the three records and then the table are renamed to their final
   names, replacing whatever an earlier run left.
   Any verdict MISS: nothing is promoted. A table and records an earlier run left under
   the final names are renamed to `.stale`, so no consumer can pair them with this run or
   read a table no current run vouches for. This run's files stay `.pending`.

Returns (promoted, scan_sha256, staled, pending) with the paths moved and left behind. A
`.stale` file from an earlier failure is not touched by a later promotion and is never
read as a scan result; the one read of it is the recall sample column, see
`default_recall_reference`.
"""
function publish_scan(; out_csv::AbstractString, csv::AbstractVector{UInt8},
                        threshold_path::AbstractString, threshold_pairs,
                        recall_path::AbstractString, recall_pairs,
                        acceptance_path::AbstractString, acceptance::AbstractDataFrame,
                        run_id::AbstractString)
    isempty(run_id) && error("publish_scan: an empty run_id identifies nothing")
    write_bytes(out_csv * PENDING_EXT, csv)
    sha = bytes2hex(sha256(read(out_csv * PENDING_EXT)))
    identity = ["run_id" => String(run_id), "scan_sha256" => sha, "derived_at" => String(run_id)]
    write_keyvals(threshold_path * PENDING_EXT, vcat(collect(threshold_pairs), identity))
    write_keyvals(recall_path * PENDING_EXT, vcat(collect(recall_pairs), identity))
    table = copy(acceptance)
    table.run_id = fill(String(run_id), nrow(table))
    table.scan_sha256 = fill(sha, nrow(table))
    write_acceptance(acceptance_path * PENDING_EXT, table)

    records = [threshold_path, recall_path, acceptance_path]
    if acceptance_blocks(table)
        staled = String[]
        for p in vcat(records, [out_csv])
            isfile(p) || continue
            mv(p, p * STALE_EXT; force = true)
            push!(staled, p * STALE_EXT)
        end
        return (promoted = false, scan_sha256 = sha, staled = staled,
                pending = [p * PENDING_EXT for p in vcat(records, [out_csv])])
    end
    # The table last: a run killed between two renames leaves new records beside an older
    # table, and the hash those records carry is not that table's, so a consumer that
    # hashes the table it reads refuses the pair.
    for p in vcat(records, [out_csv])
        mv(p * PENDING_EXT, p; force = true)
    end
    return (promoted = true, scan_sha256 = sha, staled = String[], pending = String[])
end

# ---------------------------------------------------------------------------
# Local file helpers. Offline, but they read and write files; tested on temp files.
# ---------------------------------------------------------------------------

"""
    parse_entry_tsv(files) -> DataFrame

Read UniProt TSV pages (accession, organism_id, length, xref_interpro) into
(accession, taxid, len, interpro). Every page repeats the `Entry` header line, so it is
skipped per file. Rows are returned as read, duplicates included: the caller compares the
row count with X-Total-Results BEFORE deduplicating, so a page fetched twice or a row
that failed to parse shows up as a count mismatch instead of vanishing.
"""
function parse_entry_tsv(files)
    df = DataFrame(accession = String[], taxid = Int[], len = Int[], interpro = String[])
    for f in files
        first_line = true
        for line in eachline(f)
            if first_line
                first_line = false
                startswith(line, "Entry\t") && continue
            end
            isempty(strip(line)) && continue
            p = split(line, '\t')
            length(p) >= 4 || continue
            tid = tryparse(Int, strip(p[2]))
            tid === nothing && continue
            push!(df, (String(strip(p[1])), tid,
                       something(tryparse(Int, strip(p[3])), 0), String(p[4])))
        end
    end
    return df
end

"""
    fasta_accessions(path) -> Set{String}

The accessions a FASTA file holds, read from its headers. Empty for a missing file.
"""
function fasta_accessions(path::AbstractString)
    s = Set{String}()
    isfile(path) || return s
    for line in eachline(path)
        startswith(line, '>') && push!(s, uniprot_accession(SubString(line, 2)))
    end
    return s
end

"""
    filter_fasta(files, keep, out_path) -> Int

Stream FASTA records out of `files`, keeping those whose accession is in `keep`, each
accession ONCE. The search space is the union of two downloads that overlap by
construction (most family entries are superfamily entries too), so the same record
arrives in both page sets. Written twice, hmmsearch would report it twice and every
per-taxon hit count would double for exactly the taxa that matter most.

Written through a `.part` file and renamed, because later stages decide "already done"
by what is on disk, and a truncated FASTA left by a killed run would otherwise be
scanned as if it were the whole candidate set.
"""
function filter_fasta(files, keep::Set{String}, out_path::AbstractString)
    n = 0
    seen = Set{String}()
    open(out_path * ".part", "w") do io
        for f in files
            emit = false
            for line in eachline(f)
                if startswith(line, '>')
                    acc = uniprot_accession(SubString(line, 2))
                    emit = acc in keep && !(acc in seen)
                    if emit
                        push!(seen, acc)
                        n += 1
                    end
                end
                emit && println(io, line)
            end
        end
    end
    mv(out_path * ".part", out_path; force = true)
    return n
end

"""
    sync_subset_fasta(out_path, files, accs; label) -> (n = Int, rebuilt = Bool)

Make `out_path` hold exactly the records for `accs`, reusing it only when its accession
set already equals `accs`. A subset FASTA left over from a different candidate set is the
failure this prevents: the hmmsearch that follows would scan yesterday's sequences and
its table would look exactly like a fresh one. `rebuilt` tells the caller that anything
derived from the old file is stale too. A requested accession that the pages do not
contain is an error, because a silently smaller FASTA is a silently smaller search.
"""
function sync_subset_fasta(out_path::AbstractString, files, accs::Set{String};
                           label::AbstractString)
    have = fasta_accessions(out_path)
    isfile(out_path) && have == accs && return (n = length(have), rebuilt = false)
    isfile(out_path) && @warn "$label: cached FASTA does not hold the requested accession set, rebuilding" n_cached = length(have) n_requested = length(accs)
    n = filter_fasta(files, accs, out_path)
    n == length(accs) ||
        error("$label: $(length(accs) - n) requested accession(s) have no record in the downloaded FASTA pages")
    return (n = n, rebuilt = true)
end

"""
    verify_download(label, rows, fasta_files, expected_total) -> Int

Check one paged download against the total the server announced for its query: every
TSV row parsed, no accession twice, one FASTA record per entry, and the FASTA holding the
same accessions as the TSV. Any mismatch is an error. A page lost to a dropped connection
or fetched twice cannot be told apart from a real change in UniProtKB by anything later
in the pipeline; it would just lower or inflate the flagged count.
"""
function verify_download(label::AbstractString, rows::AbstractDataFrame, fasta_files,
                         expected_total::Integer)
    fa = Set{String}()
    n_fa = 0
    for f in fasta_files, line in eachline(f)
        startswith(line, '>') || continue
        n_fa += 1
        push!(fa, uniprot_accession(SubString(line, 2)))
    end
    tsv = Set{String}(rows.accession)
    problems = String[]
    nrow(rows) == expected_total ||
        push!(problems, "TSV rows $(nrow(rows)) != X-Total-Results $expected_total")
    length(tsv) == nrow(rows) ||
        push!(problems, "TSV holds $(nrow(rows) - length(tsv)) duplicate accession row(s)")
    n_fa == expected_total ||
        push!(problems, "FASTA records $n_fa != X-Total-Results $expected_total")
    fa == tsv ||
        push!(problems, "FASTA and TSV accession sets differ ($(length(setdiff(tsv, fa))) TSV-only, $(length(setdiff(fa, tsv))) FASTA-only)")
    isempty(problems) ||
        error("$label download incomplete or inconsistent: " * join(problems, "; ") *
              ". Delete its page directories and re-run.")
    return expected_total
end

"""
    csv_recall_sample(path; required) -> Vector{Int}

The taxa an earlier scan table marks `in_prefilter_recall_sample`, sorted. A missing file
is an empty sample when `required` is false (there is no earlier run to extend) and an
error when true (the file was named on the command line). A table without that column, or
one that marks no taxon, is an error either way: it is not a scan table whose sample can
be extended, and treating it as "no earlier sample" would switch the check off unnoticed.
"""
function csv_recall_sample(path::AbstractString; required::Bool)
    if !isfile(path)
        required && error("recall reference $(basename(path)) does not exist")
        return Int[]
    end
    t = CSV.read(path, DataFrame)
    ("taxid" in names(t) && "in_prefilter_recall_sample" in names(t)) ||
        error("recall reference $(basename(path)) has no taxid / in_prefilter_recall_sample columns")
    s = sort(Int.(t.taxid[coalesce.(t.in_prefilter_recall_sample, false)]))
    isempty(s) && error("recall reference $(basename(path)) marks no taxon as sampled")
    return s
end

"""
    default_recall_reference(out_csv) -> String

The earlier scan table a run's recall sample is checked against when no reference is named
on the command line: the published table, or, when a failing run has displaced it, the
`.stale` table that run left — the last one a run published. Without that fallback a
single MISS would switch the check off for the run after it, which would then find no
table, report "nothing to compare" and publish a sample nobody compared. The `.stale`
table is read for its sample column only, never as a scan result. When neither exists the
published path is returned, and `csv_recall_sample` reads that as no earlier sample.
"""
function default_recall_reference(out_csv::AbstractString)
    isfile(out_csv) && return String(out_csv)
    isfile(out_csv * STALE_EXT) && return String(out_csv * STALE_EXT)
    return String(out_csv)
end

proteome_file(prot_dir, t) = joinpath(prot_dir, "tax_$(t).fasta")
release_record_file(prot_dir, t) = proteome_file(prot_dir, t) * ".release"

"""
    sync_recall_fasta(out_path, prot_dir, sample_taxa) -> (n = Int, rebuilt = Bool)

Concatenate the sampled taxa's proteome downloads into one FASTA for the recall search,
reusing an existing concatenation only while its manifest (`<out_path>.manifest`) names
exactly this sample and the byte size of every file in it. The failure this prevents:
reusing the concatenation whenever it exists lets a changed --recall-n or frame score the
recall on the previous sample's sequences while the TaxID map is built from the new
sample's files — only taxa in both samples can come out flagged, the measured sample
shrinks in silence, and the record still gives the requested size. `rebuilt` tells the
caller to drop the search table built on the displaced file.

The manifest is written after the FASTA is renamed into place, so a run killed between
the two leaves a FASTA without a matching manifest, which the next run rebuilds.
"""
function sync_recall_fasta(out_path::AbstractString, prot_dir::AbstractString,
                           sample_taxa::AbstractVector{<:Integer})
    files = [proteome_file(prot_dir, t) for t in sample_taxa]
    for (t, f) in zip(sample_taxa, files)
        isfile(f) || error("recall: no proteome download for sampled taxon $t")
    end
    want = join(("$t:$(filesize(f))" for (t, f) in zip(sample_taxa, files)), ";")
    manifest = out_path * ".manifest"
    isfile(out_path) && get(read_keyvals(manifest), "sources", "") == want &&
        return (n = length(files), rebuilt = false)
    isfile(out_path) && @warn "recall: the cached proteome concatenation was not built from this sample, rebuilding" n_sampled = length(files)
    open(out_path * ".part", "w") do io
        for f in files
            b = read(f)
            write(io, b)
            # Two files glued without a newline would merge a sequence line and a header.
            isempty(b) || last(b) == UInt8('\n') || write(io, '\n')
        end
    end
    mv(out_path * ".part", out_path; force = true)
    write_keyvals(manifest, ["n_taxa" => string(length(files)), "sources" => want])
    return (n = length(files), rebuilt = true)
end

"""
    proteome_taxon_map(prot_dir, sample_taxa) -> (acc2tax, n_entries, n_other_taxid, n_no_taxid)

Map every sampled proteome entry to the TaxID its own OX= token names, and count the
entries whose OX differs from the taxon they were downloaded under. The count is what
turns "a proteome may hold entries of another TaxID" from an assumption into a measured
number. A header with no OX= token has no TaxID of its own, so it falls back to the
downloaded taxon — the one place the download decides — and `n_no_taxid` counts those
headers, so a proteome whose headers lost the token cannot pass as one that has none of
another organism.
"""
function proteome_taxon_map(prot_dir::AbstractString, sample_taxa)
    acc2tax = Dict{String,Int}()
    n, other, no_ox = 0, 0, 0
    for t in sample_taxa, line in eachline(proteome_file(prot_dir, t))
        startswith(line, '>') || continue
        ox = header_taxid(line)
        n += 1
        ox === nothing && (no_ox += 1)
        ox === nothing || ox == t || (other += 1)
        acc2tax[uniprot_accession(SubString(line, 2))] = ox === nothing ? Int(t) : ox
    end
    return (acc2tax = acc2tax, n_entries = n, n_other_taxid = other, n_no_taxid = no_ox)
end

"""
    proteome_release_audit(prot_dir, sample_taxa, release) -> NamedTuple

Which UniProtKB release each sampled proteome download came from, read from the
`.release` record written beside it, against `release`, the search space's: the counts
matching, recorded-but-other, and unrecorded, plus the distinct releases seen. A recorded
release other than the search space's means the recall compares two snapshots, and a
"prefilter miss" could then be an entry the newer release added; the acceptance stage
gates on it. A download without a record predates the records and cannot be checked
either way, so it is counted and disclosed, never counted as matching.

`n_no_reference_proteome` counts the recorded downloads whose record names no proteome:
UniProt listed no reference proteome for the TaxID, so the download is every UniProtKB
entry of that TaxID instead. A record without the field cannot vouch for a reference
proteome either, so it is counted the same way.
"""
function proteome_release_audit(prot_dir::AbstractString, sample_taxa, release::AbstractString)
    matching, other, unrecorded, no_ref = 0, 0, 0, 0
    seen = Set{String}()
    for t in sample_taxa
        rec = read_keyvals(release_record_file(prot_dir, t))
        r = get(rec, "X-UniProt-Release", "")
        if isempty(r)
            unrecorded += 1
        else
            push!(seen, r)
            !isempty(release) && r == release ? (matching += 1) : (other += 1)
            isempty(get(rec, "proteome", "")) && (no_ref += 1)
        end
    end
    return (n_matching = matching, n_other = other, n_unrecorded = unrecorded,
            n_recorded = matching + other, n_no_reference_proteome = no_ref,
            releases = sort(collect(seen)))
end

# ---------------------------------------------------------------------------
# Network helpers.
# ---------------------------------------------------------------------------

# UniProt gzips responses even when identity encoding is requested; without the
# defensive gunzip every parse fails silently and the scan comes back empty. Same guard
# as src/validate_classifier.jl and scripts/prokaryote_background.jl.
function fetch(url::AbstractString; timeout = 900, tries = 5)
    for attempt in 1:tries
        try
            r = HTTP.get(url, ["Accept-Encoding" => "identity"];
                         status_exception = false, readtimeout = timeout)
            if r.status != 200
                @warn "HTTP $(r.status)" attempt
                sleep(5 * attempt)
                continue
            end
            b = r.body
            body = (length(b) >= 2 && b[1] == 0x1f && b[2] == 0x8b) ?
                   String(transcode(GzipDecompressor, b)) : String(b)
            return (body = body, resp = r)
        catch e
            @warn "fetch failed" attempt exception = e
            sleep(5 * attempt)
        end
    end
    error("giving up after $tries attempts on $url")
end

function header_value(resp, key::AbstractString)
    for (k, v) in resp.headers
        lowercase(k) == lowercase(key) && return String(v)
    end
    return ""
end

# The `rel="next"` cursor URL, or "" at the last page.
function next_link(resp)
    for (k, v) in resp.headers
        lowercase(k) == "link" || continue
        m = match(r"<([^>]+)>;\s*rel=\"next\"", String(v))
        m === nothing || return String(m.captures[1])
    end
    return ""
end

search_url(query::AbstractString, fmt::AbstractString; fields = "") =
    "https://rest.uniprot.org/uniprotkb/search?query=" * HTTP.escapeuri(query) *
    "&format=$fmt&size=$PAGE" * (isempty(fields) ? "" : "&fields=$fields")

"""
    download_pages(dir, first_url, ext) -> Vector{String}

Page through a UniProt REST search, one file per page under `dir`.

Resumability is the point. Each page body is written to `.part` and renamed only once it
is complete, and the cursor for the NEXT page is written only after that rename. A
restart therefore resumes at the first page that has no cursor file, and a half-written
page can never be mistaken for a finished one.
"""
function download_pages(dir::AbstractString, first_url::AbstractString, ext::AbstractString)
    mkpath(dir)
    page, url = 1, first_url
    while true
        f = joinpath(dir, @sprintf("page_%04d.%s", page, ext))
        c = joinpath(dir, @sprintf("page_%04d.cursor", page))
        if isfile(f) && isfile(c)
            url = read(c, String)
            isempty(url) && break
            page += 1
            continue
        end
        r = fetch(url)
        write(f * ".part", r.body)
        mv(f * ".part", f; force = true)
        nxt = next_link(r.resp)
        write(c, nxt)
        page % 20 == 0 && @printf("    page %4d\n", page)
        isempty(nxt) && break
        url, page = nxt, page + 1
    end
    return sort(filter(f -> endswith(f, "." * ext), readdir(dir; join = true)))
end

# ---------------------------------------------------------------------------
# Driver. Everything above is offline and side-effect free apart from the temp files the
# tests hand it; everything below downloads, runs HMMER and writes the outputs. The guard
# lets test/test_background_scan.jl include this file without touching the network.
# ---------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__

include(joinpath(REPO, "src", "validate_classifier.jl"))
using .ValidateClassifier: win_to_wsl, run_hmmalign, reorder_references!, fetch_control_sequences
using .ValidateClassifier.ModuleHmmer: calibrate_active_site, parse_hmmalign_results

# HOW HMMER IS CALLED, ON EACH PLATFORM. HMMER has no native Windows build, so on Windows
# every invocation is prefixed with `wsl` and every path is translated to its /mnt/<drive>
# form. On Linux and macOS the binary is called directly with the paths as they are.
# `win_to_wsl` returns a non-Windows path unchanged, so one helper covers both, and the
# option line is identical either way — the two arms of this scan are compared to each
# other, so nothing about the invocation may depend on where it ran.
const USE_WSL = Sys.iswindows()
hmmer_path(p::AbstractString) = USE_WSL ? win_to_wsl(p) : String(p)
hmmer_argv(args::AbstractString...) = USE_WSL ? String["wsl", args...] : String[args...]

"""
    run_hmmsearch(fasta, tag; threshold) -> DataFrame

Run hmmsearch (through WSL on Windows, directly elsewhere) and parse the table. Every
call uses the SAME invocation — same model, same options, only the sequence file
differs — because the calibration arm and the background arm are compared to each other
and any option that changes what is reported would put the two on different scales.

`-T` is not optional here: the model carries no GA/TC/NC cutoffs, so without it
hmmsearch falls back to its own E-value default, which depends on the size of the
database being searched. Two databases of different sizes would then be scanned at two
different effective stringencies.

A cached table is reused only if HMMER's own option line says it was produced at this
`-T`. The file name does not carry the threshold, so a re-run at a different cutoff
would otherwise silently reuse the old one and report the old scores as new ones.
"""
function run_hmmsearch(fasta::AbstractString, tag::AbstractString; threshold::Real)
    tbl = joinpath(SCAN_DIR, "$(tag).tblout")
    dom = joinpath(SCAN_DIR, "$(tag).domtblout")
    if isfile(tbl)
        cached = tblout_threshold(tbl)
        if cached === nothing || abs(cached - threshold) > 1e-6
            @warn "$tag: cached table was not produced at this -T, re-running" cached requested = Float64(threshold)
            rm(tbl; force = true)
            rm(dom; force = true)
        end
    end
    if !isfile(tbl)
        thr = @sprintf("%.4f", threshold)
        argv = hmmer_argv("hmmsearch", "--cpu", string(HMM_CPU), "-T", thr, "--incT", thr,
                          "--tblout", hmmer_path(tbl * ".part"), "--domtblout", hmmer_path(dom),
                          "-o", "/dev/null", hmmer_path(HMM), hmmer_path(fasta))
        println("* ", join(argv, " "))
        run(Cmd(argv))
        mv(tbl * ".part", tbl; force = true)   # a killed run leaves no usable .tblout
    end
    return parse_tblout(tbl)
end

"""
    drop_search_outputs(tag)

Delete a hmmsearch table whose input FASTA has just been rebuilt. Keeping it would mean
scanning the new sequence set only in name.
"""
function drop_search_outputs(tag::AbstractString)
    for f in (joinpath(SCAN_DIR, "$(tag).tblout"), joinpath(SCAN_DIR, "$(tag).domtblout"))
        isfile(f) || continue
        @warn "$tag: input FASTA changed, dropping the cached search output" file = basename(f)
        rm(f; force = true)
    end
end

# The superfamily download predates the union and keeps its original file names: its
# release record is the one the downstream habitat summary reads, and re-paging it only
# to rename directories would cost time and buy nothing.
release_file(entry) = joinpath(SCAN_DIR, entry == IPR_SUPERFAMILY ? "uniprot_release.txt" :
                                             "uniprot_release_$(entry).txt")
page_dir(entry, kind) = joinpath(SCAN_DIR, entry == IPR_SUPERFAMILY ? "family_$(kind)" :
                                               "$(lowercase(entry))_$(kind)")
entry_query(entry) = "(xref:interpro-$entry) AND ((taxonomy_id:2) OR (taxonomy_id:2157))"

"""
    download_entry(entry) -> NamedTuple

Page one InterPro entry's prokaryotic UniProtKB members to disk (TSV and FASTA), record
the query, release and announced total the first time, and verify the pages against that
total. The release record is written once, at download time, so it describes the pages
on disk rather than whatever release is live when the script is re-run.
"""
function download_entry(entry::AbstractString)
    query = entry_query(entry)
    relf = release_file(entry)
    if !isfile(relf)
        probe = fetch(search_url(query, "tsv"; fields = "accession") * "&size=1")
        rel = header_value(probe.resp, "X-UniProt-Release")
        write_keyvals(relf, [
            "X-UniProt-Release" => rel,
            "X-UniProt-Release-Date" => header_value(probe.resp, "X-UniProt-Release-Date"),
            "X-Total-Results" => header_value(probe.resp, "X-Total-Results"),
            "interpro_entry" => entry,
            "query" => query,
            "fetched_at" => string(now()),
        ])
    end
    rec = read_keyvals(relf)
    get(rec, "query", query) == query ||
        error("$(basename(relf)) records a different query than this run asks for: " *
              "$(rec["query"]) vs $query")
    total = tryparse(Int, get(rec, "X-Total-Results", ""))
    total === nothing && error("$(basename(relf)) has no usable X-Total-Results")
    println("  downloading $entry annotations ...")
    tsv_files = download_pages(page_dir(entry, "tsv"), search_url(query, "tsv";
                               fields = "accession,organism_id,length,xref_interpro"), "tsv")
    println("  downloading $entry sequences ...")
    fa_files = download_pages(page_dir(entry, "fasta"), search_url(query, "fasta"), "fasta")
    rows = parse_entry_tsv(tsv_files)
    verify_download(entry, rows, fa_files, total)
    @printf("  %s: %d entries (X-Total-Results %d), UniProtKB release %s\n",
            entry, nrow(rows), total, get(rec, "X-UniProt-Release", ""))
    return (entry = entry, query = query, release = get(rec, "X-UniProt-Release", ""),
            release_date = get(rec, "X-UniProt-Release-Date", ""), total = total,
            rows = rows, fa_files = fa_files)
end

DataInputs.require_or_exit("profile_hmm")
println("=== background PQQ-ADH scan ===")

# ---- 1. the frame: reference-proteome taxa that are not census taxa ----------------
# Both inputs are named before they are parsed: a missing one is the registry's own
# refusal, naming the input and where the layout puts it, not a CSV parse error and a
# stack trace.
DataInputs.require_or_exit("census_habitat", "pool_taxa")
# Created only once every input has been found, so a refused run leaves nothing behind.
mkpath(SCAN_DIR)
census = CSV.read(CENSUS_CSV, DataFrame; types = Dict(:TaxID => String, :accession => String))
census_taxids = Set{Int}()
for t in skipmissing(census.TaxID)
    v = tryparse(Int, String(strip(t)))
    v === nothing || push!(census_taxids, v)
end
census_accs = Set{String}(String(strip(a)) for a in skipmissing(census.accession))
pool = CSV.read(POOL_CSV, DataFrame)
frame = frame_taxids(pool.taxid, census_taxids)
frame_set = Set{Int}(frame)
@printf("pool %d taxa | census %d TaxIDs, %d accessions | frame %d taxa\n",
        nrow(pool), length(census_taxids), length(census_accs), length(frame))

# The recall sample depends only on the frame, the seed and --recall-n, so it is drawn and
# checked against the earlier sample HERE, before an hour of downloads: a sample that does
# not extend the earlier one is reported at the start, not discovered at the end. The
# verdict itself is an acceptance row, so it is persisted with everything else.
const SAMPLE = recall_sample(frame, OPT.recall_n, RECALL_SEED)
const RECALL_REF_GIVEN = !isempty(OPT.recall_reference)
const RECALL_REF = RECALL_REF_GIVEN ? OPT.recall_reference : default_recall_reference(OUT_CSV)
const EARLIER_SAMPLE = csv_recall_sample(RECALL_REF; required = RECALL_REF_GIVEN)
const EXTENDS = sample_extends(EARLIER_SAMPLE, SAMPLE.draw_order)
if EXTENDS.ok === nothing
    @printf("recall sample: %d taxa (seed %d); no earlier sample to compare\n",
            length(SAMPLE.taxa), RECALL_SEED)
else
    @printf("recall sample: %d taxa (seed %d); the earlier sample of %d (%s) is %s its first %d draws\n",
            length(SAMPLE.taxa), RECALL_SEED, EXTENDS.n_earlier, basename(RECALL_REF),
            EXTENDS.ok ? "exactly" : "NOT", EXTENDS.n_earlier)
    EXTENDS.ok || @warn "the recall sample does not extend the earlier one; the acceptance row will be a MISS" n_matched = EXTENDS.n_matched n_earlier = EXTENDS.n_earlier
end

# ---- 2. the search space: union of the superfamily and family downloads ------------
# One union serves the background candidates, the calibration set AND the whole-frame
# family membership, so all three see entries retrieved the same way from one release.
sf  = download_entry(IPR_SUPERFAMILY)
fml = download_entry(IPR_FAMILY)
sf.release == fml.release && !isempty(sf.release) ||
    error("the two downloads come from different UniProtKB releases " *
          "($IPR_SUPERFAMILY: '$(sf.release)', $IPR_FAMILY: '$(fml.release)'); " *
          "delete the older one's page directories and release record and re-run")

# Superfamily rows first, so an accession present in both keeps the row the earlier,
# superfamily-only runs parsed; the two rows are the same UniProtKB record anyway.
fam = vcat(sf.rows, fml.rows)
unique!(fam, :accession)
fa_files = vcat(sf.fa_files, fml.fa_files)
fam_is_family = [occursin(IPR_FAMILY, s) for s in fam.interpro]
fam_outside   = fam_is_family .& [!occursin(IPR_SUPERFAMILY, s) for s in fam.interpro]
@printf("search space: %d entries over %d taxa (%d superfamily + %d family, %d in both)\n",
        nrow(fam), length(unique(fam.taxid)), nrow(sf.rows), nrow(fml.rows),
        nrow(sf.rows) + nrow(fml.rows) - nrow(fam))
@printf("  family entries outside the superfamily: %d over %d taxa (%d frame taxa)\n",
        count(fam_outside), length(unique(fam.taxid[fam_outside])),
        length(unique([t for t in fam.taxid[fam_outside] if t in frame_set])))

# ---- 3. whole-frame family membership, independent of hits and of the prefilter ----
fam_tax = ipr017512_by_taxon(fam, frame)
n_family_taxa = count(fam_tax.in_ipr017512_ge500)
@printf("frame taxa with a family entry: %d (%d at >= %d aa, %d with an entry the superfamily prefilter would miss)\n",
        count(fam_tax.in_ipr017512), n_family_taxa, MIN_FAMILY_LEN,
        count(fam_tax.ipr017512_outside_prefilter))

# ---- 4. calibration: what cutoff does the census itself imply? ---------------------
cal_accs  = Set{String}(a for a in fam.accession if a in census_accs)
cand_rows = fam[[t in frame_set for t in fam.taxid], :]
cand_accs = Set{String}(cand_rows.accession)
n_cand = Dict{Int,Int}()
for t in cand_rows.taxid
    n_cand[t] = get(n_cand, t, 0) + 1
end
@printf("calibration set: %d census accessions still active and in the search space\n", length(cal_accs))
@printf("candidates:      %d sequences from %d frame taxa\n", nrow(cand_rows), length(n_cand))

cal_fa  = joinpath(SCAN_DIR, "calibration.fasta")
cand_fa = joinpath(SCAN_DIR, "candidates.fasta")
sync_subset_fasta(cal_fa, fa_files, cal_accs; label = "calibration").rebuilt &&
    drop_search_outputs("calibration")
sync_subset_fasta(cand_fa, fa_files, cand_accs; label = "candidates").rebuilt &&
    drop_search_outputs("candidates")

# -T 0 reports every calibration sequence that scores at all: a reporting cutoff chosen
# before the scores are seen would define its own answer. (--max switches HMMER's heuristic
# filters off; re-running with it and comparing the row count and the low tail shows
# whether the filters discard any of the calibration set, and the arms stay on one scale
# only if they do not.)
cal_hits = run_hmmsearch(cal_fa, "calibration"; threshold = 0.0)
cal_best = combine(groupby(cal_hits, :accession), :full_bits => maximum => :bits)
scores = cal_best.bits
thr_min = calibration_threshold(scores; rule = :min)
thr_p1  = calibration_threshold(scores; rule = :p1)
@printf("calibration: %d of %d accessions scored | min %.1f | 1st pct %.1f | median %.1f | max %.1f\n",
        length(scores), length(cal_accs), thr_min, thr_p1, median(scores), maximum(scores))
@printf("  recovery at min      %.2f %%\n", 100 * recovery(scores, thr_min, length(cal_accs)))
@printf("  recovery at 1st pct  %.2f %%\n", 100 * recovery(scores, thr_p1, length(cal_accs)))

# Which part of the census sets the cutoff? The census keeps rows whose alignment never
# reaches the active-site column ("Gap") alongside rows with a read motif, and those two
# groups do not score alike. Printing the split stops the single number above from being
# read as "this is what a PQQ-ADH scores".
cls_of = Dict{String,String}()
for r in eachrow(census)
    a = ismissing(r.accession) ? "" : String(strip(r.accession))
    isempty(a) || (cls_of[a] = ismissing(r.classification) ? "" : String(r.classification))
end
println("  calibration scores by census class:")
cls_scores = Dict{String,Vector{Float64}}()
for cls in ["Lanthanide", "Calcium", "Gap", "Unknown"]
    s = sort([b for (a, b) in zip(cal_best.accession, cal_best.bits) if get(cls_of, String(a), "") == cls])
    isempty(s) && continue
    cls_scores[cls] = s
    @printf("    %-11s n = %5d   min %7.1f   1st pct %7.1f   median %7.1f   max %7.1f\n",
            cls, length(s), first(s), quantile(s, 0.01), median(s), last(s))
end
# Kept for the sweep below: the cutoff the census implies if only its rows with a READ
# active-site motif are treated as the thing to match. It is not the operating rule —
# that would quietly redefine "carries a PQQ-ADH" as "carries one the classifier could
# read" — but the two are far enough apart that the difference has to be visible.
const THR_LN = haskey(cls_scores, "Lanthanide") ? quantile(cls_scores["Lanthanide"], 0.01) : NaN

# The 1st percentile is the operating rule: the minimum is decided by one sequence, and
# whichever census accession happens to be the most truncated fragment should not get to
# set a method parameter on its own. Whether that matters is an empirical question, not
# an assumption — the two numbers are printed above and the sweep below shows how much
# of the background sits between them. A recorded threshold given with --threshold
# replaces the rule outright; the calibration still runs, so its recovery is reported
# against whichever cutoff is in force.
const THRESHOLD_SOURCE = OPT.threshold === nothing ? "calibration" : "cli"
const THRESHOLD = OPT.threshold === nothing ? thr_p1 : Float64(OPT.threshold)
# Both searches run at the lower of the calibration minimum and the operating cutoff and
# are filtered afterwards, so every sequence at or above the cutoff is reported whichever
# of the two rules is in force, and one scan serves the sweep as well.
const SEARCH_FLOOR = min(thr_min, THRESHOLD)
const REC = recovery(scores, THRESHOLD, length(cal_accs))
@printf("threshold: %.2f bits (%s), recovering %.2f %% of the calibration set\n",
        THRESHOLD, THRESHOLD_SOURCE == "cli" ? "supplied with --threshold" :
                   "1st-percentile rule", 100 * REC)
REC >= MIN_RECOVERY || @warn "calibration recovery below the criterion" REC MIN_RECOVERY
CSV.write(joinpath(SCAN_DIR, "calibration_scores.csv"), cal_best)

# ---- 5. scan the candidates --------------------------------------------------------
# Searched at SEARCH_FLOOR and filtered afterwards, so both rules — and a sensitivity
# sweep past them — come out of one scan instead of one scan per cutoff.
all_hits = run_hmmsearch(cand_fa, "candidates"; threshold = SEARCH_FLOOR)
all_hits = leftjoin(all_hits, select(cand_rows, :accession, :taxid, :interpro), on = :accession)
dropmissing!(all_hits, :taxid)
hits = all_hits[all_hits.full_bits .>= THRESHOLD, :]
@printf("hits: %d sequences from %d taxa at or above %.2f bits\n",
        nrow(hits), length(unique(hits.taxid)), THRESHOLD)
# The part of the search space only the union reaches: frame entries of the family that
# lack the superfamily annotation, and how many of them are hits. This is the measured
# size of the hole the union closes.
const OUTSIDE_ACCS = Set{String}(a for (a, t, o) in zip(fam.accession, fam.taxid, fam_outside)
                                 if o && t in frame_set)
const N_OUTSIDE_HITS = count(a -> a in OUTSIDE_ACCS, hits.accession)
@printf("  of those, %d are among the %d frame family entries outside the superfamily\n",
        N_OUTSIDE_HITS, length(OUTSIDE_ACCS))

# How hard does the answer lean on the cutoff? The two calibration rules can land close
# together, so the rows past them are there to show where the background would
# start changing if a recorded census threshold ever replaces the empirical one.
println("sensitivity of the flagged count to the cutoff:")
for cut in sort(unique(Float64[thr_min, THRESHOLD, 50.0, 100.0, THR_LN, 200.0, 300.0]))
    isnan(cut) && continue
    sub = all_hits[all_hits.full_bits .>= cut, :]
    note = cut == THRESHOLD ? "   <- operating point" :
           cut == THR_LN    ? "   <- Ln-class 1st percentile" : ""
    @printf("  >= %7.2f bits   %5d sequences   %4d taxa%s\n", cut, nrow(sub),
            length(unique(sub.taxid)), note)
end

# ---- 6. class each hit at the calibrated active-site column ------------------------
# Reuses the production calibration/parsing path, so a hit is classed exactly the way a
# census row is: hmmalign to the model, force the Ln/Ca reference pair to the top, read
# the collision column, then read the motif at that column. The column is found in model
# coordinates (insert states stripped), so adding sequences to the alignment cannot move
# it; it is written to threshold.txt with the rest of the method record.
classes = Dict{String,String}()
active_site_col = 0
if nrow(hits) > 0
    refs_fa   = joinpath(SCAN_DIR, "refs.fasta")
    hits_only = joinpath(SCAN_DIR, "hits_only.fasta")
    hit_fa    = joinpath(SCAN_DIR, "hits_with_refs.fasta")
    afa       = joinpath(SCAN_DIR, "hits_with_refs.afa")
    # A cached reference FASTA built on a different pair would class every hit against
    # the wrong collision column while looking exactly like the right file.
    if fasta_accessions(refs_fa) != Set([LN_REF, CA_REF])
        isfile(refs_fa) && @warn "cached reference FASTA holds a different pair, refetching" pair = (LN_REF, CA_REF)
        fetch_control_sequences([LN_REF, CA_REF]; out_fasta = refs_fa)
    end
    # The alignment must hold exactly the pair plus the current hits. A stale alignment
    # missing a new hit would leave that hit "Unscored", which an Ln-type exclusion reads
    # as "not Ln" without any error.
    hit_accs = Set{String}(hits.accession)
    if fasta_accessions(hit_fa) != union(Set([LN_REF, CA_REF]), hit_accs)
        isfile(hit_fa) && @warn "cached hit FASTA does not hold the reference pair plus the current hits, rebuilding"
        filter_fasta(fa_files, hit_accs, hits_only)
        open(hit_fa * ".part", "w") do io
            write(io, read(refs_fa, String))
            write(io, read(hits_only, String))
        end
        mv(hit_fa * ".part", hit_fa; force = true)
        rm(afa; force = true)
    end
    isfile(afa) || run_hmmalign(HMM, hit_fa; wsl = USE_WSL, out_afa = afa)
    reorder_references!(afa, LN_REF, CA_REF)
    active_site_col = calibrate_active_site(afa)
    @printf("active-site column: %d\n", active_site_col)
    for r in eachrow(parse_hmmalign_results(afa, active_site_col))
        classes[String(r.accession)] = String(r.classification)
    end
end
n_unscored = count(a -> !haskey(classes, String(a)), hits.accession)

# The derived cutoff, the column the hits were classed at and the model they came from
# are one method record. It is assembled here, once the column is known, and written only
# by publish_scan, together with the table it describes — writing it now would put this
# run's cutoff on disk beside whatever table an earlier run left.
const THRESHOLD_RULE = THRESHOLD_SOURCE == "cli" ?
    "supplied on the command line" :
    "1st percentile of the census accessions' full-sequence bit scores"
const THRESHOLD_PAIRS = [
    "rule" => THRESHOLD_RULE,
    "threshold_bits" => @sprintf("%.4f", THRESHOLD),
    "threshold_source" => THRESHOLD_SOURCE,
    "threshold_min_bits" => @sprintf("%.4f", thr_min),
    "search_floor_bits" => @sprintf("%.4f", SEARCH_FLOOR),
    "n_calibration_accessions" => string(length(cal_accs)),
    "n_scored" => string(length(scores)),
    "recovery" => @sprintf("%.6f", REC),
    "min" => @sprintf("%.2f", thr_min),
    "median" => @sprintf("%.2f", median(scores)),
    "max" => @sprintf("%.2f", maximum(scores)),
    "active_site_column" => string(active_site_col),
    "hmm" => basename(HMM),
    "hmm_sha256" => bytes2hex(sha256(read(HMM))),
]

# ---- 7. recall check: is the union prefilter safe? ---------------------------------
sample_taxa = SAMPLE.taxa
sample_set = Set{Int}(sample_taxa)
prot_dir = joinpath(SCAN_DIR, "proteomes")
mkpath(prot_dir)
recall_fa = joinpath(SCAN_DIR, "recall_proteomes.fasta")
for (i, t) in enumerate(sample_taxa)
    f = proteome_file(prot_dir, t)
    isfile(f) && continue
    up = fetch("https://rest.uniprot.org/proteomes/stream?query=" *
               HTTP.escapeuri("(organism_id:$t) AND (reference:true)") *
               "&format=tsv&fields=upid").body
    lines = filter(!isempty, strip.(split(up, '\n')))
    upid = length(lines) >= 2 ? String(strip(lines[2])) : ""
    q = isempty(upid) ? "(organism_id:$t)" : "(proteome:$upid)"
    r = fetch("https://rest.uniprot.org/uniprotkb/stream?query=" *
              HTTP.escapeuri(q) * "&format=fasta")
    write(f * ".part", r.body)
    mv(f * ".part", f; force = true)
    # Which release the proteome came from, recorded at download time: the file is reused
    # on its name from now on, so this is the only moment the release can be known.
    write_keyvals(release_record_file(prot_dir, t), [
        "X-UniProt-Release" => header_value(r.resp, "X-UniProt-Release"),
        "X-UniProt-Release-Date" => header_value(r.resp, "X-UniProt-Release-Date"),
        "proteome" => upid,
        "query" => q,
        "n_records" => string(count(==('>'), r.body)),
        "fetched_at" => string(now()),
    ])
    @printf("  [%3d/%3d] taxid %-9d %-12s %6d records\n",
            i, length(sample_taxa), t, upid, count(==('>'), r.body))
end
sync_recall_fasta(recall_fa, prot_dir, sample_taxa).rebuilt && drop_search_outputs("recall")
const RELEASES = proteome_release_audit(prot_dir, sample_taxa, sf.release)
@printf("recall proteome releases: %d from %s, %d from another release, %d downloaded before releases were recorded\n",
        RELEASES.n_matching, sf.release, RELEASES.n_other, RELEASES.n_unrecorded)
@printf("  %d of the %d recorded downloads had no reference proteome and hold every entry of the TaxID\n",
        RELEASES.n_no_reference_proteome, RELEASES.n_recorded)

recall_all = run_hmmsearch(recall_fa, "recall"; threshold = SEARCH_FLOOR)
recall_hits = recall_all[recall_all.full_bits .>= THRESHOLD, :]
# Every full-proteome hit counts for the organism its entry names (OX=), and only when
# that organism is a sampled taxon — the same exact-TaxID rule the prefilter's hits
# follow. Entries of another TaxID and the hits they carry are counted, not assumed away.
const TAXMAP = proteome_taxon_map(prot_dir, sample_taxa)
const FULL = recall_flagged(recall_hits.accession, TAXMAP.acc2tax, sample_set)
pre_flagged = intersect(Set{Int}(hits.taxid), sample_set)
rs = recall_summary(n_sampled = length(sample_taxa), full_flagged = FULL.flagged,
                    pre_flagged = pre_flagged,
                    n_never_searched = count(t -> get(n_cand, t, 0) == 0, sample_taxa),
                    bar = MIN_RECALL, seed = RECALL_SEED)
@printf("recall: %d of %d sampled taxa flagged by the full proteome scan; %d of those also by the prefilter (%s; one-sided 95 %% lower bound %s)\n",
        rs.n_full_flagged, rs.n_sampled, rs.n_recovered,
        isnan(rs.rate) ? "undefined" : @sprintf("%.1f %%", 100 * rs.rate),
        isnan(rs.lower_bound_95) ? "undefined" : @sprintf("%.4f", rs.lower_bound_95))
# Nothing here compares the bound with the bar. The gate below is on the rate, and the
# statement the bound decides is decided by its consumer alone, so this script holds no
# second copy of that rule to drift from the first.
@printf("  %d sampled taxa had no candidate entry, so only the full scan searched them\n",
        rs.n_never_searched_in_sample)
@printf("  %d of %d proteome entries name a TaxID other than the taxon downloaded; %d of %d hit sequences set aside for it\n",
        TAXMAP.n_other_taxid, TAXMAP.n_entries, FULL.n_outside, nrow(recall_hits))
@printf("  %d proteome entries name no TaxID and are attributed to the taxon downloaded\n",
        TAXMAP.n_no_taxid)
isempty(rs.missed) || @printf("  missed by the prefilter: %s\n", join(rs.missed, ", "))
extra_pre = setdiff(pre_flagged, FULL.flagged)
isempty(extra_pre) || @printf("  flagged by the prefilter only: %s\n", join(sort(collect(extra_pre)), ", "))

# ---- 8. one row per frame taxon ----------------------------------------------------
by_tax = Dict{Int,NamedTuple}()
for g in groupby(hits, :taxid)
    i = argmax(g.full_bits)
    cls = String[get(classes, String(a), "Unscored") for a in g.accession]
    by_tax[g.taxid[1]] = (n = nrow(g), acc = String(g.accession[i]),
                          bits = g.full_bits[i], ev = g.full_evalue[i],
                          cls = join(sort(unique(cls)), ";"),
                          n_ln_or_ca = count(in(LN_OR_CA), cls),
                          hit_family = any(occursin(IPR_FAMILY, s) for s in g.interpro))
end
fam_by_tax = Dict(r.taxid => r for r in eachrow(fam_tax))
out = DataFrame(taxid = Int[], n_candidates = Int[], n_hits = Int[], best_acc = String[],
                best_bits = Union{Missing,Float64}[], best_evalue = Union{Missing,Float64}[],
                hit_classes = String[], in_ipr017512 = Bool[],
                in_prefilter_recall_sample = Bool[],
                n_ipr017512_entries = Int[], n_ipr017512_ge500 = Int[],
                in_ipr017512_ge500 = Bool[], ipr017512_outside_prefilter = Bool[],
                hit_ipr017512 = Bool[], n_hits_ln_or_ca = Int[])
for t in frame
    h = get(by_tax, t, nothing)
    f = fam_by_tax[t]
    push!(out, (t, get(n_cand, t, 0),
                h === nothing ? 0 : h.n,
                h === nothing ? "" : h.acc,
                h === nothing ? missing : h.bits,
                h === nothing ? missing : h.ev,
                h === nothing ? "" : h.cls,
                f.in_ipr017512,
                t in sample_set,
                f.n_ipr017512_entries, f.n_ipr017512_ge500,
                f.in_ipr017512_ge500, f.ipr017512_outside_prefilter,
                h === nothing ? false : h.hit_family,
                h === nothing ? 0 : h.n_ln_or_ca))
end

# ---- 9. acceptance: the questions this step exists to answer -----------------------
# Every gated row carries the numbers its verdict is derived from (accept! derives it), so
# the table cannot say PASS beside numbers that fail and a consumer can check it. A
# printed criterion that nothing reads and that exits 0 is not a criterion; it is a
# comment.
println("\n--- acceptance ---")
accept = acceptance_table()
flagged_taxids = Set{Int}(out.taxid[out.n_hits .> 0])
nflag = length(flagged_taxids)

# Every frame taxon with a full-length member of the methanol/ethanol family must come
# out flagged. The denominator is the whole-frame roll-up, not the prefilter's candidate
# rows: drawn from the search space, it could only contain taxa the search space reached
# and would pass exactly when the search space has a hole.
family_missed = ipr017512_unflagged(fam_tax, flagged_taxids)
accept!(accept, "$IPR_FAMILY frame taxa flagged (>= $MIN_FAMILY_LEN aa, whole-frame membership)";
        value = n_family_taxa - length(family_missed), denominator = n_family_taxa,
        bar = n_family_taxa,
        note = "bar computed by this run: every frame taxon with a family entry at the length floor" *
               (isempty(family_missed) ? "" : "; unflagged: " * join(family_missed, ";")))
accept!(accept, "calibration recovery"; value = REC, denominator = length(cal_accs),
        bar = MIN_RECOVERY,
        note = "fixed bar; share of the census accessions in the search space at or above the cutoff, unscored ones counted as missed")
# The recall: gated on the rate, the lower bound disclosed beside it (recall_acceptance!).
recall_acceptance!(accept, rs)
# The sample extends the earlier run's only if the earlier taxa are its first draws. The
# row's two numbers re-derive the verdict: n_matched >= n_earlier is sample_extends' ok.
accept!(accept, "recall sample extends the earlier sample";
        value = EXTENDS.ok === nothing ? missing : EXTENDS.n_matched,
        denominator = EXTENDS.ok === nothing ? missing : EXTENDS.n_earlier,
        bar = EXTENDS.ok === nothing ? missing : EXTENDS.n_earlier,
        note = EXTENDS.ok === nothing ?
            "no earlier scan table on disk to compare; this run's sample is $(length(sample_taxa)) taxa, seed $RECALL_SEED" :
            "bar computed by this run: all $(EXTENDS.n_earlier) taxa sampled in $(basename(RECALL_REF)) must be the first draws of this run's $(length(sample_taxa)), seed $RECALL_SEED")
# A recorded release other than the search space's fails; a download from before the
# records existed can be neither, so it gets its own disclosed row instead of being folded
# into either side. With no recorded download at all there is nothing to gate on.
accept!(accept, "recall proteome downloads recorded from the search-space release";
        value = RELEASES.n_recorded > 0 ? RELEASES.n_matching : missing,
        denominator = RELEASES.n_recorded > 0 ? RELEASES.n_recorded : missing,
        bar = RELEASES.n_recorded > 0 ? RELEASES.n_recorded : missing,
        note = "bar computed by this run: every download whose release was recorded; search space $(sf.release), recorded $(isempty(RELEASES.releases) ? "none" : join(RELEASES.releases, ";"))")
accept!(accept, "recall proteome downloads without a recorded release";
        value = RELEASES.n_unrecorded, denominator = length(sample_taxa),
        note = "no pre-registered bar, disclosed; downloaded before the release was recorded beside each proteome")
accept!(accept, "recall proteome downloads without a reference proteome";
        value = RELEASES.n_recorded > 0 ? RELEASES.n_no_reference_proteome : missing,
        denominator = RELEASES.n_recorded > 0 ? RELEASES.n_recorded : missing,
        note = "no pre-registered bar, disclosed; every UniProtKB entry of the TaxID was downloaded instead; counted over the downloads whose record exists, the others cannot say")
# The floor on the flagged count is computed, not remembered: every frame taxon carrying
# a full-length family entry is a carrier by annotation alone, so no complete flagged set
# can be smaller than their number. A floor copied from an earlier run would describe
# that run's data, not this one's.
accept!(accept, "flagged background taxa"; value = nflag, denominator = nrow(out),
        bar = n_family_taxa,
        note = "bar computed by this run: the number of frame taxa with a family entry at the length floor, which every flagged set must contain")
accept!(accept, "hit sequences classed at the active-site column";
        value = nrow(hits) - n_unscored, denominator = nrow(hits), bar = nrow(hits),
        note = "bar computed by this run: every hit; an unclassed hit reads as not-Lanthanide downstream without an error")
# WHAT THIS CONTROL ASSERTS, SPELLED OUT, because a bare hit count beside an organism's name
# reads as a claim about the organism. It tests that the scan FINDS SEQUENCES for a taxon
# carrying the enzyme family. It does not test that the organism is methylotrophic, and the
# two come apart here: the per-taxon notes on POSITIVE_CONTROLS above record that 217069
# has published evidence AGAINST growth on methanol as a sole carbon and energy source.
# Read those notes rather than this line - they also carry the limit that matters, which
# is that the 2004 assay was calcium-supplemented with no lanthanides at a date when
# lanthanide-dependent dehydrogenases were unknown, so it could not have detected one.
# Sequence presence is what is measured and it holds; what that evidence rules out is only
# the rationale "this taxon is a control because it grows on methanol", which is not the
# one the list was selected by (public annotation, as the notes above say).
for (t, name) in POSITIVE_CONTROLS
    r = out[out.taxid .== t, :]
    accept!(accept, "positive control $name ($t)";
            value = nrow(r) == 1 ? r.n_hits[1] : missing, bar = 1,
            note = nrow(r) == 1 ?
                   "fixed bar; hit sequences of the control taxon. SEQUENCE PRESENCE ONLY " *
                   "- not that the organism is methylotrophic; see the per-taxon notes on " *
                   "POSITIVE_CONTROLS for the published evidence and its limits" :
                   "fixed bar; the control taxon is absent from the frame")
end

# Disclosed, not gated. The flagged set is bounded from one side by every criterion
# above, so the quantities that would show over-exclusion, and the ones that size the
# union's reach, are written down here instead of being tested against a bar invented
# after the scan.
spec = specificity_disclosure(out.n_hits, out.n_hits_ln_or_ca)
accept!(accept, "flagged taxa with no Lanthanide- or Calcium-classed hit";
        value = spec.n_without_ln_or_ca, denominator = spec.n_flagged,
        note = "no pre-registered bar, disclosed")
accept!(accept, "$IPR_FAMILY entries outside the $IPR_SUPERFAMILY prefilter, frame taxa";
        value = count(fam_tax.ipr017512_outside_prefilter), denominator = nrow(out),
        note = "no pre-registered bar, disclosed; $(length(OUTSIDE_ACCS)) such entries in the frame, $N_OUTSIDE_HITS of them hits")
accept!(accept, "full-proteome hits set aside by the exact-TaxID rule";
        value = FULL.n_outside, denominator = nrow(recall_hits),
        note = "no pre-registered bar, disclosed; $(TAXMAP.n_other_taxid) of $(TAXMAP.n_entries) sampled proteome entries name a TaxID other than the taxon downloaded, $(TAXMAP.n_no_taxid) name none and are attributed to it")

for r in eachrow(accept)
    @printf("%-9s %-66s value %-10s of %-6s bar %s\n", r.verdict, r.check, r.value,
            isempty(r.denominator) ? "-" : r.denominator, isempty(r.bar) ? "-" : r.bar)
end

# One stamp for the run, taken as it publishes: run_id and every record's derived_at.
const RUN_ID = string(now())
pub = publish_scan(; out_csv = OUT_CSV, csv = csv_bytes(out),
                     threshold_path = THRESH_TXT, threshold_pairs = THRESHOLD_PAIRS,
                     recall_path = RECALL_TXT,
                     recall_pairs = recall_keyvals(rs; threshold_bits = THRESHOLD,
                                                   uniprot_release = sf.release),
                     acceptance_path = ACCEPT_TSV, acceptance = accept, run_id = RUN_ID)
@printf("\nrun_id %s  scan_sha256 %s\n", RUN_ID, pub.scan_sha256)
if !pub.promoted
    println(stderr, "\nACCEPTANCE MISS: nothing was promoted. This run's table and records are left as ",
            join(basename.(pub.pending), ", "), ".")
    isempty(pub.staled) ||
        println(stderr, "Renamed so no consumer reads them beside this failing run: ",
                join(basename.(pub.staled), ", "), ".")
    exit(1)
end
@printf("promoted %s (%d rows, %d flagged) with %s, %s and %s\n", OUT_CSV, nrow(out), nflag,
        basename(THRESH_TXT), basename(RECALL_TXT), basename(ACCEPT_TSV))

# ---- 10. closing report -------------------------------------------------------------
has_ln = [occursin("Lanthanide", c) for c in out.hit_classes]
println("\n--- summary ---")
@printf("frame taxa                                  %5d\n", nrow(out))
@printf("  with any hit                              %5d\n", nflag)
@printf("  with a Lanthanide-classed hit             %5d\n", count(has_ln))
@printf("  with a Lanthanide- or Calcium-classed hit %5d\n", count(out.n_hits_ln_or_ca .> 0))
@printf("  never searched (no candidate entry)       %5d\n", count(out.n_candidates .== 0))
@printf("family (%s) frame taxa                 %5d, of which flagged %d\n",
        IPR_FAMILY, n_family_taxa, n_family_taxa - length(family_missed))
sf_frame = Set{Int}(t for t in sf.rows.taxid if t in frame_set)
n_sf_unreached = count(t -> !(t in sf_frame), fam_tax.taxid[fam_tax.in_ipr017512])
@printf("  of those, never reachable by the superfamily prefilter alone %d\n", n_sf_unreached)
kept, removed = exclude_carriers(DataFrame(taxid = frame), out)
@printf("background after exclusion: %d taxa kept, %d removed\n", nrow(kept), nrow(removed))
println("recall record:     ", RECALL_TXT)
println("acceptance record: ", ACCEPT_TSV)

end # driver
