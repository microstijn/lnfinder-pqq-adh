#=
    Offline self-test for the UniProt/UniParc accession → TaxID mapper. No network: every
    request goes through the module's single injected `fetch` seam — the dependency-
    injection pattern ModuleEnvironmentTools uses for `path_of` — and is answered here from
    inline fixtures.

    The failure this suite exists to prevent: a rescue path for deleted entries that reads
    `uni_entry.organisms[1].taxonId`, when a UniParc id-mapping result states its taxonomy
    as `commonTaxons[{topLevel, commonTaxon, commonTaxonId}]` and carries no `organisms`
    array at all. That branch never fires: every deleted accession comes back without a
    taxid, the run continues without a word, and rows with no taxonomy are
    indistinguishable from organisms that genuinely have none. Nothing crashes; the answer
    is quietly wrong — which is why the rescue has a floor (`check_rescue`) that throws
    instead of returning a shrunken map.

    The fixtures below mirror the KEY SHAPES of the live payloads (UniProtKB search,
    UniParc upis, UniParc id-mapping, and /uniparc/<upi>/databases). Every value in them is
    synthetic: invented accessions, UniParc ids and taxids, so the suite asserts on
    structure and never on data.

    Run:  julia --project=. test/test_accession_mapper.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
using JSON3
using DataFrames
include(joinpath(@__DIR__, "..", "src", "ModuleAccessionMapper.jl"))
using .AccessionMapper
const AM = AccessionMapper

# ---------------------------------------------------------------------------
# Fixtures — live key shapes, synthetic values
# ---------------------------------------------------------------------------
# UniParc id-mapping (from=UniParc, to=UniParc). Note commonTaxons/commonTaxonId, and note
# that there is no `organisms` key anywhere — that is the shape a taxonId rescue expects.
const UNIPARC_PAYLOAD = """
{"results":[
  {"from":"UPI0000000001","to":{
     "uniParcId":"UPI0000000001","crossReferenceCount":3,
     "commonTaxons":[{"topLevel":"cellular organisms","commonTaxon":"Synthetica",
                      "commonTaxonId":901}],
     "uniProtKBAccessions":["A0DEL001"],
     "sequence":{"value":"MKVAAA","length":6,"molWeight":600,
                 "crc64":"0000000000000001","md5":"00000000000000000000000000000001"}}},
  {"from":"UPI0000000002","to":{
     "uniParcId":"UPI0000000002","crossReferenceCount":2,
     "commonTaxons":[{"topLevel":"cellular organisms","commonTaxon":"Synthetica beta",
                      "commonTaxonId":333}],
     "uniProtKBAccessions":["A0DEL002"],
     "sequence":{"value":"MKVBBB","length":6,"molWeight":600,
                 "crc64":"0000000000000002","md5":"00000000000000000000000000000002"}}}
]}
"""

# /uniparc/UPI0000000001/databases — three cross-references to three different organisms.
# The one that IS the rescued accession is inactive, which is what a deleted entry looks
# like, so picking "the active one" would pick the wrong organism.
const DB_PAYLOAD_ONE = """
{"results":[
  {"database":"UniProtKB/TrEMBL","id":"A0OTHER01","active":true,
   "created":"2001-01-01","lastUpdated":"2002-02-02",
   "organism":{"scientificName":"Aliena prima","taxonId":777}},
  {"database":"UniProtKB/TrEMBL","id":"A0DEL001","active":false,
   "created":"2003-03-03","lastUpdated":"2004-04-04",
   "organism":{"scientificName":"Synthetica exempla","taxonId":222}},
  {"database":"UniProtKB/Swiss-Prot","id":"A0OTHER02","active":true,
   "created":"2005-05-05","lastUpdated":"2006-06-06",
   "organism":{"scientificName":"Aliena secunda","taxonId":888}}
]}
"""

# /uniparc/UPI0000000002/databases — the rescued accession is not among the
# cross-references, so the common taxon has to answer for it.
const DB_PAYLOAD_TWO = """
{"results":[
  {"database":"UniProtKB/TrEMBL","id":"A0OTHER03","active":true,
   "created":"2007-07-07","lastUpdated":"2008-08-08",
   "organism":{"scientificName":"Aliena tertia","taxonId":999}}
]}
"""

# Small single-entry UniParc payloads for the commonTaxons cases.
uniparc_entry(upi, taxons) = """
{"results":[{"from":"$upi","to":{"uniParcId":"$upi","crossReferenceCount":1,
  "commonTaxons":$taxons,"uniProtKBAccessions":["A0TEST001"],
  "sequence":{"value":"MKV","length":3,"molWeight":300,
              "crc64":"0000000000000003","md5":"00000000000000000000000000000003"}}}]}
"""

const CELLULAR_ONE = """[{"topLevel":"cellular organisms","commonTaxon":"Synthetica","commonTaxonId":901}]"""
const METAGENOME_AND_CELLULAR = """[
  {"topLevel":"unclassified entries","commonTaxon":"metagenomes","commonTaxonId":990001},
  {"topLevel":"cellular organisms","commonTaxon":"Synthetica","commonTaxonId":901}]"""
const TWO_CELLULAR = """[
  {"topLevel":"cellular organisms","commonTaxon":"Synthetica","commonTaxonId":901},
  {"topLevel":"cellular organisms","commonTaxon":"Aliena","commonTaxonId":902}]"""
const NO_TAXONS = "[]"

# The same two UniParc ids, but with candidates no rule can decide between.
const AMBIGUOUS_UNIPARC_PAYLOAD = """
{"results":[
  {"from":"UPI0000000001","to":{"uniParcId":"UPI0000000001","crossReferenceCount":3,
     "commonTaxons":$TWO_CELLULAR,"uniProtKBAccessions":["A0DEL001"]}},
  {"from":"UPI0000000002","to":{"uniParcId":"UPI0000000002","crossReferenceCount":2,
     "commonTaxons":$TWO_CELLULAR,"uniProtKBAccessions":["A0DEL002"]}}
]}
"""

# The pre-2026 UniParc shape, kept only so an archived payload still parses.
const LEGACY_PAYLOAD = """
{"results":[{"from":"UPI0000000050","to":{
  "uniParcId":"UPI0000000050",
  "organisms":[{"taxonId":555,"name":"Synthetica legacy"}],
  "uniProtKBAccessions":["A0LEGACY1"]}}]}
"""

# ---------------------------------------------------------------------------
# T2 — commonTaxons
# ---------------------------------------------------------------------------
@testset "T2 — commonTaxons drive the rescue, in every shape the service returns" begin
    pick(taxons) = choose_taxid(nothing,
                                parse_uniparc_common(JSON3.read(uniparc_entry("UPI0000000101", taxons)))["UPI0000000101"])

    @testset "a single cellular common taxon answers" begin
        hits = parse_uniparc_common(JSON3.read(uniparc_entry("UPI0000000101", CELLULAR_ONE)))
        @test length(hits["UPI0000000101"]) == 1
        @test hits["UPI0000000101"][1].taxid == 901
        @test hits["UPI0000000101"][1].toplevel == "cellular organisms"
        c = pick(CELLULAR_ONE)
        @test c.taxid == 901 && c.source == "uniparc_common_taxon"
        @test ismissing(c.alt_cellular_taxid)
    end

    @testset "a metagenome bucket alongside a cellular taxon does not win" begin
        c = pick(METAGENOME_AND_CELLULAR)
        @test c.taxid == 901
        @test c.source == "uniparc_common_taxon"
    end

    @testset "two cellular taxa are ambiguous: no common pick" begin
        c = pick(TWO_CELLULAR)
        @test ismissing(c.taxid) && c.source == "unresolved"
    end

    @testset "an empty commonTaxons list is unresolved, not an error" begin
        hits = parse_uniparc_common(JSON3.read(uniparc_entry("UPI0000000101", NO_TAXONS)))
        @test hits["UPI0000000101"] == AM.TaxonHit[]
        c = pick(NO_TAXONS)
        @test ismissing(c.taxid) && c.source == "unresolved"
    end

    @testset "the legacy `organisms` shape is still parsed" begin
        hits = parse_uniparc_common(JSON3.read(LEGACY_PAYLOAD))
        @test hits["UPI0000000050"][1].taxid == 555
        # No top-level branch is stated in that shape, so it must not be mistaken for a
        # metagenome bucket.
        @test hits["UPI0000000050"][1].toplevel == ""
        c = choose_taxid(nothing, hits["UPI0000000050"])
        @test c.taxid == 555 && c.source == "uniparc_common_taxon"
    end

    @testset "the uniParcCrossReferences shape is accepted as well" begin
        payload = """
        {"results":[{"from":"UPI0000000060","to":{"uniParcId":"UPI0000000060",
          "uniParcCrossReferences":[{"database":"UniProtKB/TrEMBL","id":"A0XREF001",
            "organism":{"scientificName":"Synthetica xref","taxonId":666}}]}}]}
        """
        hits = parse_uniparc_common(JSON3.read(payload))
        @test hits["UPI0000000060"][1].taxid == 666
    end
end

# ---------------------------------------------------------------------------
# T3 — the entry's own cross-reference
# ---------------------------------------------------------------------------
@testset "T3 — own_xref picks the accession's own organism, active or not" begin
    db = JSON3.read(DB_PAYLOAD_ONE)
    hit = own_xref(db, "A0DEL001")

    @test hit !== nothing
    @test hit.taxid == 222                    # not 777 (first listed) and not 888
    @test hit.active == false                 # a deleted accession is inactive by definition
    @test hit.created == "2003-03-03"
    @test hit.last_updated == "2004-04-04"
    @test hit.organism == "Synthetica exempla"

    # Nothing to claim when the entry holds no cross-reference for the accession: the
    # caller must fall through to the common taxon rather than take a neighbour's organism.
    @test own_xref(JSON3.read(DB_PAYLOAD_TWO), "A0DEL002") === nothing
    @test own_xref(JSON3.read(DB_PAYLOAD_ONE), "A0ABSENT9") === nothing

    # Where the same accession appears twice, the live cross-reference is preferred.
    twice = """
    {"results":[
      {"database":"UniProtKB/TrEMBL","id":"A0DUP0001","active":false,
       "created":"2001-01-01","lastUpdated":"2002-02-02",
       "organism":{"scientificName":"Synthetica stale","taxonId":10}},
      {"database":"UniProtKB/Swiss-Prot","id":"A0DUP0001","active":true,
       "created":"2003-03-03","lastUpdated":"2004-04-04",
       "organism":{"scientificName":"Synthetica live","taxonId":20}}]}
    """
    @test own_xref(JSON3.read(twice), "A0DUP0001").taxid == 20
end

# ---------------------------------------------------------------------------
# T4 — the own-cross-reference precedence rule
# ---------------------------------------------------------------------------
@testset "T4 — own cross-reference wins; a metagenome taxid is kept, not swapped" begin
    common = parse_uniparc_common(JSON3.read(uniparc_entry("UPI0000000101",
                                                           METAGENOME_AND_CELLULAR)))["UPI0000000101"]
    mg = (taxid = 990001, active = false, created = "2009-09-09",
          last_updated = "2010-10-10", organism = "soil metagenome")

    c = choose_taxid(mg, common)
    @test c.taxid == 990001                   # kept: it is what the source record says
    @test c.source == "uniparc_own_xref"
    @test c.alt_cellular_taxid == 901         # …with the cellular alternative recorded

    # The same call without the metagenome context leaves the alternative empty rather than
    # inventing one: a cellular organism's own taxid needs no alternative.
    cell = (taxid = 222, active = false, created = "", last_updated = "",
            organism = "Synthetica exempla")
    c2 = choose_taxid(cell, common)
    @test c2.taxid == 222 && ismissing(c2.alt_cellular_taxid)

    # The own cross-reference outranks a common taxon that disagrees with it.
    c3 = choose_taxid(cell, parse_uniparc_common(JSON3.read(uniparc_entry("UPI0000000101",
                                                                         CELLULAR_ONE)))["UPI0000000101"])
    @test c3.taxid == 222 && c3.source == "uniparc_own_xref"

    # No cellular candidate at all: a single non-cellular one answers, and says so.
    only_mg = parse_uniparc_common(JSON3.read(uniparc_entry("UPI0000000101",
        """[{"topLevel":"unclassified entries","commonTaxon":"metagenomes","commonTaxonId":990001}]""")))["UPI0000000101"]
    c4 = choose_taxid(nothing, only_mg)
    @test c4.taxid == 990001 && c4.source == "uniparc_common_noncellular"
end

# ---------------------------------------------------------------------------
# T5 — the rescue floor
# ---------------------------------------------------------------------------
@testset "T5 — check_rescue fails loudly instead of returning a shrunken map" begin
    @test_throws ErrorException check_rescue(2, 1; min_frac = 0.98)
    @test check_rescue(2, 2; min_frac = 0.98)
    @test check_rescue(0, 0; min_frac = 0.98)          # nothing needed rescuing
    @test_throws ErrorException check_rescue(100, 0)   # the shape-change failure mode
    @test check_rescue(100, 99; min_frac = 0.9)

    # The message has to name the counts, or the operator cannot tell a shape change from
    # a handful of genuinely unmappable entries.
    msg = try
        check_rescue(50, 3)
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("3", msg) && occursin("50", msg)
end

# ---------------------------------------------------------------------------
# Regression — the payload shape a `taxonId`-only reader rescues nothing from
# ---------------------------------------------------------------------------
@testset "REGRESSION — commonTaxons-only payloads rescue every UPI" begin
    payload = """
    {"results":[
      {"from":"UPI0000000011","to":{"uniParcId":"UPI0000000011","crossReferenceCount":1,
        "commonTaxons":[{"topLevel":"cellular organisms","commonTaxon":"Synthetica",
                         "commonTaxonId":901}],"uniProtKBAccessions":["A0OLD0001"]}},
      {"from":"UPI0000000012","to":{"uniParcId":"UPI0000000012","crossReferenceCount":1,
        "commonTaxons":[{"topLevel":"cellular organisms","commonTaxon":"Aliena",
                         "commonTaxonId":902}],"uniProtKBAccessions":["A0OLD0002"]}},
      {"from":"UPI0000000013","to":{"uniParcId":"UPI0000000013","crossReferenceCount":1,
        "commonTaxons":[{"topLevel":"cellular organisms","commonTaxon":"Tertia",
                         "commonTaxonId":903}],"uniProtKBAccessions":["A0OLD0003"]}}
    ]}
    """
    json = JSON3.read(payload)
    common = parse_uniparc_common(json)

    # The `organisms` key a taxonId-only reader looks for is absent from every entry.
    @test all(!haskey(item.to, :organisms) for item in json.results)
    @test length(common) == 3

    rescued = count(upi -> !ismissing(choose_taxid(nothing, common[upi]).taxid), keys(common))
    @test rescued == 3
    @test sort([choose_taxid(nothing, common[u]).taxid for u in keys(common)]) == [901, 902, 903]

    # And the dead access pattern must not come back into the module.
    src = read(joinpath(@__DIR__, "..", "src", "ModuleAccessionMapper.jl"), String)
    @test !occursin("organisms[1]", src)
end

# ===================================================================
# The three-pass stateless route
# ===================================================================
# No id-mapping job is ever submitted, and the fixture fetch below throws if one is:
# robots.txt disallows /idmapping/ and every /*/stream/ path, and a job can answer short
# with failedIds empty, has that answer cached against its id set, and fails outright
# under sustained volume (see the header of ModuleAccessionMapper.jl). What follows
# exercises the route offline: /uniprotkb/search for live AND deleted entries,
# /uniparc/upis for the deleted ones, and a per-accession cross-reference only where a
# UniParc entry is shared between organisms.

const SEARCH_BODY = """{"results":[
 {"primaryAccession":"A0LIVE1","entryType":"UniProtKB unreviewed (TrEMBL)","organism":{"taxonId":111}},
 {"primaryAccession":"A0LIVE2","entryType":"UniProtKB reviewed (Swiss-Prot)","organism":{"taxonId":222}},
 {"primaryAccession":"A0DEAD1","entryType":"Inactive","inactiveReason":{"inactiveReasonType":"DELETED","deletedReason":"Not part of a reference proteome"},"extraAttributes":{"uniParcId":"UPI0000000001"}},
 {"primaryAccession":"A0DEAD2","entryType":"Inactive","inactiveReason":{"inactiveReasonType":"DELETED"},"extraAttributes":{"uniParcId":"UPI0000000002"}}
]}"""

const UPIS_BODY = """{"results":[
 {"uniParcId":"UPI0000000001","uniProtKBAccessions":["A0DEAD1.1"],"organisms":[{"taxonId":333}]},
 {"uniParcId":"UPI0000000002","uniProtKBAccessions":["A0DEAD2"],"organisms":[{"taxonId":444},{"taxonId":555}]}
]}"""

const XREF_BODY = """{"results":[
 {"database":"UniProtKB/TrEMBL","id":"A0DEAD2","active":false,"organism":{"taxonId":444},
  "created":"2015-01-01","lastUpdated":"2020-01-01"}
]}"""

const IDS = ["A0LIVE1", "A0LIVE2", "A0DEAD1", "A0DEAD2"]

function make_fetch(; total = "4", log = String[], search_body = SEARCH_BODY,
                    upis_body = UPIS_BODY)
    function f(method::Symbol, url::AbstractString; body = nothing,
               headers = Dict{String, String}())
        push!(log, String(url))
        method === :GET || error("the route must never POST — got $method to $url")
        occursin("/idmapping/", url) && error("robots-disallowed path requested: $url")
        occursin("/stream/", url) && error("robots-disallowed stream path requested: $url")
        if occursin("/uniprotkb/search", url)
            return (status = 200,
                    headers = Dict("x-total-results" => total,
                                   "x-uniprot-release" => "2099_01"),
                    body = search_body)
        elseif occursin("/uniparc/upis", url)
            return (status = 200, headers = Dict("x-uniprot-release" => "2099_01"),
                    body = upis_body)
        elseif occursin("/databases", url)
            return (status = 200, headers = Dict{String, String}(), body = XREF_BODY)
        end
        error("no fixture for $url")
    end
    return f
end

@testset "three-pass route resolves live and deleted accessions" begin
    log = String[]
    df = resolve_taxids(IDS; fetch = make_fetch(log = log), verbose = false, sleeper = _ -> nothing)

    @test nrow(df) == 4
    byacc = Dict(String(r.accession) => r for r in eachrow(df))
    @test byacc["A0LIVE1"].TaxID == 111
    @test byacc["A0LIVE1"].taxid_source == "uniprotkb_active"
    @test byacc["A0LIVE2"].TaxID == 222
    @test byacc["A0DEAD1"].TaxID == 333
    @test byacc["A0DEAD1"].taxid_source == "uniparc_single_organism"
    @test byacc["A0DEAD1"].entry_status == "inactive"
    # the shared UniParc entry must NOT be answered from organisms[]; only pass 3 decides
    @test byacc["A0DEAD2"].TaxID == 444
    @test byacc["A0DEAD2"].taxid_source == "uniparc_own_xref"
    @test byacc["A0DEAD2"].uniprot_release == "2099_01"

    @test any(u -> occursin("/uniprotkb/search", u), log)
    @test any(u -> occursin("/uniparc/upis", u), log)
    @test any(u -> occursin("id=A0DEAD2", u), log)
    @test !any(u -> occursin("/idmapping/", u), log)
end

@testset "a short pass-1 block is retried on the same ids, then aborts" begin
    short = replace(SEARCH_BODY, r",\n \{\"primaryAccession\":\"A0DEAD2\".*?\}\n" => "\n")
    log = String[]
    err = nothing
    try
        resolve_taxids(IDS; fetch = make_fetch(log = log, search_body = short),
                       verbose = false, sleeper = _ -> nothing)
    catch e
        err = e
    end
    @test err !== nothing
    msg = sprint(showerror, err)
    @test occursin("short twice", msg)
    @test occursin("A0DEAD2", msg)
    @test occursin("is a loss", msg)
    # asked twice with the SAME ids, never split into smaller blocks
    searches = filter(u -> occursin("/uniprotkb/search", u), log)
    @test length(searches) == 2
    @test searches[1] == searches[2]
end

@testset "politeness: 429 is not retried through, 5xx is" begin
    calls = Ref(0)
    on429 = (m, u; kw...) -> (calls[] += 1;
             (status = 429, headers = Dict{String, String}(), body = "slow down"))
    err = nothing
    try
        AM.polite_get("https://rest.uniprot.org/x"; fetch = on429, pace_s = 0,
                      sleeper = _ -> nothing)
    catch e
        err = e
    end
    @test err !== nothing
    @test occursin("not to try harder", sprint(showerror, err))
    @test calls[] == length(AM.BACKOFF_429_S) + 1   # bounded, not endless

    n = Ref(0)
    flaky = (m, u; kw...) -> (n[] += 1;
             n[] < 3 ? (status = 503, headers = Dict{String, String}(), body = "") :
                       (status = 200, headers = Dict{String, String}(), body = """{"results":[]}"""))
    (json, _) = AM.polite_get("https://rest.uniprot.org/y"; fetch = flaky, pace_s = 0,
                              sleeper = _ -> nothing)
    @test n[] == 3
end

# URL construction only; the ledger has its own testset below.
@testset "the block cap and the query shape are enforced locally" begin
    @test_throws ErrorException AM.search_url(["A$i" for i in 1:(AM.SEARCH_BLOCK + 1)])
    @test_throws ErrorException AM.upis_url(["UPI$i" for i in 1:(AM.UPIS_BLOCK + 1)])
    u = AM.search_url(["A0LIVE1", "A0DEAD1"])
    @test occursin("accession%3AA0LIVE1+OR+accession%3AA0DEAD1", u) ||
          occursin("accession%3AA0LIVE1%20OR%20accession%3AA0DEAD1", u)
    @test occursin("size=500", u)      # or the response silently pages at 25
end

@testset "a paged cross-reference listing is followed by its Link header" begin
    # Paginated cross-reference listings are followed by the Link header, not guessed at.
    @test AM._next_link("""<https://rest.uniprot.org/uniparc/UPI0000000001/databases?cursor=x>; rel="next\"""") ==
          "https://rest.uniprot.org/uniparc/UPI0000000001/databases?cursor=x"
    @test AM._next_link("") == ""
end

@testset "the accession cache makes a re-run cost nothing" begin
    dir = mktempdir()
    cache = joinpath(dir, "taxid_cache.csv")
    log1 = String[]
    resolve_taxids(IDS; fetch = make_fetch(log = log1), taxid_cache = cache,
                   verbose = false, sleeper = _ -> nothing)
    @test isfile(cache)

    log2 = String[]
    df2 = resolve_taxids(IDS; fetch = make_fetch(log = log2), taxid_cache = cache,
                         verbose = false, sleeper = _ -> nothing)
    @test isempty(log2)                              # zero HTTP on the second run
    byacc = Dict(String(r.accession) => r for r in eachrow(df2))
    @test byacc["A0LIVE1"].TaxID == 111
    @test byacc["A0DEAD2"].TaxID == 444
end

# An inactive entry with no UniParc pointer is the shape a renamed or moved key produces,
# and it must trip the rescue floor rather than pass as 0 needed / 0 resolved.
@testset "an inactive entry with no UniParc pointer trips the floor, not 0/0" begin
    pointerless_body = """{"results":[
     {"primaryAccession":"A0DEL001","entryType":"Inactive","inactiveReason":{"inactiveReasonType":"DELETED","deletedReason":"Redundant sequence"}},
     {"primaryAccession":"A0DEL002","entryType":"Inactive","inactiveReason":{"inactiveReasonType":"DELETED","deletedReason":"Redundant sequence"}}
    ]}"""
    function pointerless(method::Symbol, url::AbstractString; body = nothing,
                         headers = Dict{String, String}())
        occursin("/uniprotkb/search", url) &&
            return (status = 200,
                    headers = Dict("x-total-results" => "2", "x-uniprot-release" => "2099_01"),
                    body = pointerless_body)
        error("the rescue asked for $url although no UniParc id was ever returned")
    end
    @test_throws ErrorException resolve_taxids(["A0DEL001", "A0DEL002"];
                                               fetch = pointerless, verbose = false,
                                               sleeper = _ -> nothing)
end

# =======================================================================================
# REGRESSIONS - guards whose failure would be silent.
# Each drives the path its guard protects. A route can pass its own suite while returning a
# fraction of the accessions it was given, so "the tests are green" is not evidence here;
# what the tests REACH is.
# =======================================================================================

@testset "REGRESSION - a cache row that is not an established resolution is refused" begin
    dir = mktempdir()
    cache = joinpath(dir, "poisoned.csv")
    open(cache, "w") do io
        println(io, join(AM.TAXID_CACHE_COLS, ","))
        println(io, "A0GOOD1,111,uniprotkb_active,,true,2099_01,2099-01-01T00:00:00")
        println(io, "A0POISON,,unresolved,,false,2099_01,2099-01-01T00:00:00")
        println(io, "A0BADSRC,222,made_up_source,,true,2099_01,2099-01-01T00:00:00")
    end
    loaded = AM.load_taxid_cache(cache)
    @test haskey(loaded, "A0GOOD1")
    # An accession with no TaxID is an ABSENCE of a result, not a result. Loading it would
    # make `todo` skip it forever and no guard would ever see it again.
    @test !haskey(loaded, "A0POISON")
    @test !haskey(loaded, "A0BADSRC")
    @test !("unresolved" in AM.TAXID_SOURCES)
end

@testset "REGRESSION - a run that fails a guard leaves no cache behind" begin
    # The cache must be written only AFTER every guard has passed. Written before them, an
    # aborting run still leaves its rows; the next run then finds `todo` empty, skips every
    # pass and every guard, and publishes the short answer the abort existed to prevent.
    dir = mktempdir()
    cache = joinpath(dir, "taxid_cache.csv")
    pointerless_body = """{"results":[
     {"primaryAccession":"A0DEL001","entryType":"Inactive","inactiveReason":{"inactiveReasonType":"DELETED","deletedReason":"Redundant sequence"}},
     {"primaryAccession":"A0DEL002","entryType":"Inactive","inactiveReason":{"inactiveReasonType":"DELETED","deletedReason":"Redundant sequence"}}
    ]}"""
    function pointerless(method::Symbol, url::AbstractString; body = nothing,
                         headers = Dict{String, String}())
        occursin("/uniprotkb/search", url) &&
            return (status = 200,
                    headers = Dict("x-total-results" => "2", "x-uniprot-release" => "2099_01"),
                    body = pointerless_body)
        error("no fixture for $url")
    end
    @test_throws ErrorException resolve_taxids(["A0DEL001", "A0DEL002"];
                                               fetch = pointerless, taxid_cache = cache,
                                               verbose = false, sleeper = _ -> nothing)
    @test !isfile(cache)
end

@testset "REGRESSION - pass 1 checks membership, not cardinality" begin
    # A response carrying an accession nobody asked for holds the COUNT at the block size
    # while one of OUR ids is absent. A size check calls that block complete.
    swapped = """{"results":[
     {"primaryAccession":"A0LIVE1","organism":{"taxonId":111}},
     {"primaryAccession":"A0NEVER1","organism":{"taxonId":999}}
    ]}"""
    calls = String[]
    function substituting(method::Symbol, url::AbstractString; body = nothing,
                          headers = Dict{String, String}())
        push!(calls, String(url))
        occursin("/uniprotkb/search", url) &&
            return (status = 200,
                    headers = Dict("x-total-results" => "2", "x-uniprot-release" => "2099_01"),
                    body = swapped)
        error("no fixture for $url")
    end
    @test_throws ErrorException resolve_taxids(["A0LIVE1", "A0LIVE2"];
                                               fetch = substituting, verbose = false,
                                               sleeper = _ -> nothing)
    @test length(calls) == 2      # asked again for the SAME ids before giving up
end

@testset "REGRESSION - every block is fetched, not only the first" begin
    # The failure this module guards against is a request that answers for part of its
    # input. With a single-block fixture a suite cannot tell "fetched every block" from
    # "fetched the first and stopped", so this one deliberately spans two.
    n = AM.SEARCH_BLOCK + 50
    ids = ["A0BULK" * lpad(i, 4, '0') for i in 1:n]
    seen_blocks = Ref(0)
    function bulk(method::Symbol, url::AbstractString; body = nothing,
                  headers = Dict{String, String}())
        if occursin("/uniprotkb/search", url)
            asked = [m.match for m in eachmatch(r"A0BULK[0-9]{4}", url)]
            seen_blocks[] += 1
            entries = join(["""{"primaryAccession":"$a","organism":{"taxonId":$(parse(Int, a[7:end]))}}"""
                            for a in asked], ",")
            return (status = 200,
                    headers = Dict("x-total-results" => string(length(asked)),
                                   "x-uniprot-release" => "2099_01"),
                    body = """{"results":[$entries]}""")
        end
        error("no fixture for $url")
    end
    df = resolve_taxids(ids; fetch = bulk, verbose = false, sleeper = _ -> nothing)
    @test seen_blocks[] == 2
    @test nrow(df) == n
    @test all(!ismissing, df.TaxID)      # a dropped second block leaves 50 unresolved
    @test df.TaxID[end] == n
end

@testset "REGRESSION - the pass-3 cross-reference cache is actually consulted" begin
    # The lookup key must be the key the WRITER wrote. A Tuple looked up in a String-keyed
    # Dict never hits, so the documented resumability is silently dead and every shared
    # UniParc entry is re-fetched, one serial request at a time, against a service this
    # module exists to stop overloading. A round-trip test of the cache file cannot see
    # that: the writer and the reader can both be fine while the CALL SITE is not. So this
    # drives the real path twice and asserts the second run asks for no cross-reference at
    # all.
    dir = mktempdir()
    xref = joinpath(dir, "xref_cache.csv")

    log1 = String[]
    df1 = resolve_taxids(IDS; fetch = make_fetch(log = log1), xref_cache = xref,
                         verbose = false, sleeper = _ -> nothing)
    @test any(u -> occursin("/databases", u), log1)     # run 1 pays for it
    @test isfile(xref)

    # No taxid_cache here on purpose: that layer would short-circuit the whole run and hide
    # whether pass 3 consulted its own cache.
    log2 = String[]
    df2 = resolve_taxids(IDS; fetch = make_fetch(log = log2), xref_cache = xref,
                         verbose = false, sleeper = _ -> nothing)
    @test !any(u -> occursin("/databases", u), log2)    # run 2 must not
    @test Dict(String(r.accession) => r.TaxID for r in eachrow(df2))["A0DEAD2"] ==
          Dict(String(r.accession) => r.TaxID for r in eachrow(df1))["A0DEAD2"]
end

# ---------------------------------------------------------------------------------------
# CACHE-WARM RE-RUNS. Caching only established resolutions and measuring the rescue floor
# are each correct on their own, and they interact: a healthy re-run must not be aborted by
# the failures it legitimately left behind, a failing run must leave no cache even when
# some of its rows resolved, and a warm re-run must keep the release its rows were
# resolved under.
# ---------------------------------------------------------------------------------------

# One live accession, one deleted-and-rescuable, one deleted with no UniParc pointer. The
# third can never resolve, which is the point: it is what a healthy run legitimately leaves
# behind, and what the NEXT run must not be punished for.
const MIXED_SEARCH = """{"results":[
 {"primaryAccession":"A0ACT1","organism":{"taxonId":11}},
 {"primaryAccession":"A0RESC1","entryType":"Inactive","inactiveReason":{"inactiveReasonType":"DELETED"},"extraAttributes":{"uniParcId":"UPI_R1"}},
 {"primaryAccession":"A0FAIL1","entryType":"Inactive","inactiveReason":{"inactiveReasonType":"DELETED"}}
]}"""
const MIXED_UPIS = """{"results":[
 {"uniParcId":"UPI_R1","uniProtKBAccessions":["A0RESC1"],"organisms":[{"taxonId":22}]}
]}"""

function mixed_fetch(; log = String[])
    function f(method::Symbol, url::AbstractString; body = nothing,
               headers = Dict{String, String}())
        push!(log, String(url))
        if occursin("/uniprotkb/search", url)
            # Answer only for the accessions actually asked for, so a second run that asks
            # for a single leftover gets a one-entry response and not the whole fixture.
            asked = [m.match for m in eachmatch(r"A0[A-Z]+1", url)]
            keep = filter(e -> any(a -> occursin(a, e), asked),
                          split(MIXED_SEARCH, "\n")[2:(end - 1)])
            n = length(keep)
            return (status = 200,
                    headers = Dict("x-total-results" => string(n),
                                   "x-uniprot-release" => "2099_01"),
                    body = "{" * "\"results\":[" * join(rstrip.(keep, [',']), ",") * "]}")
        elseif occursin("/uniparc/upis", url)
            return (status = 200, headers = Dict("x-uniprot-release" => "2099_01"),
                    body = MIXED_UPIS)
        end
        error("no fixture for $url")
    end
    return f
end

@testset "REGRESSION - a healthy re-run is not aborted by its own leftovers" begin
    # Only resolutions are cached, so after run 1 the cache holds the successes and `todo`
    # on run 2 holds exactly the failures. A rescue floor measured over `todo` alone is then
    # measured over a population that is 100% known failures: 0/k by construction. The floor
    # has to count the cached rescues on BOTH sides so it keeps describing the census.
    dir = mktempdir()
    cache = joinpath(dir, "taxid_cache.csv")
    ids = ["A0ACT1", "A0RESC1", "A0FAIL1"]

    df1 = resolve_taxids(ids; fetch = mixed_fetch(), taxid_cache = cache,
                         min_rescue_frac = 0.5, verbose = false, sleeper = _ -> nothing)
    by1 = Dict(String(r.accession) => r for r in eachrow(df1))
    @test by1["A0ACT1"].TaxID == 11
    @test by1["A0RESC1"].TaxID == 22
    @test ismissing(by1["A0FAIL1"].TaxID)        # legitimately unresolved
    @test isfile(cache)

    # Identical input, identical service. This must NOT throw.
    df2 = resolve_taxids(ids; fetch = mixed_fetch(), taxid_cache = cache,
                         min_rescue_frac = 0.5, verbose = false, sleeper = _ -> nothing)
    by2 = Dict(String(r.accession) => r for r in eachrow(df2))
    @test by2["A0ACT1"].TaxID == 11
    @test by2["A0RESC1"].TaxID == 22
    @test ismissing(by2["A0FAIL1"].TaxID)
end

@testset "REGRESSION - a failing run leaves no cache even when some rows resolved" begin
    # A fixture in which every accession is unresolved makes this test vacuous: `fresh_rows`
    # is empty and no cache is written whatever the ordering. This one resolves one
    # accession and THEN fails the floor, so the file appears if and only if the write is
    # ahead of the guards.
    dir = mktempdir()
    cache = joinpath(dir, "taxid_cache.csv")
    body = """{"results":[
     {"primaryAccession":"A0ACT1","organism":{"taxonId":11}},
     {"primaryAccession":"A0FAIL1","entryType":"Inactive","inactiveReason":{"inactiveReasonType":"DELETED"}},
     {"primaryAccession":"A0FAIL2","entryType":"Inactive","inactiveReason":{"inactiveReasonType":"DELETED"}}
    ]}"""
    g = (method, url; kwargs...) -> begin
        occursin("/uniprotkb/search", String(url)) &&
            return (status = 200,
                    headers = Dict("x-total-results" => "3", "x-uniprot-release" => "2099_01"),
                    body = body)
        error("no fixture for $url")
    end
    @test_throws ErrorException resolve_taxids(["A0ACT1", "A0FAIL1", "A0FAIL2"];
                                               fetch = g, taxid_cache = cache,
                                               verbose = false, sleeper = _ -> nothing)
    @test !isfile(cache)      # A0ACT1 resolved, so fresh_rows was NOT empty
end

@testset "REGRESSION - a cache-warm re-run keeps the release it was resolved under" begin
    # The release is a property of the ROW, not of the run. Stamped per run, a fully cached
    # re-run (where pass 1 never executes and its release is "") blanks the column for every
    # row, and the census manifest then records an empty taxonomy provenance.
    dir = mktempdir()
    cache = joinpath(dir, "taxid_cache.csv")
    df1 = resolve_taxids(IDS; fetch = make_fetch(), taxid_cache = cache,
                         verbose = false, sleeper = _ -> nothing)
    @test all(r -> r.uniprot_release == "2099_01", eachrow(df1))

    df2 = resolve_taxids(IDS; fetch = make_fetch(), taxid_cache = cache,
                         verbose = false, sleeper = _ -> nothing)
    @test all(r -> r.uniprot_release == "2099_01", eachrow(df2))
end

# ---------------------------------------------------------------------------------------
# WHAT A RUN REMEMBERS. A short cross-reference answer cached as an absence, and a cache
# that drops the description columns, both leave every count unchanged, so a test that
# checks counts cannot see either. These two check what is written.
# ---------------------------------------------------------------------------------------

# A cross-reference page that says it holds more records than it returned. This is NOT an
# error status, so nothing upstream retries it: the service answers 200 with a short body,
# the paging loop runs out, and the accession looks as though it simply has no own
# cross-reference.
const SHORT_XREF_BODY = """{"results":[
 {"database":"UniProtKB/TrEMBL","id":"SOMEONE_ELSE","active":false,"organism":{"taxonId":999}}
]}"""

@testset "REGRESSION - a SHORT cross-reference answer is not cached as an absence" begin
    # Written to the xref cache unconditionally, `nothing` turns a transient short page into
    # a permanent "this accession has no cross-reference", and no later run ever asks again
    # - the census loses a taxon silently and for good.
    dir = mktempdir()
    xref = joinpath(dir, "xref_cache.csv")

    short_fetch(log) = function (method::Symbol, url::AbstractString; body = nothing,
                                 headers = Dict{String, String}())
        push!(log, String(url))
        if occursin("/uniprotkb/search", url)
            return (status = 200,
                    headers = Dict("x-total-results" => "4",
                                   "x-uniprot-release" => "2099_01"),
                    body = SEARCH_BODY)
        elseif occursin("/uniparc/upis", url)
            return (status = 200, headers = Dict("x-uniprot-release" => "2099_01"),
                    body = UPIS_BODY)
        elseif occursin("/databases", url)
            # Declares five, returns one, and that one is a different accession.
            return (status = 200, headers = Dict("x-total-results" => "5"),
                    body = SHORT_XREF_BODY)
        end
        error("no fixture for $url")
    end

    # min_rescue_frac = 0 on purpose: a short answer legitimately leaves the accession
    # unresolved, so the rescue floor would abort the run before the cache is even written.
    # That floor is the subject of T5; what is under test here is what the run REMEMBERS.
    log1 = String[]
    resolve_taxids(IDS; fetch = short_fetch(log1), xref_cache = xref, min_rescue_frac = 0.0,
                   verbose = false, sleeper = _ -> nothing)
    @test any(u -> occursin("/databases", u), log1)

    # THE ASSERTION THAT MATTERS, and it runs unconditionally. Written as `if isfile` it
    # would SKIP when the file is absent, which is the very outcome it is meant to confirm -
    # a test that disappears exactly when it passes proves nothing on the day it should fail.
    cached_text = isfile(xref) ? read(xref, String) : ""
    @test !occursin("A0DEAD2", cached_text)

    # ... and the next run must therefore ask again rather than inherit the gap.
    log2 = String[]
    resolve_taxids(IDS; fetch = short_fetch(log2), xref_cache = xref, min_rescue_frac = 0.0,
                   verbose = false, sleeper = _ -> nothing)
    @test any(u -> occursin("/databases", u), log2)
end

@testset "REGRESSION - the accession cache carries what the row SAYS, not only its taxid" begin
    # A cache that holds only the resolution drops `inactive_reason`, `xref_created` and
    # `xref_last_updated` on the way in. Pass 1 is not re-run for a cached accession, so a
    # warm re-run would rebuild the same taxa with three empty columns and report the same
    # census.
    dir = mktempdir()
    cache = joinpath(dir, "taxid_cache.csv")
    AM.append_taxid_cache(cache, [(accession = "A0DEAD2", taxid = 444,
                                   source = "uniparc_own_xref", uniparc_id = "UPI0000000202",
                                   active = false, release = "2099_01",
                                   fetched_at = "2099-01-01T00:00:00",
                                   inactive_reason = "deleted: merged, with a comma",
                                   xref_created = "2015-01-01",
                                   xref_last_updated = "2020-01-01")])
    back = AM.load_taxid_cache(cache)
    @test haskey(back, "A0DEAD2")
    @test back["A0DEAD2"].created == "2015-01-01"
    @test back["A0DEAD2"].updated == "2020-01-01"
    # The comma must not have split the row into a new column.
    @test back["A0DEAD2"].reason == "deleted: merged; with a comma"

    # A cache without those columns is still a cache of valid resolutions: refusing it over
    # a schema difference would throw away paid-for work.
    narrow = joinpath(dir, "narrow_schema.csv")
    open(narrow, "w") do io
        println(io, "accession,taxid,source,uniparc_id,active,uniprot_release,fetched_at")
        println(io, "A0DEAD2,444,uniparc_own_xref,UPI0000000202,false,2099_01,2099-01-01")
    end
    narrow_back = AM.load_taxid_cache(narrow)
    @test haskey(narrow_back, "A0DEAD2")
    @test narrow_back["A0DEAD2"].taxid == 444
    @test narrow_back["A0DEAD2"].reason === missing
    @test narrow_back["A0DEAD2"].created === missing
end

@testset "the ledger is checked against the TABLE, not against itself" begin
    # A sum of the counters cannot fail: the assembly loop runs once per accession and every
    # branch increments exactly one counter, so the total restates the loop's bound, and on
    # its own it would report a balanced ledger over a table that had lost half its rows.
    # What is asserted here is the part that can disagree - the counters against the table.
    dupes = ["A0LIVE1", "A0DEAD2", "A0LIVE1", "A0DEAD2", "A0LIVE2", "A0DEAD1"]
    df = resolve_taxids(dupes; fetch = make_fetch(), verbose = false, sleeper = _ -> nothing)

    # One row per DISTINCT accession, however many times it was asked for. A duplicate row
    # here would multiply every downstream join on the census silently.
    @test nrow(df) == length(unique(dupes))
    @test length(unique(df.accession)) == nrow(df)
    @test Set(String.(df.accession)) == Set(dupes)

    # An unresolved accession is exactly a row with no TaxID, counted two different ways.
    @test count(ismissing, df.TaxID) == count(r -> r.taxid_source == "unresolved", eachrow(df))

    # Every resolved row names the route that resolved it; "unresolved" is not a route.
    resolved = [r for r in eachrow(df) if !ismissing(r.TaxID)]
    @test !isempty(resolved)
    @test all(r -> r.taxid_source in AM.TAXID_SOURCES, resolved)
end

@testset "choose_taxid is not called by the pipeline, and the census columns it would feed stay empty" begin
    # The metagenome safeguard in `choose_taxid` does not act on the census, and the repair
    # is NOT to connect it: pass 2 resolves through single-organism UniParc entries, so the
    # two sources this function returns are not in the live vocabulary and a row carrying
    # one is refused on the way back in. Pinned here so that a reader who finds an unused
    # exported function does not helpfully wire it up.
    @test !("uniparc_common_taxon" in AM.TAXID_SOURCES)
    @test !("uniparc_common_noncellular" in AM.TAXID_SOURCES)

    # Those really are the sources it returns, so the check above is about this function and
    # not about two strings nothing produces.
    hits = AM.parse_uniparc_common(JSON3.read(
        """{"results":[{"uniParcId":"UPI0000000909","commonTaxons":[
             {"topLevel":"cellular organisms","commonTaxon":"Bacillus","commonTaxonId":1409}]}]}"""))
    @test choose_taxid(nothing, hits["UPI0000000909"]).source == "uniparc_common_taxon"

    # And the two columns it would feed are empty on every row of a real run, by
    # construction.
    df = resolve_taxids(IDS; fetch = make_fetch(), verbose = false, sleeper = _ -> nothing)
    @test all(r -> r.common_taxon_ids == "", eachrow(df))
    @test all(r -> ismissing(r.alt_cellular_taxid), eachrow(df))
end

# Every id that was asked for, PLUS one that was not. A secondary accession resolving to a
# different primary looks exactly like this.
const SURPLUS_SEARCH_BODY = """{"results":[
 {"primaryAccession":"A0LIVE1","entryType":"UniProtKB unreviewed (TrEMBL)","organism":{"taxonId":111}},
 {"primaryAccession":"A0LIVE2","entryType":"UniProtKB reviewed (Swiss-Prot)","organism":{"taxonId":222}},
 {"primaryAccession":"A0DEAD1","entryType":"Inactive","inactiveReason":{"inactiveReasonType":"DELETED","deletedReason":"Not part of a reference proteome"},"extraAttributes":{"uniParcId":"UPI0000000001"}},
 {"primaryAccession":"A0DEAD2","entryType":"Inactive","inactiveReason":{"inactiveReasonType":"DELETED"},"extraAttributes":{"uniParcId":"UPI0000000002"}},
 {"primaryAccession":"A0NEVER1","entryType":"UniProtKB unreviewed (TrEMBL)","organism":{"taxonId":777}}
]}"""

@testset "REGRESSION - a block that returns MORE than it was asked says so" begin
    # Nothing is missing here. Seen only as a disagreement with x-total-results, a surplus
    # reads as "0 of 4 missing" followed by an abort as "came back short twice" with an
    # empty list of missing ids - a message describing the opposite of what happened. The
    # unrequested ids are named from the response itself.
    err = nothing
    try
        resolve_taxids(IDS;
                       fetch = make_fetch(total = "5", search_body = SURPLUS_SEARCH_BODY),
                       verbose = false, sleeper = _ -> nothing)
    catch e
        err = e
    end
    @test err isa ErrorException
    msg = sprint(showerror, err)
    @test occursin("not asked for", msg)
    @test occursin("A0NEVER1", msg)          # named, not merely counted
    @test !occursin("came back short", msg)  # and not described as a loss
end

# The same pass-2 payload with UPI0000000002 absent: the block asks for two UniParc ids and
# the answer carries one. Synthetic values, like every fixture above.
const UPIS_BODY_SHORT = """{"results":[
 {"uniParcId":"UPI0000000001","uniProtKBAccessions":["A0DEAD1.1"],"organisms":[{"taxonId":333}]}
]}"""

@testset "REGRESSION - a short pass-2 block names the UniParc ids that did not come back" begin
    # Pass 1 has its own short-answer cases above; this is the same failure one pass later,
    # on /uniparc/upis. A UniParc id that was submitted and did not come back is a loss, and
    # a loss that is only counted reads as a smaller census: every submitted id is checked by
    # MEMBERSHIP and the absent ones are named. Accession counts are the wrong unit here,
    # because several accessions can legitimately share one entry.
    log = String[]
    err = nothing
    try
        resolve_taxids(IDS; fetch = make_fetch(log = log, upis_body = UPIS_BODY_SHORT),
                       verbose = false, sleeper = _ -> nothing)
    catch e
        err = e
    end
    @test err isa ErrorException
    msg = sprint(showerror, err)
    @test occursin("Pass 2 lost", msg)
    @test occursin("1 of 2 UniParc ids", msg)
    @test occursin("UPI0000000002", msg)     # named, not merely counted
    # One block, asked once: a short pass-2 answer is neither retried nor split.
    upis_calls = filter(u -> occursin("/uniparc/upis", u), log)
    @test length(upis_calls) == 1
    @test occursin("UPI0000000002", upis_calls[1])
    # And the run stops here: the per-accession lookups of pass 3 are never reached.
    @test !any(u -> occursin("/databases", u), log)
end

@testset "REGRESSION - the 429 back-off budget is per request, not per run" begin
    # One request that keeps being told to stop: it waits the whole ladder, then the fourth
    # 429 throws and the run ends.
    waits = Float64[]
    calls = Ref(0)
    always429 = (m, u; kw...) -> (calls[] += 1;
                 (status = 429, headers = Dict{String, String}(), body = "slow down"))
    err = nothing
    try
        AM.polite_get("https://rest.uniprot.org/a"; fetch = always429, pace_s = 0,
                      sleeper = s -> push!(waits, Float64(s)))
    catch e
        err = e
    end
    @test err isa ErrorException
    @test calls[] == length(AM.BACKOFF_429_S) + 1
    @test waits == collect(Float64.(AM.BACKOFF_429_S))
    @test occursin("for one request", sprint(showerror, err))

    # Six 429s over two requests, three each. Neither request spends its own budget, so both
    # answer and the run carries on: the counter lives in the call and not in the run, which
    # is what the module's 429 paragraph says.
    n429 = Ref(0)
    seen = Dict{String, Int}()
    spread = function (m, u; kw...)
        k = String(u)
        c = get(seen, k, 0) + 1
        seen[k] = c
        if c <= length(AM.BACKOFF_429_S)
            n429[] += 1
            return (status = 429, headers = Dict{String, String}(), body = "slow down")
        end
        return (status = 200, headers = Dict{String, String}(),
                body = """{"results":[]}""")
    end
    for u in ("https://rest.uniprot.org/b", "https://rest.uniprot.org/c")
        (json, _) = AM.polite_get(u; fetch = spread, pace_s = 0, sleeper = _ -> nothing)
        @test json !== nothing
    end
    @test n429[] == 2 * length(AM.BACKOFF_429_S)

    # Retry-After is honoured but floored at the sustained pace: a header of 0 must not buy
    # an immediate re-ask, which is the one answer a 429 rules out.
    floored = Float64[]
    seen429 = Ref(0)
    zero_retry = (m, u; kw...) -> (seen429[] += 1;
                  seen429[] == 1 ?
                  (status = 429, headers = Dict("retry-after" => "0"), body = "slow down") :
                  (status = 200, headers = Dict{String, String}(),
                   body = """{"results":[]}"""))
    AM.polite_get("https://rest.uniprot.org/d"; fetch = zero_retry, pace_s = 0,
                  sleeper = s -> push!(floored, Float64(s)))
    @test floored == [Float64(AM.PACE_S)]
end
