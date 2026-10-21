module AccessionMapper

export parse_fasta_headers, fetch_taxids_async, resolve_taxids
export search_url, upis_url, parse_uniprotkb_search, parse_uniparc_upis
export parse_uniparc_common, own_xref, choose_taxid
export check_rescue

using HTTP
using JSON3
using DataFrames
using CodecZlib
using CSV

const REST_ROOT = "https://rest.uniprot.org"

#=
    WHY THIS MODULE DOES NOT USE UNIPROT'S ID-MAPPING SERVICE.

    The obvious route to a taxon per accession is an id-mapping job: POST a set of
    accessions, poll the job, stream the result. That route fails in three ways, and the
    third is the one that matters:

      1. A job can SUCCEED and answer short: fewer results than ids submitted, with
         `failedIds` empty and the service's own `x-total-results` agreeing with the short
         number. It is not a fixed size cap and not a race — the short count is stable over
         repeated reads — and the dropped ids are ordinary live entries that map when asked
         for alone.
      2. Asking again cannot fix that, because a job is keyed on its id set: the same ids
         return the same jobId and the same cached answer. Splitting the set is the only
         recovery; a chunk that answers short can come back whole as two halves.
      3. Under a sustained run of back-to-back jobs the service starts FAILING them with an
         internal 503, while jobs of one to three ids still answer. The failures are caused
         by the client's own volume, and splitting — which turns one job into two — makes it
         worse: every level of halving still fails.

    And the decisive point: https://rest.uniprot.org/robots.txt disallows `/idmapping/`
    outright, and every `/*/stream/` path. The route below needs neither, creates no
    server-side state, and cannot have a wrong answer cached against it. The failure modes
    above were measured directly against the live service.
=#

# Pass 1: UniProtKB search, which returns ACTIVE and DELETED entries in one response —
# actives with organism.taxonId, deleted ones with entryType "Inactive" and the
# extraAttributes.uniParcId that pass 2 needs. `/uniprotkb/accessions` is NOT usable here:
# it omits deleted entries silently, which is the ambiguity this route exists to remove.
#
# 100 is the service's own cap: a query with more terms answers HTTP 400 "Too many OR
# conditions in query. Maximum allowed is 100." The default page size is 25, so size= must
# be set or the response pages behind your back.
const SEARCH_BLOCK = 100
const SEARCH_PAGE = 500
const PASS1_FIELDS = "accession,organism_id"

# Pass 2: UniParc, 500 UPIs per request (measured at 500; not probed above it, do not).
# Always send fields= — an entry is ~246 B with it and ~3.6 kB without.
const UPIS_BLOCK = 500
const PASS2_FIELDS = "upi,accession,organism_id"

# Politeness. All of this is OUR restraint, not an advertised limit: UniProt's responses,
# as observed, carry no X-RateLimit-*, no RateLimit-Policy, no Retry-After and no Server
# header, and UniProt publishes no requests-per-second figure. The only documented cause of
# a 429 anywhere in their help corpus is "too many parallel requests", so the rule that
# actually matters is the first one.
const PACE_S = 1.0              # sustained, between request starts; never go above 3 req/s
const PASS3_PACE_S = 0.3        # pass 3's calls are 30-100 ms; a shorter pause is defensible
const USER_AGENT = "LnFinder/0.1"
const HTTP_TRIES_5XX = 5        # UniProt's own published client behaviour for 5xx
const BACKOFF_429_S = (60, 120, 300)   # per request; a 429 is not meant to be retried through

_string_or_missing(x) = (x === nothing || ismissing(x)) ? missing : String(x)

# A commonTaxon is only usable as a taxonomic answer when it sits in one of these
# top-level branches. Anything else ("unclassified entries" and friends) is a
# metagenome-style bucket: it is kept, but flagged, never silently promoted.
const CELLULAR_TOP_LEVELS = ("cellular organisms", "Viruses")

# One taxon candidate, whatever shape it was parsed from. `toplevel` is "" when the
# payload did not say (the legacy `organisms` shape carries no top-level branch).
const TaxonHit = NamedTuple{(:taxid, :name, :toplevel), Tuple{Int, String, String}}

# ------------------------------------------------------------------
# 1. FASTA Parsing
# ------------------------------------------------------------------
function parse_fasta_headers(fasta_path::String)
    headers = String[]
    accessions = String[]

    println("📖 Reading FASTA file: $fasta_path")
    if !isfile(fasta_path)
        error("File not found: $fasta_path")
    end

    open(fasta_path, "r") do io
        for line in eachline(io)
            if startswith(line, ">")
                raw_header = strip(line[2:end])
                acc = ""

                # Handle >sp|ACC|ID, >UniRef_ACC, and >ACC_ID formats
                if occursin('|', raw_header)
                    parts = split(raw_header, '|')
                    acc = length(parts) >= 2 ? String(parts[2]) : String(parts[1])
                elseif startswith(raw_header, "UniRef")
                    parts = split(raw_header, '_')
                    acc = length(parts) >= 2 ? String(parts[2]) : raw_header
                else
                    m = match(r"^([A-Z0-9]+)", raw_header)
                    acc = isnothing(m) ? raw_header : String(m.captures[1])
                end

                acc = split(acc, ' ')[1] # Remove trailing descriptions
                push!(headers, raw_header)
                push!(accessions, acc)
            end
        end
    end

    u_len = length(unique(accessions))
    println("✅ Found $(length(accessions)) sequences ($u_len unique).")
    return DataFrame(original_header = headers, accession = accessions)
end

# ------------------------------------------------------------------
# 2. The single HTTP seam
# ------------------------------------------------------------------
#=
    Everything that talks to the network goes through one injected callable, the way
    EnvironmentTools injects `path_of` into aggregate_habitats. Tests pass a fetch that
    answers from inline fixtures, so the parsers below are exercised on the live JSON
    shapes without a socket.

    Contract:  fetch(method::Symbol, url; body=nothing, headers=Dict()) ->
               (status::Int, headers::Dict{String,String}, body::String)

      * method is :GET or :POST; a :POST `body` is a Dict of form fields.
      * header names in the returned Dict are lower-cased, so callers can look up
        "x-uniprot-release" or "link" without guessing the server's capitalisation.
      * redirects are NOT followed and HTTP statuses never throw: what a status means is
        decided by the caller (`polite_get`), which refuses a redirect rather than follow
        it somewhere nobody asked for.
=#
"""
    default_fetch(method, url; body=nothing, headers=Dict()) -> (status, headers, body)

Real-network implementation of the fetch seam (see the block comment above). Gzip bodies
are decoded here, so every caller sees plain text.
"""
function default_fetch(method::Symbol, url::AbstractString;
                       body = nothing, headers = Dict{String, String}())
    hdrs = [String(k) => String(v) for (k, v) in pairs(headers)]
    resp = if method === :GET
        HTTP.get(url, hdrs; redirect = false, status_exception = false)
    elseif method === :POST
        form = HTTP.Form(body === nothing ? Dict{String, String}() : body)
        HTTP.post(url, hdrs, form; redirect = false, status_exception = false)
    else
        error("default_fetch: unsupported method $method")
    end
    return (status = resp.status,
            headers = Dict(lowercase(String(k)) => String(v) for (k, v) in resp.headers),
            body = _decode_body(resp.body))
end

# The service answers `compressed=true` with gzip, but HTTP.jl may already have unwrapped
# it. Sniff the magic bytes instead of trusting Content-Encoding either way.
function _decode_body(raw)
    bytes = Vector{UInt8}(raw)
    if length(bytes) >= 2 && bytes[1] == 0x1f && bytes[2] == 0x8b
        return String(transcode(GzipDecompressor, bytes))
    end
    return String(bytes)
end

"""
    header_value(headers, name; default="") -> String

Case-insensitive lookup in a fetch response's header Dict.
"""
header_value(headers, name::AbstractString; default::AbstractString = "") =
    String(get(headers, lowercase(name), default))

# ------------------------------------------------------------------
# 3. The polite request layer
# ------------------------------------------------------------------
#=
    One serial GET at a time, paced, with an error policy that differs per status because
    the statuses mean different things:

      429  we are being told to stop. Never retry immediately and NEVER split the work —
           splitting doubles the load on a service that has just said it is busy. Honour
           Retry-After if present, else back off 60/120/300: three back-offs PER REQUEST,
           and the fourth 429 on one request throws, which stops the run. The counter is
           scoped to the request on purpose. The accession cache is written only once every
           guard has passed, so a run-wide tally would let isolated 429s hours apart throw
           away all of pass 1 and pass 2, while each request on its own already waits the
           whole ladder before it gives up. UniProt's own reference client deliberately
           omits 429 from its retry list: a 429 is not meant to be retried through.
      5xx  transient. Same request, up to HTTP_TRIES_5XX, backoff factor 0.25 — this is
           UniProt's own published client behaviour. A run that hits several of these is
           seeing the saturation signature described at the top of this file; the right
           response is to stop and come back later, not to push harder.
      400  not transient. The block size and accession syntax are validated locally
           precisely so this cannot be answered by retrying.
=#
"""
    polite_get(url; fetch=default_fetch, pace_s=PACE_S, tries=HTTP_TRIES_5XX)
        -> (json, headers)

A single GET, paced and retried per the policy above. Throws on 400, on a fourth 429
within this one request, and on a 5xx that outlives its retries. The 429 back-off budget is
per call: it does not carry from one request to the next.
"""
function polite_get(url::AbstractString; fetch = default_fetch, pace_s::Real = PACE_S,
                    tries::Int = HTTP_TRIES_5XX, sleeper = sleep)
    hdrs = Dict("User-Agent" => USER_AGENT, "Accept" => "application/json")
    n429 = 0        # this request's back-off budget; see the 429 paragraph above
    for attempt in 1:tries
        (status, headers, body) = fetch(:GET, url; headers = hdrs)

        if status == 200
            pace_s > 0 && sleeper(pace_s)
            return (JSON3.read(body), headers)

        elseif status == 429
            n429 += 1
            n429 > length(BACKOFF_429_S) &&
                error("UniProt answered 429 (too many requests) $(n429) times for one " *
                      "request: $url. Stopping: a 429 is a request to go away, not to try " *
                      "harder. Resume later; cross-references fetched so far are " *
                      "checkpointed.")
            retry_after = header_value(headers, "Retry-After")
            asked_s = isempty(retry_after) ? BACKOFF_429_S[n429] :
                      something(tryparse(Float64, retry_after), BACKOFF_429_S[n429])
            # Retry-After is honoured, but never below the sustained pace. A header of 0,
            # or of a fraction of a second, would turn "stop" into an immediate re-ask —
            # the one response a 429 rules out — and it arrives from the service that is
            # already saying it has too much work.
            wait_s = max(asked_s, PACE_S)
            @warn "UniProt answered 429; backing off" seconds = wait_s url
            sleeper(wait_s)

        elseif status >= 500
            attempt == tries &&
                error("UniProt answered HTTP $status on all $tries attempts for $url. " *
                      "Several of these in one run is the saturation signature: stop and " *
                      "resume later rather than pushing harder.")
            sleeper(0.25 * 2.0^(attempt - 1))

        else
            error("UniProt answered HTTP $status for $url — not a transient status, so it " *
                  "is not retried. Body: $(first(String(body), 300))")
        end
    end
    error("polite_get: exhausted $tries attempts for $url")
end

# ------------------------------------------------------------------
# 4. Pass 1 — accessions to (taxid | UPI), via UniProtKB search
# ------------------------------------------------------------------
"""
    search_url(accessions; fields=PASS1_FIELDS, page=SEARCH_PAGE) -> String

The pass-1 query for one block. Pure.

`accession:A OR accession:B ...`, at most `SEARCH_BLOCK` terms — 100 is the service's cap
and it is a clause cap, not a URL-length one, so the compact `accession:(A OR B)` form does
not buy more. Parameters are emitted in a fixed order so that a repeated or resumed run
inside the 12 h `Cache-Control: public, max-age=43200` window replays CDN hits rather than
making the server work again.
"""
function search_url(accessions::AbstractVector{<:AbstractString};
                    fields::AbstractString = PASS1_FIELDS, page::Int = SEARCH_PAGE)
    length(accessions) <= SEARCH_BLOCK ||
        error("search_url: $(length(accessions)) accessions exceeds the service's " *
              "$SEARCH_BLOCK-clause cap; split the block before calling.")
    query = join(("accession:" * String(a) for a in accessions), " OR ")
    return string(REST_ROOT, "/uniprotkb/search",
                  "?query=", HTTP.escapeuri(query),
                  "&fields=", fields,
                  "&format=json",
                  "&size=", page)
end

"""
    parse_uniprotkb_search(json) -> (active, inactive, no_upi)

Partition one pass-1 payload. Pure.

  * `active::Dict{String,Int}` — accession => `organism.taxonId` for live entries.
  * `inactive::DataFrame[accession, uniparc_id, reason_type, deleted_reason]` — one row per
    `entryType == "Inactive"` entry. The UniParc id is what pass 2 consumes.
  * `no_upi::Vector{String}` — inactive entries with no UniParc pointer. A legitimate
    non-resolution rather than a loss, reported rather than guessed at.

`entryType` decides the branch before `organism` is read, so an inactive entry can never be
taken for a live one on the strength of a left-over organism field.
"""
function parse_uniprotkb_search(json)
    active = Dict{String, Int}()
    inactive = DataFrame(accession = String[], uniparc_id = Union{String, Missing}[],
                         reason_type = Union{String, Missing}[],
                         deleted_reason = Union{String, Missing}[])
    no_upi = String[]

    for item in _result_items(json)
        acc = _string_or_missing(get(item, :primaryAccession, nothing))
        acc === missing && continue
        etype = String(get(item, :entryType, ""))

        if startswith(etype, "Inactive")
            reason = get(item, :inactiveReason, nothing)
            rtype = reason === nothing ? missing :
                    _string_or_missing(get(reason, :inactiveReasonType, nothing))
            dreason = reason === nothing ? missing :
                      _string_or_missing(get(reason, :deletedReason, nothing))
            extra = get(item, :extraAttributes, nothing)
            upi = extra === nothing ? missing :
                  _string_or_missing(get(extra, :uniParcId, nothing))
            upi === missing && push!(no_upi, String(acc))
            push!(inactive, (String(acc), upi, rtype, dreason))
        else
            org = get(item, :organism, nothing)
            if org !== nothing && haskey(org, :taxonId)
                active[String(acc)] = Int(org.taxonId)
            end
        end
    end
    return (active = active, inactive = inactive, no_upi = no_upi)
end

"""
    pass1(accessions; fetch, verbose, pace_s, sleeper)
        -> (active, inactive, no_upi, release)

Resolve a whole accession list through UniProtKB search, one `SEARCH_BLOCK` at a time.

The per-block guard is the point. `x-total-results` must equal the number submitted AND the
returned accession set must equal the submitted set. Either failing means ids were lost —
on this route absence is unambiguous, because a deleted entry comes back as a first-class
`Inactive` result rather than by being missing. A short block is retried ONCE with the same
ids (these are stateless cacheable GETs; there is no poisoned job cache to escape and
splitting would only add load), and a block still short after that aborts the run.
"""
function pass1(accessions::Vector{String}; fetch = default_fetch, verbose::Bool = true,
               pace_s::Real = PACE_S, sleeper = sleep)
    active = Dict{String, Int}()
    inactive = DataFrame(accession = String[], uniparc_id = Union{String, Missing}[],
                         reason_type = Union{String, Missing}[],
                         deleted_reason = Union{String, Missing}[])
    no_upi = String[]
    release = ""

    blocks = [accessions[i:min(i + SEARCH_BLOCK - 1, end)]
              for i in 1:SEARCH_BLOCK:length(accessions)]
    verbose && println("\n🔎 PASS 1: $(length(accessions)) accessions in " *
                       "$(length(blocks)) blocks of $SEARCH_BLOCK (UniProtKB search)")

    for (n, block) in enumerate(blocks)
        part, total = nothing, nothing
        for attempt in 1:2
            (json, headers) = polite_get(search_url(block); fetch = fetch, pace_s = pace_s,
                                         sleeper = sleeper)
            rel = header_value(headers, "X-UniProt-Release")
            isempty(rel) || (release = rel)
            total = tryparse(Int, header_value(headers, "x-total-results"))
            part = parse_uniprotkb_search(json)

            seen = union(Set(keys(part.active)), Set(String.(part.inactive.accession)))
            # MEMBERSHIP, NOT CARDINALITY. Counting is not accounting: a response carrying
            # an accession we did not ask for (a secondary that resolved to a different
            # primary) holds `length(seen)` at `length(block)` while one of OUR ids is
            # absent, and a size check then calls that block complete. Ask WHICH ids came
            # back, never how many.
            asked = Set(block)
            missing_accs = [a for a in block if !(a in seen)]
            # NAMED, NOT INFERRED FROM A COUNT. Seen only as a disagreement with
            # `x-total-results`, a surplus lets the run say "0 of N missing" and then abort
            # as "short" while listing no missing id at all — a message that describes the
            # opposite of what happened and sends the reader looking for a loss that is not
            # there. The extras are right here in `seen`.
            extra_accs = sort!([a for a in seen if !(a in asked)])
            came_short = !isempty(missing_accs)
            came_over  = !isempty(extra_accs) || (total !== nothing && total != length(block))
            if !came_short && !came_over
                break
            end
            if attempt == 1
                verbose && println("   block $n: " *
                                   (came_short ?
                                    "$(length(missing_accs)) of $(length(block)) missing" :
                                    "$(length(extra_accs)) unrequested accession(s) returned") *
                                   " (x-total-results=$total); asking again for the " *
                                   "same ids after a pause.")
                sleeper(max(pace_s, 2.0))
            elseif came_short
                error("Pass 1 block $n came back short twice: $(length(missing_accs)) of " *
                      "$(length(block)) accessions neither active nor inactive in the " *
                      "response (x-total-results=$total). On this route an absent accession " *
                      "is a loss, not a deleted entry — deleted entries come back as " *
                      "Inactive. First missing: $(join(first(missing_accs, 5), ", ")).")
            else
                error("Pass 1 block $n returned accessions that were not asked for, twice: " *
                      "$(length(extra_accs)) unrequested against $(length(block)) requested " *
                      "(x-total-results=$total). Every id this block asked for did come " *
                      "back, so nothing is lost — but a response that carries ids of its " *
                      "own is not the response to this query, and taking taxonomy from it " *
                      "would attribute another accession's organism to one of ours. " *
                      "First unrequested: $(join(first(extra_accs, 5), ", ")).")
            end
        end

        merge!(active, part.active)
        append!(inactive, part.inactive)
        append!(no_upi, part.no_upi)
        verbose && n % 25 == 0 &&
            println("   block $n/$(length(blocks)): $(length(active)) active, " *
                    "$(nrow(inactive)) inactive so far")
    end

    verbose && println("   pass 1 done: $(length(active)) active, $(nrow(inactive)) inactive, " *
                       "$(length(no_upi)) inactive without a UniParc id")
    return (active = active, inactive = inactive, no_upi = no_upi, release = release)
end

# ------------------------------------------------------------------
# 5. Pass 2 — UPIs to taxid, via UniParc
# ------------------------------------------------------------------
"""
    upis_url(upis; fields=PASS2_FIELDS) -> String

The pass-2 request for one block of UniParc ids. Pure. 500 per request is measured; do not
raise it without measuring. `fields=` is not optional: an entry is ~246 B with it against
~3.6 kB without.
"""
function upis_url(upis::AbstractVector{<:AbstractString}; fields::AbstractString = PASS2_FIELDS)
    length(upis) <= UPIS_BLOCK ||
        error("upis_url: $(length(upis)) UPIs exceeds the measured $UPIS_BLOCK per request.")
    return string(REST_ROOT, "/uniparc/upis",
                  "?upis=", join(upis, ","),
                  "&fields=", fields)
end

# UniParc reports the accessions of an entry with a version suffix on some and not others
# ("A0A0Q5E2N2.1"), so every comparison strips at the dot.
_bare_accession(a) = String(first(split(String(a), '.')))

"""
    parse_uniparc_upis(json) -> Dict{String, NamedTuple}

`upi => (accessions, taxids)` for one pass-2 payload. Pure.

`taxids` is every organism on the entry. A single-element `taxids` is the answer for every
accession on that entry; more than one means the UniParc record is shared and only the
per-cross-reference lookup of pass 3 can say which taxon belongs to which accession.
"""
function parse_uniparc_upis(json)
    out = Dict{String, NamedTuple{(:accessions, :taxids), Tuple{Vector{String}, Vector{Int}}}}()
    for item in _result_items(json)
        upi = _string_or_missing(get(item, :uniParcId, nothing))
        upi === missing && continue
        accs = String[]
        for a in get(item, :uniProtKBAccessions, ())
            push!(accs, _bare_accession(a))
        end
        taxids = Int[]
        for o in get(item, :organisms, ())
            haskey(o, :taxonId) && push!(taxids, Int(o.taxonId))
        end
        out[String(upi)] = (accessions = accs, taxids = unique(taxids))
    end
    return out
end
# ------------------------------------------------------------------
# 6. Pure parsers
# ------------------------------------------------------------------
# Every payload this module reads wraps its rows in `results`: the UniProtKB search, the
# UniParc upis and id-mapping payloads, and /uniparc/<upi>/databases. A bare array is
# accepted too, so a payload can be either.
function _result_items(json)
    json === nothing && return ()
    json isa AbstractVector && return json
    return haskey(json, :results) ? json.results : ()
end

"""
    parse_uniparc_common(json) -> Dict{String, Vector{TaxonHit}}

Taxon candidates per UniParc id, from a UniParc id-mapping payload (a plain UniParc search
payload is accepted too). Pure.

A live entry states its taxonomy as `commonTaxons[{topLevel, commonTaxon, commonTaxonId}]`,
and the key is `commonTaxonId`. Reading `taxonId` off an `organisms` array — which current
payloads do not carry at all — makes the rescue branch dead code: it never fires, and every
deleted entry comes back without a taxon. The two older shapes,
`organisms[].taxonId` and `uniParcCrossReferences[].organism.taxonId`, are still parsed so
an archived payload keeps working; neither carries a top-level branch, so their hits get
`toplevel == ""`.

Candidates are de-duplicated by taxid, first mention winning.
"""
function parse_uniparc_common(json)
    out = Dict{String, Vector{TaxonHit}}()
    for item in _result_items(json)
        entry = haskey(item, :to) ? item.to : item
        upi = haskey(item, :from) ? String(item.from) :
              (haskey(entry, :uniParcId) ? String(entry.uniParcId) : "")
        isempty(upi) && continue

        hits = TaxonHit[]
        if haskey(entry, :commonTaxons)
            for c in entry.commonTaxons
                haskey(c, :commonTaxonId) || continue
                push!(hits, (taxid = Int(c.commonTaxonId),
                             name = String(get(c, :commonTaxon, "")),
                             toplevel = String(get(c, :topLevel, ""))))
            end
        end
        if haskey(entry, :organisms)
            for o in entry.organisms
                haskey(o, :taxonId) || continue
                push!(hits, (taxid = Int(o.taxonId),
                             name = String(get(o, :name, "")), toplevel = ""))
            end
        end
        if haskey(entry, :uniParcCrossReferences)
            for x in entry.uniParcCrossReferences
                (haskey(x, :organism) && haskey(x.organism, :taxonId)) || continue
                push!(hits, (taxid = Int(x.organism.taxonId),
                             name = String(get(x.organism, :scientificName, "")),
                             toplevel = ""))
            end
        end

        seen = Set{Int}()
        kept = TaxonHit[]
        for h in hits
            h.taxid in seen && continue
            push!(seen, h.taxid); push!(kept, h)
        end
        out[upi] = kept
    end
    return out
end

"""
    own_xref(db_json, accession) -> NamedTuple or nothing

The UniParc cross-reference that IS the accession being rescued, from a
`/uniparc/<upi>/databases` payload. Pure.

A UniParc entry is a sequence, not an organism: one entry routinely covers several species.
The entry's own cross-reference is therefore the only record that states which organism
*this* accession came from, which is why it outranks the common taxon. An
obsolete accession is inactive there by definition, so `active` is reported, never used to
filter. Returns `(taxid, active, created, last_updated, organism)`, or `nothing` when the
payload holds no cross-reference for the accession.
"""
function own_xref(db_json, accession::AbstractString)
    best = nothing
    for x in _result_items(db_json)
        (haskey(x, :id) && String(x.id) == accession) || continue
        (haskey(x, :organism) && haskey(x.organism, :taxonId)) || continue
        hit = (taxid = Int(x.organism.taxonId),
               active = Bool(get(x, :active, false)),
               created = String(get(x, :created, "")),
               last_updated = String(get(x, :lastUpdated, "")),
               organism = String(get(x.organism, :scientificName, "")))
        if best === nothing || (hit.active && !best.active)
            best = hit
        end
    end
    return best
end

"""
    is_cellular_toplevel(toplevel) -> Bool

Whether a commonTaxon's top-level branch is a usable organismal answer. An empty string
means the payload did not say (the legacy shapes never do) and counts as usable: absence of
evidence, not evidence of a metagenome bucket.
"""
is_cellular_toplevel(toplevel::AbstractString) =
    isempty(toplevel) || toplevel in CELLULAR_TOP_LEVELS

# Non-cellular own organism: either the common list says so for that exact taxid, or the
# organism name is an NCBI metagenome label ("<source> metagenome").
function _own_is_noncellular(own, common)
    for c in common
        c.taxid == own.taxid && return !is_cellular_toplevel(c.toplevel)
    end
    return occursin("metagenome", lowercase(String(get(own, :organism, ""))))
end

"""
    choose_taxid(own, common) -> (taxid, source, alt_cellular_taxid)

Pick one taxid for a rescued accession. Pure.

Precedence:
 1. `own` — the accession's own UniParc cross-reference — wins outright (`uniparc_own_xref`);
 2. otherwise the single cellular (or viral) commonTaxon (`uniparc_common_taxon`);
 3. otherwise, when no cellular candidate exists at all, a single non-cellular one
    (`uniparc_common_noncellular`);
 4. otherwise `missing` / `unresolved` — an ambiguous entry stays unanswered rather than
    being handed an arbitrary first-listed organism.

A metagenome taxid is KEPT rather than swapped, because it is what the source record says.
Where a cellular commonTaxon also exists it is reported beside it in `alt_cellular_taxid`,
so a downstream analysis can choose without the original being lost.

!!! warning "NOT CALLED BY THE PIPELINE"

    Nothing in the pipeline calls this, and wiring it in would break the census rather than
    improve it. Pass 2 does not resolve through commonTaxon. It resolves through UniParc
    entries that name a single organism, which is why the live vocabulary is
    `uniparc_single_organism` and why `TAXID_SOURCES` does not contain either of the two
    sources returned below. A row carrying `uniparc_common_taxon` or
    `uniparc_common_noncellular` is REFUSED by `load_taxid_cache` and fails the ledger's
    source check, so connecting this function to the assembly would not enable a safeguard -
    it would produce rows the rest of the module throws away.

    It is kept, with its tests, because it is the written record of the metagenome policy:
    a metagenome taxid is kept rather than swapped. That policy holds. What does NOT
    hold is any impression that the policy is being ENFORCED on the census: the
    `common_taxon_ids` and `alt_cellular_taxid` columns are structurally empty on every row
    (see the assembly), and a reader who takes them for a live safeguard is reading a column
    nothing writes to. `TAXID_SOURCES` pins that in a test.
"""
function choose_taxid(own, common)
    cellular = unique(h.taxid for h in common if is_cellular_toplevel(h.toplevel))
    other = unique(h.taxid for h in common if !is_cellular_toplevel(h.toplevel))

    if own !== nothing
        alt = (length(cellular) == 1 && first(cellular) != own.taxid &&
               _own_is_noncellular(own, common)) ? first(cellular) : missing
        return (taxid = own.taxid, source = "uniparc_own_xref", alt_cellular_taxid = alt)
    end
    length(cellular) == 1 &&
        return (taxid = first(cellular), source = "uniparc_common_taxon",
                alt_cellular_taxid = missing)
    (isempty(cellular) && length(other) == 1) &&
        return (taxid = first(other), source = "uniparc_common_noncellular",
                alt_cellular_taxid = missing)
    return (taxid = missing, source = "unresolved", alt_cellular_taxid = missing)
end

# ------------------------------------------------------------------
# 7. Rescue guard
# ------------------------------------------------------------------
"""
    check_rescue(n_needed, n_resolved; min_frac=0.98) -> true

Throw unless the rescue reached `min_frac` of the entries that needed one. The failure this
guards against is silent: a payload whose shape has changed rescues nothing, the run
carries on, and the missing taxa read as biology instead of as a dead parser. The caller
writes the unresolved accessions out first, so the message can be acted on.
"""
function check_rescue(n_needed::Integer, n_resolved::Integer; min_frac::Real = 0.98)
    n_needed == 0 && return true
    frac = n_resolved / n_needed
    frac >= min_frac || error(
        "UniParc rescue resolved $n_resolved of $n_needed entries " *
        "(fraction $(round(frac, digits = 4)), required $min_frac). " *
        "Check the UniParc payload shape before trusting any downstream taxonomy.")
    return true
end

# ------------------------------------------------------------------
# 8. Resumable per-UPI cross-reference cache
# ------------------------------------------------------------------
const XREF_CACHE_COLS = ["uniparc_id", "accession", "taxid", "active", "created",
                         "last_updated", "organism"]

_xref_key(upi, acc) = string(upi, "|", acc)

# A cached row with an empty taxid records a UPI that was looked up and had no own
# cross-reference, so a resumed run does not ask for it again.
function load_xref_cache(path::Union{String, Nothing})
    cache = Dict{String, Any}()
    (path === nothing || !isfile(path)) && return cache
    # The dates are stored as the strings UniProt gave, so the column types are pinned:
    # left to itself the CSV reader turns them into Dates and the round trip stops matching.
    types = Dict("uniparc_id" => String, "accession" => String, "taxid" => Int,
                 "active" => Bool, "created" => String, "last_updated" => String,
                 "organism" => String)
    for row in CSV.File(path; types = types)
        key = _xref_key(row.uniparc_id, row.accession)
        cache[key] = ismissing(row.taxid) ? nothing :
            (taxid = Int(row.taxid), active = coalesce(row.active, false),
             created = coalesce(row.created, ""),
             last_updated = coalesce(row.last_updated, ""),
             organism = coalesce(row.organism, ""))
    end
    return cache
end

function append_xref_cache(path::Union{String, Nothing}, upi, acc, hit)
    path === nothing && return nothing
    fresh = !isfile(path)
    open(path, "a") do io
        fresh && println(io, join(XREF_CACHE_COLS, ","))
        println(io, join([upi, acc,
                          hit === nothing ? "" : string(hit.taxid),
                          hit === nothing ? "" : string(hit.active),
                          hit === nothing ? "" : hit.created,
                          hit === nothing ? "" : hit.last_updated,
                          hit === nothing ? "" : replace(hit.organism, ',' => ' ')], ","))
    end
    return nothing
end

"""
    fetch_own_xref(upi, accession; fetch=default_fetch) -> (hit, complete)

Fetch `/uniparc/<upi>/databases` and pick the accession's own cross-reference. The `?id=`
filter keeps the payload to this accession's cross-reference; `Link` headers are followed,
so an answer that spans more than one page still yields it.

WHY THIS REPORTS COMPLETENESS AND NOT ONLY A HIT. A miss means one of two different things
and they may not be treated alike. Either the entry genuinely holds no cross-reference for
this accession — a fact about the data, and cacheable forever — or this particular response
came back short. A short answer here is NOT an error: the service returns HTTP 200 with
fewer records than it says it holds, so `polite_get` has nothing to retry and the paging
loop simply ends early. Cached, that becomes a permanent "no cross-reference" for an
accession that has one, and no later run ever asks again.

So the count of records actually seen is compared against the `x-total-results` the service
declares, exactly as pass 1 compares them, and `complete` says whether the miss is a
statement about the data. Only a complete miss may be written to the cache. A hit needs no
such test: having found the cross-reference, how many other rows the page held is moot.
"""
function fetch_own_xref(upi::AbstractString, accession::AbstractString;
                        fetch = default_fetch, pace_s::Real = PASS3_PACE_S, sleeper = sleep)
    # `?id=` filters on the accession itself, so the response is the one cross-reference we
    # care about. Preferred over `?dbTypes=UniProtKB/TrEMBL`, which assumes the accession
    # sits in TrEMBL rather than Swiss-Prot — an assumption nothing here can guarantee.
    url = "$REST_ROOT/uniparc/$upi/databases?format=json&size=500" *
          "&id=$(HTTP.escapeuri(accession))"
    seen = 0
    declared = nothing
    while !isempty(url)
        (json, headers) = polite_get(url; fetch = fetch, pace_s = pace_s, sleeper = sleeper)
        # Read once, from the first page: the header restates the same total on every page.
        declared === nothing &&
            (declared = tryparse(Int, header_value(headers, "x-total-results")))
        seen += length(_result_items(json))
        hit = own_xref(json, accession)
        hit === nothing || return (hit = hit, complete = true)
        url = _next_link(header_value(headers, "Link"))
    end
    # A miss the service did not account for is not a fact about the data. If the header is
    # absent the answer cannot be shown complete either, and is treated as short rather than
    # assumed whole — the cost of that is one repeated request, and the cost of assuming is a
    # wrong negative that outlives every run.
    return (hit = nothing, complete = declared !== nothing && seen == declared)
end

# `Link: <url>; rel="next"` — the only rel the results endpoints emit.
function _next_link(link_header::AbstractString)
    isempty(link_header) && return ""
    m = match(r"<([^>]+)>\s*;\s*rel=\"next\"", link_header)
    return m === nothing ? "" : String(m.captures[1])
end

# ------------------------------------------------------------------
# 9. Main entry point
# ------------------------------------------------------------------
function empty_taxid_table()
    DataFrame(accession = String[], TaxID = Union{Int, Missing}[],
              taxid_source = String[], entry_status = String[],
              inactive_reason = Union{String, Missing}[],
              uniparc_id = Union{String, Missing}[], common_taxon_ids = String[],
              alt_cellular_taxid = Union{Int, Missing}[],
              xref_created = Union{String, Missing}[],
              xref_last_updated = Union{String, Missing}[], uniprot_release = String[])
end

_join_taxids(hits) = join((string(h.taxid) for h in hits), ";")

function _reason_text(rtype, dreason)
    ismissing(rtype) && return ismissing(dreason) ? missing : String(dreason)
    return ismissing(dreason) ? String(rtype) : string(rtype, ": ", dreason)
end

# ------------------------------------------------------------------
# 9b. The accession-keyed cache — layer 1
# ------------------------------------------------------------------
#=
    The (UPI, accession) cache above serves pass 3, which is the expensive step. What it
    cannot do is spare a FRESH run any work, because it cannot be consulted until
    `accession -> UPI` has already been re-derived over the network.

    This layer is keyed on the bare accession and nothing else, so a re-run is a local join
    and costs no HTTP at all while the release is unchanged.

    `source` records which request shape produced each taxon, so one class can be
    re-derived later without re-deriving all of them. `uniprot_release` is the ageing
    handle: taxids are merged and deleted upstream, and a cache without a release stamp
    cannot be invalidated honestly.

    One trap, because it would be easy and wrong: do NOT seed this from a census's own
    filled TaxID column. That would confirm the census's numbers with the census's
    numbers. Rows of uniparc_xref_cache.csv are fetched independently and promote freely;
    a census's own values are a validation target, not an input.
=#
const TAXID_CACHE_COLS = ["accession", "taxid", "source", "uniparc_id", "active",
                          "uniprot_release", "fetched_at",
                          # CARRIED SO THAT A CACHE-WARM RE-RUN DESCRIBES THE SAME CENSUS.
                          # These three are not resolutions, they are what the resolution
                          # SAYS, and the table this cache feeds carries them. Pass 1 is not
                          # re-run for a cached accession, so without them a warm re-run
                          # rebuilds the same taxa with three empty columns and calls it the
                          # same census. A file without these columns still loads: the
                          # reader keys off the header it finds.
                          "inactive_reason", "xref_created", "xref_last_updated"]

# "unresolved" is deliberately NOT a member. A failure to resolve is not a fact about an
# accession, it is the absence of one; caching it makes a transient failure permanent and
# invisible, because `todo` then skips that accession on every later run and no guard ever
# sees it again. Only established resolutions are cacheable.
const TAXID_SOURCES = ("uniprotkb_active", "uniparc_single_organism", "uniparc_own_xref")

# The two sources that mean "pass 1 called this Inactive and a later pass rescued it". A
# cached row carrying one of these both NEEDED a rescue and RECEIVED one, which is why the
# rescue floor has to count it on both sides.
const RESCUE_SOURCES = ("uniparc_single_organism", "uniparc_own_xref")

"""
    load_taxid_cache(path) -> Dict{String, NamedTuple}

Accession => cached resolution. Missing file or `nothing` gives an empty cache.
"""
function load_taxid_cache(path::Union{String, Nothing})
    out = Dict{String, NamedTuple{(:taxid, :source, :uniparc_id, :release,
                                   :reason, :created, :updated),
                                  Tuple{Union{Int, Missing}, String,
                                        Union{String, Missing}, String,
                                        Union{String, Missing}, Union{String, Missing},
                                        Union{String, Missing}}}}()
    (path === nothing || !isfile(path)) && return out
    # WHAT THE HEADER ACTUALLY HOLDS, not what the writer emits. A cache without the three
    # description columns is still a cache of valid resolutions; refusing it would throw
    # away paid-for work over a schema difference.
    hdr = split(strip(first(eachline(path))), ',')
    _opt(row, name) = name in hdr ?
        (getproperty(row, Symbol(name)) === missing ||
         isempty(getproperty(row, Symbol(name))) ? missing :
         String(getproperty(row, Symbol(name)))) : missing
    n_refused = 0
    for row in CSV.Rows(path; types = String, reusebuffer = false)
        acc = row.accession
        (acc === missing || isempty(acc)) && continue
        tx = (row.taxid === missing || isempty(row.taxid)) ? missing :
             something(tryparse(Int, row.taxid), missing)
        src = row.source === missing ? "" : String(row.source)
        # REFUSED, NOT LOADED. A row with no TaxID, or with a source outside the allowed
        # set, is never written by this module; such a row can only come from a writer that
        # caches failures. Replaying it would reinstate a non-answer that `todo` then
        # skips and no guard can see, so it is dropped here and counted out loud.
        if tx === missing || !(src in TAXID_SOURCES)
            n_refused += 1
            continue
        end
        out[String(acc)] = (taxid = tx,
                            source = src,
                            uniparc_id = (row.uniparc_id === missing ||
                                          isempty(row.uniparc_id)) ? missing :
                                         String(row.uniparc_id),
                            release = row.uniprot_release === missing ? "" :
                                      String(row.uniprot_release),
                            reason = _opt(row, "inactive_reason"),
                            created = _opt(row, "xref_created"),
                            updated = _opt(row, "xref_last_updated"))
    end
    n_refused > 0 && @warn("Accession cache held rows that are not established " *
                           "resolutions; they were refused, not loaded. Those accessions " *
                           "will be fetched again.", path, n_refused)
    return out
end

"""
    append_taxid_cache(path, rows) -> nothing

Append resolutions, creating the file with a header when it does not exist. Append-only.

Called ONCE, at the end of a successful `resolve_taxids`, never per block: a run that fails
a guard must not leave rows behind for the next run to read. An interrupted run therefore
loses its pass-1/2 bookkeeping, but not the expensive part - `append_xref_cache` still
writes per accession as pass 3 goes, so a resumed run does not re-fetch cross-references.
"""
function append_taxid_cache(path::Union{String, Nothing}, rows::AbstractVector)
    (path === nothing || isempty(rows)) && return nothing
    mkpath(dirname(path))
    fresh = !isfile(path)
    open(path, "a") do io
        fresh && println(io, join(TAXID_CACHE_COLS, ","))
        for r in rows
            # `get` and not `r.field`: promote_xref_cache builds rows from the xref layer
            # and carries no reason. A missing value is written empty, never guessed.
            reason = get(r, :inactive_reason, missing)
            println(io, join((r.accession,
                              r.taxid === missing ? "" : string(r.taxid),
                              r.source,
                              r.uniparc_id === missing ? "" : r.uniparc_id,
                              string(r.active),
                              r.release,
                              r.fetched_at,
                              reason === missing ? "" : replace(String(reason), ',' => ';'),
                              get(r, :xref_created, missing) === missing ? "" :
                                  String(get(r, :xref_created, "")),
                              get(r, :xref_last_updated, missing) === missing ? "" :
                                  String(get(r, :xref_last_updated, ""))), ","))
        end
    end
    return nothing
end

"""
    promote_xref_cache(xref_path, taxid_path; release) -> Int

Lift the rows the (UPI, accession) cache already holds into the accession-keyed layer, so
work paid for on an earlier run is not paid for again. No HTTP. Returns the number promoted.
"""
function promote_xref_cache(xref_path::Union{String, Nothing},
                            taxid_path::Union{String, Nothing};
                            release::AbstractString = "", now_str::AbstractString = "")
    (xref_path === nothing || taxid_path === nothing || !isfile(xref_path)) && return 0
    have = load_taxid_cache(taxid_path)
    rows = []
    for (key, hit) in load_xref_cache(xref_path)
        # The xref cache is keyed by the STRING "<upi>|<accession>". Destructuring that
        # as a pair walks its CHARACTERS instead, which promotes one-character accessions
        # and throws on every no-xref row (those are stored as `nothing`, no `.taxid`).
        parts = split(key, '|'; limit = 2)
        length(parts) == 2 || continue
        upi, acc = String(parts[1]), String(parts[2])
        (acc in keys(have)) && continue
        hit === nothing && continue
        push!(rows, (accession = acc, taxid = hit.taxid, source = "uniparc_own_xref",
                     uniparc_id = upi, active = false, release = release,
                     fetched_at = now_str,
                     xref_created = get(hit, :created, missing),
                     xref_last_updated = get(hit, :last_updated, missing)))
    end
    append_taxid_cache(taxid_path, rows)
    return length(rows)
end

# ------------------------------------------------------------------
# 10. resolve_taxids — the three passes, guarded end to end
# ------------------------------------------------------------------
"""
    resolve_taxids(accessions; fetch=default_fetch, min_rescue_frac=0.98,
                   xref_cache=nothing, taxid_cache=nothing, unresolved_csv=nothing,
                   pace_s=PACE_S, pass3_pace_s=PASS3_PACE_S, verbose=true) -> DataFrame

Map UniProt accessions to NCBI taxids over three stateless passes. No id-mapping jobs, no
stream endpoints, no parallelism — see the header block for why.

  1. `/uniprotkb/search` resolves live entries outright and hands back a UniParc id for
     every deleted one, 100 accessions at a time.
  2. `/uniparc/upis` turns those UniParc ids into taxa, 500 at a time. A taxon is accepted
     only where the entry carries exactly ONE organism.
  3. Entries spanning several organisms get one `/uniparc/<upi>/databases?id=<acc>` each,
     which returns the cross-reference for that accession alone and so its own taxon.

`commonTaxonId` is never used as a fallback: in spot checks it did not name the accession's
organism, in one case naming a whole class where the truth was a strain.

Returns one row per unique accession with `accession, TaxID, taxid_source, entry_status,
inactive_reason, uniparc_id, common_taxon_ids, alt_cellular_taxid, xref_created,
xref_last_updated, uniprot_release` — one contract, whichever pass resolved the row.
"""
function resolve_taxids(accessions::AbstractVector{<:AbstractString};
                        fetch = default_fetch,
                        min_rescue_frac::Real = 0.98,
                        xref_cache::Union{String, Nothing} = nothing,
                        taxid_cache::Union{String, Nothing} = nothing,
                        unresolved_csv::Union{String, Nothing} = nothing,
                        pace_s::Real = PACE_S,
                        pass3_pace_s::Real = PASS3_PACE_S,
                        sleeper = sleep,
                        verbose::Bool = true)
    uniq = unique(String.(accessions))
    isempty(uniq) && return empty_taxid_table()
    # ISO-ish stamp without taking a dependency on Dates for one string.
    stamp = Libc.strftime("%Y-%m-%dT%H:%M:%S", time())

    # --- Layer 1: anything already resolved on an earlier run costs nothing ---
    cached = load_taxid_cache(taxid_cache)
    todo = [a for a in uniq if !haskey(cached, a)]
    verbose && !isempty(cached) &&
        println("\n💾 accession cache: $(length(uniq) - length(todo)) of $(length(uniq)) " *
                "already resolved; $(length(todo)) to fetch")

    # --- Pass 1 ---
    p1 = isempty(todo) ?
         (active = Dict{String, Int}(),
          inactive = DataFrame(accession = String[], uniparc_id = Union{String, Missing}[],
                               reason_type = Union{String, Missing}[],
                               deleted_reason = Union{String, Missing}[]),
          no_upi = String[], release = "") :
         pass1(todo; fetch = fetch, verbose = verbose, pace_s = pace_s, sleeper = sleeper)
    release = p1.release

    # --- Pass 2: UniParc, single-organism entries answered here ---
    upi_of = Dict{String, String}()
    for row in eachrow(p1.inactive)
        ismissing(row.uniparc_id) || (upi_of[String(row.accession)] = String(row.uniparc_id))
    end
    upis = unique(collect(values(upi_of)))
    single = Dict{String, Int}()       # accession => taxid, from a one-organism entry
    multi = String[]                   # accessions needing pass 3
    if !isempty(upis)
        verbose && println("\n🧬 PASS 2: $(length(upis)) UniParc ids in " *
                           "$(cld(length(upis), UPIS_BLOCK)) blocks of $UPIS_BLOCK")
        entries = Dict{String, Any}()
        for i in 1:UPIS_BLOCK:length(upis)
            block = upis[i:min(i + UPIS_BLOCK - 1, end)]
            (json, headers) = polite_get(upis_url(block); fetch = fetch, pace_s = pace_s,
                                         sleeper = sleeper)
            rel = header_value(headers, "X-UniProt-Release")
            isempty(rel) || (release = rel)
            part = parse_uniparc_upis(json)
            # Every submitted UPI must come back. Do NOT count accessions here: several
            # census accessions can legitimately share one UniParc entry, so fewer entries
            # than accessions is correct, not short.
            gone = [u for u in block if !haskey(part, u)]
            isempty(gone) ||
                error("Pass 2 lost $(length(gone)) of $(length(block)) UniParc ids: " *
                      "$(join(first(gone, 5), ", ")).")
            merge!(entries, part)
        end
        for (acc, upi) in upi_of
            e = get(entries, upi, nothing)
            e === nothing && continue
            # Require the accession to actually reappear on the entry. UniParc indexes
            # other databases' ids too, so without this an EMBL/RefSeq collision could
            # slip a wrong taxon in. No such collision has been observed, and none is
            # specifically tested.
            acc in e.accessions || continue
            if length(e.taxids) == 1
                single[acc] = e.taxids[1]
            else
                push!(multi, acc)
            end
        end
        verbose && println("   pass 2 done: $(length(single)) resolved from single-organism " *
                           "entries, $(length(multi)) need a per-accession lookup")
    end

    # --- Pass 3: one request per accession whose UniParc entry is shared ---
    own = Dict{String, Any}()
    n_xref_short = 0
    xcache = load_xref_cache(xref_cache)
    if !isempty(multi)
        verbose && println("\n🔬 PASS 3: $(length(multi)) per-accession cross-references")
        for (n, acc) in enumerate(multi)
            upi = upi_of[acc]
            key = _xref_key(upi, acc)
            if haskey(xcache, key)
                own[acc] = xcache[key]
                continue
            end
            answer = fetch_own_xref(upi, acc; fetch = fetch, pace_s = pass3_pace_s,
                                    sleeper = sleeper)
            own[acc] = answer.hit
            # ONLY A COMPLETE ANSWER MAY BE CACHED. A short response and a genuine absence
            # look identical here — both leave the accession unresolved for this run — but
            # only one of them is true tomorrow. Caching the other is how a transient gap
            # becomes permanent, silently, for every run that follows.
            if answer.complete
                append_xref_cache(xref_cache, upi, acc, answer.hit)
            else
                n_xref_short += 1
            end
            verbose && n % 250 == 0 && println("   $n/$(length(multi))")
        end
        # SAID OUT LOUD, because a silent counter is the same defect one layer up. These
        # accessions were asked about and the answer could not be shown complete, so nothing
        # was written for them: they will be asked again next run, which is the point. A
        # number that stays high across runs is the service answering short repeatedly, not
        # a census with missing cross-references.
        n_xref_short > 0 &&
            @warn("$n_xref_short of $(length(multi)) cross-reference lookups came back " *
                  "without an accountable total and were NOT cached; they are re-asked on " *
                  "the next run rather than remembered as absent")
    end

    # --- Assemble one row per accession, and a ledger that must add up ---
    by_acc = Dict(String(row.accession) => row for row in eachrow(p1.inactive))
    df = empty_taxid_table()
    fresh_rows = []
    n_active = n_single = n_own = n_unresolved = n_cached = 0

    for acc in uniq
        taxid = missing; source = "unresolved"; status = "unmapped"
        # Per row, not per run: a cached resolution carries the release it was established
        # under, and this run's release says nothing about it.
        row_release = release
        upi = get(upi_of, acc, missing)
        row = get(by_acc, acc, nothing)
        reason = row === nothing ? missing : _reason_text(row.reason_type, row.deleted_reason)
        created = missing; updated = missing

        if haskey(cached, acc)
            c = cached[acc]
            taxid = c.taxid; source = c.source; upi = c.uniparc_id
            # A cached row's release is the only release that describes it. If the row
            # carries none, the column stays EMPTY: stamping this run's release onto a
            # resolution established under an unknown one is a provenance claim nothing
            # supports, and it is the kind of quiet false stamp this module exists to stop.
            row_release = c.release
            status = c.source == "uniprotkb_active" ? "active" : "inactive"
            # WHAT THE ROW SAYS, restored with it. Pass 1 is not re-run for a cached
            # accession, so `by_acc` does not hold it and `reason` would otherwise be
            # blank on every warm run - the same row, described differently, which makes a
            # warm re-run a different census.
            c.reason === missing || (reason = c.reason)
            created = c.created
            updated = c.updated
            # load_taxid_cache refuses TaxID-less rows, so this cannot hold unless the
            # cache reached this point by some other route. Fail loudly rather than count a
            # non-answer as a cache hit - counted as a hit, it lets the ledger balance while
            # rows are missing.
            taxid === missing && error("Cache hit for $acc carries no TaxID. A cached row " *
                                       "must be an established resolution; delete the " *
                                       "cache file and re-run.")
            n_cached += 1
        elseif haskey(p1.active, acc)
            taxid = p1.active[acc]; source = "uniprotkb_active"; status = "active"
            n_active += 1
        elseif haskey(single, acc)
            taxid = single[acc]; source = "uniparc_single_organism"; status = "inactive"
            n_single += 1
        elseif haskey(own, acc)
            hit = own[acc]
            if hit !== nothing && hit.taxid !== missing
                taxid = hit.taxid; source = "uniparc_own_xref"; status = "inactive"
                created = hit.created; updated = hit.last_updated
                n_own += 1
            else
                n_unresolved += 1
            end
        else
            n_unresolved += 1
        end

        push!(df, (accession = acc, TaxID = taxid, taxid_source = source,
                   entry_status = status, inactive_reason = reason, uniparc_id = upi,
                   # STRUCTURALLY EMPTY, ON EVERY ROW, AND NOT AN OVERSIGHT TO BE FIXED
                   # BY FILLING THEM IN. They are the surface of the commonTaxon route (see
                   # `choose_taxid`), which pass 2 does not use. Populating them means
                   # resolving through that route, and its sources are not in
                   # `TAXID_SOURCES`. Anything reading these columns as a metagenome
                   # safeguard is reading a column nothing writes.
                   common_taxon_ids = "", alt_cellular_taxid = missing,
                   xref_created = created, xref_last_updated = updated,
                   uniprot_release = row_release))

        # ONLY ESTABLISHED RESOLUTIONS ARE CACHED. An unresolved accession goes to
        # `unresolved_csv` for the operator, never to the cache: cached, the next run skips
        # it at `todo`, reinstates the non-answer, and reaches the census with no guard
        # having seen it.
        (haskey(cached, acc) || taxid === missing) || push!(fresh_rows,
            (accession = acc, taxid = taxid, source = source, uniparc_id = upi,
             active = status == "active", release = row_release, fetched_at = stamp,
             inactive_reason = reason, xref_created = created,
             xref_last_updated = updated))
    end

    # --- The guard: every accession accounted for by a row the code actually saw ---
    #
    # THE SUM BELOW CANNOT FAIL, AND IS KEPT ONLY BECAUSE ITS MESSAGE NAMES THE TERMS. The
    # assembly loop runs exactly once per member of `uniq` and every branch of it increments
    # exactly one of these counters, so their total restates the loop's own bound. On its
    # own it would report a balanced ledger over a table that had lost half its rows. The
    # checks that follow compare the counters against the TABLE, which is built by a
    # different statement, and those can disagree.
    total = n_cached + n_active + n_single + n_own + n_unresolved
    if total != length(uniq)
        error("Ledger does not balance: $(length(uniq)) accessions in, but " *
              "$n_cached cached + $n_active active + $n_single uniparc-single + " *
              "$n_own uniparc-own-xref + $n_unresolved unresolved = $total. " *
              "Every term must count rows that were seen, never a subtraction.")
    end
    # One row per accession asked about. A `push!` that did not happen, or happened twice,
    # is invisible to the counters and fatal to everything downstream that joins on it.
    if nrow(df) != length(uniq)
        error("Ledger does not balance: $(length(uniq)) accessions in, but the table holds " *
              "$(nrow(df)) rows. The counters agree with themselves; they do not agree " *
              "with the table, so a row was dropped or written twice.")
    end
    if length(unique(df.accession)) != nrow(df)
        error("The table holds $(nrow(df)) rows for " *
              "$(length(unique(df.accession))) distinct accessions: an accession was " *
              "written more than once, and a join on it would multiply rows silently.")
    end
    # An unresolved accession is exactly a row with no TaxID. Counted one way by the loop
    # and another way by the table; if they part company, one of them is lying about the
    # census and the rescue floor below is measured over the wrong population.
    n_missing = count(ismissing, df.TaxID)
    if n_missing != n_unresolved
        error("The loop counted $n_unresolved unresolved accessions and the table holds " *
              "$n_missing rows with no TaxID. These must be the same accessions; a gap " *
              "here means a row carries a TaxID the ledger never credited, or the reverse.")
    end
    if nrow(df) - n_missing != n_cached + n_active + n_single + n_own
        error("The table holds $(nrow(df) - n_missing) resolved rows against " *
              "$(n_cached + n_active + n_single + n_own) credited to the four resolving " *
              "routes. A resolution the ledger cannot source is a resolution nothing checked.")
    end
    verbose && println("\n📒 Ledger: $(length(uniq)) = $n_cached cached + $n_active active + " *
                       "$n_single uniparc-single + $n_own own-xref + $n_unresolved unresolved " *
                       "(UniProt release $(isempty(release) ? "unknown" : release))")

    unresolved = [String(r.accession) for r in eachrow(df) if ismissing(r.TaxID)]
    if unresolved_csv !== nothing && !isempty(unresolved)
        mkpath(dirname(unresolved_csv))
        CSV.write(unresolved_csv, DataFrame(accession = unresolved))
        verbose && println("📝 Wrote $(length(unresolved)) unresolved accessions to $unresolved_csv")
    end

    # The rescue floor applies to passes 2+3 jointly: of what pass 1 called Inactive, at
    # least min_rescue_frac must end with a taxon. Counted from what NEEDED a rescue, not
    # from what we found a pointer for — otherwise a payload whose shape changed empties
    # the numerator and denominator together and a floor of 0/0 passes while the rescue
    # silently returns nothing.
    #
    # THE FLOOR MEASURES THE CENSUS, NOT THIS RUN'S SLICE. Only resolutions are cached, so
    # after a successful run `todo` holds exactly the accessions that FAILED. A floor taken
    # over `todo` alone is then taken over a population that is 100% known failures: it is
    # 0/k by construction and aborts a perfectly healthy re-run, on unchanged input, while
    # blaming the UniParc payload. A cached rescue-sourced row both needed a rescue and got
    # one, so it belongs in BOTH terms and the fraction stays a statement about the census.
    n_cached_rescued = count(a -> haskey(cached, a) && cached[a].source in RESCUE_SOURCES,
                             uniq)
    n_needed = nrow(p1.inactive) + n_cached_rescued
    n_rescued = n_single + n_own + n_cached_rescued
    n_needed > 0 && check_rescue(n_needed, n_rescued; min_frac = min_rescue_frac)

    # THE CACHE IS WRITTEN LAST, ONLY ONCE EVERY GUARD ABOVE HAS PASSED. Written earlier, a
    # run that aborts still leaves its rows behind; the NEXT run then reads them, finds
    # `todo` empty, skips every pass and every guard, and publishes the short answer the
    # abort existed to prevent - an aborting run converted into a passing one by re-running
    # it, which is the id-mapping service's own failure mode (a short answer cached against
    # the request) moved inside this module's cache.
    append_taxid_cache(taxid_cache, fresh_rows)

    return df
end

"""
    fetch_taxids_async(accessions; kwargs...) -> Dict{String, Int}

Back-compat wrapper over `resolve_taxids`: accession => TaxID for every accession that
resolved. The provenance columns are dropped, so prefer `resolve_taxids` in new code.
"""
function fetch_taxids_async(accessions::AbstractVector{<:AbstractString}; kwargs...)
    df = resolve_taxids(accessions; kwargs...)
    out = Dict{String, Int}()
    for row in eachrow(df)
        ismissing(row.TaxID) || (out[String(row.accession)] = Int(row.TaxID))
    end
    return out
end

end
