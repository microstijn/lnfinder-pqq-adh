# TOOLING — what the code in this repository does, and how

This is the method documentation for LnFinder: the algorithms, the order they run in, the
parameters and what each one means, and the external tools and services involved. It
describes the **method**, not any result obtained with it. No count, effect size, accuracy
or fitted threshold from a run appears anywhere below; where a default is calibrated from
data, this document says how it is calibrated and where the calibrated value is recorded,
not what it came out as.

For what can and cannot be re-run from this repository, see `README.md`. That question is
answered there once, honestly, and is not repeated here.

---

## 0. Software environment

| | |
|---|---|
| Language | Julia, floor declared in `Project.toml` as `[compat] julia = "1.12"` |
| Dependencies | declared in `Project.toml`, resolved exactly in `Manifest.toml` |
| Install | `julia --project=. -e 'using Pkg; Pkg.instantiate()'` |
| Tests | `julia --project=. test/runtests.jl` — offline, no external binary |

`Manifest.toml` ships on purpose: it pins the exact package closure, which is the
reproducibility claim a code release exists to make. It was resolved on Julia 1.12, and on
an older minor version `instantiate` fails on the standard-library set rather than on
anything in this repository. A reader who wants a different Julia should delete
`Manifest.toml` and re-resolve, knowing that the package versions then float.

Nothing in this repository requires a Julia version manager. Every invocation is written
`julia --project=. …`, which works with Julia installed by any route; the test suite fails
the build if a version-manager-only form (`julia +<channel> …`) appears anywhere in
`src/`, `scripts/` or `test/`, or in the prose and CI files a reader copies commands from
(`README.md`, `docs/TOOLING.md`, `docs/DATA.md`, `dat/models/PROVENANCE.md` and the CI
workflow).

---

## 1. The cofactor call

### 1.1 The equivalence being applied

PQQ-dependent alcohol dehydrogenases of the methanol/ethanol family occur in two forms that
differ in the metal held in the active site: a lanthanide ion or a calcium ion. Lanthanide
ions are larger and take a higher coordination number, and the enzymes that bind them
supply the extra coordination with an additional aspartate in the active site. That
difference is legible in sequence:

| active-site residues at the calibrated column | call |
|---|---|
| `D-x-D` (Asp, any, Asp) | lanthanide |
| `D-x-[ATS]` (Asp, any, Ala/Thr/Ser) | calcium |
| a gap at any of the three positions | gap — no call |
| anything else | unknown — no call |

A sequence whose model string ends before the column is reported as truncated, which is no
call either. At the calibrated column the rule is applied residue by residue by
`src/ModuleHmmer.jl`: `calibrate_active_site` finds the column (§1.3) and
`parse_hmmalign_results` reads every aligned sequence there. The same two patterns are
also written in PROSITE notation and matched with `BioSequences.jl` in
`src/ModuleMotifSearch.jl` (`classify_pqq_adh`), which searches a whole sequence and so
answers a different question: whether the motif occurs anywhere, not what sits at the
calibrated column.

The rule itself is not this project's; it is the published lanthanide-switch equivalence
(Voutsinos et al. 2025, doi:10.1093/ismejo/wraf057). The structural basis the calibration
relies on — an additional aspartate two positions downstream of a conserved one, where the
calcium enzymes carry Ala, Ser or Thr — is cited in `src/ModuleHmmer.jl` to Good et al. 2020
(doi:10.1074/jbc.RA120.013227). What this repository adds is (a) reading it at a
**calibrated model column** rather than at a sequence position, (b) a validation harness
that measures whether the rule recovers the cofactor on enzymes whose cofactor is known,
and (c) the tooling built around the call.

### 1.2 Model coordinates, not sequence coordinates

A residue's position differs between sequences, so the call is read at a column of the
profile HMM rather than at an offset in a sequence. Sequences are aligned to the profile
with HMMER `hmmalign` in `afa` output format; in that format, upper-case characters and
`-` are match-state (model) columns and lower-case characters are insertions relative to
the model. `src/ModuleHmmer.jl` strips insertions, which makes the remaining index the
model coordinate, identical across every aligned sequence.

The profile itself is `dat/models/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm`, third-party
and documented in `dat/models/PROVENANCE.md`. It carries no GA/TC/NC cutoff, so it defines
no inclusion threshold of its own: a threshold is a property of the search command, not of
the profile. The background search calibrates its threshold and records it with the run
(§4.2).

### 1.3 Calibrating the active-site column

`calibrate_active_site` takes the first two sequences of the alignment — one enzyme known
to bind a lanthanide, one known to bind calcium — and returns the first model column at
which the first reads `D-x-D` and the second reads `D-x-[ATS]`. Which two records those are
is the caller's decision, taken in one of two ways and never read from a data file. The
census driver (`src/scriptLn.jl`) calls `calibrate_active_site(paths.afa)`, so it calibrates
on whatever the first two records of its alignment are, and the census FASTA carries the
pair as its first two records. The validation harness is told the two accessions instead,
as the `ln_ref_acc` and `ca_ref_acc` keywords of `validate` and `validate_from_afa`
(`src/validate_classifier.jl`); `reorder_references!` then moves those two records to the
top before calibrating, so record order cannot decide what the run calibrated on (§1.4).
The background scan takes that same route with the pair fixed in its own source (`LN_REF`,
`CA_REF`, §11). Either way a run states what it calibrated on.

This is a **calibrated parameter, not a constant**. Its value is a property of the profile
and the two reference sequences used, it is computed at the start of every run, and it is
written into that run's output. It is deliberately not hard-coded anywhere in this
repository.

The procedure has a known fragility, and the code reports it rather than hiding it:
if more than one model column satisfies both reference patterns, the first is not
necessarily the catalytic one, and every downstream call inherits the choice.
`diagnose_calibration` (`src/validate_classifier.jl`) counts those "collision" columns and
reports them, so a caller can see whether the choice was unique before trusting anything
that rests on it.

### 1.4 Control-set validation

The call is validated against `dat/control_set/known_cofactor_mdhs.tsv`, a curated table in
which every row is a factual claim about a published enzyme: accession, gene, organism,
domain, the metal, a one-line statement of the evidence, the citation it comes from, any
PDB entries holding that metal, and the provenance of the row. Its header documents the
schema in full; the parts that govern the procedure are these:

**Roles.** Rows marked `reference` are calibration references. The two a run calibrates on
are passed to it by accession (`ln_ref_acc`, `ca_ref_acc`) and locate the active-site
column; any further reference row takes no part in the calibration. No reference row enters
a held-out stratum, so none can flatter the held-out result: they are reported apart, in a
`reference` stratum of their own and in the evaluation over all rows. Rows marked `test` are
held out, and only those enter a held-out accuracy.

**Evidence tiers**, strongest first:

| tier | meaning |
|---|---|
| `ref` | a reference row; belongs to no held-out stratum |
| `A` | direct determination for this protein — the metal resolved in a crystal structure of the enzyme itself, or measured in (or required by) the purified enzyme |
| `B` | indirect but gene specific — expression, deletion mutant or growth evidence ties this gene to a Ca- or Ln-dependent activity, but no metal was measured in, reconstituted into or shown to be required by the purified enzyme, or the enzyme was not purified at all |
| `C` | homology only — the cofactor is inherited from close relatives and was never measured for this protein |
| `X` | contested — reported with its call and excluded from every accuracy |

The tiers exist so that rows of different evidential weight are scored apart instead of
averaged together, and so that a downstream analysis that must not be told the answer by
homology can exclude the homology tier explicitly (see §9.1).

**The run.** `ValidateClassifier.validate` loads the table, fetches each sequence from its
database of record (UniProtKB by default; a row's `source` column can route it to another,
and a deleted UniProtKB entry is fetched from UniParc under the UniParc id UniProt reports
for it, which the row's evidence cell also names), aligns them to the profile with
`hmmalign`, checks that the alignment's match-column
count equals the profile's `LENG` (an alignment in foreign coordinates cannot be scored and
is refused rather than scored wrongly), moves the two reference rows to the top, runs the
calibration diagnostic, calibrates, reads the call at that column for every row, and scores
it. `validate_from_afa` is the same tail starting from an alignment made elsewhere, which
takes HMMER off the critical path entirely.

**What is reported.** Sensitivity, specificity and accuracy over all scorable rows and over
the held-out rows alone, plus a per-stratum summary with exact (binomial) confidence
intervals. A gap, unknown or truncated call counts as a misclassification of the row's
known cofactor, and a contested row is listed with the call it received and never scored.
Every run writes its own fetch manifest — which accession resolved to what, from
which database, with the release header the service reported — beside the predictions, so a
later re-run can be compared against it rather than against memory.

**Near-identical held-out rows.** A control set drawn from the literature contains close
relatives, and counting them as independent overstates the evidence. The harness clusters
the held-out sequences by single linkage at a pairwise identity over the match columns —
0.9 in `validate` and `validate_from_afa`; `identity_clusters` takes another as `threshold`
— and reports both the row count and the number of clusters, so the effective size of the
held-out set is visible next to its nominal size.

**What this proves and does not prove.** It measures whether the residue rule, read at the
calibrated column, recovers the cofactor that experiment reports on labelled proteins. It
does not establish the cofactor of any uncharacterised sequence; that needs structure and
ultimately measurement.

---

## 2. The census pipeline

`src/scriptLn.jl` is the top-level driver: one linear pass from an alignment of search hits
to a census table with taxonomy and habitat attached. It reads every large input by name
from the data root (`$LNFINDER_DATA`, see [DATA.md](DATA.md)), where a variable of their
own can move two of them (`LNFINDER_AFA`, `LNFINDER_TAXDUMP`), and refuses, naming the
input, where the layout expects it and which variable moves it, if one is missing. It
writes under `LNFINDER_CENSUS_OUT`, which has no default.

**The starting point is the census hit set.** The census begins from the sequences in
`census/firstTry.fa` (`hits_fasta`), the hits of an `hmmsearch` with the shipped profile run
outside this repository, and from their alignment to the profile (`model_alignment`, which
`hmmalign` regenerates from them; [DATA.md](DATA.md) §5).

The driver's per-residue and per-sequence tables come from `src/ModuleSequenceTools.jl`, a small
library of two functions: `count_residues_per_position`, which tabulates residue counts by
model column across an aligned set, and `export_tree_fasta`, which writes the per-class FASTA
a tree is built from elsewhere. It has no suite of its own; the census driver is its only
caller.

Order of operations:

1. **Calibrate** the active-site column on the alignment (§1.3).
2. **Classify** every aligned hit at that column (§1.1).
3. **Dereplicate** to cluster centroids. Clustering is done outside this repository with
   USEARCH `-cluster_fast` at 90 % identity; the census reads the `-uc` file. The identity is
   `DEDUP_ID` (0.90) in the script, beside `USEARCH_CMD`, which documents the clustering
   command. That command is never run from here, so the two are kept in step by hand; the manifest records both. Optionally, `EXPECTED_CENTROIDS` makes a rebuild
   abort unless the clustering yields exactly the recorded centroid count.
4. **Resolve taxonomy**: each centroid's accession is resolved to a TaxID over UniProt, and
   its lineage is read from an unpacked NCBI taxonomy dump (§3.1).
5. **Attach habitat** from the habitat layer (§3.2). `RUN_HABITAT=false` skips it, for
   reproduction runs that only need the census table.
6. **Write the census and its manifest**, then compute statistics.

Two ordering rules are load-bearing and are the reason this file exists in the shape it
does:

- The full census and its TOML manifest are written **before** anything is filtered. Rows
  without a TaxID, a kingdom or a habitat belong to the census; dropping them before the
  write makes exclusions invisible and unauditable afterwards.
- Nothing mutates the census table after that write. Every analysis subset is built with
  `filter`, so the table on disk and the table in memory stay equal.

The manifest records, for every input, its path and its SHA-256 (or an explicit marker when
a file was absent, so a missing input records itself instead of aborting a run whose census
is already written), plus the run date, the calibrated active-site column (§1.3), the
dereplication settings, the UniProt release the TaxIDs were resolved against, the taxonomy
dump's own provenance (`taxdump_provenance` reads the dump's release files and digests) and
the versions of the direct package dependencies.
A census is never overwritten: the driver refuses to start if the output file already
exists.

---

## 3. Taxonomy and habitat

### 3.1 Taxonomy (`src/ModuleTaxonomyTools.jl`)

The NCBI taxonomy dump is loaded into parent/rank/name maps plus the merged and deleted id
lists. Two facts about the dump are handled explicitly because both silently corrupt a
kingdom column otherwise:

- NCBI reorganised the top of the taxonomy in March 2025: `superkingdom` became `domain`, a
  real `kingdom` rank was inserted below it, and viruses were moved under an acellular
  root. Code keyed on `superkingdom` alone returns an empty kingdom for every row of a
  later dump; code keyed on rank `kingdom` picks up the newly inserted rank instead of
  Bacteria or Archaea. Both vintages are handled, and NCBI's own `kingdom` rank is
  deliberately not used.
- A TaxID that has been merged is followed to its successor, and a deleted one is reported
  as deleted rather than dropped.

`append_taxonomy!` refuses to return a census whose kingdom column is mostly empty, so a
dump-vintage mismatch stops the run instead of travelling downstream as missing data.

The TaxIDs themselves come from UniProt. Accessions are resolved over three stateless
UniProt requests (`src/ModuleAccessionMapper.jl`, `resolve_taxids`), issued one at a time
and paced:

1. a UniProtKB search, 100 accessions per request (`SEARCH_BLOCK`, the service's own cap on
   one query), gives live entries their organism directly and hands back a UniParc id for
   every deleted one;
2. UniParc answers a deleted entry whose record names a single organism, and only where the
   accession reappears on that record;
3. a record shared between organisms is resolved through the accession's own
   cross-reference on it.

UniParc's `commonTaxonId` is never used as a fallback. This route uses neither UniProt's
id-mapping jobs nor its stream endpoints: the service's robots.txt names both, and a job
can answer short without reporting a failure. Two other steps do use a stream route, and
what robots.txt says about the paths they ask for is below. Every search block is checked for the exact set
of accessions asked for: a block that answers short is asked once more for the same
accessions, and one
that is still short, or that names an accession nobody asked for, stops the run. Every
accession is accounted for in a ledger that is checked against the output table, and the
run stops when UniParc resolves fewer than 98% of the deleted entries (`min_rescue_frac`).
Each row records the UniProt release the service reported.

An accession-keyed cache (`accession_taxid_cache.csv` in the run's output directory) holds
only established resolutions and is written after every check has passed, so a re-run into
the same directory costs no requests for the accessions already resolved. The third
request's answers are checkpointed per UniParc record and accession in
`uniparc_xref_cache.csv` beside it, and only an answer shown to be complete is kept, so a
short answer is asked again on the next run rather than remembered as an absence. The
accessions still unresolved are written to `unresolved_accessions.csv` before the 98% floor
is checked, so a run that stops there leaves them behind to investigate.

The request policy differs by status because the statuses mean different things: a 5xx is
retried, up to `HTTP_TRIES_5XX` attempts in all, with a short growing wait; a 429 is waited
out (the server's `Retry-After` when it sends one, never below `PACE_S`, else
`BACKOFF_429_S`) at most three times per request, and a fourth 429 on one request stops the
run, because a 429 is a request to go away rather than to try harder; any other status, a
redirect included, stops the run at once. Requests are spaced by `PACE_S`, and by
`PASS3_PACE_S` in the third request.

**Where a stream route is used, and what robots.txt says.** Two pipeline steps use one. Step 0
(`scripts/prokaryote_background.jl`) asks `https://rest.uniprot.org/proteomes/stream?query=`
for the reference-proteome taxa of a domain, and step 2 (`scripts/background_pqq_scan.jl`) asks
that same route for a sampled taxon's reference proteome and
`https://rest.uniprot.org/uniprotkb/stream?query=` for that proteome's sequences, one taxon at
a time, writing each answer to a file it reuses rather than asking again. Read on 2026-09-23,
`https://rest.uniprot.org/robots.txt` disallows `/docs/`, `/idmapping/`, and a
`/<dataset>/stream/` path for each dataset it serves — `/proteomes/stream/` and
`/uniprotkb/stream/` among them. Those patterns end in a slash and the two routes above carry a
query string in its place, so the requests fall outside the literal prefixes while asking of
the same endpoints; the restraint the two steps apply — one request at a time, each answer
kept on disk — is their own. No route anywhere in this repository submits an id-mapping job.

### 3.2 Habitat (`src/ModuleEnvironmentTools.jl`)

Habitat comes from OmniCrobe, which text-mines habitat relations from published work and
expresses each habitat as an OntoBiotope (OBT) term. Two consequences follow, and both
shape the method:

- **A raw habitat tally measures research effort, not ecology.** A model organism
  accumulates relations by orders of magnitude faster than an unstudied one. Every quantity
  downstream is therefore a comparison against a denominator computed from the same source
  (§4), never a raw count.
- **A term is a node in an ontology, not a category.** Each term is climbed to its OBT root
  path, and the path — not the term — decides the category. A run writes the climb to a
  cache under the data root (`habitat/obt_path_cache.csv`, see [DATA.md](DATA.md)), and every
  later run on that root reads it, so the ontology is walked once. The data set does not carry
  that cache: on a root without it, the climb goes back to the network term by term.

Mapping rules:

| step | rule |
|---|---|
| habitat **type** | the term's OBT path is tested against a table of branch ids, in order, first match wins. The specific branches (aquatic, terrestrial, geological, host-associated, anthropogenic, …) are tested before the general ones, because several general branches contain the specific ones as children |
| habitat **chemistry** | a second, non-exclusive set of condition ids (e.g. thermal, saline). A taxon may carry several, or none; these are flags, not a partition |
| **isolation sources** | terms that record where a strain was *obtained* rather than where it lives. They are commoner in a background drawn from culture collections, so leaving them in manufactures a fake depletion. Under the primary rule they are dropped term by term, which keeps a taxon that also has genuine habitat terms; the legacy rule drops the whole taxon |
| **generic terms** | terms too general to carry information are identified and excluded from the winner decision |
| **distinct-term cap** | at most `MAX_DISTINCT` (50) distinct terms per taxon enter the tally, so one exhaustively annotated organism cannot dominate |
| **ties** | when two categories tie on the deciding tally, the taxon is split between them with weight 1/k rather than assigned to one. Assigning the winner by dictionary order instead would make the answer depend on hash iteration order |

The branch table is built from the OBT paths actually returned for the terms in the data,
and each entry carries the ontology's own label for the id it names, so a reader can check
a mapping against the source without a network call. Several entries carry a comment saying
why a branch is split by child rather than labelled at the parent; those comments are part
of the method and should be read with this table.

### 3.3 The habitat call table (`scripts/habitat_reaggregate.jl`, pipeline step 1)

A habitat number taken from a live re-fetch makes a call depend on what the service answers
that afternoon, on whether the fetch succeeds at all (a failed fetch enters the analysis as
the ecological statement "this organism has no known habitat"), and on dictionary iteration
order where two categories tie. Two runs of the same pipeline are then not comparable.

The remedy is to freeze the **inputs**: a run stores the raw response bodies per taxon under
`habitat/omnicrobe_cache/relations/` in the data root and the ontology climb in the OBT cache,
and every rule variant is computed from those same bytes in one pass. For whoever holds that
cache, a re-run with the network unplugged reproduces the table exactly. The bodies and the OBT cache
are third-party response content and are not redistributed with the data set, which carries
this step's output — the call table the later steps read — together with the cache's
manifest. A run that starts without the cache therefore fetches, and what it fetches is a new
measurement rather than these bytes again. The cache carries a manifest stamping each body
with where it came from and when, and with the service and ontology releases, because a cache
with no release stamp cannot be told from one taken years later. The two releases are
constants in the script (`OMNICROBE_RELEASE`, `OBT_RELEASE`), not values read from the
service's responses.

A taxon missing from the cache is fetched live, one taxon at a time. A transient failure —
status 429, 500, 502, 503 or 504, or a connection that drops before any status comes back —
is asked again up to three times, waiting 2, 4 and 8 seconds, and a fetch that succeeds is
followed by a 0.2-second pause. Any other status is not asked again. A fetch that fails
either way is recorded as a failed fetch — never as a taxon with no habitat — and a method
gate below refuses to publish over it.

The output is one row per (configuration, taxon), carrying the winning category, the raw
distinct-term and mention tallies per type, the 1/k tie weights, the non-exclusive chemistry
flags and the isolation/generic bookkeeping — everything a downstream rule variant needs, so
no variant costs a refetch.

**Configurations.** `primary` is the rule this release computes with. `legacy` reproduces the
rule the calls carried by the census and background-pool tables were computed under, and
exists only so those calls can be re-derived and diffed field by field. The remaining
configurations each move exactly one axis
off `primary`; they are the levels the sensitivity grid (§7) sweeps.

**Gates before publication.** Three kinds of gate stand between a run and its outputs.

- **Method gates**, whose expected value is a requirement of the method: no ontology path
  may be missing after the fill, and no relation fetch may have failed. They are written in
  the script and are checked on every run.
- **Measured gates**, whose expected value is a *measurement of the dataset*: coverage and
  isolation-term recovery under `primary`. Such a value is never typed into the code. It
  arrives in a `--gate-targets` file, one row per gate, each with the note that says where
  its number came from, and a run given no such file is refused unless `--no-gates` says
  so. Under `--no-gates` each measured gate that is not also a reproduction gate prints its
  measurement as NOT CHECKED.
- **Reproduction gates**, four of them: the `legacy` configuration re-derives the calls the
  census and the background pool carry in their own habitat columns — habitat type and
  chemistry, on the census's lanthanide arm and on the background — and compares them taxon
  by taxon. Reproducing those calls is a question about the census they were made on. Asked
  of another census it is a different test, one an honest run could fail and whose pass
  would read as evidence that two censuses agree. So the four are evaluated only when the
  `--gate-targets` file declares, by its SHA-256 digest, the census they were written for
  (one row `census,,<sha256>,,<note>`: no gate named, no tolerance, a note required; the run
  prints the digest of the census it read as `Census sha256:`),
  and the census the run reads has that digest. Otherwise — `--no-gates`, a targets file
  that declares no census, or a census with another digest — each of the four is declined:
  its measurement is printed as NOT APPLICABLE with the reason, and it still counts as
  reached for the check that every declared gate was evaluated, so declining a gate cannot
  hide a deleted one.

A run's closing line says which gates it did not make. Step 1 on the deposited census runs
with `--no-gates` ([DATA.md](DATA.md) §4): its measured and reproduction gates are then not
made, and the method gates are the ones it is held to. `--census <path>` runs the step on
another census and the report names the file it read. Taxa with no phylum are counted in
the report for each analysed arm and not gated here: the rule they bear on is a rule about
the contrast, and step 3 enforces it where the arms are built (§6.1).

Outputs are written into a `_staging` directory beside their destinations and moved into
place only if every gate passed, so a failed run leaves its output inspectable in staging
and leaves the destination holding the last good run.

---

## 4. The comparison arms

### 4.1 The denominator (`scripts/prokaryote_background.jl`, step 0)

"Where do the carriers live?" cannot be answered by counting where they live; it needs a
background of comparable organisms put through the identical pipeline. The pool is UniProt
**reference proteomes** for Bacteria and Archaea, chosen because the census hits are
themselves proteins from sequenced organisms: the pool shares the census's strongest bias —
only cultured, sequenced, named organisms are represented — so that bias largely cancels in
the contrast. A pool drawn from all of taxonomy would not share it and would manufacture
enrichment. Reference proteomes are deduplicated by UniProt to roughly one per species
group, so a heavily sequenced species counts once, which is the sampling-frequency
correction wanted.

Each taxon's lineage is fetched alongside its identifier so that a phylum-matched background
can be built downstream (§6). Above a pool-size cap the taxa are subsampled with a recorded
seed; the cap, the subsample size and the seed are named constants in the script
(`MAX_TAXA`, `SUBSAMPLE`, `SEED`), so the same pool can be rebuilt. The run is resumable:
taxa already present in the output are skipped.

**Step 3 records that provenance on every run**, in its summary table under
`background_pool_*` and `background_frame_source_*`: the seed, read from the generator's own
source and cited there by file and line rather than typed into step 3; the subsample size
and the cap, read from the same source; the generator's digest; the file name, digest and
taxon count of the pool and of the frame table; and whether the pool on disk holds exactly
the subsample size, which is what says the subsample was taken. A value that cannot be read
is recorded as absent, with a state that says why (`background_pool_provenance_state`): a
default there would be step 3 inventing the seed of a subsample it did not take.

### 4.2 The frame, and the carriers hiding in it (`scripts/background_pqq_scan.jl`, step 2)

The background is the **sampling frame minus the census**: reference-proteome prokaryotes
that do not carry the enzyme family. Built by taxonomy alone, that subtraction leaves a
known contaminant: a taxon can be outside the census either because it has no such enzyme
or because the census search did not see the one it has (a proteome released later, an
entry whose TaxID moved, a sequence below the census threshold). Those are carriers sitting
in the comparator arm, and they bias every contrast toward the null — invisibly, because
nothing else in the pipeline asks whether a background taxon has the gene. This step asks,
with the census profile itself.

**Targeted, not full-proteome.** Downloading every background proteome is several gigabytes
for a question about one domain family, so the search space is narrowed first, to the
UniProtKB entries InterPro places in either of two entries:

| InterPro entry | what it is |
|---|---|
| IPR011047 | quinoprotein alcohol dehydrogenase-like superfamily (homologous superfamily) |
| IPR017512 | PQQ-dependent dehydrogenase, methanol/ethanol family (family) |

The union is taken, by accession, locally. The superfamily alone is not enough: InterPro
annotates some family members without the superfamily entry, and under a superfamily-only
prefilter such a taxon never reaches the search at all and comes out with zero hits however
permissive the cutoff. Whether the union's entries then reach the cutoff is measured, not
assumed. Each download keeps its own release record and the run stops if the two halves came
from different UniProtKB releases: a union of two snapshots is not a snapshot.

**The search.** `hmmsearch` with the shipped profile, at a bit-score cutoff, writing a
`--tblout` table that the script parses. The launcher is chosen from the platform: on
Windows the invocation goes through WSL and the paths are translated, because HMMER has no
native Windows build; on Linux and macOS `hmmsearch` is called directly. The option line is
the same either way. See the platform note in `README.md`.

**The cutoff is calibrated, not chosen.** A hand-picked bit score would be unreproducible,
so the cutoff is set empirically from the distribution of the census accessions' own scores against
the same profile, by a named rule (`calibration_threshold(scores; rule = :p1)`, the first
percentile). The calibrated value is written to a record file beside the scan output and
reported by the run; it is not a constant in this repository. A caller who wants a different
cutoff passes `--threshold` and the run records that it did so.

**The prefilter is an assumption, so it is measured.** A sample of frame taxa has its full
reference proteome downloaded and searched at the same cutoff; of the taxa the full search
flags, the share the prefilter also flags is the measured recall. The sample is the first
*n* draws of one seeded permutation of the frame, so a larger sample keeps every taxon a
smaller one drew — and the run checks that against the sample column of an earlier scan
table rather than taking it on faith. A sample in which the full search flags nothing is
recorded as a miss, not a pass: it measured no recall at all. A rate is not a bound, and the
run records both the observed rate and a Clopper–Pearson lower bound on it, because the
question "does the prefilter find what the full search finds" and the question "does the
remaining background hold no missed carrier" are different questions.

**Acceptance.** The recovery and recall bars the step is held to (`MIN_RECOVERY`,
`MIN_RECALL`) and the confidence level of the bound (`RECALL_ALPHA`) are fixed *before* a
run, as the analysis plan's own parameters, and are constants in the script. They are rules,
not measurements: they say what would count as acceptable, and the run reports what was
observed against them. A run that misses promotes nothing.

**The positive controls, and what a PASS is evidence of.** Acceptance also requires five named
frame taxa to come out flagged, one row each against a bar of one (`POSITIVE_CONTROLS`).
**They were not drawn from the literature, and this release claims no published cofactor
measurement for any of them.** Each was selected from public UniProt annotation: it is a frame
taxon — in the reference-proteome pool, absent from the census — that has a UniProtKB entry of
at least `MIN_FAMILY_LEN` residues annotated to InterPro IPR017512, the methanol/ethanol
family, under that exact TaxID. The five entries are `A0ABM9ICZ5` (TaxID 591154),
`A0ABX7PVQ1` (511746), `A0ABU0I1Z8` (418708), `A0ABU0LS59` (217069) and `A0ABQ2ZAD6`
(980347), named in the script beside the taxa they select, so the choice can be reproduced
from public records alone.

Two things follow, and both decide how a PASS should be read. The list was fixed before this
search was first run, so it cannot have been read off the scan's own output — and it is not a
set that output would produce: these five neither top nor order the scan's hit counts or bit
scores. And because the criterion is annotation rather than measurement, the check asks whether
this prefilter and this profile **recover what public annotation already claims**. That is a
wiring check on the search, not a measurement of what the search discovers, and not evidence
about these organisms' cofactors: all five of those annotations are themselves inferred from
homology rather than measured. One of the five species is the subject of a published
lanthanide-in-the-active-site determination, but for a different strain, and a different TaxID,
than the one used here — see the note above the constant, and
`dat/control_set/known_cofactor_mdhs.tsv`, which carries that evidence against the strain's own
accession. The note above the constant also states, taxon by taxon and with DOIs, what the
published literature does and does not establish for each of the five, including, for one,
published evidence against growth on methanol from an assay that could not have detected a
lanthanide-dependent enzyme. The acceptance row of a control present in the frame says what
its hit count measures: sequence presence only, not that the organism is methylotrophic.

---

## 5. Where the numbers come from, and where they do not

Three conventions run through every analysis script and are worth stating once, because they
are the reason the code looks the way it does.

**No data literal in a script.** Counts, shares, coverages, releases and test statistics are
computed or read at run time, written to a summary table, and read back out of that table
for the closing report. A number typed into a comment goes stale the moment an input changes
and then quietly contradicts the output beside it.

**A bound fixed before the run is not a measurement, and arrives differently.** Bounds the
analysis plan set before the run are supplied on the command line, never defaulted, and are
reported beside the quantity they bound together with the comparison and whether it was met.
Omitting one omits its row; nothing passes silently because a bound was not supplied.

**Estimators, rules and parsers live in one place.** The estimators are in
`src/ModuleHabitatStats.jl`, the outcome predicates in `src/ModuleOutcomeRules.jl`, and each
persisted record format is read by its own producer's reader. A second implementation of any
of them is the one that drifts, and a grid computed by drifted code measures the drift rather
than the rule.

### 5.1 "The registered specification" and "the analysis plan"

This document and the scripts use both phrases often, and a reader should know what they point
at. **They name an external record, and this release does not carry it.** There is no plan
document in this repository to go looking for, and its absence is a fact about the release
rather than an omission to be filed as a bug.

**What "registered" does and does not claim.** The record is a dated, checksummed set of rules
that map each outcome to a colour, to a caption phrase and to a sentence of the accompanying
text. Its timestamp precedes every output of the re-analysis these steps perform, so the
numbers reported here did not choose their own wording. It was **not** written before the data
were first looked at, and it does not claim to have been: its own disclosure records that an
earlier run's estimates for the categories it names, and that run's sensitivity findings, were
already known when it was written. So it is registered against the re-analysis, not against the
original analysis — and the categories it speaks for were **selected from those earlier
estimates** rather than fixed in advance, which is why they are a per-run input here and not a
constant (`docs/DATA.md` §4). Wherever this document or a script says "the registered
specification", "the registered primary" or "the analysis plan", that is what is being pointed
at, and nothing stronger.

What the phrases surface as, in the code that is here:

- **the primary configuration.** Each script has one default option set, which it prints at the
  start of a run and records in its own output. "Off the registered specification" means "from
  that configuration, one choice changed".
- **the criteria supplied per run** — `--gate-targets` (step 1), `--anchor-targets` (step 4),
  `--anchor-reference` (steps 3 and 4), `--categories` (steps 3 and 5), `--plan-*` (step 3). The
  two targets files and the category selection carry a note on every row saying where its
  number or selection came from. The anchor targets also name the files they were
  transcribed from, with their digests, which step 4 checks (§7); a step-1 targets file
  names, by digest, the census its reproduction gates were written for (§3.3); and the
  steps that read `--categories` record that file's digest. They are deliberately not
  constants in the code, for the reason §5 gives. The anchor targets for the deposited data
  ship in this repository, as `dat/criteria/anchor_targets.csv` with their two source files
  (`impl_4.csv`, `impl_w1.csv`) beside it.
- **the criteria described in the script itself.** Where a script says "see THE RECORDED
  CRITERIA below" (step 3), that section of the same file states which bounds the plan fixed
  and how a run reports them; their values are the `--plan-*` bounds above, which this
  release does not carry (`docs/DATA.md` §4).

---

## 6. Enrichment estimators (`src/ModuleHabitatStats.jl`, `scripts/habitat_enrichment.jl`, step 3)

### 6.1 Two estimators, always reported together

| estimator | definition |
|---|---|
| **crude** | (a/n1) / (b/n2), as a log2 risk ratio, with a Katz Wald interval and a Fisher exact test. When a or b is zero, and only then, a Haldane correction (+0.5 on the hit cells, +1 on the denominators) is applied — to this estimator only. Fisher's test is taken on the hit totals rounded to counts; the estimate and its interval use them unrounded |
| **standardised** | the background's share of the category reweighted to the foreground's phylum distribution: p0_std = Σ_k w_k (c_k / n0_k) with w_k = n1_k / n1. No correction is applied: with no foreground hit, or a standardised background share of zero, the estimate is undecidable (NaN) rather than corrected |

The gap between them is the phylogenetic contribution. Habitat annotation tracks phylum, so
a crude ratio can report enrichment that is nothing but the phylogenetic range of the
protein family; reweighting makes that contribution legible instead of leaving it to be
argued about. Both intervals are taken 1.96 standard errors either side of the estimate,
the 95% level.

Further rules of the estimator layer:

- **Fractional hits.** A hit vector may be boolean (a unit is in the category or not) or
  real (a unit split across tied categories, §3.2). Both are accepted, and a boolean vector
  and its 0/1 real equivalent give bitwise identical results. The delta-method variance
  treats the stratum weights as fixed, the units as independent and both totals as
  binomial; for a fractional vector each total is a sum of bounded values rather than a
  count of Bernoulli draws, whose variance is no larger, so under that model the formula is
  an **upper bound** and the resulting intervals are conservative. The form is left that
  way deliberately: adjusting it would need a model for how the weights were split, and an
  interval that is too wide is the safer error. The bound concerns the form of the
  variance, not how well its plug-in value is estimated, and it does **not** hold for a
  stratum whose background cell is empty or thin (next point).
- **The unit must be the same on both sides.** Counting sequences counts how often an
  organism was sequenced, not how many organisms carry the enzyme, so the primary unit is
  the unique taxon. Per-sequence and per-genus series are computed as secondary so that the
  size of that bias stays visible, and the genus series collapses **both** arms to genus so
  that the two counts are in the same currency.
- **Thin cells are disclosed, not hidden.** A row whose combined hit total a + b is below
  `MIN_CELL` (5) is reported but never tested. A phylum thinner than `MIN_FG_K` (10) in the
  foreground or `MIN_BG_K` (20) in the background is relabelled into a single
  "Other phyla" stratum in **both** arms before any counting, so one thin stratum cannot
  carry a large weight against a near-empty background cell; a stratum with foreground
  units but no background units is dropped, and the number of strata that actually
  contributed is reported. That collapse bounds the weight, not the emptiness: a
  contributing stratum whose background cell is empty adds nothing to the variance, so it
  is the stratum that knows least about the background and contributes no uncertainty at
  all, and a thin cell is better only in degree. Row totals cannot reveal that, so every
  row also carries how many contributing strata have an empty background cell and the
  smallest background cell among them — unrounded, because a cell built from fractional
  weights can hold a fraction of a unit.
- **A taxon with no phylum is in no arm.** Phylum is the stratum, so a taxon without one
  cannot enter the standardised contrast, and bucketing it under "Other phyla" would state a
  lineage where there is none. NCBI's catch-all placeholders (`bacterium`, `uncultured
  bacterium`, `uncultured archaeon`) are such taxa, and their habitat relations aggregate
  unrelated deposits rather than describe one organism. A taxonomy that cannot resolve those
  taxa drops them by accident and one that resolves every TaxID admits them, so the
  exclusion is a rule rather than a side effect of how complete the taxonomy is. It is
  applied where the arms are assembled, to every arm step 3 builds — the census foreground,
  the taxa a move-in level adds to it, the comparator arms and the background — and the run
  stops, with a count per arm, if an assembled arm still holds such a taxon. A phylum-less
  frame taxon stays in the frame the scan's coverage is compared against; it leaves only the
  arms.
- **Multiple testing.** Benjamini–Hochberg, within each contrast × unit panel: one
  contrast, in one unit, over its habitat categories or over its chemistry flags — never
  across panels or across a run. It is applied to the crude estimator's Fisher exact
  p-values of the rows that are tested (a + b ≥ `MIN_CELL`); an untested row carries no
  q-value. No colour and no outcome code reads a q-value.

**The anchor, and what `anchor_status` says.** Asking whether a category is new needs an
earlier result to be new against, so step 3 computes one in the same run: the same builder on
the earlier settings — the `legacy` configuration, its own chemistry column, no background
exclusion and nobody moved in. Whether the partition may use that anchor's category set is
decided in a fixed order, each rung withholding it for a reason that makes the later rungs
meaningless: `undecidable_config_absent`, the call table carries no such configuration;
`undecidable_frame_mismatch`, the two specifications sampled different frames — the frame
each drew from, not the arm each kept, so a move-in level is not a mismatch for moving taxa;
`undecidable_anchor_intervals`, a row of the anchor panel has no interval, so its category set
is not the earlier result; `computed_unverified`, no `--anchor-reference` was supplied, so the
anchor stands and nothing checked it; `undecidable_reference_missing` and
`undecidable_reference_unreadable`, a reference was named and is absent, or is not a rows
table; and then, on the two category sets, `reproduces_reference` or
`differs_from_reference`. Only `computed_unverified` and `reproduces_reference` let the
partition use the anchor at all. `anchor_verified` in `habitat/habitat_arm_summary.csv` is the
narrower statement that `anchor_status` is `reproduces_reference`, so a run given no reference
reports `computed_unverified` with `anchor_verified` false, and that pair reads as unchecked,
not as contradicted. A supplied reference is also compared row by row — rows matched, rows
agreeing, rows disagreeing and the largest cell difference — because a category set can
agree while the numbers behind it move; that comparison is a disclosure and decides nothing.

**The two move-in pool rules, and which one a run reports.** One background-exclusion level
adds instead of only removing: `hmmsearch_movein` takes the carriers out of the background
frame and gives the lanthanide-type ones to the foreground. It is a sensitivity level, never
the headline. Which taxa it moves is a rule, and there are two of them.
`census_lanthanide_hits`, the rule the recorded level names, is frame taxa carrying at least
one sequence the census itself classified lanthanide-type, read from the census's
pre-deduplication hit table, which is not a repository input and so has to be supplied
(`--movein-pool`). `scan_lanthanide_hits` is
frame taxa whose scan row carries a lanthanide-classed hit — this run's own model call on
the family members the background scan fetched, a different measurement, and the rule a run
falls back to when no hit table is supplied. A run records which rule it applied
(`movein_pool_rule`), whether that rule was supplied or fallen back to
(`movein_pool_rule_source`), whether it is therefore the registered one or descriptive
(`movein_pool_rule_role`), the rule's own note, both rule names, and the pool's source, digest,
row counts and size; whenever it can compute them it also records the other rule's pool size
and the overlap between the two, since how far apart two rules sit is the only comparison one
run can make. A fallback warns as well as records. Only this level emits those keys, because a
run of any other mode moves nobody and the keys would describe something it did not do, and a
pool that cannot be answered at all stops the run rather than turning it into the mode it
would fall back on. A moved taxon joins the foreground and leaves the frame in the same
expression, so no taxon is counted on both sides.

### 6.2 Outcome rules (`src/ModuleOutcomeRules.jl`)

The predicates that turn an estimate into a decision are a module of their own, with unit
tests, rather than expressions inline next to the code that prints their answer. Their
design principles:

- **Three-valued interval tests.** Excluding zero, including zero and *undecidable* are three
  different states. An estimate that was never computed is undecidable, not "includes zero";
  folding the two together reports a decision about a number that does not exist.
- **Anchor-relative, never constant-relative.** Whether a category is new can only be asked
  against the earlier result it would be new relative to. The earlier (anchor) set is an
  argument; there is no hard-coded list to fall back on.
- **A tier that needs a grid may not be claimed without one.** The strong robustness tiers
  are defined by what the one-at-a-time grid does, so they cannot be read off a primary row.
  Without a grid the answer is undecidable, except for the one tier that needs no grid.
- **Neutral codes, never sentences.** Every function returns symbols and label lists. Wording
  belongs to whatever writes prose, not to this repository.
- **The module's own constants are rule parameters, not measurements.** The sharing bands,
  the confidence level and the thresholds that grade a ratio are the recorded rule's own
  parameters; they are declared as such in the module and are not quantities measured from
  any dataset.

---

## 7. Sensitivity grid and the robustness flag (`scripts/habitat_sensitivity.jl`, step 4)

Every estimate rests on rule choices a reader cannot see in an interval: which habitat wins
when a taxon's terms tie, whether a chemistry flag needs a share of the evidence or a single
mention, which text a keyword fallback reads, which background counts as carrying no such
enzyme. A result that holds only under its own tie rule is not a result. This step varies
each choice **one at a time** off the registered specification and reports what survives.

The grid is not computed by new code. Every level reaches the estimators through the same
aggregation rules, the same arm builder and the same predicates the primary run used. Two
seams are not arguments of the estimator entry point and are crossed explicitly, each with
a check:

- **Type rule.** The arm builder reads the call table's weight column verbatim. A type
  variant is computed by rewriting that column from the same table's own term and mention
  tallies through the aggregation rule, then calling the unchanged entry point. Before it is
  used, recomputing the *registered* rule from those tallies must reproduce the committed
  column exactly — otherwise the rewrite is not reading the tallies the table was built from.
- **Chemistry rule.** The flag route cannot express a fractional weight, so a chemistry
  variant is carried by the weight column instead, and under the registered rule the two
  routes must agree cell for cell.

Weights are not rounded on the way in: the published call table rounds its weight column for
reading, and a three-decimal round trip changes a hit total, so a rewritten column is carried
at full precision and only the *check* against the published column uses the published
format.

**What decides the flag.** The flag reads the primary interval and the one-at-a-time variants
on the two units whose denominator is a taxon. Unit checks — the symmetric genus unit, the
per-sequence unit, and the subset of taxa whose call rests on several distinct terms — are
reported for the same grid and explicitly marked as *not* feeding the flag, because a unit
check changes what is counted rather than which rule counted it. Anchors and disclosed
background levels sit outside the variant set for the same reason.

**Anchors, and what is checked while they are not evaluated.** Two anchors outside the grid
reproduce earlier runs of this analysis: the all-legacy specification
(`anchor:all_legacy`), compared row by row against an earlier run's rows table
(`--anchor-reference`), and four levels that reproduce an earlier sensitivity run's
chemistry rules on the legacy arms. Both reproduce rows of a census other than the one
a run reads. Asked of a run reading a different census they are a different test, and a
pass would read as evidence that the two censuses give the same rows. `ANCHOR_LEVELS_RETIRED`
in the script therefore keeps them out of the grid: their definitions stay in the code, and
flipping that one constant builds them again for a run on the census they describe. While it
holds, the row-by-row comparison covers no rows and the run records `anchor_levels_retired`
and an `anchor_reference_status` of not applicable. The targets file declares
`reference_rows,0` to match; the run reads that declaration but, with the levels retired,
does not compare it. `--anchor-reference` is then not required, and a reference given anyway
is recorded and not opened (provenance key `anchor_reference_named_not_opened`, its file
name only; no digest).

The wiring check rests instead on the **anchor targets** (`--anchor-targets`; the file this
repository ships for the deposited data is `dat/criteria/anchor_targets.csv`). A target is
one quantity the grid writes for one row — `a`, `b`, `n1`, `n2`, `crude`, `std`, `lo` or
`hi` — on the registered specification or on one of its one-at-a-time levels, with a
tolerance, transcribed from the output of an independent implementation of the
specification rather than from this script's output (`README.md` describes those
implementations and what the specification left to them). Before anything is published,
step 4 requires that:

- every target is reproduced within its tolerance;
- every level named in a published breaker list is pinned by a target, because a wiring
  error confined to one level's arms would rewrite that list while every count in the run
  agreed with itself;
- the file's own declared counts — of targets, and of levels, whose names must equal the
  set the script enumerates — match what it carries;
- each file the targets name as their source, where it sits beside the targets file, has
  the digest the targets declare. The repository ships both sources there. A source that is
  absent is reported rather than refused; one that is present with other bytes stops the
  run, because the targets were then transcribed from a version of it this run cannot see.

Running without the targets has to be typed (`--no-anchors`, which declines the targets
and the reference together), and is then recorded in the output as a decision
(`anchor_checks_declined`) rather than left as a gap a reader of the output cannot see.
Step 5 refuses to draw from a grid whose anchor checks were declined or whose targets were
not all met.

---

## 8. The enrichment panel (`scripts/habitat_figure.jl`, step 5)

The panel is described here as a drawing procedure; what it draws is data the run reads.

Two panels share one row per category. The upper panel shows each arm's share of taxa in the
category, the carrier arm drawn as a solid bar inside the background arm's outline on a
common baseline, so that an overshoot is read as a length against a fixed origin; it runs no
test and therefore carries no test's colour. The lower panel shows log2 enrichment against
the background with three marks per row on three sub-rows: the crude estimate as an open ring,
the phylum-standardised estimate as a filled dot on its interval, and above both a thin tick
for the same standardised estimate computed for the comparator arm. A connector joins ring and
dot, so the contribution of phylum composition is drawn rather than asserted; the adjustment
has no fixed direction, moving some rows toward zero and others away from it, and the
connector is what tells them apart. Chemistry shares the lower panel's axis as a second block
below a rule — the same estimator in the same units — and its flags are not exclusive of each
other, which the block label says.

Rules the drawing obeys:

- **Colour is the robustness flag and nothing else.** A dot is coloured only where the
  sensitivity grid's flag gives the row a sign that the primary specification and every
  one-at-a-time variant reproduce. An interval that excludes zero under the primary
  specification alone does not colour a row. A flag that contradicts its own interval stops
  the run; it is never drawn grey instead.
- **Marks may not overlap.** The tick has a sub-row of its own, because a tick drawn at the
  dot's height disappears under the dot and shows one measurement where the run made two.
  The offsets and marker sizes are named constants whose clearances the test suite checks as
  arithmetic, and the panel height is derived from a row pitch, so a clearance that holds in
  the suite holds on the page.
- **Which categories are drawn** is decided by a minimum weighted count, computed on the
  fractional foreground count and blind to direction, so it cannot select on the outcome.
  Held-back categories are reported with their counts.
- **Nothing is published unless the inputs describe one analysis.** The rows, the arm
  summary, the flags, the grid and the grid's provenance come from separate runs; every
  identity between them is a named gate, all failures are reported together, and one failure
  stops the run before anything is staged. A flag is recomputed from the grid rows it
  summarises rather than trusted. Outputs are staged side by side and moved into place
  together, and only if the staged set is exactly the declared one.
- **Provenance names files, not locations.** The published provenance records each input's
  file name, whether its path was typed or defaulted, and its digest. A value that looks like
  an absolute path stops the run, because a local path adds nothing the digest does not and
  travels with the file.

**What the declared output set is.** Six files, all written into the figure directory
(`--figdir`, else `$LNFINDER_FIGDIR`; the step never creates one):

| file | what it is |
|---|---|
| `habitat_enrichment.pdf`, `.svg`, `.png` | the panel, the same drawing in three formats |
| `habitat_enrichment_source_data.csv` | one row per mark the panel draws or holds back, with the role that says which |
| `habitat_enrichment_caption_values.csv` | the numbers a caption would otherwise be typed from, exported as `key,value` |
| `habitat_figure_provenance.csv` | each input's file name, whether its path was typed or defaulted, and its digest |

Staging is a directory of the step's own (`_staging_habitat_figure`) beside the destinations,
so a failed run of this step never clears what another step's failed run left for inspection
in the same directory. Nothing is moved into place unless the staged set is exactly those six
names — no more, no fewer.

**Why a caption-values file exists at all.** A caption typed by hand goes stale the first time
a rerun changes a count, and nothing catches it. So every number a caption would state is
derived here from the same inputs the panel is drawn from, exported under a fixed key, and the
key set is declared in advance: a key the run cannot fill stops the run, and a key nothing asks
for is refused. The export is a mechanism; the file is an output of a run and is not shipped.

**Two exports are about a named category**, and neither touches what is drawn or coloured: the
interval exported across every specification in the grid, and the comparator tick divided by
its own standardised estimate to give the sharing ratio. **Which** categories those are is not
in this script and is not a default of it. They arrive in the `--categories` criteria file —
the same file step 3 was given, read with step 3's own reader — which the run refuses by name
when it is absent. `caption_keys` takes the selection as a required argument, and
`caption_values`, `derive` and `run_figure` take it with the inputs rather than as a defaulted
keyword, so no code path can supply a category this script chose. Each run records the
selection and that file's digest in its provenance, and the exported key names are derived from
the selected names, so a run on other categories exports keys for those categories. Where those
selections came from is stated in `docs/DATA.md` §4: not from this software, and not fixed in
advance of the estimates.

---

## 9. Specificity determinants: sequence and structure

This is the library half of the repository: it looks for what, beyond the active-site
residue, distinguishes the two classes. It ships as modules with unit tests; the
command-line drivers that ran it are not part of this release (see `README.md`).

### 9.1 Non-circular labels (`src/build_labels.jl`)

A scan for class-discriminating positions must not be told the answer by the feature it is
scanning for. Labels therefore never come from the active-site residue. They come from the
validated control set, from a structure-derived pocket table, or from both, merged in a
declared precedence order (structure first by default). Disagreements between sources are
**reported, not silently resolved**. `reconcile_with_alignment` then reports how many
alignment sequences actually receive a label and which labels map to no sequence.

Non-circularity has a second edge: when the control set carries evidence tiers, only the
tiers in `DEFAULT_LABEL_TIERS` (`ref`, `A`, `B`) become labels. Rows whose cofactor was
assigned by similarity to another protein, and rows whose call is contested, would otherwise
feed the scan the very resemblance it is meant to discover.

### 9.2 The determinant scan (`src/sdp_scan.jl`)

For each model column, the reweighted residue distribution among calcium-class sequences is
compared with the distribution among lanthanide-class sequences, and the column is scored on
how strongly it separates them. Two scores are computed: Jensen-Shannon divergence (log base
2, so it lands in [0, 1]) and mutual information in bits between residue identity and class.

Three guards, each against a specific way this analysis lies:

| guard | what it does |
|---|---|
| **circularity** | the active-site column is *masked* — excluded from the ranking — but still scored and reported, as a positive control. A correct run shows the masked column with a high raw score; the interesting output is the high-scoring columns that are not it |
| **phylogeny** | near-identical over-represented sequences inflate apparent conservation, so sequences are reweighted by the standard scheme: weight = 1 / (number of neighbours within `theta` identity). `theta` is a parameter of the scan |
| **variability is not discrimination** | a column can be hypervariable and yet identical across classes, so the null is a **label permutation**: the class labels are shuffled over the same sequences many times and the observed separation is expressed as a z-score and an empirical rate against that null |

The module reports a first-order residual diagnostic beside the ranking, which asks whether a
single-position model already accounts for the class separation or whether coupled positions
must be modelled (§9.4). Its vocabulary is a property of the method; what any dataset makes it
say is a result and is not stated here.

### 9.3 Structure-derived labels (`src/pocket_geometry.jl`, `src/export_pockets.jl`, `src/ingest_pockets.jl`)

The chemistry gives a structural criterion independent of sequence: a lanthanide ion is
larger and takes a higher coordination number than calcium, and the enzymes that bind it
donate more carboxylate oxygens to the metal. A pocket is therefore labelled by counting
coordinating carboxylate oxygens (Asp OD1/OD2, Glu OE1/OE2) and the overall coordination
number. Two regimes, and they are not equivalent:

- **Holo** structures, with the metal present: the label is ground truth read from the metal
  element itself, with no dependence on any sequence feature.
- **Apo** structures and predicted models, with no metal: the label is *inferred* from pocket
  geometry around a reference point, carries lower confidence and is flagged as
  geometry-based.

The coordination threshold is **calibrated** on structures whose metal is known —
`calibrate_threshold` does this, the same philosophy as calibrating the sequence classifier
on known-cofactor proteins, in structural space. The recipe is four things together, it
depends on the origin, so the four travel together — but they do not all come from the same
place:

| element | where it comes from |
|---|---|
| the distance cutoff | an **input**, chosen for the origin: a metal site and a residue centroid are not the same neighbourhood |
| which count is read (`n_carboxylate_O` at a metal origin, `n_acidic_res` at a residue origin) | an **input**, likewise chosen for the origin |
| the threshold on that count | the **output** of `calibrate_threshold`, which `calibration_report` (metal origin) and `apo_calibration_report` (catalytic-residue origin, metal stripped) wrap. Ties break towards the lower threshold, so it is a function of the calibration table alone |
| the width of the no-call band | the **caller's**, read off the class gap the reports print — nothing in this repository computes it |

So each report produces exactly one of the four. It is *handed* the cutoff and the count,
reads the known pockets at them, and echoes them back beside the threshold so a run's record
carries all three together; it reports the class gap and the leave-one-out accuracy the call
rests on, and the no-call width is then the reader's choice against that gap — 0 when the
classes land on adjacent values, 1 when they separate by more. None of the four is a default
or a hard-coded value anywhere in `src/`: every function in `pocket_geometry.jl` that reads a
pocket or calibrates on one takes what it needs of the recipe as a **required** argument, and
so do `label_pockets` and `ingest`. A directory of pockets cannot be labelled at a recipe
nobody named, a calibration table cannot be built at a radius nobody stated, and two runs
calibrated differently cannot be confused for each other. The two pocket suites do pass
values, because a required argument has to be given something: each is that fixture's own
arithmetic on structures the suite builds, and the test says how it follows (§13).

Two caveats are documented in the code and repeated here because they bound what the labels
are worth. First, the geometry signal *is* essentially the extra aspartate — the structural
face of the sequence motif — so these labels are only non-circular for a downstream scan if
the holo ground-truth subset is used and/or the catalytic column is masked in the scan; both
are supported. Second, a predicted model's pocket is apo and may not be in the bound
side-chain rotamer, so a geometry call is a hypothesis, not a measurement.

`export_pockets.jl` reconstructs per-residue pocket PDBs from the alignment and the
structures, using the identical definition of the catalytic residue — the residue aligned to
the calibrated active-site column — and writing each pocket translated so that the centroid
sits at the origin, which is the convention the labeller expects. It refuses a residue whose
identity does not match the alignment rather than exporting a mismatched pocket. Note the
origin caveat it carries: a holo threshold calibrated with the origin at the metal does not
transfer unchanged to an origin at the catalytic residue centroid, so the threshold must be
recalibrated there, at a cutoff and on a count chosen for that origin, before apo labels flow
downstream.

`ingest_pockets.jl` is the batch path: a directory of pocket structures becomes a label
table, which `build_labels` turns into the label dictionary the scan consumes, reconciled
against the alignment.

A minimal fixed-column PDB parser (`ATOM`/`HETATM`) is built into the module rather than
taken from a dependency, so the structural path stays self-contained and offline-testable.

### 9.4 Coupled determinants (`src/contrastive_potts.jl`, `src/plmdca_adapter.jl`)

A per-column scan sees first-order determinants only. It cannot see pairs or networks of
positions that jointly determine specificity while no single position is individually
predictive. A label-aware Potts model recovers those through its pairwise couplings.

The scaffold is dependency-injected: the coupling fit arrives as a function argument, so the
module itself stays installable and offline-testable, and the real backend is a thin separate
adapter. What is implemented:

| function | what it does |
|---|---|
| `class_meff` | per-class effective sequence count. It **gates** the coupling fit: couplings need many effective sequences, and one class is usually the limiting one. Do not fit couplings the data cannot support |
| `field_contrast` | the single-site (field) comparison between classes, in log-odds units — the generative-model view of the per-column signal, computable without a coupling fit |
| `export_class_msas` | per-class alignments, for an external tool |
| `coupling_contrast` | gated on the effective counts: runs the injected fit per class and returns the ranked contrast between the two coupling matrices |

`src/plmdca_adapter.jl` wraps `PlmDCA.jl` into that interface. It reconciles the conventions
explicitly — sequence-major versus position-major matrix orientation, residue encoding with
the gap as the last symbol, weights renormalised to sum to one, and the backend's
APC-corrected pair ranking densified into a symmetric matrix — because a silent convention
mismatch here produces a plausible-looking ranking that means nothing.

---

## 10. Phylogenetic utilities (`src/ModulePhyloTools.jl`)

Newick parsing, pruning, midpoint rooting and patristic distances, for trees built outside
this repository (FastTree output). Two design decisions are method-relevant:

- **Everything is iterative.** The trees involved are large and ladder-like: the tool that
  builds them emits long unbalanced combs, so the nesting depth of a Newick string is a large
  fraction of the leaf count rather than its logarithm. A recursive parser, writer, traversal
  or prune overflows the stack on such a file, and the crash then looks like a corrupt input
  rather than a code bug. Every walk carries its own explicit stack; no function calls itself,
  and the test suite pins that with a deep caterpillar on a small stack.
- **Internal node labels survive round trips.** The tree builder writes its local support
  values in the internal-label slot, and a support value describes the branch below its node,
  so it must survive parse → prune → root → write unchanged. It is kept verbatim as a string.

---

## 11. Parameters, and how each one is set

Only parameters that shape the method are listed. Anything calibrated is named with its
calibration procedure and the file its run records it in, never with a value.

Where a library function carries a default, the default is the **chemistry or the convention**
— a coordination rule, an interval level, a residue set — and never a value measured from the
data. Anything measured is calibrated by the run and written into that run's own output, and
it is the run's record, not the library default, that says what a result was computed at. The
active-site model column has no default anywhere: every entry point that reads a column takes
it as a required argument, so no call can silently inherit one profile's calibration. The
pocket coordination recipe of §9.3 is required the same way, and not only where a pocket is
labelled: every function in `src/pocket_geometry.jl` that reads a pocket or calibrates on one
takes what it needs of the recipe as a required argument, so no table, threshold or report is
produced at a radius or on a count nobody named.

| parameter | where | meaning | how it is set |
|---|---|---|---|
| active-site model column | `calibrate_active_site` | the column at which the cofactor call is read | **calibrated** per run from two reference sequences; written to the run's output |
| collision-column count | `diagnose_calibration` | how many columns satisfy both reference patterns | measured per run; a count above one means the calibration is not unique |
| control-set tiers admitted | `include_tiers`, `DEFAULT_LABEL_TIERS` | which evidence tiers may be scored, and which may become scan labels | fixed rule: homology-only and contested tiers never become labels |
| identity threshold for held-out clustering | `validate`, `identity_clusters` | when two held-out controls count as one | 0.9 pairwise identity over the match columns in the validation entry points; `identity_clusters` takes another as `threshold`. Reported with the result |
| `DEDUP_ID` (0.90) | `src/scriptLn.jl` | the identity the census clustering was run at | fixed to match the recorded external clustering command |
| `EXPECTED_CENTROIDS` | environment | abort unless the clustering yields exactly this many centroids | supplied only when rebuilding a recorded run |
| `SEARCH_BLOCK`, `UPIS_BLOCK` | `src/ModuleAccessionMapper.jl` | accessions per UniProtKB search, UniParc ids per UniParc request | fixed: 100 is the service's own cap on one query; 500 |
| `PACE_S`, `PASS3_PACE_S`, `HTTP_TRIES_5XX`, `BACKOFF_429_S` | `src/ModuleAccessionMapper.jl` | spacing of the TaxID requests and what a 5xx or a 429 costs | fixed restraint of this client, not a limit the service advertises: 1 s between requests (0.3 s in the third request), 5 attempts on a 5xx, waits of 60, 120 and 300 s (or the server's `Retry-After`, never below `PACE_S`) on a 429 — three per request, and a fourth on one request stops the run (§3.1) |
| rescue floor (`min_rescue_frac`) | `resolve_taxids` | the share of deleted entries UniParc must resolve | 0.98, a keyword default; below it the run stops (§3.1) |
| `MAX_DISTINCT` (50) | habitat aggregation | cap on distinct habitat terms per taxon | fixed rule, held at the library default so the legacy configuration reproduces |
| tie weight | habitat aggregation | 1/k split across k tied categories | fixed rule |
| configuration name | step 1 | which rule set the call table row belongs to | one primary, one legacy, and one for each grid level that needs its own call table (`generic_none`, `iso_taxon`, `text_first`, `text_object`) |
| relation-fetch retries, wait, backoff, pause | step 1 | how a transient OmniCrobe failure is handled | fixed in the script: up to 3 further attempts on status 429, 500, 502, 503 or 504 or a dropped connection, waiting 2 s and doubling, and 0.2 s after each fetch that succeeds (§3.3). The library's own default is one request and no retry |
| census the reproduction gates apply to | step 1 | which census the four reproduction gates were written for | declared by SHA-256 digest in the `--gate-targets` file, never defaulted; without a declaration that matches the census read, the four are declined as not applicable (§3.3) |
| pool cap, subsample size, seed | step 0 | the background pool and its reproducible subsample | fixed constants (`MAX_TAXA`, `SUBSAMPLE`, `SEED`); the seed is recorded so the pool can be rebuilt, and step 3 records all three, read from the generator's source, on every run (§4.1) |
| search cutoff (bits) | step 2 | inclusion cutoff of the background search | **calibrated** at run time from the census accessions' own score distribution by a named percentile rule; written to a record file beside the output. `--threshold` overrides and the run records that it did |
| recall sample size, permutation seed | step 2 | the measured-recall sample | the size defaults to `RECALL_N_DEFAULT` and `--recall-n` overrides it; the seed (`RECALL_SEED`) is a fixed constant. The sample is the first *n* draws of one seeded permutation, and an extension is checked against the earlier sample |
| `MIN_RECOVERY`, `MIN_RECALL`, `RECALL_ALPHA` | step 2 | acceptance bars and the confidence level of the recall bound | fixed **before** the run as plan parameters; rules, not measurements |
| `LN_REF`, `CA_REF` | step 2 | the reference pair the step's own alignment is calibrated on (§1.3) | fixed in the script to the pair the census calibrated on, so the step reads its hits at the same column the census read; record order is load-bearing and the script says so |
| `POSITIVE_CONTROLS` | step 2 | five named taxa that must come out flagged, or the run publishes nothing | **selected from public UniProt annotation, not from the literature**: frame taxa with a full-length UniProtKB entry annotated to the methanol/ethanol family under that exact TaxID, the five entries named beside them in the script. Fixed **before** the search they check, so the check asks whether the search recovers what public annotation already claims — a wiring check on the search, not a measurement of it, and **not evidence about any of these organisms' cofactors**. Each is reported with its own hit count against a bar of one. §4.2 |
| `MIN_FAMILY_LEN` | step 2 | below which a family member is a fragment rather than a carrier | fixed rule; the column names carry the floor so a table states its own convention |
| background exclusion mode | step 3 | how carriers are removed from the background arm | one registered primary; every other mode is descriptive and the run records which was in force |
| `MIN_CELL`, `MIN_BG_K`, `MIN_FG_K` | estimators | below which a row's combined hit total a + b is not tested, and below which a phylum stratum is collapsed | fixed reporting rules |
| interval level | estimators | 95% | fixed: 1.96 standard errors, a literal constant in both estimators |
| `ANCHOR_LEVELS_RETIRED` | step 4 | whether the anchor levels, which reproduce rows of another census, are built | fixed in the script to true: the levels stay defined and are not evaluated, and the run records the setting (§7) |
| plan bounds (`--plan-*`) | step 3 | bounds the analysis plan set before this step ran — written from an earlier, approximate version of the same exclusion, so expectations for the run rather than bounds set blind | supplied on the command line, never defaulted; disclosed, never gating. `docs/DATA.md` §4 |
| category selection (`--categories`) | steps 3 and 5 | which category each recorded outcome rule speaks for, which categories the grid counts cover, and which phylum's shares are quoted | the categories are **selected from the data**, not fixed in advance: chosen from this analysis's own earlier estimates, with those estimates in hand. The phylum is **determined, not chosen** — the heaviest standardisation stratum, with the figure step refusing to assert the weight identities unless the run's own heaviest stratum is still that phylum. Supplied per run in a criteria file with no default anywhere, refused by name when absent, and shipped with a required note per row recording where each selection came from. The summary-key families that carry a selected name are derived from it, so a run on other categories emits keys for those categories. `docs/DATA.md` §4 |
| thinnest-stratum list length | step 5 | how many drawn rows the thinnest-background-stratum list names | fixed listing convention, and **outcome-independent by construction rather than by date**: the list is extended through a tie at its last place and the smallest count among the rows it leaves out is exported beside it, so the list plus that bound is a statement about every drawn row. It was written during figure construction, after the results were known |
| coordination recipe: cutoff, feature, threshold, no-call margin | `call_cofactor`, `label_pockets`, `ingest`, and every calibration entry point of `pocket_geometry.jl` | the origin-dependent recipe a pocket's Ca/Ln call is read at (§9.3) | mixed, and §9.3's table says which is which: the cutoff and the feature are **inputs** chosen for the origin, the threshold is **calibrated** on structures whose metal is known, and the no-call margin is the caller's, read off the class gap the calibration reports print — nothing here computes it. No default anywhere: each function requires what it uses, so no call can silently inherit another run's calibration |
| `theta` | `sdp_scan.jl` | identity radius of the sequence-reweighting scheme | parameter of the scan, reported with the run |
| permutation count | `sdp_scan.jl` | size of the label-permutation null | parameter of the scan |
| effective-count gate | `contrastive_potts.jl` | minimum per-class effective sequence count before couplings are fitted | fixed rule; a fit below it is refused, not attempted |

---

## 12. External tools and services

### 12.1 Binaries

| tool | version | how it is used | invoked by this repository? |
|---|---|---|---|
| HMMER `hmmalign` | HMMER 3 | aligns sequences to the profile in `afa` format; this is what makes the model coordinate well defined | **yes**, with a configurable binary name; a missing binary raises |
| HMMER `hmmsearch` | HMMER 3 | the background carrier search, writing `--tblout` | **yes** — directly on Linux and macOS, through WSL on Windows; see the platform note in `README.md` |
| HMMER `hmmsearch`, the census search | HMMER 3 | produced the hit set the census begins from (§2) | **no** — run outside this repository; the census starts from its archived output |
| USEARCH | 11 | clusters the census sequences; the pipeline reads the `-uc` output | **no** — run outside this repository; the command is documented in `src/scriptLn.jl` (`USEARCH_CMD`) and in the run manifest |
| FastTree | — | produces the Newick trees the phylogenetic module parses | **no** — output consumed |
| structure prediction / structural alignment | — | produces the predicted structures the pocket path reads | **no** — output consumed |

The profile itself was built with HMMER 3.0 `hmmbuild` (its header says so) and is read by
HMMER 3.x. No HMMER version is pinned by this repository, and no run manifest records one:
the `--tblout` tables step 2 writes carry the version in their own trailer. Record the version
you used.

### 12.2 Web services

| service | endpoints | used for |
|---|---|---|
| UniProt REST | `rest.uniprot.org` — UniProtKB search and entries, UniParc entries and cross-references, reference-proteome streams | control-set sequence fetch, accession-to-TaxID resolution, the background pool, the targeted prefilter |
| OmniCrobe | `omnicrobe.migale.inrae.fr` — relation search and OBT term lookup | the habitat relations and the ontology climb |
| NCBI E-utilities | `eutils.ncbi.nlm.nih.gov` — `efetch` | the one control row whose sequence is fetched from the NCBI protein record |

Four things follow, and all four are handled in code rather than left to the reader:

- **Requests name the client, and nothing about who runs it.** No call sends an address, an
  account or an API key, and none sends NCBI's `tool` or `email` parameters. Two request
  paths identify the software in their User-Agent header: OmniCrobe relation requests send
  `LnFinder/1.0 (research; serial)` (the default of `fetch_relations`), and the TaxID
  resolution's UniProt requests send `LnFinder/0.1` (`USER_AGENT` in
  `src/ModuleAccessionMapper.jl`). The other calls set no User-Agent of their own. NCBI
  rate-limits unidentified clients, so a caller doing many NCBI requests should space them.
- **Where a failed request is asked again, the rule is written down.** Step 1's relation
  fetches ask a transient failure again up to three times with a growing wait, and pause
  after each fetch that succeeds (§3.3). The TaxID resolution retries a 5xx, waits out a
  429 a bounded number of times, and stops on any other status (§3.1). In both, a status
  that is not transient is an answer about the request rather than about the service's
  load, and is not asked again. Step 2's downloads make up to five attempts in all,
  waiting 5 s longer after each failure.
- **UniProt may gzip a response even when identity encoding is requested.** Every fetch
  therefore checks for the gzip magic bytes and decompresses defensively. Without that guard
  the parse fails silently and the result comes back empty rather than erroring.
- **These services are not versioned like a package.** They change, and a re-run later will
  not reproduce an earlier run byte for byte. That is why the control-set fetch, the TaxID
  resolution and the background scan record the release header UniProt reported, why the
  habitat layer is replayed from a stored cache rather than re-fetched, and why release
  strings are stamped into the cache manifest. A reproduction should replay a cache; a
  fresh fetch is a new measurement, not a reproduction.

---

## 13. Tests

`julia --project=. test/runtests.jl` runs every suite in `test/`; the runner and the
directory are kept as equal sets, because a suite that ships but is not run is a suite nobody
notices breaking. The suite is offline by construction: no network, no external binary, and no
input the repository does not carry. Every fixture is either built by the suite that uses it —
sometimes out of real public identifiers, as the taxonomy dumps are — or a response recorded
from a public service or a local tool and copied verbatim, so that the parsers are held to the shapes those services and tools really return rather than
to invented ones. The recorded fixtures are of three kinds: OntoBiotope ontology paths together
with one OmniCrobe relation record; UniProt response shapes, whose structure is real and whose
every value is invented; and one `hmmsearch --tblout` table for three public Swiss-Prot
entries scored against the profile HMM this repository ships. **No fixture carries a value from
any run of the analysis this software was built for.** Every value any fixture carries is
either public — an ontology path, a Swiss-Prot accession, a bit score recomputable from the
profile HMM that ships here — or invented for the case under test, in which case the expected
result follows arithmetically from the fixture and the suite says how. Where a suite pins a
calibration, it pins one performed on the structures or the table that suite builds, at a
recipe that fixture implies. The shipped data files the suite touches are these. The control
set and the recorded validation predictions
(`dat/control_set/control_validation_predictions.csv`): it only checks that they are present
and spelled exactly as the code or the documents name them. The profile HMM: the same check,
and the classifier suite also reads its header, for the match-state count (`LENG`) the
alignment guard falls back to when a caller names no profile. The anchor targets and the two
source files they name: the same check, and the sources' bytes against the SHA-256 digests
the targets declare. The profile's provenance note, `dat/models/PROVENANCE.md`, it reads as
prose, like `README.md`.

Some shipped pieces have no suite of their own: the census driver (§2), whose
environment-variable gates and never-overwrite-a-census guard are not unit-tested; step 0
(§4.1); `src/ModuleSequenceTools.jl`, which the census driver calls (§2);
`src/ModuleMotifSearch.jl`, which no test calls; and `scripts/verify_data.jl`, which the
suite checks for presence and never runs. `README.md` names the same five, at more length.

One suite is structural rather than numerical: it resolves every `include` target in `src/`,
`scripts/` and `test/` against a directory listing rather than against a file-existence
check, because a case-insensitive filesystem answers "yes" for the wrong spelling and hides a
load failure that is fatal everywhere else. It also refuses two files whose names differ only
by case and any invocation that depends on a version manager. The CI workflow runs the whole
suite on a case-sensitive and a case-insensitive platform for the same reason.
