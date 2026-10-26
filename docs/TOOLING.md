# TOOLING: what the code does

Method documentation for LnFinder: the algorithms, the order they run in, the parameters, and
the external tools and services. No count, effect size, accuracy or fitted threshold from a
run appears below; where a default is calibrated from data, this document says how it is
calibrated and where the calibrated value is recorded.

For what can and cannot be re-run from this repository, see `README.md`.

---

## 0. Software environment

| | |
|---|---|
| Language | Julia, floor declared in `Project.toml` as `[compat] julia = "1.12"` |
| Dependencies | declared in `Project.toml`, resolved exactly in `Manifest.toml` |
| Install | `julia --project=. -e 'using Pkg; Pkg.instantiate()'` |
| Tests | `julia --project=. test/runtests.jl` — offline, no external binary |

`Manifest.toml` ships on purpose: it pins the exact package closure. It was resolved on Julia
1.12, and on an older minor version `instantiate` fails on the standard-library set, not on
anything here. To use a different Julia, delete `Manifest.toml` and re-resolve; the package
versions then float.

Every invocation is written `julia --project=. …`, and the test suite fails the build if a
version-manager-only invocation appears in `src/`, `scripts/`, `test/`, `README.md`,
`docs/TOOLING.md`, `docs/DATA.md`, `dat/models/PROVENANCE.md` or the CI workflow.

---

## 1. The cofactor call

### 1.1 The equivalence being applied

PQQ-dependent alcohol dehydrogenases of the methanol/ethanol family hold either a lanthanide
or a calcium ion, and the difference is legible in sequence:

| active-site residues at the calibrated column | call |
|---|---|
| `D-x-D` (Asp, any, Asp) | lanthanide |
| `D-x-[ATS]` (Asp, any, Ala/Thr/Ser) | calcium |
| a gap at any of the three positions | gap — no call |
| anything else | unknown — no call |

A sequence whose model string ends before the column is truncated, also no call.
`src/ModuleHmmer.jl` applies the rule residue by residue: `calibrate_active_site` finds the
column (§1.3), `parse_hmmalign_results` reads every aligned sequence there. The same two
patterns, in PROSITE notation, are matched with `BioSequences.jl` in
`src/ModuleMotifSearch.jl` (`classify_pqq_adh`), which searches a whole sequence and answers a
different question: whether the motif occurs anywhere.

The rule is the published lanthanide-switch equivalence (Voutsinos et al. 2025,
doi:10.1093/ismejo/wraf057). Its structural basis, the additional aspartate two positions
downstream of a conserved one where calcium enzymes carry Ala, Ser or Thr, is cited in
`src/ModuleHmmer.jl` to Good et al. 2020 (doi:10.1074/jbc.RA120.013227). This repository adds
reading the rule at a calibrated model column instead of a sequence position, a harness
measuring whether it recovers the cofactor on enzymes whose cofactor is known, and the tooling
around the call.

### 1.2 Model coordinates

A residue's position differs between sequences, so the call is read at a profile-HMM column.
`hmmalign` writes `afa`, where upper-case characters and `-` are match-state columns and
lower-case characters are insertions; `src/ModuleHmmer.jl` strips the insertions, making the
remaining index the model coordinate, identical across every aligned sequence.

The profile is `dat/models/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm`, third-party and
documented in `dat/models/PROVENANCE.md`. It carries no GA/TC/NC cutoff, so it sets no
inclusion threshold; a threshold belongs to the search command, and the background search
calibrates and records its own (§4.2).

### 1.3 Calibrating the active-site column

`calibrate_active_site` takes the first two sequences of the alignment, one enzyme known to
bind a lanthanide and one known to bind calcium, and returns the first model column where the
first reads `D-x-D` and the second `D-x-[ATS]`. Which two records those are is the caller's
decision, never read from a data file: the census driver (`src/scriptLn.jl`) calls
`calibrate_active_site(paths.afa)` on a FASTA that carries the pair first, while `validate`
and `validate_from_afa` (`src/validate_classifier.jl`) are told accessions, as `ln_ref_acc`
and `ca_ref_acc`, and `reorder_references!` moves those records to the top first, so record
order cannot decide what a run calibrated on. The background scan takes that route with the
pair fixed in its own source (`LN_REF`, `CA_REF`, §9). Either way a run states what it
calibrated on.

The column is a property of the profile and those two sequences, computed at the start of
every run, written into its output, hard-coded nowhere. If more than one column satisfies both
patterns, the first is not necessarily the catalytic one and every downstream call inherits
the choice; `diagnose_calibration` counts those collision columns.

### 1.4 Control-set validation

The call is validated against `dat/control_set/known_cofactor_mdhs.tsv`, whose every row is a
factual claim about a published enzyme: accession, gene, organism, domain, the metal, a
one-line statement of the evidence, the citation, any PDB entries holding that metal, and the
row's provenance. Its header documents the schema in full.

Rows marked `reference` are calibration references, passed by accession (`ln_ref_acc`,
`ca_ref_acc`); none enters a held-out stratum, so none can flatter the held-out result, and
they are reported in a `reference` stratum of their own and in the evaluation over all rows.
Only rows marked `test` enter a held-out accuracy. Evidence tiers, strongest first:

| tier | meaning |
|---|---|
| `ref` | a reference row; belongs to no held-out stratum |
| `A` | direct determination for this protein: the metal resolved in a crystal structure of the enzyme itself, or measured in or required by the purified enzyme |
| `B` | indirect but gene specific: expression, deletion mutant or growth evidence ties this gene to a Ca- or Ln-dependent activity, but no metal was measured in, reconstituted into or shown to be required by the purified enzyme |
| `C` | homology only: the cofactor is inherited from close relatives and was never measured for this protein |
| `X` | contested; reported with its call and excluded from every accuracy |

`ValidateClassifier.validate` loads the table, fetches each sequence from its database of
record (UniProtKB by default; a row's `source` column can route it elsewhere, and a deleted
entry is fetched from UniParc under the UniParc id the row's evidence cell also names), aligns
with `hmmalign`, refuses rather than scores an alignment whose match-column count differs from
the profile's `LENG`, moves the two reference rows to the top, runs the calibration
diagnostic, calibrates, then reads and scores the call at that column for every row.
`validate_from_afa` is the same tail from an alignment made elsewhere, which takes HMMER off
the critical path.

Reported: sensitivity, specificity and accuracy over all scorable rows and over the held-out
rows alone, a per-stratum summary, and exact (binomial) confidence intervals. A gap, unknown
or truncated call counts as a misclassification of the row's known cofactor; a contested row
is listed with its call and never scored. Each run writes a fetch manifest beside the
predictions: which accession resolved to what, from which database, with the release header
the service reported. Close relatives are not independent, so held-out sequences are clustered
by single linkage on pairwise identity over the match columns, 0.9 in both entry points
(`identity_clusters` takes another as `threshold`), and the cluster count is reported beside
the row count.

This measures whether the residue rule recovers the cofactor experiment reports on labelled
proteins. It establishes the cofactor of no uncharacterised sequence; that needs structure and
measurement.

---

## 2. The census pipeline

`src/scriptLn.jl` is the top-level driver: one linear pass from an alignment of search hits to
a census table with taxonomy and habitat attached. It reads every large input by name from the
data root (`$LNFINDER_DATA`, see [DATA.md](DATA.md)), where `LNFINDER_AFA` and
`LNFINDER_TAXDUMP` can move two of them; a missing input stops the run, naming the input,
where the layout expects it and which variable moves it. It writes under
`LNFINDER_CENSUS_OUT`, which has no default.

The census begins from the sequences in `census/firstTry.fa` (`hits_fasta`), the hits of an
`hmmsearch` with the shipped profile run outside this repository, and from their alignment to
the profile (`model_alignment`, which `hmmalign` regenerates from them; [DATA.md](DATA.md) §5).
Its per-residue and per-sequence tables come from `src/ModuleSequenceTools.jl`:
`count_residues_per_position` tabulates residue counts by model column across an aligned set,
`export_tree_fasta` writes the per-class FASTA a tree is built from elsewhere. That module has
no suite of its own; the census driver is its only caller.

1. Calibrate the active-site column on the alignment (§1.3).
2. Classify every aligned hit at that column (§1.1).
3. Dereplicate to cluster centroids. Clustering runs outside this repository with USEARCH
   `-cluster_fast` at 90 % identity; the census reads the `-uc` file. The identity is
   `DEDUP_ID` (0.90), beside `USEARCH_CMD`, which documents the clustering command. That
   command is never run from here, so the two are kept in step by hand and the manifest
   records both. `EXPECTED_CENTROIDS` optionally aborts a rebuild unless the clustering yields
   exactly the recorded centroid count.
4. Resolve taxonomy: each centroid's accession to a TaxID over UniProt, its lineage from an
   unpacked NCBI taxonomy dump (§3.1).
5. Attach habitat (§3.2). `RUN_HABITAT=false` skips it, for runs that need only the census
   table.
6. Write the census and its manifest, then compute statistics.

Two ordering rules are load-bearing. The full census and its TOML manifest are written before
anything is filtered, since rows without a TaxID, a kingdom or a habitat belong to the census
and dropping them before the write makes exclusions unauditable. And nothing mutates the table
after that write, every analysis subset being built with `filter`.

The manifest records each input's path and SHA-256, or an explicit marker when a file was
absent, so a missing input records itself instead of aborting a run whose census is already
written. It also records the run date, the calibrated column (§1.3), the dereplication
settings, the UniProt release the TaxIDs were resolved against, the taxonomy dump's own
provenance (`taxdump_provenance` reads the dump's release files and digests) and the versions
of the direct dependencies. A census is never overwritten: the driver refuses to start if the
output file exists.

---

## 3. Annotation layers

### 3.1 Taxonomy (`src/ModuleTaxonomyTools.jl`)

The NCBI taxonomy dump is loaded into parent/rank/name maps plus the merged and deleted id
lists. Two dump facts are handled explicitly, because both silently corrupt a kingdom column.
NCBI reorganised the top of the taxonomy in March 2025: `superkingdom` became `domain`, a real
`kingdom` rank was inserted below it, and viruses moved under an acellular root. Both vintages
are handled, and NCBI's own `kingdom` rank is deliberately not used. And a merged TaxID is
followed to its successor, a deleted one reported as deleted rather than dropped.
`append_taxonomy!` refuses to return a census whose kingdom column is mostly empty, so a
dump-vintage mismatch stops the run.

The TaxIDs come from UniProt, over three stateless requests (`src/ModuleAccessionMapper.jl`,
`resolve_taxids`), issued one at a time and paced: a UniProtKB search, 100 accessions per
request (`SEARCH_BLOCK`, the service's own cap on one query), which gives live entries their
organism and hands back a UniParc id for every deleted one; UniParc, for a deleted entry whose
record names a single organism and where the accession reappears on that record; and, for a
record shared between organisms, the accession's own cross-reference on it. UniParc's
`commonTaxonId` is never a fallback, and the route uses neither UniProt's id-mapping jobs nor
its stream endpoints, both named in robots.txt.

Every search block is checked against the accessions asked for; a short block is asked once
more, and one still short, or naming an accession nobody asked for, stops the run. A ledger
accounts for every accession against the output table, and the run stops when UniParc resolves
fewer than 98% of the deleted entries (`min_rescue_frac`). Each row records the UniProt
release the service reported. `accession_taxid_cache.csv`, in the run's output directory,
holds only established resolutions and is written after every check has passed; the third
request is checkpointed per UniParc record and accession in `uniparc_xref_cache.csv` beside
it, keeping only answers shown to be complete; and accessions still unresolved go to
`unresolved_accessions.csv` before the 98% floor is checked. Status decides the retry policy
(§9): a 5xx is retried; a 429 is waited out at most three times per request, a fourth stopping
the run; any other status, a redirect included, stops the run at once.

Two steps use a stream route. Step 0 (`scripts/prokaryote_background.jl`) asks
`https://rest.uniprot.org/proteomes/stream?query=` for a domain's reference-proteome taxa, and
step 2 (`scripts/background_pqq_scan.jl`) asks it for a sampled taxon's reference proteome and
`https://rest.uniprot.org/uniprotkb/stream?query=` for that proteome's sequences, one taxon at
a time, writing each answer to a file it reuses. Read on 2026-09-23,
`https://rest.uniprot.org/robots.txt` disallows `/docs/`, `/idmapping/`, and a
`/<dataset>/stream/` path per dataset, `/proteomes/stream/` and `/uniprotkb/stream/` among
them; those patterns end in a slash and these two routes carry a query string in its place, so
the requests fall outside the literal prefixes while asking of the same endpoints. No route
here submits an id-mapping job.

### 3.2 Habitat (`src/ModuleEnvironmentTools.jl`)

Habitat comes from OmniCrobe, which text-mines habitat relations from published work and
expresses each habitat as an OntoBiotope (OBT) term. Two consequences shape the method. A raw
habitat tally measures research effort, so every quantity downstream is a comparison against a
denominator computed from the same source (§4). And a term is a node in an ontology: it is
climbed to its OBT root path, and the path decides the category. A run writes the climb to
`habitat/obt_path_cache.csv` under the data root (see [DATA.md](DATA.md)) and every later run
on that root reads it; the data set does not carry that cache, and on a root without it the
climb goes back to the network term by term.

Mapping rules:

| step | rule |
|---|---|
| habitat type | the term's OBT path is tested against a table of branch ids, in order, first match wins. The specific branches (aquatic, terrestrial, geological, host-associated, anthropogenic, …) are tested first |
| habitat chemistry | a second, non-exclusive set of condition ids (e.g. thermal, saline). A taxon may carry several, or none; these are flags and do not partition the taxa |
| isolation sources | terms recording where a strain was *obtained* instead of where it lives. They are commoner in a background drawn from culture collections, so leaving them in manufactures a fake depletion. The primary rule drops them term by term, keeping a taxon that also has genuine habitat terms; the legacy rule drops the whole taxon |
| generic terms | terms too general to carry information are excluded from the winner decision |
| distinct-term cap | at most `MAX_DISTINCT` (50) distinct terms per taxon enter the tally |
| ties | on a tied deciding tally the taxon is split between the categories with weight 1/k |

The branch table is built from the OBT paths actually returned for the terms in the data, and
each entry carries the ontology's own label for the id it names, so a mapping can be checked
against the source without a network call. Several entries carry a comment saying why a branch
is split by child instead of labelled at the parent; those comments are part of the method.

### 3.3 The habitat call table (`scripts/habitat_reaggregate.jl`, pipeline step 1)

A habitat number from a live re-fetch depends on what the service answers that afternoon, on
whether the fetch succeeds at all (a failed fetch enters the analysis as the ecological
statement "this organism has no known habitat"), and on dictionary iteration order where two
categories tie. The remedy is to freeze the inputs: a run stores the raw response bodies per
taxon under `habitat/omnicrobe_cache/relations/` in the data root and the ontology climb in
the OBT cache, then computes every rule variant from those same bytes in one pass, so for
whoever holds that cache a re-run with the network unplugged reproduces the table exactly.

The bodies and the OBT cache are third-party response content and are not redistributed: the
data set carries this step's output, the call table the later steps read, plus the cache's
manifest. A run that starts without the cache fetches, and what it fetches is a new
measurement. The manifest stamps each body with where it came from and when, and with the
service and ontology releases, which are constants in the script (`OMNICROBE_RELEASE`,
`OBT_RELEASE`), not values read from the responses. A taxon missing from the cache is fetched
live, one at a time, under the retry rule in §9; a failed fetch is recorded as a failed fetch,
never as a taxon with no habitat.

The output is one row per (configuration, taxon): the winning category, the raw distinct-term
and mention tallies per type, the 1/k tie weights, the non-exclusive chemistry flags and the
isolation/generic bookkeeping, so no downstream rule variant costs a refetch. `primary` is the
rule this release computes with. `legacy` reproduces the rule the census and background-pool
tables carry their calls under, so those calls can be re-derived and diffed field by field.
The rest each move exactly one axis off `primary`; they are the levels the sensitivity grid
(§7) sweeps.

Three kinds of gate stand between a run and its outputs. Method gates, whose expected value is
a requirement of the method, are checked on every run: no ontology path missing after the
fill, no failed relation fetch. Measured gates, whose expected value is a *measurement of the
dataset*, are coverage and isolation-term recovery under `primary`; such a value is never
typed into the code but arrives in a `--gate-targets` file, one row per gate with a note
saying where its number came from, and a run given no such file is refused unless `--no-gates`
says so, under which each measured gate that is not also a reproduction gate prints NOT
CHECKED.

The four reproduction gates have the `legacy` configuration re-derive the calls the census and
the background pool carry in their own habitat columns, habitat type and chemistry, on the
census's lanthanide arm and on the background, and compare them taxon by taxon. That is a
question about the census those calls were made on; asked of another census it is a different
test, one an honest run could fail and whose pass would read as evidence that two censuses
agree. So the four are evaluated only when the `--gate-targets` file declares, by SHA-256
digest, the census they were written for (one row `census,,<sha256>,,<note>`: no gate named,
no tolerance, a note required; the run prints the digest it read as `Census sha256:`) and the
census read has that digest. Otherwise each is declined, printing NOT APPLICABLE with the
reason, and still counts as reached for the check that every declared gate was evaluated, so
declining one cannot hide a deleted one.

A run's closing line says which gates it did not make. Step 1 on the deposited census runs
with `--no-gates` ([DATA.md](DATA.md) §4), so only the method gates are made.
`--census <path>` runs the step on another census and the report names the file read. Taxa with no
phylum are counted per analysed arm in the report and not gated here; step 3 enforces that
rule where the arms are built (§6.1). Outputs go to a `_staging` directory beside their
destinations and move into place only if every gate passed.

---

## 4. The comparison arms

### 4.1 The denominator (`scripts/prokaryote_background.jl`, step 0)

"Where do the carriers live?" needs a background of comparable organisms put through the
identical pipeline. The pool is UniProt reference proteomes for Bacteria and Archaea. Since
the census hits are themselves proteins from sequenced organisms, the pool shares the census's
strongest bias, that only cultured, sequenced, named organisms are represented, and that bias
largely cancels in the contrast. UniProt deduplicates reference proteomes to roughly one per species group,
the sampling-frequency correction wanted. Each taxon's lineage is fetched alongside its
identifier, so a phylum-matched background can be built downstream (§6). Above a pool-size cap
the taxa are subsampled with a recorded seed; cap, subsample size and seed are named constants
(`MAX_TAXA`, `SUBSAMPLE`, `SEED`). The run is resumable: taxa already in the output are
skipped.

Step 3 records that provenance on every run, under `background_pool_*` and
`background_frame_source_*` in its summary table: the seed, the subsample size and the cap,
read from the generator's own source and cited there by file and line; the generator's digest;
the file name, digest and taxon count of the pool and of the frame table; and whether the pool
on disk holds exactly the subsample size, which is what says the subsample was taken. A value
that cannot be read is recorded as absent with a state saying why
(`background_pool_provenance_state`); a default there would be step 3 inventing the seed of a
subsample it did not take.

### 4.2 Carriers hiding in the frame (`scripts/background_pqq_scan.jl`, step 2)

The background is the sampling frame minus the census: reference-proteome prokaryotes that do
not carry the enzyme family. Subtracting by taxonomy alone leaves a known contaminant, since a
taxon can be outside the census for having no such enzyme or because the census search did not
see the one it has (a proteome released later, an entry whose TaxID moved, a sequence below
the census threshold). Those carriers sit in the comparator arm and bias every contrast toward
the null, and nothing else in the pipeline asks whether a background taxon has the gene. This
step asks, with the census profile itself.

Fetching every background proteome is several gigabytes for a question about one domain
family, so the search space is narrowed first, to the UniProtKB entries InterPro places in
either of two entries:

| InterPro entry | what it is |
|---|---|
| IPR011047 | quinoprotein alcohol dehydrogenase-like superfamily (homologous superfamily) |
| IPR017512 | PQQ-dependent dehydrogenase, methanol/ethanol family (family) |

The union is taken by accession, locally. The superfamily alone is not enough: InterPro
annotates some family members without the superfamily entry, and under a superfamily-only
prefilter such a taxon never reaches the search and comes out with zero hits however
permissive the cutoff. Whether the union's entries reach the cutoff is measured. Each download
keeps its own release record, and the run stops if the two halves came from different
UniProtKB releases.

The search is `hmmsearch` with the shipped profile at a bit-score cutoff, writing a `--tblout`
table the script parses. On Windows the invocation goes through WSL with the paths translated,
HMMER having no native Windows build; on Linux and macOS `hmmsearch` is called directly, with
the same option line. See the platform note in `README.md`. The cutoff is calibrated at run
time, a hand-picked bit score being unreproducible: it is set from the distribution of the
census accessions' own scores against the same profile by a named rule
(`calibration_threshold(scores; rule = :p1)`, the first percentile), written to a record file
beside the scan output and reported by the run. It is no constant of this repository.
`--threshold` sets another, and the run records that it did.

The prefilter is an assumption, so it is measured. A sample of frame taxa has its full
reference proteome downloaded and searched at the same cutoff; of the taxa the full search
flags, the share the prefilter also flags is the measured recall. The sample is the first *n*
draws of one seeded permutation of the frame, so a larger sample keeps every taxon a smaller
one drew, and the run checks that against the sample column of an earlier scan table. A sample
in which the full search flags nothing is recorded as a miss: it measured no recall at all.
The run records the observed rate and a Clopper–Pearson lower bound on it, because "does the
prefilter find what the full search finds" and "does the remaining background hold no missed
carrier" are different questions. The bars it is held to (`MIN_RECOVERY`, `MIN_RECALL`) and
the confidence level of the bound (`RECALL_ALPHA`) are constants fixed before a run, as plan
parameters rather than measurements. A run that misses promotes nothing.

Acceptance also requires five named frame taxa to come out flagged, one row each against a bar
of one (`POSITIVE_CONTROLS`). They were not drawn from the literature, and this release claims
no published cofactor measurement for any of them. Each was selected from public UniProt
annotation: a frame taxon, in the reference-proteome pool and absent from the census, with a
UniProtKB entry of at least `MIN_FAMILY_LEN` residues annotated to InterPro IPR017512, the
methanol/ethanol family, under that exact TaxID. The five entries are `A0ABM9ICZ5` (TaxID
591154), `A0ABX7PVQ1` (511746), `A0ABU0I1Z8` (418708), `A0ABU0LS59` (217069) and `A0ABQ2ZAD6`
(980347), named in the script beside the taxa they select.

Two things decide how a PASS reads. The list was fixed before this search was first run, so it
cannot have been read off the scan's own output, and these five neither top nor order its hit
counts or bit scores. And because the criterion is annotation, the check asks whether this
prefilter and this profile recover what public annotation already claims: a wiring check on
the search, and **not evidence about these organisms' cofactors**, those five annotations being
themselves inferred from homology. One of the five species is the subject of a published
lanthanide-in-the-active-site determination, but for a different strain, and a different
TaxID, than the one used here; see the note above the constant, and
`dat/control_set/known_cofactor_mdhs.tsv`, which carries that evidence against the strain's
own accession. That note also states, taxon by taxon and with DOIs, what the literature does
and does not establish for each of the five, including, for one, published evidence against
growth on methanol from an assay that could not have detected a lanthanide-dependent enzyme.
The acceptance row of a control present in the frame says what its hit count measures:
sequence presence only, not that the organism is methylotrophic.

---

## 5. Where the numbers come from

Three conventions run through every analysis script.

No data literal in a script. Counts, shares, coverages, releases and test statistics are
computed or read at run time, written to a summary table, and read back out of it for the
closing report.

Bounds the analysis plan fixed before the run are supplied on the command line, never
defaulted, and are reported beside the quantity they bound with the comparison and whether it
was met. Omitting one omits its row; nothing passes silently because a bound was not supplied.

Estimators, rules and parsers live in one place: the estimators in
`src/ModuleHabitatStats.jl`, the outcome predicates in `src/ModuleOutcomeRules.jl`, each
persisted record format read by its own producer's reader.

### 5.1 What "the registered specification" points at

This document and the scripts use that phrase, and "the analysis plan", often. They name an
external record, and this release does not carry it; there is no plan document in this
repository to go looking for.

The record is a dated, checksummed set of rules mapping each outcome to a colour, to a caption
phrase and to a sentence of the accompanying text. Its timestamp precedes every output of the
re-analysis these steps perform, so the numbers reported here did not choose their own
wording. It was not written before the data were first looked at, and does not claim to have
been: its own disclosure records that an earlier run's estimates for the categories it names,
and that run's sensitivity findings, were already known when it was written. It is registered
against the re-analysis, and the categories it speaks for were selected from those earlier
estimates, which is why they are a per-run input here (`docs/DATA.md` §4). That is what "the
registered specification", "the registered primary" and "the analysis plan" point at, and
nothing stronger. In the code that is here, they surface as:

- The primary configuration: one default option set per script, printed at the start of a run
  and recorded in its own output. "Off the registered specification" means "from that
  configuration, one choice changed".
- The criteria supplied per run: `--gate-targets` (step 1), `--anchor-targets` (step 4),
  `--anchor-reference` (steps 3 and 4), `--categories` (steps 3 and 5), `--plan-*` (step 3).
  The two targets files and the category selection carry a note on every row saying where its
  number or selection came from; the anchor targets also name their source files with digests,
  which step 4 checks (§7); a step-1 targets file names, by digest, the census its
  reproduction gates were written for (§3.3); and the steps that read `--categories` record
  that file's digest. The anchor targets for the deposited data ship here as
  `dat/criteria/anchor_targets.csv`, with their two source files (`impl_4.csv`, `impl_w1.csv`)
  beside it.
- The criteria described in the script itself. Where step 3 says "see THE RECORDED CRITERIA
  below", that section states which bounds the plan fixed and how a run reports them; their
  values are the `--plan-*` bounds above, which this release does not carry
  (`docs/DATA.md` §4).

---

## 6. Enrichment estimators (`src/ModuleHabitatStats.jl`, `scripts/habitat_enrichment.jl`, step 3)

### 6.1 Two estimators

| estimator | definition |
|---|---|
| crude | (a/n1) / (b/n2), as a log2 risk ratio, with a Katz Wald interval and a Fisher exact test. When a or b is zero, and only then, a Haldane correction (+0.5 on the hit cells, +1 on the denominators) is applied to this estimator only. Fisher's test is taken on the hit totals rounded to counts; the estimate and its interval use them unrounded |
| standardised | the background's share of the category reweighted to the foreground's phylum distribution: p0_std = Σ_k w_k (c_k / n0_k) with w_k = n1_k / n1. No correction: with no foreground hit, or a standardised background share of zero, the estimate is undecidable (NaN) |

The two are reported together, because the gap between them is the phylogenetic contribution:
habitat annotation tracks phylum, so a crude ratio can report enrichment that is only the
phylogenetic range of the protein family. Both intervals are 1.96 standard errors either side
of the estimate, the 95% level.

Further rules of the estimator layer:

- Fractional hits. A hit vector may be boolean, or real where a unit is split across tied
  categories (§3.2); a boolean vector and its 0/1 real equivalent give bitwise identical
  results. For a fractional vector the delta-method variance is an upper bound, so the
  intervals are conservative. That bound is about the form of the variance, not how well its
  plug-in value is estimated, and it does not hold for a stratum whose background cell is
  empty or thin.
- The same unit on both sides. Counting sequences counts how often an organism was sequenced,
  so the primary unit is the unique taxon. Per-sequence and per-genus series are secondary,
  keeping the size of that bias visible, and the genus series collapses both arms to genus.
- Thin cells are disclosed. A row whose combined hit total a + b is below `MIN_CELL` (5) is
  reported, never tested. A phylum thinner than `MIN_FG_K` (10) in the foreground or
  `MIN_BG_K` (20) in the background is relabelled into a single "Other phyla" stratum in both
  arms before any counting; a stratum with foreground units but no background units is
  dropped, and the number of contributing strata is reported. That collapse bounds the weight,
  not the emptiness, so every row also carries how many contributing strata have an empty
  background cell and the smallest background cell among them, unrounded.
- A taxon with no phylum is in no arm, phylum being the stratum. NCBI's catch-all placeholders
  (`bacterium`, `uncultured bacterium`, `uncultured archaeon`) are such taxa. The exclusion is
  a rule, applied to every arm step 3 builds: the census foreground, the taxa a move-in level
  adds to it, the comparator arms and the background. The run stops, with a count per arm, if
  an assembled arm still holds one. A phylum-less frame taxon stays in the frame the scan's
  coverage is compared against; it leaves only the arms.
- Multiple testing. Benjamini–Hochberg within each contrast × unit panel: one contrast, in one
  unit, over its habitat categories or over its chemistry flags, never across panels or across
  a run. It is applied to the crude estimator's Fisher exact p-values of the tested rows
  (a + b ≥ `MIN_CELL`); an untested row carries no q-value, and no colour or outcome code
  reads one.

The anchor. Asking whether a category is new needs an earlier result to be new against, so
step 3 computes one in the same run: the same builder on the `legacy` configuration, its own
chemistry column, no background exclusion, nobody moved in. Whether the partition may use that
anchor's category set is decided in a fixed order of statuses, each rung withholding it for a
reason that makes the later rungs meaningless: `undecidable_config_absent`,
`undecidable_frame_mismatch` (the frame each specification drew from, not the arm each kept,
so a move-in level is no mismatch), `undecidable_anchor_intervals`, `computed_unverified`,
`undecidable_reference_missing`, `undecidable_reference_unreadable`, then
`reproduces_reference` or `differs_from_reference`. Only `computed_unverified` and
`reproduces_reference` let the partition use the anchor. `anchor_verified` in
`habitat/habitat_arm_summary.csv` is the narrower statement that `anchor_status` is
`reproduces_reference`, so a run given no reference reports `computed_unverified` with
`anchor_verified` false, which reads as unchecked. A supplied reference is also compared row
by row (rows matched, agreeing, disagreeing, largest cell difference); that comparison decides
nothing.

Move-in pools. One background-exclusion level adds instead of only removing:
`hmmsearch_movein` takes the carriers out of the background frame and gives the
lanthanide-type ones to the foreground. It is a sensitivity level, never the headline, and two
rules can decide which taxa move. `census_lanthanide_hits`, the one the recorded level names,
is frame taxa carrying at least one sequence the census itself classified lanthanide-type,
read from the census's pre-deduplication hit table, which is not a repository input and has to
be supplied (`--movein-pool`). `scan_lanthanide_hits` is frame taxa whose scan row carries a
lanthanide-classed hit: this run's own model call on the family members the background scan
fetched, a different measurement, and the rule a run falls back to when no hit table is
supplied. A run records the rule applied (`movein_pool_rule`), whether it was supplied or
fallen back to (`movein_pool_rule_source`), whether it is therefore the registered one or
descriptive (`movein_pool_rule_role`), the rule's own note, both rule names, and the pool's
source, digest, row counts and size; where it can, it also records the other rule's pool size
and the overlap. A fallback warns as well as records. Only this level emits those keys, and a
pool that cannot be answered stops the run instead of turning it into the mode it would fall
back on. No taxon is counted on both sides.

### 6.2 Outcome rules (`src/ModuleOutcomeRules.jl`)

The predicates that turn an estimate into a decision are a module of their own, with unit
tests. Their design principles:

- Three-valued interval tests: excluding zero, including zero and *undecidable*. An estimate
  that was never computed is undecidable.
- Anchor-relative. The earlier (anchor) set is an argument; there is no hard-coded list.
- A tier that needs a grid may not be claimed without one. The strong robustness tiers are
  defined by what the one-at-a-time grid does, so without a grid the answer is undecidable,
  except for the one tier that needs no grid.
- Neutral codes: every function returns symbols and label lists.
- The module's own constants are rule parameters. The sharing bands, the confidence level and
  the thresholds that grade a ratio are the recorded rule's parameters, declared as such, and
  are not measured from any dataset.

---

## 7. The sensitivity grid (`scripts/habitat_sensitivity.jl`, step 4)

Every estimate rests on rule choices an interval does not show: which habitat wins when a
taxon's terms tie, whether a chemistry flag needs a share of the evidence or a single mention,
which text a keyword fallback reads, which background counts as carrying no such enzyme. This
step varies each choice one at a time off the registered specification and reports what
survives.

No new code computes the grid: every level reaches the estimators through the same aggregation
rules, arm builder and predicates the primary run used. Two seams are not arguments of the
estimator entry point and are crossed explicitly. For the type rule, the arm builder reads the
call table's weight column verbatim, so a type variant rewrites that column from the same
table's own term and mention tallies and then calls the unchanged entry point; before it is
used, recomputing the *registered* rule from those tallies must reproduce the committed column
exactly. For the chemistry rule, the flag route cannot express a fractional weight, so a
chemistry variant is carried by the weight column instead, and under the registered rule the
two routes must agree cell for cell. Weights are not rounded on the way in, a three-decimal
round trip changing a hit total, so a rewritten column is carried at full precision and only
the *check* against the published column uses the published format.

The flag reads the primary interval and the one-at-a-time variants on the two units whose
denominator is a taxon. Unit checks (the symmetric genus unit, the per-sequence unit, and the
subset of taxa whose call rests on several distinct terms) are reported for the same grid and
marked as *not* feeding the flag, a unit check changing what is counted instead of which rule
counted it. Anchors and disclosed background levels sit outside the variant set for the same
reason.

Two anchors outside the grid reproduce earlier runs of this analysis: the all-legacy
specification (`anchor:all_legacy`), compared row by row against an earlier run's rows table
(`--anchor-reference`), and four levels reproducing an earlier sensitivity run's chemistry
rules on the legacy arms. Both reproduce rows of another census, so they are a different test,
for the reason §3.3 gives of its reproduction gates. `ANCHOR_LEVELS_RETIRED` keeps them out of
the grid: their definitions stay in the code, and flipping that one constant builds them again
for a run on the census they describe. While it holds, the run records
`anchor_levels_retired` and an `anchor_reference_status` of not applicable, the targets file
declares `reference_rows,0` to match, `--anchor-reference` is not required, and a reference
given anyway is recorded and not opened (provenance key
`anchor_reference_named_not_opened`, its file name only; no digest).

The wiring check rests instead on the anchor targets (`--anchor-targets`; the file this
repository ships for the deposited data is `dat/criteria/anchor_targets.csv`). A target is one
quantity the grid writes for one row (`a`, `b`, `n1`, `n2`, `crude`, `std`, `lo` or `hi`) on
the registered specification or one of its one-at-a-time levels, with a tolerance, transcribed
from the output of an independent implementation of the specification (`README.md` describes
those implementations and what the specification left to them). Before anything is published,
step 4 requires that:

- every target is reproduced within its tolerance;
- every level named in a published breaker list is pinned by a target;
- the file's own declared counts, of targets and of levels, whose names must equal the set the
  script enumerates, match what it carries;
- each file the targets name as their source, where it sits beside the targets file, has the
  digest the targets declare. The repository ships both sources there. An absent source is
  reported; one present with other bytes stops the run, the targets having been transcribed
  from a version this run cannot see.

Running without the targets has to be typed (`--no-anchors`, which declines the targets and
the reference together) and is recorded as a decision (`anchor_checks_declined`). Step 5
refuses to draw from a grid whose anchor checks were declined or whose targets were not all
met.

---

## 8. The enrichment panel (`scripts/habitat_figure.jl`, step 5)

The panel is a drawing procedure; what it draws is data the run reads.

Two panels share one row per category. The upper panel shows each arm's share of taxa in the
category, the carrier arm a solid bar inside the background arm's outline on a common
baseline; it runs no test and carries no test's colour. The lower panel shows log2 enrichment
against the background, three marks per row on three sub-rows: the crude estimate as an open
ring, the phylum-standardised estimate as a filled dot on its interval, and above both a thin
tick for the same standardised estimate computed for the comparator arm. A connector joins ring and dot, drawing the
contribution of phylum composition, which has no fixed direction. Chemistry shares the lower
panel's axis as a second block below a rule, the same estimator in the same units, and its
flags are not exclusive of each other, which the block label says.

Rules the drawing obeys:

- Colour is the robustness flag and nothing else. A dot is coloured only where the sensitivity
  grid's flag gives the row a sign that the primary specification and every one-at-a-time
  variant reproduce; an interval that excludes zero under the primary specification alone does
  not colour a row. A flag that contradicts its own interval stops the run.
- Marks may not overlap, so the tick has a sub-row of its own. The offsets and marker sizes
  are named constants whose clearances the test suite checks as arithmetic, and the panel
  height derives from a row pitch.
- Which categories are drawn is decided by a minimum weighted count, computed on the
  fractional foreground count and blind to direction, so it cannot select on the outcome.
  Held-back categories are reported with their counts.
- Nothing is published unless the inputs describe one analysis. The rows, the arm summary, the
  flags, the grid and the grid's provenance come from separate runs; every identity between
  them is a named gate, all failures are reported together, and one failure stops the run
  before anything is staged. A flag is recomputed from the grid rows it summarises, and
  outputs are staged side by side and moved into place together.
- Provenance names files, not locations: each input's file name, whether its path was typed or
  defaulted, and its digest. A value that looks like an absolute path stops the run.

The declared output set is six files, all written into the figure directory (`--figdir`, else
`$LNFINDER_FIGDIR`; the step never creates one):

| file | what it is |
|---|---|
| `habitat_enrichment.pdf`, `.svg`, `.png` | the panel, the same drawing in three formats |
| `habitat_enrichment_source_data.csv` | one row per mark the panel draws or holds back, with the role that says which |
| `habitat_enrichment_caption_values.csv` | the numbers a caption would otherwise be typed from, exported as `key,value` |
| `habitat_figure_provenance.csv` | each input's file name, whether its path was typed or defaulted, and its digest |

Staging is a directory of the step's own (`_staging_habitat_figure`) beside the destinations,
so a failed run of this step never clears what another step's failed run left for inspection
there. Nothing moves into place unless the staged set is exactly those six names.

The caption-values file follows §5's rule against data literals: every number a caption would
state is derived here from the same inputs the panel is drawn from and exported under a fixed
key, the key set declared in advance, so a key the run cannot fill stops the run and a key
nothing asks for is refused. The file is an output of a run and is not shipped.

Two exports are about a named category and touch neither what is drawn nor what is coloured:
the interval exported across every specification in the grid, and the comparator tick divided
by its own standardised estimate to give the sharing ratio. Which categories those are is not
in this script and is not a default of it: they arrive in the `--categories` criteria file,
the same file step 3 was given, read with step 3's own reader, and the run refuses by name
when it is absent. `caption_keys` takes the selection as a required argument, and
`caption_values`, `derive` and `run_figure` take it with the inputs instead of as a defaulted
keyword, so no code path can supply a category this script chose. Each run records the
selection and that file's digest in its provenance. Where the selections came from is stated
in `docs/DATA.md` §4: not from this software, and not fixed in advance of the estimates.

---

## 9. Parameters

Only parameters that shape the method are listed, and anything calibrated is named with its
calibration procedure and the file its run records it in, never with a value. A library
default is the chemistry or the convention, never a value measured from the data. The
active-site model column has no default anywhere: every entry point that reads a column takes
it as a required argument, so no call can silently inherit one profile's calibration.

| parameter | where | meaning | how it is set |
|---|---|---|---|
| active-site model column | `calibrate_active_site` | the column the cofactor call is read at | calibrated per run from two reference sequences; written to the run's output |
| collision-column count | `diagnose_calibration` | columns satisfying both reference patterns | measured per run; above one, the calibration is not unique |
| control-set tiers admitted | `include_tiers`, `DEFAULT_LABEL_TIERS` | which tiers may be scored, and which may become scan labels | fixed rule: homology-only and contested tiers never become labels |
| identity threshold for held-out clustering | `validate`, `identity_clusters` | when two held-out controls count as one | 0.9 pairwise identity over the match columns in the validation entry points; `identity_clusters` takes another as `threshold`. Reported with the result |
| `DEDUP_ID` (0.90) | `src/scriptLn.jl` | the identity the census clustering was run at | fixed to match the recorded external clustering command |
| `EXPECTED_CENTROIDS` | environment | abort unless the clustering yields exactly this many centroids | supplied only when rebuilding a recorded run |
| `SEARCH_BLOCK`, `UPIS_BLOCK` | `src/ModuleAccessionMapper.jl` | accessions per UniProtKB search, UniParc ids per UniParc request | fixed: 100, the service's own cap on one query; 500 |
| `PACE_S`, `PASS3_PACE_S`, `HTTP_TRIES_5XX`, `BACKOFF_429_S` | `src/ModuleAccessionMapper.jl` | request spacing, and what a 5xx or a 429 costs | this client's own restraint, not a limit the service advertises: 1 s between requests (0.3 s in the third), 5 attempts on a 5xx, waits of 60, 120 and 300 s on a 429 (or the server's `Retry-After`, never below `PACE_S`), three per request, a fourth stopping the run (§3.1) |
| rescue floor (`min_rescue_frac`) | `resolve_taxids` | share of deleted entries UniParc must resolve | 0.98, a keyword default; below it the run stops (§3.1) |
| `MAX_DISTINCT` (50) | habitat aggregation | cap on distinct habitat terms per taxon | fixed rule, held at the library default so `legacy` reproduces |
| tie weight | habitat aggregation | 1/k split across k tied categories | fixed rule |
| configuration name | step 1 | which rule set the call table row belongs to | one primary, one legacy, and one per grid level needing its own call table (`generic_none`, `iso_taxon`, `text_first`, `text_object`) |
| relation-fetch retries, wait, backoff, pause | step 1 | how a transient OmniCrobe failure is handled | fixed in the script: up to 3 further attempts on status 429, 500, 502, 503 or 504 or a dropped connection, waiting 2 s and doubling, 0.2 s after each fetch that succeeds. Any other status is not retried. The library's own default is one request and no retry |
| census the reproduction gates apply to | step 1 | which census the four gates were written for | declared by SHA-256 digest in `--gate-targets`, never defaulted; without a matching declaration the four are declined as not applicable (§3.3) |
| pool cap, subsample size, seed | step 0 | the background pool and its reproducible subsample | fixed constants (`MAX_TAXA`, `SUBSAMPLE`, `SEED`); step 3 records all three, read from the generator's source, on every run (§4.1) |
| search cutoff (bits) | step 2 | inclusion cutoff of the background search | calibrated at run time from the census accessions' own score distribution by a named percentile rule; written to a record file beside the output. `--threshold` overrides and the run records that it did |
| recall sample size, permutation seed | step 2 | the measured-recall sample | size defaults to `RECALL_N_DEFAULT`, `--recall-n` overrides; the seed (`RECALL_SEED`) is fixed. The sample is the first *n* draws of one seeded permutation, and an extension is checked against the earlier sample |
| `MIN_RECOVERY`, `MIN_RECALL`, `RECALL_ALPHA` | step 2 | acceptance bars, and the confidence level of the recall bound | fixed *before* the run as plan parameters; rules, not measurements |
| `LN_REF`, `CA_REF` | step 2 | the reference pair the step's own alignment is calibrated on (§1.3) | fixed to the pair the census calibrated on, so the step reads its hits at the same column; record order is load-bearing and the script says so |
| `POSITIVE_CONTROLS` | step 2 | five named taxa that must come out flagged, or the run publishes nothing | selected from public UniProt annotation and fixed before the search they check, so a PASS is a wiring check and **not evidence about any of these organisms' cofactors**. Each is reported with its own hit count against a bar of one. §4.2 |
| `MIN_FAMILY_LEN` | step 2 | below which a family member is a fragment | fixed rule; the column names carry the floor, so a table states its own convention |
| background exclusion mode | step 3 | how carriers are removed from the background arm | one registered primary; every other mode is descriptive, and the run records which was in force |
| `MIN_CELL`, `MIN_BG_K`, `MIN_FG_K` | estimators | below which a row's hit total a + b is not tested, and below which a phylum stratum is collapsed | fixed reporting rules |
| interval level | estimators | 95% | fixed: 1.96 standard errors, a literal constant in both estimators |
| `ANCHOR_LEVELS_RETIRED` | step 4 | whether the anchor levels, which reproduce rows of another census, are built | fixed to true: the levels stay defined, are not evaluated, and the run records the setting (§7) |
| plan bounds (`--plan-*`) | step 3 | bounds the plan set before this step ran, written from an earlier, approximate version of the same exclusion, so expectations rather than bounds set blind | supplied on the command line, never defaulted; disclosed, never gating. `docs/DATA.md` §4 |
| category selection (`--categories`) | steps 3 and 5 | which category each recorded outcome rule speaks for, which categories the grid counts cover, which phylum's shares are quoted | the categories are selected from the data, from this analysis's own earlier estimates; the phylum is determined, as the heaviest standardisation stratum, and the figure step refuses to assert the weight identities unless the run's own heaviest stratum is still that phylum. Supplied per run in a criteria file with no default anywhere, refused by name when absent, with a required note per row. `docs/DATA.md` §4 |
| thinnest-stratum list length | step 5 | how many drawn rows the thinnest-background-stratum list names | fixed listing convention, outcome-independent by construction: the list is extended through a tie at its last place and the smallest count among the rows it leaves out is exported beside it, so the list plus that bound is a statement about every drawn row. Written during figure construction, after the results were known |

---

## 10. External tools and services

### 10.1 Binaries

| tool | version | how it is used | invoked by this repository? |
|---|---|---|---|
| HMMER `hmmalign` | HMMER 3 | aligns sequences to the profile in `afa` format, which makes the model coordinate well defined | yes, with a configurable binary name; a missing binary raises |
| HMMER `hmmsearch` | HMMER 3 | the background carrier search, writing `--tblout` | yes: directly on Linux and macOS, through WSL on Windows (§4.2) |
| HMMER `hmmsearch`, the census search | HMMER 3 | produced the hit set the census begins from (§2) | no: run outside this repository; the census starts from its archived output |
| USEARCH | 11 | clusters the census sequences; the pipeline reads the `-uc` output | no: run outside this repository; the command is documented in `src/scriptLn.jl` (`USEARCH_CMD`) and in the run manifest |
| FastTree | — | produces the Newick trees the phylogenetic module parses | no: output consumed |
| structure prediction / structural alignment | — | produces the predicted structures the pocket path reads | no: output consumed |

The profile was built with HMMER 3.0 `hmmbuild` (its header says so) and is read by HMMER 3.x.
No HMMER version is pinned here and no run manifest records one: the `--tblout` tables step 2
writes carry the version in their own trailer. Record the version you used.

### 10.2 Web services

| service | endpoints | used for |
|---|---|---|
| UniProt REST | `rest.uniprot.org` — UniProtKB search and entries, UniParc entries and cross-references, reference-proteome streams | control-set sequence fetch, accession-to-TaxID resolution, the background pool, the targeted prefilter |
| OmniCrobe | `omnicrobe.migale.inrae.fr` — relation search and OBT term lookup | the habitat relations and the ontology climb |
| NCBI E-utilities | `eutils.ncbi.nlm.nih.gov` — `efetch` | the one control row whose sequence is fetched from the NCBI protein record |

Four things are handled in code:

- Requests name the client, and nothing about who runs it. No call sends an address, an
  account or an API key, and none sends NCBI's `tool` or `email` parameters. Two paths
  identify the software in their User-Agent header: OmniCrobe relation requests send
  `LnFinder/1.0 (research; serial)` (the default of `fetch_relations`), and the TaxID
  resolution's UniProt requests send `LnFinder/0.1` (`USER_AGENT` in
  `src/ModuleAccessionMapper.jl`). The other calls set none. NCBI rate-limits unidentified
  clients, so a caller doing many NCBI requests should space them.
- Where a failed request is asked again, the rule is written down: step 1's relation fetches
  (§3.3), the TaxID resolution (§3.1), and step 2's downloads, which make up to five attempts
  in all, waiting 5 s longer after each failure. In all three, a status that is not transient
  is an answer about the request, not about the service's load, and is not asked again.
- UniProt may gzip a response even when identity encoding is requested, so every fetch checks
  for the gzip magic bytes and decompresses defensively.
- These services are not versioned like a package. They change, and a re-run later will not
  reproduce an earlier run byte for byte. That is why the control-set fetch, the TaxID
  resolution and the background scan record the release header UniProt reported, why the
  habitat layer is replayed from a stored cache, and why release strings are stamped into the
  cache manifest. A reproduction should replay a cache (§3.3).

---

## 11. Tests

`julia --project=. test/runtests.jl` runs every suite in `test/`, and the runner and the
directory are kept as equal sets. The suite is offline by construction: no network, no
external binary, no input the repository does not carry. Every fixture is either built by the
suite that uses it, sometimes out of real public identifiers as the taxonomy dumps are, or a
response recorded verbatim from a public service or a local tool. The recorded fixtures are of
three kinds: OntoBiotope ontology paths with one OmniCrobe relation record; UniProt response
shapes, whose structure is real and whose every value is invented; and one
`hmmsearch --tblout` table for three public Swiss-Prot entries scored against the profile HMM
this repository ships. No fixture carries a value from any run of the analysis this software was
built for. Every value a fixture carries is either public (an ontology path, a Swiss-Prot
accession, a bit score recomputable from the profile HMM that ships here) or invented for the
case under test, where the expected result follows arithmetically from the fixture and the
suite says how. Where a suite pins a calibration, it pins one performed on the structures or
the table that suite builds.

The shipped data files the suite touches: the control set and the recorded validation
predictions (`dat/control_set/control_validation_predictions.csv`), checked only for presence
and for being spelled exactly as the code and the documents name them; the profile HMM, the
same check plus its header, which the classifier suite reads for the match-state count
(`LENG`) the alignment guard falls back to when a caller names no profile; the anchor targets
and the two source files they name, the same check plus the sources' bytes against the SHA-256
digests the targets declare; and the profile's provenance note, `dat/models/PROVENANCE.md`,
read as prose like `README.md`.

Some shipped pieces have no suite of their own: the census driver (§2), whose
environment-variable gates and never-overwrite-a-census guard are not unit-tested; step 0
(§4.1); `src/ModuleSequenceTools.jl`, which the census driver calls (§2);
`src/ModuleMotifSearch.jl`, which no test calls; and `scripts/verify_data.jl`, checked for
presence and never run. `README.md` names the same five, at more length.

One suite is structural: it resolves every `include` target in `src/`, `scripts/` and `test/`
against a directory listing instead of a file-existence check, which a case-insensitive
filesystem answers for the wrong spelling. It also refuses two files whose names differ only
by case, and any invocation that depends on a version manager. CI runs the whole suite on a
case-sensitive and a case-insensitive platform for the same reason.
