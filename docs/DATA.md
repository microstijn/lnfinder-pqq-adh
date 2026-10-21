# The data this repository does not carry

LnFinder's code is here; the files it reads are not. They are large, most of them are
derived, and several of them are results of the analysis the software was written for. They
are archived as one data set, and this document is the contract between that data set and
this code: **one data root, one layout underneath it, one name per input.**

Read this before trying to run anything past the test suite.
[`README.md`](../README.md) says which stages can be run at all; this file says what each
stage needs, where it expects to find it, and what to do when you do not have it.

---

## 1. The contract

Everything the analysis code reads that this repository does not ship lives under a single
directory — the **data root** — in a fixed layout.

| | |
|---|---|
| **The root** | the directory named by `$LNFINDER_DATA`. When that is not set, the root is `dat/` inside the clone, which is where the shipped inputs already live and which `.gitignore` already covers for everything else. |
| **The layout** | fixed, and declared in exactly one place in the code: [`src/ModuleDataInputs.jl`](../src/ModuleDataInputs.jl). The tree is in §3 below. |
| **Per-input overrides** | two, and only two: `$LNFINDER_AFA` and `$LNFINDER_TAXDUMP`. They exist because those two files are the ones a reader *makes* or *fetches* rather than unpacks, and both are large enough to be kept outside a project tree. Everything else resolves under the root and has no variable of its own. |
| **When something is missing** | the step stops before it does any work and prints a refusal naming the input, where the layout puts it and which variable moves it. It does not print a stack trace and it does not print a path off anybody's machine. That is for the inputs the step resolved from the root; a path you typed yourself is reported as the path you typed, because the layout is not what is wrong with it. |

So: unpack the data set at one place, point one variable at it, and every path the data set
carries is satisfied at once — every path but the criteria files, which are deliberately
passed on the command line for the reason §4 gives. What the data set does not carry, and
what each step does without it, is §5.

```
# everything under one root, anywhere you like
export LNFINDER_DATA=/scratch/lnfinder-data   # PowerShell: $env:LNFINDER_DATA = "C:\lnfinder-data"
julia scripts/verify_data.jl                  # check it before running anything
```

A relative value is made absolute when it is read, so one variable means one directory
whatever a step was started from.

Two more variables name where a run **writes**, and are not inputs:
`$LNFINDER_CENSUS_OUT` (the directory a census run writes its table and manifest into) and
`$LNFINDER_FIGDIR` (the figure directory of habitat steps 4 and 5). A third,
`$LNFINDER_RELATION_SEED`, primes step 1's relation cache from directories of
already-fetched response bodies.

That is six `LNFINDER_*` variables — the root, the two per-input overrides and those
three. **Two more switches are read by the census driver and carry no prefix**:
`RUN_HABITAT` (`false` turns off its habitat layer, which is its only per-taxon network
route besides accession resolution) and `EXPECTED_CENTROIDS` (a whole number it checks its
centroid count against, so a clustering that silently changed size stops the run). Eight
names in all, every one of them declared in `src/ModuleDataInputs.jl`. Two tests fail the
build if a script reads a ninth: one scans for the `LNFINDER_` prefix, and one scans for
the *access* — every `ENV["…"]`, `get(ENV, "…", …)` and `haskey(ENV, "…")` in `src/` and
`scripts/`, whatever the name is — because a scan for the prefix is exactly what could not
see the two switches.

---

## 2. How to obtain it, and when

**The data set is not public yet, and there is no DOI to give you.** It contains results —
see §6 — so it is released with the accompanying publication and not before. When it is
released, this line is replaced by its archive identifier:

> **Data set DOI: `<to be assigned at release — this placeholder is filled in when the
> archive is published, and the same identifier is added to CITATION.cff>`**

Until then, nothing here promises an archive exists. What this document does promise is
that the *layout* the archive will have is the layout the code already reads, so this
section is the only part of this file that changes on the day it appears.

### Verifying what you got

[`docs/data_manifest.sha256`](data_manifest.sha256) carries the sha256 digest of the
archived files, recorded when they were assembled. It is plain `sha256sum` format, so there
are two independent ways to check a data root:

```
julia scripts/verify_data.jl                       # also checks the contract itself
julia scripts/verify_data.jl --root <path> --manifest <path>

cd <data root> && sha256sum -c .../docs/data_manifest.sha256    # shasum -a 256 -c on macOS
```

A data root without the optional files makes `sha256sum -c` report each of them as FAILED;
`--ignore-missing` skips them, and with them any missing required file, which
`verify_data.jl` does not.

`verify_data.jl` needs nothing but Julia — not even an instantiated environment — so it can
be run while you are still setting up. It hashes every file the manifest lists, names what
is missing, names what differs, additionally reports any declared input the root does not
satisfy, and **exits non-zero** if anything is wrong, so it can stand in front of a long run
in a script.

The archive carries its own `MANIFEST.sha256`, which `verify_data.jl` prefers when it finds
one at the root. **That copy is the complete one.** On its own, though, it would establish
only that the archive agrees with itself — bytes and record travel together — so whenever
some other record is used, the repository's copy is **also** consulted for the paths it
covers, and a disagreement between the two is a failure in its own right. That is the
reason this file is shipped at all: it is the second, independent location.

The repository's copy leaves one thing out on purpose:

- **`structures/structure_accessions.txt`**, because it is one identifier per member of the
  census, so a list of them is the census's membership, and an integrity record is not the
  place to publish one. The twelve experimental structures further down **do** have digests,
  and that is a deliberate distinction rather than an inconsistency: twelve published
  crystal entries, which the accompanying publication names in any case, are not a
  census-wide identifier list.

Everything else the archive carries has a digest here, `criteria/category_selection.csv`
included: its digest is that of the copy the archive is assembled from, so the repository's
record and the archive's own name the same bytes for it. The habitat relation cache's stored
bodies and the ontology climb beside them are not in the archive at all (§5), so an integrity
record of the archive holds no digest of them; the archive does carry that cache's own
`manifest.csv`, and its digest is in the repository copy.

**One consequence worth stating, because it is easy to get the order wrong.** A step-4 run
records `anchor_targets_sha256` in its provenance file, and step 5 records the digest of
that provenance file in its own. Both are digests of whichever copy the run read. So the
copy of `dat/criteria/anchor_targets.csv` this repository tracks has to be the copy the
published step-4 run read: any edit to it, a reworded note included, is followed by a
step-4 and a step-5 run before their provenance is published. Editing afterwards would leave
a published provenance pointing at bytes nobody can obtain. The two files beside it are
pinned the other way round: the targets file declares their digests, and step 4 refuses to
publish when a file it finds beside the targets is not those bytes (§4).

The same holds for `category_selection.csv`, and twice over: step 3 records its digest as
`selection_sha256` in the arm summary and step 5 records it as `input_categories_sha256` in
the figure provenance, so the deposited copy has to be the copy both of those runs read.

---

## 3. The layout

```
<data root>/                          $LNFINDER_DATA, else dat/ in the clone
  MANIFEST.sha256                     the archive's own checksums
  habitat/                            the habitat comparison, steps 0 to 5
    census_habitat.csv
    prokaryote_taxa.csv
    prokaryote_background.csv
    obt_path_cache.csv                NOT in the archive — written by steps 0 and 1
    habitat_taxon_calls.csv
    background_pqq_scan.csv
    background_scan/                  the scan's own records
      recall.txt  threshold.txt  acceptance.tsv
      uniprot_release.txt  uniprot_release_IPR017512.txt
      refs.fasta  recall_proteomes.fasta.manifest  calibration_scores.csv
    omnicrobe_cache/                  the habitat relation cache — see §5
      manifest.csv                    what the cache holds; the archive carries this one
      relations/<taxid>.json          NOT in the archive — one stored response body per taxon
    habitat_enrichment.csv            written by step 3, not archived
    habitat_arm_summary.csv           written by step 3, not archived
    habitat_enrichment_anchor.csv     written by step 3, not archived; no stage reads it
  figure/                             the environment figure's three tables
    prevalence.csv  prevalence_phylum_average.csv  thermal_options.csv
  criteria/                           what a run is held to; passed on the command line
    category_selection.csv
  census/                             the census route
    firstTry.fa
    clusters.uc
    clusters_50.uc
    final_deduplicated_data.csv
  determinant/
    pocket_labels_real.csv
  structures/
    structure_accessions.txt          identifiers; the structures are NOT archived here
  derived/                            NOT in the archive — you make or fetch these
    aligned_to_model.afa              $LNFINDER_AFA moves it
    taxdump/                          $LNFINDER_TAXDUMP moves it
      nodes.dmp  names.dmp  merged.dmp  delnodes.dmp
    structures/
      crystals/<pdbid>.pdb            lower-case four-character identifier
```

In the clone, `dat/criteria/` also holds three files this repository tracks —
`anchor_targets.csv` and the two files its cells were transcribed from (§4). They resolve from
the repository whatever `$LNFINDER_DATA` says, so an archive unpacked at `dat/` puts
`category_selection.csv` beside them and a root elsewhere leaves them where they are. The
archive carries no file of those three names, so unpacking it in place overwrites none of them.

Under `habitat/omnicrobe_cache/` the archive carries `manifest.csv` and nothing else: the
stored response bodies under `relations/` are not redistributed, and neither is
`obt_path_cache.csv` beside it (§5). A reader who holds those caches puts them where the tree
above does; a run without them fetches, which §5 and §7 spell out.

Total, as assembled: **about 88 MB** unpacked, of which `census/firstTry.fa` and
`habitat/habitat_taxon_calls.csv` are the bulk. Plan disk for that order of size — the
per-file figures are the archive's own, and `du` on the unpacked root reports them once you
have it.

---

## 4. Every file: what it is, who reads it

The `name` column is the name the code uses. It is the only handle you need: the layout can
be read off `src/ModuleDataInputs.jl`, and every refusal quotes both.

### Shipped with the code — not under the data root

| name | path | size | what it is |
|---|---|---|---|
| `profile_hmm` | `dat/models/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm` | 0.26 MB | the profile the cofactor call is calibrated against. Third-party and unlicensed upstream — see `dat/models/PROVENANCE.md`. Read by the census driver, habitat step 2 and control-set validation |
| `control_set` | `dat/control_set/known_cofactor_mdhs.tsv` | 0.01 MB | the curated set of enzymes whose cofactor is known from the literature, one row per published enzyme with its citation. Read by control-set validation |
| `anchor_targets` | `dat/criteria/anchor_targets.csv` | 1.1 MB | the transcribed anchor cells step 4 checks its grid against — the registered specification and each of the fourteen registered levels, on three habitat and three chemistry categories — with a tolerance and a note per cell, the names of the registered levels, and the name and sha256 of the two files the cells were transcribed from. `scripts/habitat_sensitivity.jl --anchor-targets dat/criteria/anchor_targets.csv`. It carries result values (§6) |
| `anchor_source_primary` | `dat/criteria/impl_4.csv` | small | the independent implementation the registered specification's cells are transcribed from: an `anchor,unit,category,quantity,value` table computed from the call table and the scan table by code written from the specification without reading this repository's grid script. Step 4 checks it against the digest the targets file declares, and refuses to publish when the two differ |
| `anchor_source_variants` | `dat/criteria/impl_w1.csv` | 0.03 MB | the same, for the fourteen registered levels |
| `figure_tables` | `dat/figure/prevalence.csv`, `dat/figure/prevalence_phylum_average.csv`, `dat/figure/thermal_options.csv` | 0.004 MB | the three tables the environment figure draws from. DERIVED, not primary: `scripts/prevalence.py` and `scripts/prevalence_phylum_average.py` rebuild them byte for byte from `habitat/habitat_taxon_calls.csv` and `habitat/background_pqq_scan.csv` in the data set. They are tracked so that `scripts/fold_change_figure.jl` renders the published figure from a bare clone, with no data set present |

**The three files under `dat/criteria/` are stored byte for byte.** `.gitattributes` marks them
`-text`, because the targets file pins its two sources by sha256 and step 4 digests the bytes
on disk: a checkout that rewrote their line endings would fail step 4's `target sources` gate
on files nobody edited. Step 4 reports a declared source it cannot find beside the targets
file without refusing, so `test/test_data_inputs.jl` is what checks that both sources are
present, beside the targets, and are the bytes the targets name.

**What the anchor targets check, and what they do not.** Every target carries a note saying
what the specification behind it did and did not determine. The recorded analysis plan does
not determine every rule the specification states; for three of them the specification
takes the rule from this repository's implementation: the thresholds of the thin-stratum collapse, which every level
shares; the rule definitions of four levels (`T:alphabetical`, `T:unresolved`, `C:share`,
`C:winner_fractional`); and one edge case of a fifth (`I:iso_taxon`). Where that is so, a
target checks that the grid is wired to its data and computes the rule it states; it does not
independently check the choice of rule.

**One more file is shipped, and it is a result rather than an input.**
`dat/control_set/control_validation_predictions.csv` (small) is the classifier's call on each
control of the shipped control set — the `control_predictions.csv` that control-set
validation writes, with the columns `original_header,accession,classification,active_site`
and one row per control. No code reads it, so it has no name in the registry; it is tracked
as the record of that validation (§6).

### `habitat/` — the habitat comparison

| name | path | size | what it is, and which stage reads it |
|---|---|---|---|
| `census_habitat` | `habitat/census_habitat.csv` | 8.5 MB | the census table carrying the habitat layer, one row per sequence. **Steps 1, 2, 3 and 4 all read it, and nothing in this repository writes it.** Its identifier, cofactor-call, active-site, `pos2_aa` and taxonomy columns (`TaxID` to `species`) are those of `census/final_deduplicated_data.csv`, row for row. Its seventeen habitat-layer columns — `report_group`, `isolation`, `chemistry_share` and the `habitat_*`, `published_*` and `pillar_*` families — are a per-row annotation joined to it on `id` by a step that is not part of this release; only step 1 reads them, for its replication report, and the habitat comparison takes every call from `habitat_calls`. Those columns are computed from relations the microbe–habitat service returns, and the habitat names in them are the habitat ontology's own term labels (§8) |
| `pool_taxa` | `habitat/prokaryote_taxa.csv` | 0.34 MB | the reference-proteome taxon pool the comparison arm is drawn from. Steps 1 and 2 read it, and step 3 records its digest and taxon count in its provenance when it is there; step 0 writes it, from a live query whose answer moves, which is why the archived copy matters |
| `pool_background` | `habitat/prokaryote_background.csv` | 0.62 MB | the habitat annotation of that pool. Step 1 reads it, and step 3 records its digest and taxon count in its provenance when it is there; step 0 writes it |
| `obt_path_cache` | `habitat/obt_path_cache.csv` | — | the habitat-ontology climb, cached, so the ontology is walked once. Steps 0 and 1 read it, and both write it back. **Not in the archive** (§5): it is a verbatim extract of a third-party ontology, keyed by ontology term id, and it is not redistributed here |
| `relation_cache` | `habitat/omnicrobe_cache` | — | the habitat relation cache: the directory the stored response bodies and their manifest live in. Step 1. The archive carries the manifest and not the bodies (§5), so **this is the item whose absence makes a step-1 run a new measurement rather than a replay** — the service behind it is text-mined and unversioned, so a fresh fetch is a different data set, not the same one |
| `relation_bodies` | `habitat/omnicrobe_cache/relations` | — | the stored bodies themselves, one `<taxid>.json` per taxon. **Not in the archive** (§5). Step 1 **requires them by name** unless it is given seed directories: a cache miss here is a live fetch, so requiring the directory that holds them instead would let a cache that exists with nothing underneath it pass every check and then go to the network on a run that reports itself as replaying from disk |
| `relation_manifest` | `habitat/omnicrobe_cache/manifest.csv` | 0.69 MB | what that cache holds, where each body came from and when: one row per taxon with a relation count, a source, a fetch stamp, and the service and ontology releases the run records. It records the fetch and carries none of the responses, and it is in the archive. Steps 1 and 3. Step 1 **rewrites it in place** at the end of a run, wherever `--out` sends the call table |
| `habitat_calls` | `habitat/habitat_taxon_calls.csv` | 24.4 MB | one habitat call per taxon per rule variant. Steps 3, 4 and 5 read it; step 1 writes it. **It is the archive's entry point to steps 3, 4 and 5**, which run from it with no network; step 1, which writes it, cannot be replayed from the archive (§5). It is computed from relations the microbe–habitat service returns, and its habitat names are the habitat ontology's own term labels (§8) |
| `background_scan_dir` | `habitat/background_scan` | 0.14 MB | the scan's own records, listed under §3. Step 3 reads four of them; **without them it records the evidence as unknown and withholds the background licence rather than refusing** — unknown withholds the licence instead of erroring or passing. The rest record which sequences and which release the threshold was set on, and exist so that step 2 can be audited without re-running it |
| `background_scan_table` | `habitat/background_pqq_scan.csv` | 0.26 MB | which taxa of the comparison arm carry the enzyme at the calibrated threshold. Steps 3 and 4 read it; step 2 writes it, which costs a network crawl, HMMER and several hundred megabytes of downloads. Its frame is the pool less the census members the scan's own run saw; step 3 compares that with the frame of the census it reads, names every taxon that is in one and not the other, and withholds the complete-removal statement when they differ |
| `scan_recall` | `habitat/background_scan/recall.txt` | small | the measured recall of the scan's prefilter with the counts it was measured on. Step 3 |
| `scan_threshold` | `habitat/background_scan/threshold.txt` | small | the calibrated score threshold and how it was set. Step 3 |
| `scan_acceptance` | `habitat/background_scan/acceptance.tsv` | small | the scan's acceptance rows, one per check it had to pass. Step 3 |
| `uniprot_release` | `habitat/background_scan/uniprot_release.txt` | small | the sequence-database release the comparison arm was built from. Step 3 |
| `enrichment_table` | `habitat/habitat_enrichment.csv` | — | **written by step 3, not archived.** The estimates, one row per contrast. Step 5 reads it. Step 4 does **not**: it includes step 3 and recomputes |
| `arm_summary` | `habitat/habitat_arm_summary.csv` | — | **written by step 3, not archived.** What each arm was made of, and every provenance key the panel echoes. Step 5 |

### `criteria/` — what a run is held to

The criteria are passed on the command line, deliberately: their values are measurements of
this data set, or selections made from it, so they are not constants in the code. A step
given none of them refuses to run rather than quietly skipping the check or falling back
to a value the software picked. The archive carries one of them; `anchor_targets.csv` is
tracked in this repository (above).

| name | path | size | what it is |
|---|---|---|---|
| `category_selection` | `criteria/category_selection.csv` | small | **which categories the analysis is about.** `--categories <path>` to steps 3 and 5 |

**Two criteria options take no file from this data set, and say so in their output.**

- **Step 1's `--gate-targets`** names a file of expected values for its measured replication
  gates. No such file is carried for this census, because the expected values on record are
  measurements of another one, so the call table here is built with `--no-gates`, and step
  1's report states gate by gate that those gates were not checked (§7). Step 1's four
  reproduction gates — the ones that compare a run with the call table and the habitat
  columns of the census they were written for — are evaluated only when a `--gate-targets`
  file declares that census by its sha256 and the run reads a census with that digest.
  Under `--no-gates`, with no such declaration, or on another census, each of the four is
  reported as not applicable with its reason, rather than checked against a census it was
  not written for.
- **Step 4's `--anchor-reference`** names an earlier run's rows for the all-legacy anchor
  level. That level reproduces rows of a census other than the one this data set carries,
  so it is not evaluated (`ANCHOR_LEVELS_RETIRED` in `scripts/habitat_sensitivity.jl`), and
  while that holds step 4 needs `--anchor-targets` only: `--anchor-reference` is not
  required, a path given anyway is recorded and not opened, `anchor_targets.csv` declares
  `reference_rows,0`, and the provenance records `anchor_reference_status` as not applicable
  with no row compared.

#### `category_selection.csv`, and why it is not in the code

Three of the recorded outcome rules speak for one habitat category each — one whose
depletion a robustness tier decides, one whose genus clause and category-drop are decided,
one whose comparator-to-estimate sharing ratio is decided — the figure step tabulates one
or more categories' intervals across every registered specification, and both steps quote
the shares of one phylum beside the standardisation. **Which** categories those are was
chosen from this data set's own estimates, with the estimates in hand — after those estimates
existed, and not fixed in advance of them. The habitat classification scheme the names are
drawn from is a mapping of an external ontology and predates every estimate
(`src/ModuleEnvironmentTools.jl`); selecting four names out of it does not, and this release
does not present the second as though it were the first.

**The phylum is the one entry that is determined rather than chosen**, and the file says so in
its own note: it is the heaviest standardisation stratum of the analysed arms, so a run on
other data has another one. It is named in the file because the summary's share keys are named
after it and a reader has to be able to derive that name — and the figure step publishes the
identities tying those two shares to the standardisation weights **only** while the run's own
heaviest stratum is still that phylum, and refuses to assert them otherwise. It was written
down after the shares had been computed, so it is not pre-specified either; it is simply not a
degree of freedom.

So the selection is a per-run input with **no default anywhere in the code**, and a run
that does not supply it refuses, naming the option. The file is a four-column table —
`kind,role,value,note` — with one `select` row per role, one `expect,selections,<n>` row
stating its own size, and a required `note` on every row saying where that selection came
from. Roles left unfilled, roles the step does not know, a role filled twice, a missing
note, a count that disagrees with the rows, or two selections whose names reduce to one
summary-key prefix are each a refusal rather than a quieter success.

The names also reach the outputs. The summary's rule-keyed families and the phylum-share
keys are named after the selected categories, and those names are **derived** from the
file by one rule — lower case, non-alphanumeric characters removed — so a run given other
categories emits keys for those categories instead of filing other numbers under names the
software chose. Each run records the criteria file's digest and every selected name beside
its numbers (`selection_*` in step 3's summary, `input_categories_sha256` and
`selection_*` in step 5's provenance), so two runs can be told apart by what they were
asked about before a single estimate is compared.

**Two criteria of step 3 are not files, and are not in the archive.**
`scripts/habitat_enrichment.jl` also takes `--plan-bg-size-below` and
`--plan-bg-phylum-pct-max`: two upper bounds written down in the analysis plan
before the step ran and supplied on the command line, which is why the script gives them
no default — a bound it invented would be the script's opinion, and a bound it defaulted
to would be indistinguishable from one that was checked. *Before the step ran* is the whole
of the claim: they were written from an earlier, approximate version of the same exclusion,
so they are the plan's expectation for this run and not bounds set blind — which is also why
nothing gates on them and a violated one is reported rather than fatal. **They are not part of this
release**, because the plan they were written in is not part of it either (see *"The
registered specification" and "the analysis plan"* in `README.md`). The consequence for a
reader is exact and small: a run without them is identical on every computed quantity, and
its summary carries twelve fewer provenance keys — the `criterion_bg_size_*` and
`criterion_bg_<phylum>_pct_*` families, which record whether each bound was met (the second
family is named after the phylum the criteria file selects). Nothing else in the run moves.

### `census/` — the census route

| name | path | size | what it is |
|---|---|---|---|
| `hits_fasta` | `census/firstTry.fa` | 33.0 MB | the unaligned sequences the census was built from, as retrieved. Read by `src/scriptLn.jl`. It is archived, so the census is rebuilt from these exact sequences |
| `uc_clusters` | `census/clusters.uc` | 4.1 MB | the clustering of those sequences that defines the census centroids. Read by `src/scriptLn.jl`. Archived because regenerating it needs USEARCH 11, a proprietary binary that is not and cannot be redistributed. The clustering command is documented in `src/scriptLn.jl` (`USEARCH_CMD`) |
| `uc_clusters_50` | `census/clusters_50.uc` | 4.3 MB | the coarser clustering used as the sub-family unit of the phylogenetic control. **No stage in this release reads it**; it is archived because it cannot be rebuilt without that same proprietary binary |
| `census_table` | `census/final_deduplicated_data.csv` | 7.6 MB | the census as `src/scriptLn.jl` writes it: one row per centroid with its cofactor call, active-site residues, TaxID and lineage, and the route, entry status and sequence-database release each TaxID was resolved through. It was written with `RUN_HABITAT=false`, so its `report_group` reads `_Not Mapped_` on every row. `habitat/census_habitat.csv` carries its identifiers, calls and taxonomy row for row. **No stage in this release reads it**: it is here so that the habitat input can be checked against the table it came from |

### `determinant/` and `structures/`

| name | path | size | what it is |
|---|---|---|---|
| `pocket_labels` | `determinant/pocket_labels_real.csv` | 2.8 MB | the per-accession metal label set of the whole census, as the structure route produced it. **No driver in this release reads it** — the determinant and pocket modules ship as a library with unit tests and no command-line driver — but it is archived because rebuilding it means every predicted structure plus a full pocket export |
| `structure_accessions` | `structures/structure_accessions.txt` | not stated (see below) | one identifier per structure the analysis read, in sections by archive. **The structures themselves are not archived here** — see §5 |

**Why one archived file above carries no size.** A dash in the size column marks an entry whose
bytes the archive does not carry — a file it leaves out, or a directory it carries only as the
place the manifest sits. The one file it does carry without a figure is
`structures/structure_accessions.txt`, and that is the distinction §2 already draws for
digests, applied to sizes: a file that is one short identifier per line is very nearly its own
line count, and this one is one identifier per member of the census. A size would therefore
say how long a list this repository is not publishing yet — so it is not stated here. Nothing
a reader has to do needs it: the sizes in this document are here to help you plan a download,
which §3's total does, and verification is by digest — `verify_data.jl` and `sha256sum -c`
compare hashes, never sizes. The figure is readable off the archive itself, which is released
with the publication. (A table is a much weaker proxy for its own length, because its rows vary in
width; the sizes of the tables above are stated, and §6 says which of them are results.)

---

## 5. What is not in the archive, and how to get it

Five things the analysis read are not archived. Four of them are large and come from public
sources, and **three of the five have digests** in `docs/data_manifest.sha256` under
*optional*, so whatever you make or fetch can be compared against what the analysis
actually read. Two have none. Step 2's bulk download cache has none because it is hundreds of
megabytes of third-party response bodies that no later step reads, and the records that
**are** read from that step are archived in their own right. The fifth — the habitat relation
cache and the ontology climb beside it — has none for a different reason: the other four come
from archives you can ask for a release or an identifier, while the service behind this one is
re-mined and unversioned, so a re-query answers with a different data set rather than with
these bytes. It has its own subsection at the end.

### `model_alignment` — `derived/aligned_to_model.afa` (348 MB): regenerate it

One HMMER run over a file that **is** archived, against the profile that ships here:

```
mkdir -p "$LNFINDER_DATA/derived"
hmmalign --trim --outformat afa \
    dat/models/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm \
    "$LNFINDER_DATA/census/firstTry.fa" > "$LNFINDER_DATA/derived/aligned_to_model.afa"
```

`$LNFINDER_AFA` points at it if you would rather keep it elsewhere. Compare your result
against the recorded digest before using it: HMMER versions differ.

### `taxdump` — `derived/taxdump/` (449 MB unpacked): re-fetch it

NCBI Taxonomy, from the dated archive at
<https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/>. Unpack it so that
`nodes.dmp`, `names.dmp`, `merged.dmp` and `delnodes.dmp` sit directly under
`derived/taxdump/`; the loader opens those four and ignores the rest. `$LNFINDER_TAXDUMP`
points at it elsewhere — taxonomy dumps are commonly shared between projects.

**This one will not match its recorded digest, and that is expected.** The archive holds
snapshots at month boundaries, and the dump the analysis used was taken between two of
them. Take the nearest snapshot at or before **2025-02-05**, and treat the digest as a
record of how far your copy has drifted rather than as a test to pass.

### `structures` — `derived/structures/`: re-fetch by identifier

Not archived, and deliberately: every one of these files is a public record with a stable
identifier, and `structures/structure_accessions.txt` carries every identifier the analysis
used. Its sections say where each one comes from:

| section | archive | URL pattern |
|---|---|---|
| `afdb` | predicted structures | `https://alphafold.ebi.ac.uk/files/AF-<ACCESSION>-F1-model_v4.pdb` |
| `afdb_subsample` | a subset of the same set, listed separately because the selection is recorded nowhere else | as above |
| `rcsb` | experimental structures | `https://files.rcsb.org/download/<ID>.pdb` |

The predicted set is version **v4** of that archive (`ALPHAFOLD MONOMER V2.0`, headers dated
`01-JUN-22`). Later releases of the same archive change the files, so a re-fetch today may
not be byte-identical and some identifiers may no longer resolve; the digests of the twelve
experimental structures are recorded, so those at least can be checked. Save them as
`derived/structures/crystals/<lower-case id>.pdb` if you want `verify_data.jl` to check them
for you.

### `habitat/background_scan/` bulk downloads (about 640 MB)

Step 2's download cache and its HMMER working files are not archived. The *records* step 3
reads are (§4), and the caches themselves are re-fetchable from the sequence database at the
release pinned in `uniprot_release.txt`.

### `relation_cache`, `relation_bodies` and `obt_path_cache`: not redistributed

`habitat/omnicrobe_cache/relations/<taxid>.json` is one stored response body per taxon from
the microbe–habitat service the code queries, and `habitat/obt_path_cache.csv` is the climb
over the habitat ontology those responses are read against — an `obtid,path` table keyed by
ontology term id, carrying no taxon and no census value. Both are third-party content as
retrieved, and they are **not redistributed here**. What the archive carries is the habitat
layer's derived output instead: `habitat/habitat_taxon_calls.csv`, one call per taxon per
rule variant, and `habitat/omnicrobe_cache/manifest.csv`, the record of which taxon was
fetched, from where and when, with the service and ontology releases the run declares.

**What that costs, plainly.** Habitat steps 3, 4 and 5 run from the archive with no network
(§7), because they read the call table. Step 1, which writes the call table, does not: without
the cache it fetches, and the service behind it is text-mined, unversioned and periodically
re-mined, so what comes back is a different data set rather than these bytes again. A step-1
run made that way is a **new measurement, not a reproduction** of the call table in the
archive. The call table is therefore the reproducible starting point, not an intermediate a
reader regenerates.

**What step 1 does when the caches are not there.** It stops before it does any work and
refuses, naming the input, where the layout puts it and which variable moves it, and saying
that without it the step fetches instead of replaying. The refusal also names the one route
that is not a refusal: prime the cache from directories of already-fetched `<taxid>.json`
bodies, given as bare arguments to step 1 or as a semicolon-separated list in
`$LNFINDER_RELATION_SEED`. A run given seed directories is left alone, because bootstrapping
the cache is the one case where not having it yet is the point — and seeding is not a replay
either: any taxon the seed does not cover is fetched live, and any ontology term the OBT cache
does not hold is looked up live.

Cite the service and the ontology in their own right (§8). The tables the archive carries are
computed from them.

---

## 6. Which of these are results

Stated plainly, because a reader who re-runs and gets different numbers should know which
files were the answer and which were the question. **These files are outputs of the
unpublished analysis**, which is why the data set is released with the accompanying
publication and not before:

| file | what it discloses |
|---|---|
| `census/final_deduplicated_data.csv` | the census: one row per centroid with its cofactor call, active-site residues and resolved taxonomy |
| `habitat/census_habitat.csv` | the census table carrying the habitat layer as published, one row per sequence |
| `habitat/habitat_taxon_calls.csv` | per-taxon habitat calls with the per-taxon carrier counts the comparison is computed from |
| `habitat/background_pqq_scan.csv` | which background taxa carry the enzyme at the calibrated threshold |
| `habitat/background_scan/` records | the calibrated threshold and its measured recall — small files, but fitted quantities |
| `determinant/pocket_labels_real.csv` | the per-accession metal label set of the whole census |
| `criteria/category_selection.csv` | the categories the analysis's own estimates led to |

**This repository carries result values too.** `dat/criteria/anchor_targets.csv` and the two
files beside it hold the registered specification's and the registered levels' estimates on
six categories, computed by implementations independent of this code; they ship because
step 4 checks its grid against them, and a check a reader cannot run is not a check.
`dat/control_set/control_validation_predictions.csv` holds the classifier's call on every
control of the shipped control set.

The rest are not results: `habitat/prokaryote_taxa.csv` is a public taxon list,
`habitat/omnicrobe_cache/manifest.csv` is bookkeeping — which taxon was fetched, from where
and when — `census/firstTry.fa` is public sequence as retrieved, and `census/clusters.uc` and
`census/clusters_50.uc` are clusterings of those sequences that carry no cofactor claim. The
two caches the archive does not carry, `habitat/omnicrobe_cache/relations` and
`habitat/obt_path_cache.csv`, are third-party content as retrieved rather than results of
this analysis (§5).

---

## 7. Which stage needs which file

| stage | reads | writes |
|---|---|---|
| census — `src/scriptLn.jl` | `model_alignment`, `hits_fasta`, `taxdump`, `uc_clusters`, `profile_hmm`; and a sequence database over the network, always | `$LNFINDER_CENSUS_OUT/census_rebuild_<date>/` |
| habitat step 0 — `scripts/prokaryote_background.jl` | the network; on a rerun, its own three outputs (it is resumable) | `pool_taxa`, `pool_background`, `obt_path_cache`, of which the archive carries the first two |
| habitat step 1 — `scripts/habitat_reaggregate.jl` | `census_habitat` (or the file `--census` names), `pool_taxa`, `pool_background`, `relation_bodies`, `relation_manifest`, `obt_path_cache` — **two of those three cache inputs, `relation_bodies` and `obt_path_cache`, are not in the archive (§5), so a run that has only the archive fetches rather than replays** — and `--no-gates` on the command line (§4); a reader's own `--gate-targets` file makes the four reproduction gates only if it carries one `census,,<sha256>,,<note>` row matching the digest of the census read | `habitat_calls` (or `--out`), and `relation_manifest` and `obt_path_cache`, **always in place under the data root** |
| habitat step 2 — `scripts/background_pqq_scan.jl` | `profile_hmm`, `census_habitat`, `pool_taxa`; plus HMMER and the network | `background_scan_table`, `background_scan_dir` |
| habitat step 3 — `scripts/habitat_enrichment.jl` | `habitat_calls`, `background_scan_table`, `census_habitat`, `scan_recall`, `scan_threshold`, `scan_acceptance`, `uniprot_release`, `relation_manifest`; `pool_taxa` and `pool_background` for its provenance, when they are there; **`category_selection` on the command line as `--categories`, required, no default**; and two plan bounds on the command line, which this release does not carry (§4) | `enrichment_table`, `arm_summary`, and the legacy configuration's rows beside `--out` (`habitat_enrichment_anchor.csv`), which no stage reads |
| habitat step 4 — `scripts/habitat_sensitivity.jl` | `habitat_calls`, `background_scan_table`, `census_habitat`, and `anchor_targets` on the command line, with `anchor_source_primary` and `anchor_source_variants` read from beside it; no `--anchor-reference` (§4). **Not the step-3 tables**: it includes step 3 and recomputes, so it runs whether or not they are there | the grid, its robustness table, its panel and its provenance, under `--figdir`, else `$LNFINDER_FIGDIR`, else `figs/` in the clone |
| habitat step 5 — `scripts/habitat_figure.jl` | `enrichment_table`, `arm_summary`, `habitat_calls`, step 4's outputs from the figure directory, and **`category_selection` as `--categories`** — the same file step 3 was given, so the two steps cannot disagree about which categories the run is about | the panel and its five companion files, under `--figdir`, else `$LNFINDER_FIGDIR` |
| control-set validation | `control_set`, `profile_hmm`; plus HMMER `hmmalign` and the network | its own outputs, where you tell it |

**Step 1 writes into the data root, so give it `--out`.** Its call table goes to
`habitat/habitat_taxon_calls.csv` under the root unless `--out <path>` says otherwise — and
that file is in the archive, its digest is in `docs/data_manifest.sha256`, and it is where
steps 3, 4 and 5 start. A run with no `--out` replaces it with the output of a fresh fetch,
which is a different data set (§5), and no re-run brings the archived bytes back. Send it
somewhere outside the data root:

```
julia --project=. scripts/habitat_reaggregate.jl \
    --no-gates \
    --out <a directory of your own>/habitat_taxon_calls.csv
```

`--out` moves the call table and nothing else. The relation cache's `manifest.csv` is
**rewritten in place under the data root either way**, because it is the cache's own record
and the cache is where it lives — and that file is in the archive, with a digest. The OBT
cache is written back under the data root too, on every run. So run step 1 against a copy of
the data root, or keep your own copy of the two files it writes in place: `--out` protects the
call table and nothing else.

`--census <path>` runs step 1 on a census other than the deposited one. The four
reproduction gates follow the same rule there (§4): they are declined as not applicable,
and the run says so, unless a `--gate-targets` file declares that census by its sha256.

### What this buys, honestly

- **With the archive in place, habitat steps 3, 4 and 5 run with no network and no
  external binary.** They read `habitat_calls`, and nothing in them reaches the service the
  habitat layer is computed from.
- **Step 1 is not one of those three.** It has **two** network paths, not one — the relation
  fetch and the ontology climb — and the archive answers neither, because the stored bodies
  and `obt_path_cache.csv` are not redistributed (§5). Step 1 requires `relation_bodies`,
  `relation_manifest` and `obt_path_cache` by name unless it is given seed directories, so a
  missing cache is a refusal rather than a quiet crawl, and a run that fetches or seeds is a
  new measurement rather than a replay.
- **Steps 0 and 2 are not offline, and do not need to be**: their outputs are archived, but
  for the ontology cache step 0 writes (§5). Step 2 additionally needs HMMER, and on Windows
  WSL.
- **The census route is not offline even with the archive.** `src/scriptLn.jl` resolves
  accessions against a live service before any other switch applies, so a re-run will not
  reproduce the census byte for byte.
- **Nothing here makes the upstream data collection re-runnable.** The database search
  behind `census/firstTry.fa`, the proteome pull behind step 0 and the pull behind step 2
  all went to live services whose contents have moved since.

---

## 8. Third-party content, and what to cite

| item | origin | terms |
|---|---|---|
| `habitat/habitat_taxon_calls.csv`, the habitat-layer columns of `habitat/census_habitat.csv` and of `habitat/prokaryote_background.csv`, and `habitat/omnicrobe_cache/manifest.csv` | this project's own tables, computed from relations returned by the microbe–habitat service the code queries, over the habitat ontology those responses are climbed against; the habitat names in them are that ontology's term labels. The manifest is the fetch record for the cache — a relation count per taxon and none of the retrieved content | the tables are the authors'. The service and the ontology are **cited, not redistributed** — cite each in its own right |
| `habitat/omnicrobe_cache/relations`, `habitat/obt_path_cache.csv` | response bodies from that service, as retrieved, and the climb over that ontology. **Not redistributed here** (§5) | the service's own terms and the ontology's own licence; cite each in its own right |
| `census/firstTry.fa`, `habitat/background_scan/refs.fasta` | public sequence records, the second being the two calibration reference records as retrieved | the sequence database's own licence (CC BY 4.0 at the time of writing); cite the release |
| `census/clusters.uc`, `census/clusters_50.uc` | output of USEARCH 11 | the output is the authors'. **The binary is not, and is not redistributed** |
| `determinant/pocket_labels_real.csv` | derived from a public structure-prediction archive | that archive's licence (CC BY 4.0 at the time of writing) |
| `structures/structure_accessions.txt` | identifiers only, in sections by archive: the predicted-structure archive and the experimental-structure archive (§5 says which section is which) | each archive's own terms; the list itself is the authors' |
| `dat/models/…hmm` | third party, redistributed here unmodified | **no licence upstream** — see `dat/models/PROVENANCE.md` and `LICENSE` |

Citing LnFinder does not cite any of these. If you use the data, cite each in its own right.
