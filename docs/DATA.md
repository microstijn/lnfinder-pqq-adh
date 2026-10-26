# The data this repository does not carry

LnFinder's code is here; the files it reads are not. They are large, mostly derived, and
several are results of the analysis the software was written for. They are archived as one
data set, and this document is the contract between it and this code: one data root, one
layout underneath it, one name per input. [`README.md`](../README.md) says which stages can
be run at all; this file says what each stage reads and where it looks.

---

## 1. The contract

Everything the analysis code reads that this repository does not ship lives under one
directory, the data root, in a fixed layout.

| | |
|---|---|
| The root | the directory named by `$LNFINDER_DATA`. Unset, it is `dat/` in the clone, where the shipped inputs live and which `.gitignore` covers for the rest. |
| The layout | declared in one place in the code, [`src/ModuleDataInputs.jl`](../src/ModuleDataInputs.jl). The tree is §3. |
| Per-input overrides | `$LNFINDER_AFA` and `$LNFINDER_TAXDUMP`, for the two files you make or fetch and may want outside a project tree. Everything else resolves under the root. |
| When something is missing | the step stops before doing any work and names the input, where the layout puts it and which variable moves it. No stack trace, and no path off anybody's machine. A path you typed yourself is reported as you typed it. |

Unpack the data set at one place, point one variable at it, and every path it carries is
satisfied at once. The criteria files are the exception, passed on the command line for the
reason §4 gives; §5 covers what the data set leaves out.

```
# everything under one root, anywhere you like
export LNFINDER_DATA=/scratch/lnfinder-data   # PowerShell: $env:LNFINDER_DATA = "C:\lnfinder-data"
julia scripts/verify_data.jl                  # check it before running anything
```

A relative value is made absolute when it is read.

Two variables name where a run writes: `$LNFINDER_CENSUS_OUT` (a census run's table and
manifest) and `$LNFINDER_FIGDIR` (the figure directory of habitat steps 4 and 5).
`$LNFINDER_RELATION_SEED` primes step 1's relation cache from directories of already-fetched
response bodies. That is six `LNFINDER_*` variables; the census driver reads two more without
the prefix, `RUN_HABITAT` (`false` turns off its habitat layer, its only per-taxon network
route besides accession resolution) and `EXPECTED_CENTROIDS` (a whole number it checks its
centroid count against, so a clustering that silently changed size stops the run). All eight
are declared in `src/ModuleDataInputs.jl`, and two tests fail the build if a script reads a
ninth: one scans for the `LNFINDER_` prefix, the other for every `ENV["…"]`,
`get(ENV, "…", …)` and `haskey(ENV, "…")` in `src/` and `scripts/`.

---

## 2. How to obtain it

Part of the data set is published: the habitat route's inputs and the figure's summary
tables, the `habitat/` and `figure/` trees of §3.

> Data set DOI: [`10.5281/zenodo.23667326`](https://doi.org/10.5281/zenodo.23667326), under
> CC BY 4.0. That is the concept DOI, so it resolves to whichever version is current.

Unpack that archive at the data root and every required entry of
[`docs/data_manifest.sha256`](data_manifest.sha256) is satisfied.

The `census/`, `criteria/`, `determinant/` and `structures/` trees are not in it, nor are
three of step 2's bulk working files. Several are results of the analysis this software was
written for (§6), so they are released with the accompanying publication, and
`src/ModuleDataInputs.jl` marks them `:unreleased`: a step that needs one names it and says
why it is absent, and `verify_data.jl` does not require it.
[`docs/REPRODUCING.md`](REPRODUCING.md) says which steps the published archive is enough for.

### Verifying what you got

[`docs/data_manifest.sha256`](data_manifest.sha256) carries the sha256 digest of the archived
files as assembled, in plain `sha256sum` format, so a data root can be checked two ways:

```
julia scripts/verify_data.jl                       # also checks the contract itself
julia scripts/verify_data.jl --root <path> --manifest <path>

cd <data root> && sha256sum -c .../docs/data_manifest.sha256    # shasum -a 256 -c on macOS
```

Without the optional files, `sha256sum -c` reports each of them as FAILED;
`--ignore-missing` skips them, and with them any missing required file, which
`verify_data.jl` does not. `verify_data.jl` needs nothing but Julia, not even an instantiated
environment: it hashes every file the manifest lists, names what is missing and what differs,
reports any declared input the root does not satisfy, and **exits non-zero** if anything is
wrong.

The archive carries its own `MANIFEST.sha256`, the complete record, which `verify_data.jl`
prefers at the root; the repository's copy is still consulted for the paths it covers, since
that record alone shows only that the archive agrees with itself. The repository's copy omits
the relation cache's thirteen thousand stored response bodies, digesting the cache's
`manifest.csv` and the ontology path cache instead, and `structures/structure_accessions.txt`,
one identifier per member of the census, which the published archive omits too. The twelve
experimental structures further down do have digests, and files the archive does not carry yet
keep theirs in the `optional` section, `criteria/category_selection.csv` among them, at the
digest of the copy the archive will be assembled from.

Order matters. Step 4 records `anchor_targets_sha256` in its provenance and step 5 the digest
of that provenance file, each of whichever copy the run read, so the tracked copy of
`dat/criteria/anchor_targets.csv` has to be the copy the published step-4 run read: any edit
to it, a reworded note included, is followed by a step-4 and a step-5 run before their
provenance is published. The two files beside it are pinned the other way round, the targets
file declaring their digests (§4). `category_selection.csv` is pinned twice, as
`selection_sha256` in step 3's arm summary and `input_categories_sha256` in step 5's figure
provenance.

---

## 3. The layout

```
<data root>/                          $LNFINDER_DATA, else dat/ in the clone
  MANIFEST.sha256                     the archive's own checksums
  habitat/                            the habitat comparison, steps 0 to 5
    census_habitat_rerun_v2.csv
    prokaryote_taxa.csv
    prokaryote_background.csv
    obt_path_cache.csv                written by steps 0 and 1; see §5
    habitat_taxon_calls.csv
    background_pqq_scan.csv
    background_scan/
      recall.txt  threshold.txt  acceptance.tsv
      uniprot_release.txt  uniprot_release_IPR017512.txt
      refs.fasta  recall_proteomes.fasta.manifest  calibration_scores.csv
    omnicrobe_cache/                  see §5
      manifest.csv
      relations/<taxid>.json
    habitat_enrichment.csv            written by step 3, not archived
    habitat_arm_summary.csv           written by step 3, not archived
    habitat_enrichment_anchor.csv     written by step 3, not archived; no stage reads it
  figure/
    prevalence.csv  prevalence_phylum_average.csv  thermal_options.csv
  criteria/                           passed on the command line
    category_selection.csv
  census/
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

`dat/criteria/` in the clone also holds three tracked files, `anchor_targets.csv` and the two
its cells were transcribed from (§4). They resolve from the repository whatever
`$LNFINDER_DATA` says, and the archive carries no file of those three names, so an archive
unpacked at `dat/` puts `category_selection.csv` beside them and overwrites none of them.

The archive carries `habitat/omnicrobe_cache/manifest.csv`, the stored bodies under
`relations/` and `obt_path_cache.csv` beside the cache (§5). It is 26 MB compressed, about
200 MB unpacked, of which `habitat/habitat_taxon_calls.csv` and the thirteen thousand bodies
under `habitat/omnicrobe_cache/relations/` are the bulk. The largest file not in it yet is
`census/firstTry.fa`; §4 lists the rest.

---

## 4. Every file: what it is, who reads it

The `name` column is the name the code uses; the layout is in `src/ModuleDataInputs.jl`, and
every refusal quotes both.

### Shipped with the code

| name | path | size | what it is |
|---|---|---|---|
| `profile_hmm` | `dat/models/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm` | 0.26 MB | the profile the cofactor call is calibrated against. Third-party, unlicensed upstream; see `dat/models/PROVENANCE.md`. Read by the census driver, habitat step 2 and control-set validation |
| `control_set` | `dat/control_set/known_cofactor_mdhs.tsv` | 0.01 MB | enzymes whose cofactor is known from the literature, one row per published enzyme with its citation. Read by control-set validation |
| `anchor_targets` | `dat/criteria/anchor_targets.csv` | 1.1 MB | the transcribed anchor cells step 4 checks its grid against: the registered specification and each of the fourteen registered levels, on three habitat and three chemistry categories, a tolerance and a note per cell, the levels' names, and the name and sha256 of the two files the cells came from. `scripts/habitat_sensitivity.jl --anchor-targets dat/criteria/anchor_targets.csv`. It carries result values (§6) |
| `anchor_source_primary` | `dat/criteria/impl_4.csv` | small | the independent implementation those cells come from: an `anchor,unit,category,quantity,value` table computed from the call table and the scan table by code written from the specification without reading this repository's grid script. Step 4 refuses to publish when its bytes are not the declared digest |
| `anchor_source_variants` | `dat/criteria/impl_w1.csv` | 0.03 MB | the same, for the fourteen registered levels |
| `figure_tables` | `dat/figure/prevalence.csv`, `dat/figure/prevalence_phylum_average.csv`, `dat/figure/thermal_options.csv` | 0.004 MB | the three tables the environment figure draws from, derived: `scripts/prevalence.py` and `scripts/prevalence_phylum_average.py` rebuild them byte for byte from `habitat/habitat_taxon_calls.csv` and `habitat/background_pqq_scan.csv`. Tracked so `scripts/fold_change_figure.jl` renders the figure from a bare clone |

These travel with the code and do not move when the data root moves. The three under
`dat/criteria/` are stored byte for byte, marked `-text` in `.gitattributes`: step 4 digests
the bytes on disk, so a checkout that rewrote their line endings would fail its
`target sources` gate on files nobody edited. Step 4 reports a missing source without
refusing, so `test/test_data_inputs.jl` checks that both sit beside the targets with the
declared bytes.

Every anchor target carries a note saying what the specification behind it did and did not
determine. The recorded analysis plan does not determine every rule: for three the
specification takes the rule from this repository's implementation, the thin-stratum collapse
thresholds shared by every level, four levels' rule definitions (`T:alphabetical`,
`T:unresolved`, `C:share`, `C:winner_fractional`) and one edge case of a fifth
(`I:iso_taxon`). There a target checks that the grid computes the rule it states, and not the
choice of rule.

One more shipped file is a result: `dat/control_set/control_validation_predictions.csv`
(small), the `control_predictions.csv` control-set validation writes, with the classifier's
call on each control in the columns `original_header,accession,classification,active_site`.
No code reads it; it is tracked as the record of that validation (§6).

### `habitat/` — the habitat comparison

| name | path | size | what it is, and which stage reads it |
|---|---|---|---|
| `census_habitat` | `habitat/census_habitat_rerun_v2.csv` | 8.5 MB | the census table with the habitat layer, one row per sequence; steps 1 to 4 read it, nothing here writes it. Its identifier, cofactor-call, active-site, `pos2_aa` and taxonomy columns (`TaxID` to `species`) are those of `census/final_deduplicated_data.csv`, row for row; its seventeen habitat-layer columns (`report_group`, `isolation`, `chemistry_share` and the `habitat_*`, `published_*` and `pillar_*` families) are joined on `id` outside this release and read only by step 1's replication report. Habitat names: §8 |
| `pool_taxa` | `habitat/prokaryote_taxa.csv` | 0.34 MB | the reference-proteome taxon pool the arm is drawn from. Steps 1 and 2, step 3's provenance; step 0 writes it, from a live query whose answer moves |
| `pool_background` | `habitat/prokaryote_background.csv` | 0.62 MB | the habitat annotation of that pool. Step 1, step 3's provenance; step 0 writes it |
| `obt_path_cache` | `habitat/obt_path_cache.csv` | — | the habitat-ontology climb, cached. Steps 0 and 1 read and rewrite it. In the archive (§5): a verbatim third-party extract keyed by ontology term id, no taxon and no census value (§8) |
| `relation_cache` | `habitat/omnicrobe_cache` | — | the relation cache: the directory the stored bodies and their manifest live in. Step 1. In the archive (§5) |
| `relation_bodies` | `habitat/omnicrobe_cache/relations` | — | the stored bodies, one `<taxid>.json` per taxon. In the archive (§5), unmodified. Step 1 requires them by name, not just their directory, unless given seed directories |
| `relation_manifest` | `habitat/omnicrobe_cache/manifest.csv` | 0.69 MB | what the cache holds: one row per taxon with a relation count, a source, a fetch stamp and the service and ontology releases the run records, and none of the responses. Steps 1 and 3; step 1 rewrites it in place, wherever `--out` sends the call table |
| `habitat_calls` | `habitat/habitat_taxon_calls.csv` | 24.4 MB | one call per taxon per rule variant, and the archive's entry point to steps 3, 4 and 5. They read it; step 1 writes it, and cannot itself be replayed from the archive (§5). Habitat names: §8 |
| `background_scan_dir` | `habitat/background_scan` | 0.14 MB | the scan's own records, listed in §3. Step 3 reads four of them; without them it records the evidence as unknown and withholds the background licence, which is neither an error nor a pass. The rest record the sequences and release behind the threshold |
| `background_scan_table` | `habitat/background_pqq_scan.csv` | 0.26 MB | which taxa of the arm carry the enzyme at the calibrated threshold. Steps 3 and 4 read it; step 2 writes it, at the cost of a network crawl, HMMER and several hundred megabytes of downloads. Its frame is the pool less the census members the scan's run saw; step 3 names every taxon in one frame and not the other, and withholds the complete-removal statement when they differ |
| `scan_recall` | `habitat/background_scan/recall.txt` | small | the scan prefilter's measured recall, with the counts behind it. Step 3 |
| `scan_threshold` | `habitat/background_scan/threshold.txt` | small | the calibrated score threshold and how it was set. Step 3 |
| `scan_acceptance` | `habitat/background_scan/acceptance.tsv` | small | the scan's acceptance rows, one per check. Step 3 |
| `uniprot_release` | `habitat/background_scan/uniprot_release.txt` | small | the sequence-database release the arm was built from. Step 3 |
| `enrichment_table` | `habitat/habitat_enrichment.csv` | — | written by step 3, not archived. The estimates, one row per contrast. Step 5; step 4 includes step 3 and recomputes |
| `arm_summary` | `habitat/habitat_arm_summary.csv` | — | written by step 3, not archived. What each arm was made of, and every provenance key the panel echoes. Step 5 |

### `criteria/` — what a run is held to

The criteria are passed on the command line because their values are measurements of this
data set, or selections made from it. A step given none of them refuses to run. The archive
carries one; `anchor_targets.csv` is tracked here (above).

| name | path | size | what it is |
|---|---|---|---|
| `category_selection` | `criteria/category_selection.csv` | small | which categories the analysis is about. `--categories <path>` to steps 3 and 5 |

Two criteria options take no file from this data set, and say so in their output.

- Step 1's `--gate-targets` names a file of expected values for its measured replication
  gates. None is carried for this census, the expected values on record being measurements of
  another one, so the call table here is built with `--no-gates` and step 1's report states
  gate by gate that those gates were not checked (§7). Its four reproduction gates, which
  compare a run with the call table and the habitat columns of the census they were written
  for, need a `--gate-targets` file declaring that census by its sha256; otherwise each is
  reported as not applicable with its reason.
- Step 4's `--anchor-reference` names an earlier run's rows for the all-legacy anchor level.
  That level reproduces rows of another census, so it is not evaluated
  (`ANCHOR_LEVELS_RETIRED` in `scripts/habitat_sensitivity.jl`) and step 4 needs
  `--anchor-targets` only: a path given anyway is recorded and not opened,
  `anchor_targets.csv` declares `reference_rows,0`, and the provenance records
  `anchor_reference_status` as not applicable.

#### Why `category_selection.csv` is not in the code

Three of the recorded outcome rules speak for one habitat category each, the figure step
tabulates categories' intervals across every registered specification, and both steps quote
the shares of one phylum beside the standardisation. Which categories those are was chosen
from this data set's own estimates, with the estimates in hand. The classification scheme the
names come from maps an external ontology and predates every estimate
(`src/ModuleEnvironmentTools.jl`); selecting four names out of it does not.

The phylum is the one determined entry: the heaviest standardisation stratum of the analysed
arms, so a run on other data has another one. The summary's share keys are named after it,
and the figure step asserts the identities tying those two shares to the standardisation
weights only while the run's own heaviest stratum is still that phylum. It was written down
after the shares were computed, so it is not pre-specified either.

The selection has no default anywhere in the code, and a run that does not supply it refuses,
naming the option. The file is a four-column table, `kind,role,value,note`, with one `select`
row per role, one `expect,selections,<n>` row stating its own size, and a required `note` per
row saying where that selection came from. Roles left unfilled, roles the step does not know,
a role filled twice, a missing note, a count that disagrees with the rows, or two selections
whose names reduce to one summary-key prefix are each a refusal.

The summary's rule-keyed families and phylum-share keys are derived from the selected names by
lower-casing and removing non-alphanumeric characters, so a run given other categories emits
keys for those. Each run records the file's digest and every selected name beside its numbers:
`selection_*` in step 3's summary, `input_categories_sha256` and `selection_*` in step 5's
provenance.

Two criteria of step 3 are not files, and are not in the archive.
`scripts/habitat_enrichment.jl` takes `--plan-bg-size-below` and `--plan-bg-phylum-pct-max`,
two upper bounds written down in the analysis plan before the step ran, supplied on the
command line and never defaulted. They came from an earlier, approximate version of the same
exclusion, so nothing gates on them and a violated one is reported without being fatal. They
are not part of this release, because the plan is not either (see *"The registered
specification" and "the analysis plan"* in `README.md`). A run without them is identical on
every computed quantity, with twelve fewer
provenance keys: the `criterion_bg_size_*` and `criterion_bg_<phylum>_pct_*` families, which
record whether each bound was met (the second named after the phylum the file selects).

### `census/` — the census route

| name | path | size | what it is |
|---|---|---|---|
| `hits_fasta` | `census/firstTry.fa` | 33.0 MB | the unaligned sequences the census was built from, as retrieved. Read by `src/scriptLn.jl`; archived, so the census is rebuilt from these exact sequences |
| `uc_clusters` | `census/clusters.uc` | 4.1 MB | the clustering of those sequences that defines the census centroids. Read by `src/scriptLn.jl`. Archived because regenerating it needs USEARCH 11, a proprietary binary that cannot be redistributed; the command is in `src/scriptLn.jl` (`USEARCH_CMD`) |
| `uc_clusters_50` | `census/clusters_50.uc` | 4.3 MB | the coarser clustering used as the sub-family unit of the phylogenetic control. No stage in this release reads it; archived because it too needs that binary |
| `census_table` | `census/final_deduplicated_data.csv` | 7.6 MB | the census as `src/scriptLn.jl` writes it: one row per centroid with its cofactor call, active-site residues, TaxID and lineage, and the route, entry status and sequence-database release each TaxID was resolved through. Written with `RUN_HABITAT=false`, so its `report_group` reads `_Not Mapped_` on every row. No stage reads it; it is here so the habitat input, which carries its identifiers, calls and taxonomy row for row, can be checked against its source |

### `determinant/` and `structures/`

| name | path | size | what it is |
|---|---|---|---|
| `pocket_labels` | `determinant/pocket_labels_real.csv` | 2.8 MB | the per-accession metal label set of the whole census, as the structure route produced it. No driver in this release reads it, the determinant and pocket modules shipping as a library with unit tests; it is archived because rebuilding it would mean every predicted structure plus a full pocket export |
| `structure_accessions` | `structures/structure_accessions.txt` | not stated (see below) | one identifier per structure the analysis read, in sections by archive. The structures themselves are not archived here; see §5 |

A dash in the size column marks an entry whose bytes the archive does not carry.
`structures/structure_accessions.txt` is carried without a figure, for the reason §2 gives
for its digest: one identifier per line is nearly its own line count. Sizes here are for
planning a download, which §3's total covers; verification is by digest, and
`verify_data.jl` and `sha256sum -c` compare hashes, never sizes.

---

## 5. What the archive leaves out

Four things the analysis read are not in the archive: the model alignment, the taxonomy dump,
the structures and step 2's bulk downloads. The first three have digests in
`docs/data_manifest.sha256` under *optional*, so what you make or fetch can be compared
against what the analysis read. Step 2's bulk download cache has no digest, being hundreds of
megabytes of third-party response bodies no later step reads, while the records step 3 reads
from that step are archived in their own right. The habitat relation cache and the ontology
climb beside it are a fifth case: in the archive, with no digest here.

### `model_alignment` — `derived/aligned_to_model.afa` (348 MB): regenerate it

One HMMER run over a file that is archived, against the profile that ships here:

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
<https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/>. Unpack it so that `nodes.dmp`,
`names.dmp`, `merged.dmp` and `delnodes.dmp` sit directly under `derived/taxdump/`; the loader
opens those four and ignores the rest. `$LNFINDER_TAXDUMP` points at it elsewhere. It will not
match its recorded digest: the archive holds snapshots at month boundaries, and the dump the
analysis used was taken between two of them. Take the nearest snapshot at or before
**2025-02-05**, and treat the digest as a record of how far your copy has drifted.

### `structures` — `derived/structures/`: re-fetch by identifier

Each is a public record with a stable identifier, and `structures/structure_accessions.txt`
carries every identifier the analysis used, in sections:

| section | archive | URL pattern |
|---|---|---|
| `afdb` | predicted structures | `https://alphafold.ebi.ac.uk/files/AF-<ACCESSION>-F1-model_v4.pdb` |
| `afdb_subsample` | a subset of the same set, listed separately because the selection is recorded nowhere else | as above |
| `rcsb` | experimental structures | `https://files.rcsb.org/download/<ID>.pdb` |

The predicted set is version **v4** of that archive (`ALPHAFOLD MONOMER V2.0`, headers dated
`01-JUN-22`). Later releases change the files, so a re-fetch today may not be byte-identical
and some identifiers may no longer resolve; the digests of the twelve experimental structures
are recorded. Save them as `derived/structures/crystals/<lower-case id>.pdb` if you want
`verify_data.jl` to check them.

### `habitat/background_scan/` bulk downloads (about 640 MB)

Step 2's download cache and its HMMER working files are not archived. The *records* step 3
reads are (§4), and the caches themselves are re-fetchable from the sequence database at the
release pinned in `uniprot_release.txt`.

### `relation_cache`, `relation_bodies` and `obt_path_cache`: in the archive

`habitat/omnicrobe_cache/relations/<taxid>.json` is one stored response body per taxon from
the microbe–habitat service the code queries; `habitat/obt_path_cache.csv` is the climb over
the habitat ontology those responses are read against, an `obtid,path` table keyed by ontology
term id with no taxon and no census value. Both are third-party content as retrieved, both are
in the published archive unmodified, under the service's own CC BY terms (§8), and beside them
it carries `habitat/habitat_taxon_calls.csv` and `habitat/omnicrobe_cache/manifest.csv`, which
records which taxon was fetched, from where and when, with the service and ontology releases.

Step 1 reads the stored bodies for every taxon the cache covers, which `manifest.csv` names,
so the archived call table can be recomputed from the bytes it was computed from. The service
is text-mined, unversioned and periodically re-mined: a query made today returns a different
data set, so a step-1 run that fetches is a fresh measurement, and neither cache has a digest
here.

Without the caches, step 1 refuses before doing any work, naming the input, where the layout
puts it, which variable moves it, and that without it the step fetches. The one alternative it
names: prime the cache from directories of already-fetched `<taxid>.json` bodies, as bare
arguments or a semicolon-separated list in `$LNFINDER_RELATION_SEED`. Seeding is not a replay:
an uncovered taxon is fetched live, as is any ontology term the OBT cache does not hold. Cite
the service and the ontology in their own right (§8).

---

## 6. Which of these are results

These files are outputs of the unpublished analysis, which is why the data set is released
with the accompanying publication:

| file | what it discloses |
|---|---|
| `census/final_deduplicated_data.csv` | the census: one row per centroid with its cofactor call, active-site residues and resolved taxonomy |
| `habitat/census_habitat_rerun_v2.csv` | the census table carrying the habitat layer as published |
| `habitat/habitat_taxon_calls.csv` | per-taxon habitat calls with the carrier counts the comparison is computed from |
| `habitat/background_pqq_scan.csv` | which background taxa carry the enzyme at the calibrated threshold |
| `habitat/background_scan/` records | the calibrated threshold and its measured recall: small files, but fitted quantities |
| `determinant/pocket_labels_real.csv` | the per-accession metal label set of the whole census |
| `criteria/category_selection.csv` | the categories the analysis's own estimates led to |

This repository carries result values too. `dat/criteria/anchor_targets.csv` and the two files
beside it hold the registered specification's and the registered levels' estimates on six
categories, computed by implementations independent of this code; they ship because step 4
checks its grid against them, and a check a reader cannot run is not a check.
`dat/control_set/control_validation_predictions.csv` holds the classifier's call on every
control of the shipped control set.

The rest are not results: `habitat/prokaryote_taxa.csv` is a public taxon list,
`habitat/omnicrobe_cache/manifest.csv` is bookkeeping, `census/firstTry.fa` is public sequence
as retrieved, `census/clusters.uc` and `census/clusters_50.uc` are clusterings that carry no
cofactor claim, and the third-party caches `habitat/omnicrobe_cache/relations` and
`habitat/obt_path_cache.csv` are content as retrieved (§5).

---

## 7. Which stage needs which file

| stage | reads | writes |
|---|---|---|
| census — `src/scriptLn.jl` | `model_alignment`, `hits_fasta`, `taxdump`, `uc_clusters`, `profile_hmm`; and a sequence database over the network, always | `$LNFINDER_CENSUS_OUT/census_rebuild_<date>/` |
| habitat step 0 — `scripts/prokaryote_background.jl` | the network; on a rerun, its own three outputs (it is resumable) | `pool_taxa`, `pool_background`, `obt_path_cache`, of which the archive carries the first two |
| habitat step 1 — `scripts/habitat_reaggregate.jl` | `census_habitat` (or the file `--census` names), `pool_taxa`, `pool_background`, `relation_bodies`, `relation_manifest`, `obt_path_cache` (the three cache inputs are in the archive, §5), and `--no-gates` (§4); a reader's own `--gate-targets` file makes the four reproduction gates only if it carries one `census,,<sha256>,,<note>` row matching the digest of the census read | `habitat_calls` (or `--out`), and `relation_manifest` and `obt_path_cache`, always in place under the data root |
| habitat step 2 — `scripts/background_pqq_scan.jl` | `profile_hmm`, `census_habitat`, `pool_taxa`; plus HMMER and the network | `background_scan_table`, `background_scan_dir` |
| habitat step 3 — `scripts/habitat_enrichment.jl` | `habitat_calls`, `background_scan_table`, `census_habitat`, `scan_recall`, `scan_threshold`, `scan_acceptance`, `uniprot_release`, `relation_manifest`; `pool_taxa` and `pool_background` for its provenance; `category_selection` as `--categories`, required, no default; and two plan bounds this release does not carry (§4) | `enrichment_table`, `arm_summary`, and the legacy configuration's rows beside `--out` (`habitat_enrichment_anchor.csv`), which no stage reads |
| habitat step 4 — `scripts/habitat_sensitivity.jl` | `habitat_calls`, `background_scan_table`, `census_habitat`, and `anchor_targets` on the command line, with `anchor_source_primary` and `anchor_source_variants` beside it; no `--anchor-reference` (§4). It includes step 3 and recomputes | the grid, its robustness table, its panel and its provenance, under `--figdir`, else `$LNFINDER_FIGDIR`, else `figs/` in the clone |
| habitat step 5 — `scripts/habitat_figure.jl` | `enrichment_table`, `arm_summary`, `habitat_calls`, step 4's outputs from the figure directory, and `category_selection` as `--categories`, the same file step 3 was given | the panel and its five companion files, under `--figdir`, else `$LNFINDER_FIGDIR` |
| control-set validation | `control_set`, `profile_hmm`; plus HMMER `hmmalign` and the network | its own outputs, where you tell it |

**Step 1 writes into the data root, so give it `--out`.** Its call table goes to
`habitat/habitat_taxon_calls.csv` under the root unless `--out <path>` says otherwise, and that
file is in the archive, its digest is in `docs/data_manifest.sha256`, and it is where steps 3,
4 and 5 start. A run with no `--out` replaces it with a fresh fetch's output, a different data
set (§5), which no re-run brings back.

```
julia --project=. scripts/habitat_reaggregate.jl \
    --no-gates \
    --out <a directory of your own>/habitat_taxon_calls.csv
```

`--out` moves the call table and nothing else: the relation cache's `manifest.csv` is
rewritten in place under the data root either way, and so is the OBT cache, on every run. Both
are in the archive, so run step 1 against a copy of the root or keep your own copy of those
two files. `--census <path>` runs step 1 on a census other than the deposited one; there the
four reproduction gates are declined as not applicable, and the run says so, unless a
`--gate-targets` file declares that census by its sha256 (§4).

### Which steps need a network

- With the archive in place, habitat steps 3, 4 and 5 run with no network and no external
  binary; nothing in them reaches the service the habitat layer comes from.
- Step 1's two network paths, the relation fetch and the ontology climb, are both answered by
  the archive (§5), bounded by what the cache covers: `manifest.csv` says which taxa those
  are, and a taxon outside it is fetched, which makes that run a fresh measurement.
- Steps 0 and 2 are not offline, and need not be: their outputs are archived, but for the
  ontology cache step 0 writes (§5). Step 2 also needs HMMER, and on Windows WSL.
- The census route is not offline even with the archive: `src/scriptLn.jl` resolves accessions
  against a live service before any other switch applies, so a re-run will not reproduce the
  census byte for byte.
- Nothing here makes the upstream data collection re-runnable: the database search behind
  `census/firstTry.fa`, the proteome pull behind step 0 and the pull behind step 2 all went to
  live services whose contents have moved since.

---

## 8. Third-party content

| item | origin | terms |
|---|---|---|
| `habitat/habitat_taxon_calls.csv`, the habitat-layer columns of `habitat/census_habitat_rerun_v2.csv` and of `habitat/prokaryote_background.csv`, and `habitat/omnicrobe_cache/manifest.csv` | this project's own tables, computed from relations the microbe–habitat service returns, over the habitat ontology those responses are climbed against; the habitat names in them are that ontology's term labels. The manifest is the cache's fetch record, with a relation count per taxon and none of the retrieved content | the tables are the authors'. The service's own stored responses and the ontology extract beside them are redistributed, in the archive and unmodified, under the service's and the ontology's own terms; cite each in its own right |
| `habitat/omnicrobe_cache/relations`, `habitat/obt_path_cache.csv` | response bodies from that service, as retrieved, and the climb over that ontology. Not redistributed here (§5) | the service's own terms and the ontology's own licence; cite each in its own right |
| `census/firstTry.fa`, `habitat/background_scan/refs.fasta` | public sequence records, the second being the two calibration reference records as retrieved | the sequence database's own licence (CC BY 4.0 at the time of writing); cite the release |
| `census/clusters.uc`, `census/clusters_50.uc` | output of USEARCH 11 | the output is the authors'. The binary is not, and is not redistributed |
| `determinant/pocket_labels_real.csv` | derived from a public structure-prediction archive | that archive's licence (CC BY 4.0 at the time of writing) |
| `structures/structure_accessions.txt` | identifiers only, in sections by archive: the predicted-structure archive and the experimental-structure archive (§5 says which is which) | each archive's own terms; the list itself is the authors' |
| `dat/models/…hmm` | third party, redistributed here unmodified | no licence upstream; see `dat/models/PROVENANCE.md` and `LICENSE` |

Citing LnFinder does not cite any of these. If you use the data, cite each in its own right.
