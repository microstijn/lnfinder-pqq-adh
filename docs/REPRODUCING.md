# Reproducibility

The long form, moved out of [`README.md`](../README.md), which keeps the summary table.
Read [`docs/DATA.md`](DATA.md) first if you are building a data root.

---

## What is and is not reproducible here

No analysis here runs end to end from a clone alone. Every entry point needs at least one of: a
large derived input or an earlier stage's output, neither shipped; a criteria file supplied per
run; an external binary; a live web service.

The missing half is the archived data set: one directory under one environment variable,
documented file by file in [`docs/DATA.md`](DATA.md). Part of it is published and part is not
(see [Where the inputs come from](#where-the-inputs-come-from)), so tiers 2 and 3 are partly
available today.

### Tier 1: runs from a clone, offline, with nothing else

| what | command |
|---|---|
| the test suite | `julia --project=. test/runtests.jl` |
| checking a data root against its recorded checksums | `julia scripts/verify_data.jl` |
| the environment figure — `julia --project=. scripts/fold_change_figure.jl notitle` | no data set: its three inputs are the tables tracked under `dat/figure/` (`<root>/figure/` at the default root). Writes to `$LNFINDER_FIGDIR`, else beside them. Omit `notitle` for a title; `pooled` or `twophylum` draw the two thermal variants |

`verify_data.jl` needs only Julia, with no instantiated environment. With no data root it lists
every file that is missing and exits non-zero.

### Tier 2: clone plus the data set, offline, no external binary

| what | needs |
|---|---|
| the figure's tables — `python scripts/prevalence.py` then `python scripts/prevalence_phylum_average.py` | the data set, for `<root>/habitat/habitat_taxon_calls.csv` and `<root>/habitat/background_pqq_scan.csv`. They rewrite the three tables under `<root>/figure/` byte for byte, bootstrap columns included, so the tracked copies can be checked. Needed only to re-derive, before `fold_change_figure.jl` |
| habitat step 3 — `scripts/habitat_enrichment.jl` | the data set (steps 1 and 2 are in it), and `--categories` pointed at `criteria/category_selection.csv` in it: required, no default, for the reason below. The two plan bounds it also takes are not in this release; without them a run differs only by twelve provenance keys ([`docs/DATA.md`](DATA.md) §4) |
| habitat step 4 — `scripts/habitat_sensitivity.jl` | the data set, and `--anchor-targets dat/criteria/anchor_targets.csv`, which ships here. `--anchor-reference` is not needed: its anchor levels are not evaluated (`ANCHOR_LEVELS_RETIRED` in the script) and a path given is recorded unopened. Output goes to `--figdir`, else `$LNFINDER_FIGDIR`, else `figs/` in the clone; set one of the first two and pass step 5 the same directory, as step 5 has no fallback. It includes step 3, so step 3's outputs need not exist |
| habitat step 5 — `scripts/habitat_figure.jl` | the step-3 and step-4 outputs, a `--figdir` (else `$LNFINDER_FIGDIR`) that already exists, and the same `--categories` file step 3 was given |

In order, from the root of the clone, with the data set unpacked at `$LNFINDER_DATA`:

```
julia --project=. scripts/habitat_enrichment.jl \
    --categories "$LNFINDER_DATA/criteria/category_selection.csv"
julia --project=. scripts/habitat_sensitivity.jl --figdir <figure directory> \
    --anchor-targets dat/criteria/anchor_targets.csv
julia --project=. scripts/habitat_figure.jl --figdir <figure directory> \
    --categories "$LNFINDER_DATA/criteria/category_selection.csv"
```

Step 3 writes under `habitat/` in the data root, where step 5 reads its rows and arm summary by
default; none of those is archived.

This has been run end to end offline, in a fresh clone against a copy of the data set, with
outbound connections blocked and a deliberate fetch confirmed to fail first. All three
reproduced the recorded outputs value for value (the estimates table, arm summary, grid,
robustness table and the figure's two companion tables); the panel was byte-identical in PNG and
SVG, the PDF differing only in the creation date that format embeds. The data set is released
with the accompanying publication, so a reader can repeat this then. Step 1, in Tier 3, produced
the call table these three start from; a fresh run of it measures the service again instead of
reproducing the deposited table.

### Tier 3: clone plus the data set, plus a network or a binary

| what | needs | note |
|---|---|---|
| control-set validation — `ValidateClassifier.validate(; control_tsv, hmm, ln_ref_acc, ca_ref_acc, …)` | HMMER `hmmalign` on `PATH`, and network (UniProt, plus NCBI for one row) | the only entry point whose data inputs both ship here, the control set and the profile. No command-line driver: load the file in Julia and call the function. Pass one lanthanide-binding and one calcium-binding row from those whose `role` is `reference` (the recorded run used I0JWN7 and P16027); the harness reports whether the column they calibrate is unique. Outputs go where `work_dir` says, and the recorded run's call for every row ships as `dat/control_set/control_validation_predictions.csv` |
| step 0, the background pool — `julia --project=. scripts/prokaryote_background.jl` | network (UniProt reference proteomes, then the habitat service per taxon) | slow and resumable. It writes three files step 1 reads: the taxon pool and the background table, which the data set carries, and the ontology path cache, which it does not. No options; `--help` prints what it reads and writes, and any other argument is refused |
| habitat step 1 — `scripts/habitat_reaggregate.jl` | the data set, which carries the stored OmniCrobe relation responses and the OntoBiotope path cache. Also `--no-gates`, as no gate-targets file is in it, and `--out <path>` pointed **outside the data set**, or the call table overwrites the archived copy | it replays from the stored bodies, so for the taxa the cache covers it needs no network and recomputes the deposited call table from the bytes that table came from; `habitat/omnicrobe_cache/manifest.csv` lists those taxa. Anything outside the cache, a taxon or an ontology term, is fetched live from a re-mined, unversioned service, so those rows are a new measurement. Under `--no-gates` it still makes the gates the method itself fixes (no unresolved ontology path, no failed relation fetch), marks as not checked the gates whose expected values are measurements of the data, declines as not applicable the four gates that re-derive calls the archived census and background tables already carry, and closes by listing what it skipped. [`docs/DATA.md`](DATA.md) §7 has the command and the two files `--out` does not move |
| habitat step 2 — `scripts/background_pqq_scan.jl` | the data set, HMMER `hmmsearch` and `hmmalign`, and network — on Windows also WSL, since HMMER has no native Windows build | several hundred megabytes of downloads, and its outputs are in the data set, so you need not run it |
| the census — `julia --project=. src/scriptLn.jl` | the data set, plus an alignment you regenerate with HMMER and a taxonomy dump you re-fetch (both in [`docs/DATA.md`](DATA.md) §5), plus network, and `$LNFINDER_CENSUS_OUT` naming a directory outside the clone | not offline even with the data set: it resolves accessions against a live service before any other switch applies, so a re-run will not reproduce the census byte for byte. `EXPECTED_CENTROIDS=33428` stops it if the clustering yields another centroid count |

### Tier 4: not reproducible from this repository, with or without the data set

| stage | entry point | what is missing |
|---|---|---|
| determinants and structures | not in this repository | the sequence-and-structure route for what distinguishes the two classes beyond the active-site residue: neither its modules nor a driver for them is in this release, and its label table is not in the published archive either |
| phylogeny | not in this repository | no script here builds a tree, and the Newick utilities that pruned and rooted the trees built elsewhere are not in this release either |

### Where the inputs come from

Every input the analysis reads and this repository does not carry lives under one data root:
`$LNFINDER_DATA`, or `dat/` in the clone when that is unset. Unpack the archived data set there
and every path it carries is satisfied, except the criteria files, which are passed on the
command line (below). [`docs/DATA.md`](DATA.md) §5 has what the data set lacks and what each
step does without it.

| variable | what it names |
|---|---|
| `LNFINDER_DATA` | the data root. One variable moves every input |
| `LNFINDER_AFA` | the alignment of the search sequences to the profile, if kept outside the root. The one input you regenerate yourself, with one HMMER command |
| `LNFINDER_TAXDUMP` | an unpacked NCBI taxonomy dump, if kept outside the root; dumps are commonly shared between projects, hence a variable of its own |
| `LNFINDER_CENSUS_OUT` | where a census run writes. Not an input. Pick a directory **outside** the clone, so a census table and its manifest cannot land among tracked files |
| `LNFINDER_FIGDIR` | the figure directory for habitat steps 4 and 5. Not an input; `--figdir` overrides it, and step 5 never creates it |
| `LNFINDER_RELATION_SEED` | semicolon-separated directories of already-fetched habitat response bodies, priming step 1's relation cache so those taxa are read from disk. Bare arguments to step 1 do the same |

Those six are the whole `LNFINDER_*` set. The census driver reads two more without the prefix:
`RUN_HABITAT` (`false` turns off its habitat layer, its only per-taxon network route besides
accession resolution) and `EXPECTED_CENTROIDS` (checked against its centroid count). Eight in
all, declared in `src/ModuleDataInputs.jl`; two tests fail the build on a ninth, one scanning
for the `LNFINDER_` prefix and one for every environment-variable read by name.

A step that cannot find an input it resolved from the root refuses before doing any work, naming
the input, where in the layout it looked and the variable that moves it (a path you typed is
reported as you typed it). [`docs/DATA.md`](DATA.md) has what every file is, which stage reads
it, how to rebuild what the archive does not carry, and which of them are results.

```
LnFinder: the required input "census_habitat" is not there.
  what it is:  the census table carrying the habitat layer, one row per sequence
  read by:     habitat steps 1, 2, 3 and 4
  looked for:  habitat/census_habitat_rerun_v2.csv, under the data root
  data root:   $LNFINDER_DATA when that is set; dat/ in the clone when it is not
  This repository does not carry that file;
  docs/DATA.md says which stage reads it and how to obtain or rebuild it.
```

Criteria files are passed on the command line deliberately. `--gate-targets` (step 1) and
`--anchor-targets` (step 4) name files whose expected values are measurements of the data, so
the code holds no constant for them; given neither the file nor an explicit refusal
(`--no-gates`, `--no-anchors`) a step refuses to run instead of skipping the check.
`--anchor-reference` names an earlier run's rows for step 4's anchor levels, which are not
evaluated, so it is not required (Tier 2 above).

`--categories`, for steps 3 and 5, has no default for a related reason. Three recorded outcome
rules speak for one habitat category each, the figure tabulates categories across every
specification, and both steps quote one phylum's shares; which categories those are was selected
from this data set's own estimates, with the estimates in hand, so a default in the code would
publish that choice as the software's own. The external ontology mapping they come from predates
every estimate; picking names out of it does not. The phylum is the one entry the data fixes, the
heaviest standardisation stratum: step 5 refuses to assert the identities that tie the quoted
shares to the standardisation weights unless the run's own heaviest stratum is still that
phylum, and it is in the file so a reader can derive the key family named after it. The names
reach the output keys, so a run on other categories emits keys for those. See
[`docs/DATA.md`](DATA.md) §4.

The data set carries `criteria/category_selection.csv`; the anchor targets at
`dat/criteria/anchor_targets.csv` ship here with the two tables they were transcribed from.
Neither carries a gate-targets or anchor-reference file: step 1 runs on the archived census with
`--no-gates` (Tier 3 above), and step 4 compares no reference row.

### "The registered specification" and "the analysis plan"

Both phrases, used often in the code and `docs/TOOLING.md`, name an external record that this
release does not carry; there is no plan document here. It is a dated, checksummed set of rules
mapping each outcome to a colour, a caption phrase and a sentence of the accompanying text,
timestamped before every output of the re-analysis these steps perform, so the numbers reported
here did not choose their own wording. It was not written before the data were first looked at,
and does not claim to have been: it discloses that an earlier run's estimates for the categories
it names, and that run's sensitivity findings, were already known when it was written, and it
covers the re-analysis only. `docs/TOOLING.md`, under *"The registered specification" and "the
analysis plan"*, has the detail.

In the code, it surfaces as the `primary` configuration each script runs under, printed and
recorded by every script, and as the criteria files, each recording where its values came from:
the independent implementation and digest a transcribed anchor value came from, a note per gate
saying where its number came from, and a per-row note that the categories came from the
analysis's own earlier estimates. The
bounds the plan fixed are the two `--plan-*` values, which arrive on the command line and are
not part of this release; step 3's "see THE RECORDED CRITERIA below" points to the section that
reports them against the measured values ([`docs/DATA.md`](DATA.md) §4).

### The archived data set

Part of it is published. The habitat route's inputs and the figure's summary tables are
deposited and open at [doi:10.5281/zenodo.23667326](https://doi.org/10.5281/zenodo.23667326),
under CC BY 4.0. The `census/`, `criteria/`, `determinant/` and `structures/` trees are not:
several of those files are results of the analysis this software was written for, so they are
released with the accompanying publication and not before. Their layout is fixed and the code
reads it, and [`docs/data_manifest.sha256`](data_manifest.sha256) already carries a digest for
every archived file except one named exception, whose reason is in [`docs/DATA.md`](DATA.md) §2.

Two habitat inputs carry no digest, there or in the archive's own manifest, because they are not
in the archive: the stored relation responses, one per taxon, and the ontology path cache beside
them, both third-party content that this project does not redistribute. Recorded in their place
is the digest of the cache's `manifest.csv`, which says which taxa were fetched, by which route
and when, and carries none of the responses.

The shipped record is a second, independent one, so a data root can be checked against this
repository and not only against itself. With both present, `verify_data.jl` checks the
archive's, then cross-checks this one and fails if they disagree:

```
julia scripts/verify_data.jl --root <path>
```

When the archive is published, one line of `docs/DATA.md` §2 gains its identifier and
`CITATION.cff` gains the same one; nothing else about the route changes. Until then, not
everything the archive would carry can be rebuilt from public databases. `census/firstTry.fa` is
the sequence set an `hmmsearch` returned from a database that has changed since; the clusterings
need a proprietary binary that is not redistributed; the habitat relation cache is a frozen
snapshot of a re-mined service, so a fresh fetch returns different data. That is why the
snapshot had to be deposited alongside the call table step 1 computed from it (Tier 3 above).

---
