# What is reproducible here, and what is not

The long form, moved out of [`README.md`](../README.md), which keeps the summary table.
Read [`docs/DATA.md`](DATA.md) first if you are building a data root.

---

## What is and is not reproducible here

Honest summary: **no analysis in this repository runs end to end from a clone alone.** Every
analysis entry point needs at least one of: a large derived input that is not shipped, an
output of an earlier stage that is not shipped, a criteria file that is supplied per run, an
external binary, or a live web service.

There is a defined route to the missing half: **the repository plus the archived data
set**, which is one directory in a documented layout under one environment variable. What
that route does and does not buy is set out tier by tier below, and in full in
[`docs/DATA.md`](docs/DATA.md). The data set is not public yet — see
[Where the inputs come from](#where-the-inputs-come-from) — so tiers 2 and 3 describe what
will be possible with it, not something you can do today with a clone alone.

### Tier 1 — runs from a clone, offline, with nothing else

| what | command |
|---|---|
| the test suite | `julia --project=. test/runtests.jl` |
| checking a data root against its recorded checksums | `julia scripts/verify_data.jl` |
| **the environment figure** — `julia --project=. scripts/fold_change_figure.jl notitle` | needs no data set: its three inputs are the small tables tracked under `dat/figure/`, which is `<root>/figure/` when the root is left at its default. Writes to `$LNFINDER_FIGDIR`, else beside them. Omit `notitle` for the titled version; `pooled` or `twophylum` draw the two thermal variants |

`verify_data.jl` needs nothing but Julia, not even an instantiated environment. With no data
root it reports every file that is missing and exits non-zero, which is also a way to read
the list of what the next tier needs.

### Tier 2 — clone plus the data set, offline, no external binary

| what | needs |
|---|---|
| **the figure's tables** — `python scripts/prevalence.py` then `python scripts/prevalence_phylum_average.py` | the data set, for `<root>/habitat/habitat_taxon_calls.csv` and `<root>/habitat/background_pqq_scan.csv`. They rewrite the three tables under `<root>/figure/` **byte for byte**, bootstrap columns included, so the tracked copies can be checked rather than trusted. Run them before `fold_change_figure.jl` only if you want to re-derive; the figure does not need them |
| **habitat step 3** — `scripts/habitat_enrichment.jl` | the data set (steps 1 and 2 are already in it, so step 3 can be run without running them), and **`--categories` pointed at `criteria/category_selection.csv` in it** — required, with no default, because which categories the recorded rules speak for was selected from this data set's own estimates. It also takes two plan bounds on the command line that this release does not carry; a run without them differs only by twelve provenance keys ([`docs/DATA.md`](docs/DATA.md) §4) |
| **habitat step 4** — `scripts/habitat_sensitivity.jl` | the data set, and **`--anchor-targets dat/criteria/anchor_targets.csv`**: the targets ship with this repository, with the two tables they were transcribed from beside them. `--anchor-reference` is not needed: the anchor levels whose rows it would be compared with are not evaluated (`ANCHOR_LEVELS_RETIRED` in the script), and a path given is recorded without being opened. Its outputs go to `--figdir`, else `$LNFINDER_FIGDIR`, else `figs/` in the clone; set one of the first two, and give step 5 the same directory, since step 5 has no such fallback. **Not step 3's outputs**: this step includes step 3 and recomputes, so it runs whether or not they are on disk |
| **habitat step 5** — `scripts/habitat_figure.jl` | the step-3 and step-4 outputs, a `--figdir` (else `$LNFINDER_FIGDIR`) that already exists, and **the same `--categories` file step 3 was given** |

In order, from the root of the clone, with the data set unpacked at `$LNFINDER_DATA`:

```
julia --project=. scripts/habitat_enrichment.jl \
    --categories "$LNFINDER_DATA/criteria/category_selection.csv"
julia --project=. scripts/habitat_sensitivity.jl --figdir <figure directory> \
    --anchor-targets dat/criteria/anchor_targets.csv
julia --project=. scripts/habitat_figure.jl --figdir <figure directory> \
    --categories "$LNFINDER_DATA/criteria/category_selection.csv"
```

Step 3 writes its outputs under `habitat/` in the data root, where step 5 reads its rows and
its arm summary by default; none of them is an archived file.

These three are the stages whose every input the data set or this repository supplies, and
none of them reaches the network: steps 3 and 4 read the call table, the census table, the
background pool and the scan tables the data set carries, and step 5 reads what steps 3 and 4
wrote. **This has been run end to end** — offline, in a fresh clone of this repository,
against a copy of the assembled data set, with outbound connections blocked and a deliberate
fetch confirmed to fail before the run began. All three stages completed without reaching the
network, and reproduced the recorded outputs: every value of the estimates table, the arm
summary, the grid, the robustness table and the figure's two companion tables was identical,
and the rendered panel was byte-identical in PNG and SVG, the PDF differing only in the
creation date that format embeds. The data set is released with the accompanying publication,
so a reader can repeat this once it is out. **Step 1 is not among them**: it produced the call
table the three start from, it is in Tier 3, and what it produces is a new measurement rather
than a reproduction of the deposited table.

### Tier 3 — clone plus the data set, plus a network or a binary

| what | needs | note |
|---|---|---|
| **control-set validation** — `ValidateClassifier.validate(; control_tsv, hmm, ln_ref_acc, ca_ref_acc, …)` | HMMER `hmmalign` on `PATH`, and network (UniProt, plus NCBI for one row) | the only analysis entry point whose **data** inputs are both in this repository: the control set and the profile. There is no command-line driver — load the file in Julia and call the function. For the two reference accessions, pass one lanthanide-binding and one calcium-binding row from those whose `role` is `reference`; the recorded run passed I0JWN7 and P16027, and the harness reports whether the column they calibrate is unique. It writes its own outputs where `work_dir` says; this repository ships one of them from the recorded run, the call for every row, as `dat/control_set/control_validation_predictions.csv` |
| **step 0, the background pool** — `julia --project=. scripts/prokaryote_background.jl` | network (UniProt reference proteomes, then the habitat service per taxon) | takes a long time and is resumable. It writes three files step 1 reads: the taxon pool and the background table, both of which the data set carries, and the ontology path cache, which it does not. It takes **no options**: `--help` prints what it reads and writes and exits, and any other argument is refused rather than ignored |
| **habitat step 1** — `scripts/habitat_reaggregate.jl` | the data set, **plus** the stored OmniCrobe relation responses and the OntoBiotope path cache, which the data set does not carry — so network, unless you prime the cache from bodies you already hold, given as bare arguments or in `$LNFINDER_RELATION_SEED`. Also `--no-gates`, because the data set carries no gate-targets file, and **`--out <path>` pointed outside the data set**, because without it the call table is written over the archived copy | **what this step produces is a new measurement, not a reproduction of the deposited call table.** The responses are a third-party service's own output and are not redistributed with the data set, so every taxon a seed does not cover is fetched live, and so is every ontology term the path cache does not hold; the service is re-mined and unversioned, so a fetch returns different data rather than the bytes the deposited table was computed from. Given neither a seed nor a cache the step refuses, naming what it wants, rather than fetching silently. Under `--no-gates` it makes the gates the method itself fixes (no ontology path left unresolved, no failed relation fetch), reports each gate whose expected value is a measurement of the data as not checked, and declines as not applicable the four gates that re-derive the calls the archived census and background tables already carry; its closing line says which gates it did not make. See [`docs/DATA.md`](docs/DATA.md) §7 for the command and for the two files `--out` does not move |
| **habitat step 2** — `scripts/background_pqq_scan.jl` | the data set, HMMER `hmmsearch` and `hmmalign`, and network — on Windows also WSL, since HMMER has no native Windows build (see [Platform note](#platform-note-hmmer-through-wsl-on-windows)) | several hundred megabytes of downloads. Its outputs are in the data set, so with the data set you do not need to run it either |
| **the census** — `julia --project=. src/scriptLn.jl` | the data set, **plus** an alignment you regenerate with HMMER and a taxonomy dump you re-fetch (both in [`docs/DATA.md`](docs/DATA.md) §5), **plus** network, and `$LNFINDER_CENSUS_OUT` naming a directory outside the clone | this one is not offline even with the data set: it resolves accessions against a live service before any other switch applies, so a re-run will not reproduce the census byte for byte. `EXPECTED_CENTROIDS=33428` makes it stop if the clustering yields any other centroid count than the recorded census has |

### Tier 4 — not reproducible from this repository, with or without the data set

| stage | entry point | what is missing |
|---|---|---|
| determinants and structures | `src/sdp_scan.jl`, `build_labels.jl`, `pocket_geometry.jl`, `export_pockets.jl`, `ingest_pockets.jl`, `contrastive_potts.jl`, `plmdca_adapter.jl` | **no command-line driver ships in this release.** These are libraries with unit tests. The data set carries their label table, and the structures they read are public and listed by identifier, but there is no shipped program that puts the two together. You can call the modules on your own data; you cannot re-run the analysis they were written for |
| phylogeny | `src/ModulePhyloTools.jl` | it parses, prunes and roots Newick trees built elsewhere. **No script in this repository builds a tree**, so a tree is not reproducible from this repository at all |

### Where the inputs come from

Every input the analysis reads and this repository does not carry lives under **one data
root** in a documented layout. The root is the directory named by `$LNFINDER_DATA`; with
that unset it is `dat/` inside the clone. Unpack the archived data set at the root and every
path the data set carries is satisfied at once — every path but the criteria files, which
are deliberately passed on the command line for the reason given below. What the data set
does not carry, and what each step does without it, is [`docs/DATA.md`](docs/DATA.md) §5.

| variable | what it names |
|---|---|
| `LNFINDER_DATA` | the data root. One variable moves every input |
| `LNFINDER_AFA` | the alignment of the search sequences to the profile, if you keep it outside the root. It is the one input you regenerate yourself, with one HMMER command |
| `LNFINDER_TAXDUMP` | an unpacked NCBI taxonomy dump, if you keep it outside the root. Taxonomy dumps are commonly shared between projects, which is why this one has a variable of its own |
| `LNFINDER_CENSUS_OUT` | where a census run writes. Not an input. Choose a directory **outside** the clone, so a census table and its manifest cannot land among the tracked files |
| `LNFINDER_FIGDIR` | the figure directory of habitat steps 4 and 5. Not an input; `--figdir` overrides it, and step 5 never creates it |
| `LNFINDER_RELATION_SEED` | a semicolon-separated list of directories of already-fetched habitat response bodies, used to prime step 1's relation cache so that those taxa are read from disk rather than fetched. Bare arguments to step 1 do the same thing |

Those six are the whole `LNFINDER_*` set. **Two more switches, carrying no prefix, are read
by the census driver:** `RUN_HABITAT` (`false` turns off its habitat layer, which is its only
per-taxon network route besides accession resolution) and `EXPECTED_CENTROIDS` (a whole
number it checks its centroid count against, so a clustering that silently changed size stops
the run). Eight names in all, every one declared in `src/ModuleDataInputs.jl`, and two tests
fail the build if a script reads a ninth: one scans for the `LNFINDER_` prefix, and one scans
for every environment-variable read by name whatever the name is — because a scan for the
prefix is exactly what could not see those two switches.

The layout under the root, what every file in it is, which stage reads it, how to obtain or
rebuild the things the archive does not carry, and which of them are results are all in
**[`docs/DATA.md`](docs/DATA.md)**. A step that cannot find an input it resolved from the
root refuses before it does any work, naming the input, the place in the layout it looked
and the variable that moves it (a path you typed yourself is reported as the path you
typed, because the layout is not what is wrong with it):

```
LnFinder: the required input "census_habitat" is not there.
  what it is:  the census table carrying the habitat layer, one row per sequence
  read by:     habitat steps 1, 2, 3 and 4
  looked for:  habitat/census_habitat.csv, under the data root
  data root:   $LNFINDER_DATA when that is set; dat/ in the clone when it is not
  This repository does not carry that file;
  docs/DATA.md says which stage reads it and how to obtain or rebuild it.
```

Criteria files are supplied **on the command line** rather than resolved from the root, and
deliberately so. `--gate-targets` (step 1) and `--anchor-targets` (step 4) name files whose
expected values are measurements of the data, so they are not constants in the code, and a
step given neither the file nor an explicit refusal of it (`--no-gates`, `--no-anchors`)
refuses rather than quietly skipping the check. `--anchor-reference` names an earlier run's
rows for step 4's anchor levels; while those levels are not evaluated
(`ANCHOR_LEVELS_RETIRED` in the script) no row is compared, the option is not required, and a
path given is recorded without being opened.

`--categories`, for steps 3 and 5, is refused for a related but distinct reason. Three of
the recorded outcome rules speak for one habitat category each, the figure tabulates one or
more categories across every specification, and both steps quote one phylum's shares —
and **which** categories those are was selected from this data set's own estimates, with the
estimates in hand, after those estimates existed. The classification scheme they are drawn
from is an external ontology's mapping and predates every estimate; picking names out of it
does not. The phylum is the one entry that is **determined rather than chosen** — it is the
heaviest standardisation stratum, and the figure step refuses to assert the identities that
tie the quoted shares to the standardisation weights unless the run's own heaviest stratum is
still that phylum — so it is named in the file for a different reason: a reader has to be able
to derive the key family named after it.
A default in the code would publish that choice as the software's own and could not be told
apart from a selection somebody checked, so steps 3 and 5 hold no default and refuse by
name without it. The names in the file also reach the output keys, derived rather than
written down, so a run on other categories emits keys for those categories. See
[`docs/DATA.md`](docs/DATA.md) §4.

Where each of them comes from: the data set carries `criteria/category_selection.csv`, and
this repository carries the anchor targets at `dat/criteria/anchor_targets.csv`, with the two
tables they were transcribed from beside them. Neither carries a gate-targets file or an
anchor-reference file: step 1 runs on the archived census with `--no-gates` (Tier 3 above),
and step 4 compares no reference row (Tier 2 above).

### "The registered specification" and "the analysis plan"

The code and `docs/TOOLING.md` use both phrases often. **They name an external record, and
this release does not carry it.** There is no plan document in this repository to go looking
for.

**What "registered" claims, exactly.** That record is a dated, checksummed set of rules mapping
each outcome to a colour, a caption phrase and a sentence of the accompanying text. Its
timestamp precedes every output of the re-analysis these steps perform, so the numbers reported
here did not choose their own wording. It was **not** written before the data were first looked
at, and does not claim to have been: it discloses that an earlier run's estimates for the
categories it names, and that run's sensitivity findings, were already known when it was
written. It is registered against the re-analysis, not against the original analysis, and the
categories it speaks for were selected from those earlier estimates rather than fixed in
advance. `docs/TOOLING.md`, under *"The registered specification" and "the analysis plan"*,
says the same in more detail.

What they surface as, in the code you do have, is: the `primary` configuration each script
runs under (its default option set, which every script prints and records), and the criteria
files supplied per run, each of which states inside itself where its values came from — the
anchor targets name the independent implementation each value was transcribed from and that
file's digest; a gate-targets file carries a note per gate saying where its number came from;
and `--categories` carries a required note per row recording that the categories were
selected from the analysis's own earlier estimates rather than fixed in advance.
Where a script says "see THE RECORDED CRITERIA below" (step 3 does), that section of the same
file says which bounds the plan fixed and how a run reports them against the measured values;
the bounds themselves are the two `--plan-*` values, which arrive on the command line and are
not part of this release ([`docs/DATA.md`](docs/DATA.md) §4).

### The archived data set

**It is not published yet, and there is no DOI to give you.** Several of its files are
results of the analysis this software was written for, so it is released with the
accompanying publication and not before. What exists today is the other half of the work:
the layout is fixed, the code reads it, [`docs/DATA.md`](docs/DATA.md) documents it file by
file, and [`docs/data_manifest.sha256`](docs/data_manifest.sha256) already carries a digest
for every archived file **except one named exception**, whose reason is in
[`docs/DATA.md`](docs/DATA.md) §2. Two of the habitat route's inputs carry no digest at all,
in that file or in the archive's own manifest, because they are not in the archive: the stored
relation responses, one per taxon, and the ontology path cache beside them. Both are
third-party content rather than this project's output, and neither is redistributed here. What
is recorded in their place is the digest of the cache's `manifest.csv`, which says which taxa
were fetched, by which route and when, and carries none of the responses. What the shipped
record buys is the thing an archive's own manifest cannot: a **second, independent** record,
so a data root can be checked against this repository and not only against itself — and when
both records are
present, `verify_data.jl` checks the archive's, then cross-checks this one and fails if the
two disagree:

```
julia scripts/verify_data.jl --root <path>
```

When the archive is published, one line of `docs/DATA.md` §2 gains its identifier and
`CITATION.cff` gains the same one. Nothing else about the route changes — which is the point
of writing it down now.

Until then: the things the archive would carry cannot be rebuilt from public databases in
every case, and this file does not pretend otherwise. `census/firstTry.fa` is the sequence
set an `hmmsearch` returned from a database that has changed since; the clusterings need a
proprietary binary that is not redistributed; the habitat relation cache is a frozen
snapshot of a service that is re-mined, so a fresh fetch is different data rather than the
same data again. That last one is why the habitat route is the one the archive cannot make
replayable end to end: those responses are third-party content, they are not redistributed,
and re-fetching them does not give them back. What the archive deposits in their place is what
step 1 computed from them — the call table — which is where the reproducible route starts. For
everything else, the archive is the difference between reading this code and running it, and
is why this repository alone is not enough.

---
