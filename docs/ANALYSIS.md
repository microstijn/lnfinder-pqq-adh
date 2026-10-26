# How the analysis works

The method behind the code, moved out of [`README.md`](../README.md) to keep the landing
page short. Nothing here is a command to run: [`README.md`](../README.md) says what runs,
[`docs/DATA.md`](DATA.md) says what each step reads.

---

## Methods

### What this is

LnFinder was built for one analysis: a census of PQQ-dependent alcohol dehydrogenases of the
methanol/ethanol family, drawn from the retained hits of a profile-HMM search; a lanthanide
or calcium call for every census sequence, read at one calibrated column of that profile and
validated against enzymes whose cofactor is known; a taxon for every census centroid; and a
comparison of where the lanthanide-type carriers live against a background of comparable
prokaryotes, with a phylum-standardised estimator and a one-at-a-time sensitivity grid around
every rule choice. This section is the record of what that analysis computed: each method in
a few sentences, the parameter values the recorded runs used, and the values they produced.
Every measured value in it is in a file this repository or the archived data set carries, is
written by one of this repository's steps (the step is named), or is stated by the cited
paper; every parameter of this repository's own steps is in its code, except the criteria
supplied on the command line (see [Where the inputs come from](REPRODUCING.md#where-the-inputs-come-from)).
The procedure behind each subsection is set out in [`docs/TOOLING.md`](TOOLING.md),
under the section title given at its end; the data set is described in
[`docs/DATA.md`](DATA.md); the habitat and figure part of it is published, the rest is
not (see [`docs/DATA.md`](DATA.md) §2).

### Census

The census starts from `census/firstTry.fa` in the data set: **59,175 records**, of which the
first two are the calibration references and the remaining **59,173** are the hits HMMER
`hmmsearch` retained against the profile HMM this repository ships
(`dat/models/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm`). The two references are first because
the calibration takes the first two sequences of the alignment, and they are the file's only
two full UniProt `sp|` records; every other record carries a bare identifier. All 59,175 are
aligned to the same profile with `hmmalign` ([`docs/DATA.md`](DATA.md) §5 gives the
command). The census driver (`src/scriptLn.jl`) reads that alignment, calls every aligned
sequence (next subsection), and keeps the centroids of a USEARCH `-cluster_fast` clustering
at 90 % identity, `census/clusters.uc` in the data set: **33,428 centroids, one census row
each**. Rebuilding with `EXPECTED_CENTROIDS=33428`
stops the driver if the clustering yields any other count.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *The census pipeline*.

### The cofactor call

The call is read at **match state 319 of the profile**. The column is calibrated, not written
into the code: every run takes the first two sequences of its alignment — an enzyme known to
bind a lanthanide, then one known to bind calcium — and returns the first model column at
which the first reads `D-x-D` and the second `D-x-[ATS]` (`calibrate_active_site`,
`src/ModuleHmmer.jl`). The census was calibrated on **C5B120** (XoxF1, *Methylorubrum
extorquens* AM1, lanthanide) and **P16027** (MxaF, the same organism, calcium), the first two
records of `census/firstTry.fa`; the control-set validation on **I0JWN7** (XoxF,
*Methylacidiphilum fumariolicum* SolV, lanthanide) and P16027. Both calibrations give column
319, and in each it is the only column that satisfies both patterns.

The rule, exactly as `parse_hmmalign_results` applies it to model columns 319, 320 and 321 of
a sequence, insertions removed:

| model columns 319–321 | class |
|---|---|
| a gap in any of the three | `Gap` |
| Asp at 319 and Asp at 321 | `Lanthanide` |
| Asp at 319 and Ala, Thr or Ser at 321 | `Calcium` |
| anything else | `Unknown` |

Column 320 is read only for a gap, and a sequence whose model coordinates end before column
321 is `Truncated`.

The rule rests on Good et al. 2020, *J. Biol. Chem.* 295:8272–8284,
[doi:10.1074/jbc.RA120.013227](https://doi.org/10.1074/jbc.RA120.013227). In their crystal
structures of XoxF1 from *M. extorquens* AM1 the lanthanide is coordinated by Glu192, Asn276,
Asp318 and Asp320, and they report that nearly all XoxF and ExaF sequences carry that
additional aspartate two positions downstream of a conserved one, where MxaFI-type methanol
dehydrogenases carry Ala and ExaA-type ethanol dehydrogenases Ser or Thr. Asp318 and Asp320
are the `D-x-D` pair the calibration looks for.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *The equivalence being applied* and
*Calibrating the active-site column*.

### Validation

The call is validated against `dat/control_set/known_cofactor_mdhs.tsv`: **29 rows, scored by
evidence tier** — 13 tier A, 5 tier B, 7 tier C, 1 contested (tier X) and 3 reference rows —
each carrying the citation it rests on. Only tiers A and B rest on a determination made for
that protein or that gene; tier C inherits the cofactor from close relatives and is kept as a
stratum of its own for that reason; the tier definitions are quoted from the file's own
header below. The recorded run (`ValidateClassifier.validate`, Tier 3 below) fetched all 29
sequences — 25 from UniProtKB, 3 from UniParc for UniProtKB entries that have been deleted, 1
from the NCBI protein record — aligned them to the profile, calibrated on I0JWN7 and P16027,
and called every row. Its call for every row ships as
`dat/control_set/control_validation_predictions.csv`.

- **Held out: 25 of 25 correct** — lanthanide 10 of 10, calcium 15 of 15; tier A 13 of 13,
  tier B 5 of 5, tier C 7 of 7 — with an exact (Clopper–Pearson) 95 % interval of
  [0.863, 1.000].
- **Effective size: n_eff = 22.** The 25 held-out sequences form 22 single-linkage clusters
  at 90 % or more pairwise identity over the match columns (the harness's default
  threshold), so 22, not 25, is the size that belongs beside that interval.
- **The reference stratum is reported apart: 3 of 3.** The three rows with role `reference`
  are outside every held-out stratum. I0JWN7 and P16027 calibrate the run; C5B120 is the
  census's lanthanide reference, and is kept out of the held-out set for that reason.
- **One contested row is listed and never scored.** CBE67239, tier X, labelled calcium, is
  called lanthanide (`DYD`).

The evidence tiers, quoted from the file's header:

```
tier, strongest evidence first
  ref  A reference row (see above). Belongs to no held-out stratum.
  A    Direct determination for this protein: the metal resolved in a crystal structure
       of the enzyme itself, or measured in (or required by) the purified enzyme.
  B    Indirect but gene specific: expression, deletion mutant or growth evidence ties
       this gene to a Ca- or Ln-dependent activity, but no metal was measured in,
       reconstituted into or shown to be required by the purified enzyme, or the enzyme
       was not purified at all.
  C    Homology only. The cofactor is inherited from close relatives and was never
       measured for this protein. Kept as a separate stratum.
  X    Contested. Reported together with its call and excluded from every accuracy.
```

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *Control-set validation*.

### Taxonomy

Every census centroid is resolved to an NCBI TaxID through UniProt in three passes
(`resolve_taxids`, `src/ModuleAccessionMapper.jl`), issued one at a time and paced: a
UniProtKB search gives a live entry its organism and hands back the UniParc record of a
deleted one; a UniParc record that names a single organism, and lists the accession, answers
a deleted entry; and a record shared between organisms is resolved through the accession's
own cross-reference on it. UniParc's `commonTaxonId` is never used as a fallback.
**33,428 of 33,428 centroids resolved, none unresolved**: 9,130 through the UniProtKB search,
21,533 through a single-organism UniParc record and 2,765 through the accession's own
cross-reference (the `taxid_source` column of `census/final_deduplicated_data.csv` in the
data set). The lineage of each TaxID is then read from an NCBI taxonomy dump
([`docs/DATA.md`](DATA.md) §5).

**UniProt release 2026_03 is the release the taxonomy was resolved against, not the one the
hits were drawn from.** Every census row records it (`uniprot_release`), and 24,298 of the
33,428 accessions are inactive in it (`entry_status`). The release the search read is not
recorded.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *Taxonomy*.

### Habitat

Habitat comes from OmniCrobe, whose relations carry OntoBiotope (OBT) terms; each term is
climbed to its OBT root path, and the path decides the category. Step 1
(`scripts/habitat_reaggregate.jl`) computes one call per taxon per rule configuration from a
stored cache of raw relation responses, so a run that holds that cache replays the same bytes
rather than fetching again. **The data set carries the call table step 1 wrote, not the
responses it read.** Those responses are a third-party service's own output and are not
redistributed with it, which is why a reader of the data set enters the habitat
route at the call table rather than at step 1, whose inputs the data set does not carry
(see [What is and is not reproducible here](REPRODUCING.md#what-is-and-is-not-reproducible-here)). What is carried is the cache's manifest — one
row per taxon, 13,133 taxa — recording how many relations each taxon had, by which route it
was fetched and when, and none of the responses themselves. **What is recorded as the service
releases:** every row of that manifest is stamped `OmniCrobe 1.7` and
`OBT 2022-04-28`. Both strings are constants in `scripts/habitat_reaggregate.jl`
(`OMNICROBE_RELEASE`, `OBT_RELEASE`), typed into the code rather than read from a service
response, and no release date is recorded with them.

**The arms** are built by step 3 (`scripts/habitat_enrichment.jl`), and its arm summary
(`habitat_arm_summary.csv`) records every count below. **The unit is the taxon.**

- **Foreground: n1 = 2,397 taxa** — census taxa that carry at least one centroid of class
  `Lanthanide`, are Bacteria or Archaea, have a phylum, and have a usable habitat
  annotation.
- **Background: n2 = 2,740 taxa** — taxa of the UniProt reference-proteome pool (step 0)
  that are not in the census, less the 1,339 the carrier scan (step 2) flags as carrying a
  family member at the calibrated threshold, with a usable habitat annotation.
- **Taxa without a phylum are in no arm.** They are excluded when the arms are assembled,
  and step 3 stops if an assembled arm still holds one.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *Habitat*, *The habitat call table* and
*The comparison arms*.

### Estimator

Step 3 reports two estimators of the log2 risk ratio of the foreground against the
background, always together (`src/ModuleHabitatStats.jl`). Where a Fisher exact p-value
appears below it is two-sided by the **central** method — twice the smaller of the two
one-sided p-values, capped at 1 — which is what `HypothesisTests.jl` returns for a
`FisherExactTest` when no method is asked for, and is not the minimum-likelihood rule some
other packages default to. The version is pinned in `Manifest.toml` (HypothesisTests 0.11.6).

- **Crude:** log2((a/n1)/(b/n2)), with a Katz interval on the log scale and a Fisher exact
  p-value (a fractional hit total is rounded for the test only). A Haldane correction —
  +0.5 on each hit total, +1 on each denominator — is applied only when a hit total, a or b, is zero,
  and only in this estimator.
- **Phylum-standardised, the primary:** the background's share of the category re-weighted
  to the foreground's phylum mix, p0_std = Σ_k w_k (c_k / n0_k) with w_k = n1_k / n1. It
  has no zero-cell correction: where it cannot be computed it reports no estimate
  (undecidable) rather than a corrected one. Its interval is the one that decides every
  robustness class and every colour.

**Thin strata.** A phylum keeps a stratum of its own only if it has at least 10 taxa in the
foreground and at least 20 in the background (`MIN_FG_K`, `MIN_BG_K`); every other phylum is
relabelled "Other phyla" in both arms before anything is counted, and a stratum with
foreground taxa but no background taxa is dropped. Collapsing thin strata changes the
standardised Thermal estimate materially; the asymmetry of the two thresholds does not, and
the estimate is the same at each of 10 and 20 (the shipped pair), 10 and 10, 20 and 20, and
20 and 10. The shipped pair's value is step 4's, in `habitat_robustness.csv`. The other three
pairs and the uncollapsed estimate are probes of the same arms through the estimator's
`min_fg_k` and `min_bg_k` arguments (`src/ModuleHabitatStats.jl`); no file in the data set
carries them, and they are not reproduced here.

**Intervals** are the estimate ± 1.96 standard errors in both estimators, 1.96 being a
literal constant in the code.

**Multiple testing.** Benjamini–Hochberg runs **within each contrast × unit panel** — the 10
habitat categories of a habitat panel, or the 5 chemistry flags of a chemistry panel — on the
crude estimator's Fisher exact p-values, over the rows whose hit total a + b is at least 5
(`MIN_CELL`). It is not applied across the whole family of tests. The q-values belong to the
crude test and decide no colour.

**Variance.** The standardised variance is a delta-method form that treats both hit totals
as binomial and the stratum weights as fixed. For fractional (tie-split) hits it is an upper
bound under that model, so the interval is conservative there. It is not conservative where
a contributing stratum's background cell is empty — such a stratum adds nothing to the
variance — or thin, and it assumes the units are independent. Every row therefore also
carries how many contributing strata have an empty background cell, and the smallest
background cell among the contributing strata.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *Enrichment estimators* and *Two estimators,
always reported together*.

### Sensitivity and robustness

Step 4 (`scripts/habitat_sensitivity.jl`) varies each rule choice one at a time off the
registered specification: **14 variants over six axes**.

| axis | variants |
|---|---|
| type rule | `T:fractional`, `T:alphabetical`, `T:unresolved`, `T:share`, `T:any_type` |
| chemistry rule | `C:share`, `C:winner_fractional` |
| chemistry text | `X:text_first`, `X:text_object` |
| generic terms | `G:generic_none` |
| isolation | `I:iso_taxon` |
| background exclusion | `B:ipr017512`, `B:none`, `B:hmmsearch_movein` |

Two further background lines, `disclosed:B:lnca` and `disclosed:B:lntype`, are disclosed beside the grid and are
outside the set, and the unit checks (genus, per-sequence, well-supported subset) are
reported without feeding the classes.

**The rule.** A category is given a sign, and only then coloured, if the primary standardised
interval excludes zero **and every variant's interval excludes zero on the same side**. A
variant whose interval cannot be computed counts as breaking. A primary interval that
excludes zero with a variant that does not is `fragile`; a primary interval that includes
zero is `null`, whatever the variants say; and a category the primary gives no interval at
all is `undecidable`. A variant whose axis cannot reach a unit reproduces the primary there by
construction, so the robustness table also records how many variants moved each estimate
(`n_variants_that_moved_the_estimate`).

**The classes.** Step 4 writes one row per category to `habitat_robustness.csv`, carrying the
primary standardised interval, each variant's interval, the count of variants that moved the
estimate, and the resulting class in its `robust_sign` column. That file is a result: it is
not tracked in this repository and it is not in the archived data set, so its values are not
reproduced in this document. Run step 4 to obtain them.

The table holds one row per category. `robust_sign` takes five values — `+`, `-`, `fragile`,
`null` and `undecidable` — defined by the rule above. Where a category's breaking variants
include `C:share` or `C:winner_fractional`, note that those two levels have rule definitions
taken from the implementation rather than from the registered specification (see
[Reproducibility](#reproducibility)).

**The comparator inside the family.** Step 3 runs the same estimator on further contrasts
beside the lanthanide-type carriers, among them the census taxa that carry a family member
but no lanthanide-class centroid (`non-Ln PQQ-ADH carriers vs background`), which step 5
draws as a tick beside each carrier dot. The comparator is what separates a category that is
specific to the lanthanide-type carriers from one that belongs to the enzyme family as a
whole; step 5 draws both, so the comparison is read off the figure rather than asserted here.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *Sensitivity grid and the robustness flag*.

### Coverage

Habitat annotation does not cover the two arms equally. In the foreground, 65.98 % of the
candidate taxa have a usable annotation (2,397 of 3,633); in the analysed background, 87.62 %
do (2,740 of 3,127, after the 1,339 scan-flagged taxa are removed from the 4,466-taxon
frame). Both are in step 3's arm summary (`fg_coverage_pct`, `bg_analysed_coverage_pct`).
Every estimate above is therefore conditional on annotation: it compares annotated
lanthanide-type carriers with annotated background taxa, and says nothing about taxa for
which OmniCrobe holds no usable relation.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *The comparison arms*.

### Reproducibility

**The anchor targets.** `dat/criteria/anchor_targets.csv` ships with this repository and holds
**630 targets: 42 on the registered specification and 588 on its 14 variants**. They cover
six categories — Terrestrial, Host-plant and Host-animal in unit `taxa`; Alkaline,
Hypersaline and Thermal in unit `taxa (chemistry)` — with seven quantities each (n1, n2, a,
b, and the standardised estimate and its two bounds) and a tolerance per unit and quantity. **Every
step-4 run checks its grid against all of them and refuses to publish if any target is not
reproduced:** `publish!` in `scripts/habitat_sensitivity.jl` moves nothing into the figure
directory unless every gate in `PUBLISH_GATES` passes, and its `anchor targets` gate fails on
any target missed. The recorded run met all 630. Skipping the check has to be typed
(`--no-anchors`); the run then records that it was skipped, and step 5 refuses to draw from
it.

**Where the targets come from.** They are transcribed from independent plain-Python
reimplementations of the estimates, **written by AI coding agents** from a written
specification, without reading the analysis code. Each specification was implemented three
times on three toolchains: the standard library alone, pandas with numpy, and numpy. The two
halves of the file rest on separate claims.

- The **42 primary cells** are transcribed from `impl_4.csv`, the standard-library
  implementation of the registered specification. All three implementations of that
  specification agree on all 42 cells to 2.3e-13.
- The **588 variant cells** are transcribed from `impl_w1.csv`, the standard-library
  implementation of the variants. All three agree on all 588 cells to 3.3e-12, and a variant
  row is written only where at least two of the three agree to 1e-9 on every one of the 588
  cells — one implementation on its own is not enough to put a variant target in the file.

Both source tables ship beside the targets file, which declares each one's sha256. A step-4
run on the shipped file checks both digests before it publishes (its `target sources` gate).
A targets file run without its sources beside it is not refused: the run records each source
as not beside the targets file.

**What the specification does not determine.** The analysis plan does not determine every
rule the targets depend on. **Three** of them are taken from this repository's implementation
instead: the thin-stratum threshold values, 10 and 20, which every level shares, because the
plan names those thresholds without stating the collapse rule that uses them or valuing both;
the rule definitions of four variant levels — `T:alphabetical`, `T:unresolved`, `C:share` and
`C:winner_fractional`; and one edge case of a fifth, `I:iso_taxon`, which affects two taxa.
For those three the anchors check the wiring and the arithmetic, not the rule choice.
[`docs/DATA.md`](DATA.md) lists the same three beside the targets file it describes.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *Sensitivity grid and the robustness flag* and
*Where the numbers come from, and where they do not*.

---
