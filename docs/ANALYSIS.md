# How the analysis works

The method behind the code. [`README.md`](../README.md) says what runs;
[`docs/DATA.md`](DATA.md) says what each step reads.

---

## Methods

### What this is

LnFinder was built for one analysis: a census of PQQ-dependent alcohol dehydrogenases of the
methanol/ethanol family, a lanthanide or calcium call and a taxon for every census sequence,
and a comparison of where the lanthanide-type carriers live against comparable prokaryotes.

Every measured value below is in a file this repository or the archived data set carries, is
written by the step named with it, or is stated by the cited paper. Every parameter is in the
code, except the criteria supplied on the command line
([Where the inputs come from](REPRODUCING.md#where-the-inputs-come-from)). Each subsection
ends with its procedure in [`docs/TOOLING.md`](TOOLING.md). [`docs/DATA.md`](DATA.md)
describes the data set; its §2 says the habitat and figure part is published and the rest is
not.

### Census

`census/firstTry.fa` in the data set holds 59,175 records: two calibration references, then
the 59,173 hits HMMER `hmmsearch` retained against the profile HMM this repository ships
(`dat/models/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm`). The references are its only two full
UniProt `sp|` records. All 59,175 are aligned to that profile with `hmmalign`
([`docs/DATA.md`](DATA.md) §5 gives the command). The census driver (`src/scriptLn.jl`) reads
the alignment, calls every sequence, and keeps the centroids of a USEARCH `-cluster_fast`
clustering at 90 % identity (`census/clusters.uc`): 33,428 centroids, one census row each.
`EXPECTED_CENTROIDS=33428` stops the driver on any other count.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *The census pipeline*.

### The cofactor call

The call is read at match state 319. Every run calibrates that column from the first two
sequences of its alignment, a known lanthanide binder then a known calcium binder, taking the
first model column where they read `D-x-D` and `D-x-[ATS]`
(`calibrate_active_site`, `src/ModuleHmmer.jl`). The census was calibrated on C5B120 (XoxF1,
*Methylorubrum extorquens* AM1, lanthanide) and P16027 (MxaF, same organism, calcium); the
control-set validation on I0JWN7 (XoxF, *Methylacidiphilum fumariolicum* SolV, lanthanide) and
P16027. Both give column 319, in each the only column satisfying both patterns.

`parse_hmmalign_results` then reads model columns 319, 320 and 321, insertions removed:

| model columns 319–321 | class |
|---|---|
| a gap in any of the three | `Gap` |
| Asp at 319 and Asp at 321 | `Lanthanide` |
| Asp at 319 and Ala, Thr or Ser at 321 | `Calcium` |
| anything else | `Unknown` |

Column 320 is read only for a gap; a sequence whose model coordinates end before 321 is
`Truncated`.

The rule rests on Good et al. 2020, *J. Biol. Chem.* 295:8272–8284,
[doi:10.1074/jbc.RA120.013227](https://doi.org/10.1074/jbc.RA120.013227). In their XoxF1
structures from *M. extorquens* AM1 the lanthanide is coordinated by Glu192, Asn276, Asp318
and Asp320, and they report that nearly all XoxF and ExaF sequences carry that extra
aspartate two positions downstream of a conserved one, where MxaFI-type methanol
dehydrogenases carry Ala and ExaA-type ethanol dehydrogenases Ser or Thr. Asp318 and Asp320
are the `D-x-D` pair.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *The equivalence being applied* and
*Calibrating the active-site column*.

### Validation

The call is validated against `dat/control_set/known_cofactor_mdhs.tsv`: 29 rows by evidence
tier, 13 tier A, 5 tier B, 7 tier C, 1 contested (tier X), 3 reference rows, each with its
citation. The recorded run (`ValidateClassifier.validate`, Tier 3 below) fetched all 29
sequences, 25 from UniProtKB, 3 from UniParc for deleted UniProtKB entries, 1 from the NCBI
protein record, aligned them, and called every row. Its calls ship as
`dat/control_set/control_validation_predictions.csv`.

- Held out: 25 of 25 correct. Lanthanide 10 of 10, calcium 15 of 15; tier A 13 of 13, tier B
  5 of 5, tier C 7 of 7. Exact (Clopper–Pearson) 95 % interval [0.863, 1.000].
- n_eff = 22: the 25 held-out sequences form 22 single-linkage clusters at 90 % or more
  pairwise identity over the match columns (the harness's default threshold). 22 is the size
  that belongs beside that interval.
- Reference stratum, reported apart: 3 of 3. Its rows (role `reference`) are outside every
  held-out stratum. I0JWN7 and P16027 calibrate the run; C5B120 is the census's lanthanide
  reference.
- CBE67239, tier X, labelled calcium, is called lanthanide (`DYD`), listed and never scored.

The tiers, quoted from the file's header:

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

Every centroid is resolved to an NCBI TaxID through UniProt in three passes (`resolve_taxids`,
`src/ModuleAccessionMapper.jl`), issued one at a time and paced; UniParc's `commonTaxonId` is
never used as a fallback. All 33,428 centroids resolved: 9,130 by UniProtKB search, 21,533 by
a single-organism UniParc record that lists the accession, 2,765 by the accession's own
cross-reference on a record shared between organisms (`taxid_source` in the data set's
`census/final_deduplicated_data.csv`). Lineages come from an NCBI taxonomy dump
([`docs/DATA.md`](DATA.md) §5).

**UniProt release 2026_03 is the release the taxonomy was resolved against, not the release
the hits were drawn from, which is not recorded.** Every census row records it
(`uniprot_release`); 24,298 of the 33,428 accessions are inactive in it (`entry_status`).

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *Taxonomy*.

### Habitat

Habitat comes from OmniCrobe, whose relations carry OntoBiotope (OBT) terms; each term is
climbed to its OBT root path, and the path decides the category. Step 1
(`scripts/habitat_reaggregate.jl`) computes one call per taxon per rule configuration from a
stored cache of raw relation responses, so a run holding that cache replays the same bytes.

**The data set carries the call table step 1 wrote, not the responses it read**: those are a
third-party service's output, not redistributed, so a reader enters the habitat route at the
call table (see
[What is and is not reproducible here](REPRODUCING.md#what-is-and-is-not-reproducible-here)).
The cache's manifest is carried instead: one row per taxon, 13,133 taxa, each with its
relation count, fetch route and fetch time, every row stamped `OmniCrobe 1.7` and
`OBT 2022-04-28`. Both strings are constants typed into `scripts/habitat_reaggregate.jl`
(`OMNICROBE_RELEASE`, `OBT_RELEASE`); no service response confirms them, and no release date
is recorded.

Step 3 (`scripts/habitat_enrichment.jl`) builds the arms; its arm summary
(`habitat_arm_summary.csv`) records every count below. The unit is the taxon.

- Foreground, n1 = 2,397 taxa: census taxa with at least one `Lanthanide` centroid, Bacteria
  or Archaea, with a phylum and a usable habitat annotation.
- Background, n2 = 2,740 taxa: taxa of the UniProt reference-proteome pool (step 0) outside
  the census, less the 1,339 the carrier scan (step 2) flags as carrying a family member at
  the calibrated threshold, with a usable habitat annotation.
- Taxa without a phylum are in no arm; step 3 stops if an assembled arm holds one.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *Habitat*, *The habitat call table* and
*The comparison arms*.

### Estimator

Step 3 reports two estimators of the log2 risk ratio of foreground against background, always
together (`src/ModuleHabitatStats.jl`).

- Crude: log2((a/n1)/(b/n2)), with a Katz interval on the log scale and a Fisher exact
  p-value (a fractional hit total is rounded for the test only). A Haldane correction, +0.5
  per hit total and +1 per denominator, applies only when a hit total, a or b, is zero, and
  only here.
- Phylum-standardised, the primary: p0_std = Σ_k w_k (c_k / n0_k) with w_k = n1_k / n1, the
  background's share of the category re-weighted to the foreground's phylum mix. It has no
  zero-cell correction: where it cannot be computed it reports no estimate (undecidable). Its
  interval decides every robustness class and colour.

Intervals are the estimate ± 1.96 standard errors, 1.96 a literal constant in the code.
Fisher exact p-values are two-sided by the central method, twice the smaller one-sided p-value
capped at 1, not the minimum-likelihood rule; that is the `HypothesisTests.jl` default for
`FisherExactTest`, pinned in `Manifest.toml` (HypothesisTests 0.11.6).

A phylum keeps its own stratum only with at least 10 foreground and 20 background taxa
(`MIN_FG_K`, `MIN_BG_K`); every other phylum becomes "Other phyla" in both arms before
counting, and a stratum with foreground but no background taxa is dropped. Collapsing thin
strata changes the standardised Thermal estimate materially; the threshold asymmetry does
not, and the estimate is equal at 10 and 20 (the shipped pair), 10 and 10, 20 and 20, and 20
and 10. The shipped pair's value is step 4's, in `habitat_robustness.csv`. The other three
pairs and the uncollapsed estimate come from the estimator's `min_fg_k` and `min_bg_k`
arguments, are in no deposited file, and are not reproduced here.

Benjamini–Hochberg runs within each contrast × unit panel, the 10 habitat categories or the 5
chemistry flags, on the crude Fisher exact p-values, over rows with hit total a + b at least 5
(`MIN_CELL`), never across the whole family of tests. Its q-values decide no colour.

The standardised variance is a delta-method form: both hit totals binomial, stratum weights
fixed, units independent. For fractional (tie-split) hits it is an upper bound, so the
interval is conservative there. It is not conservative where a contributing stratum's
background cell is empty, since such a stratum adds nothing to the variance, or where one is
thin. Every row carries the number of contributing strata with an empty background cell and
the smallest background cell among them.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *Enrichment estimators* and *Two estimators,
always reported together*.

### Sensitivity

Step 4 (`scripts/habitat_sensitivity.jl`) varies each rule choice one at a time off the
registered specification: 14 variants over six axes.

| axis | variants |
|---|---|
| type rule | `T:fractional`, `T:alphabetical`, `T:unresolved`, `T:share`, `T:any_type` |
| chemistry rule | `C:share`, `C:winner_fractional` |
| chemistry text | `X:text_first`, `X:text_object` |
| generic terms | `G:generic_none` |
| isolation | `I:iso_taxon` |
| background exclusion | `B:ipr017512`, `B:none`, `B:hmmsearch_movein` |

`disclosed:B:lnca` and `disclosed:B:lntype` are two further background lines, disclosed beside
the grid and outside the set; the unit checks (genus, per-sequence, well-supported subset) are
reported without feeding the classes.

A category is signed, and only then coloured, when the primary standardised interval excludes
zero and every variant's does so on the same side; a variant whose interval cannot be computed
counts as breaking. A primary interval that excludes zero with a variant that does not is
`fragile`; one that includes zero is `null`, whatever the variants say; a category with no
primary interval is `undecidable`. A variant whose axis cannot reach a unit reproduces the
primary there, so the table also records how many variants moved each estimate
(`n_variants_that_moved_the_estimate`).

Step 4 writes one row per category to `habitat_robustness.csv`: the primary standardised
interval, each variant's interval, the count of variants that moved the estimate, and the
class in `robust_sign`, one of `+`, `-`, `fragile`, `null` and `undecidable`. That file is a
result, in neither this repository nor the archived data set, so its values are not reproduced
here; run step 4 to obtain them. Where a category's breaking variants include `C:share` or
`C:winner_fractional`, those levels take their rule definitions from the implementation (see
[Reproducibility](#reproducibility)).

Step 3 runs the same estimator on further contrasts, among them census taxa carrying a family
member but no lanthanide-class centroid (`non-Ln PQQ-ADH carriers vs background`), which step
5 draws as a tick beside each carrier dot. It separates a category specific to the
lanthanide-type carriers from one that is family-wide.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *Sensitivity grid and the robustness flag*.

### Coverage

Habitat annotation does not cover the two arms equally: 65.98 % of foreground candidate taxa
have a usable annotation (2,397 of 3,633), against 87.62 % of the analysed background (2,740
of 3,127, after the 1,339 scan-flagged taxa are removed from the 4,466-taxon frame). Both are
in step 3's arm summary (`fg_coverage_pct`, `bg_analysed_coverage_pct`). Every estimate above
is therefore conditional on annotation, and says nothing about taxa for which OmniCrobe holds
no usable relation.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *The comparison arms*.

### Reproducibility

`dat/criteria/anchor_targets.csv` ships with this repository and holds 630 targets: 42 on the
registered specification, 588 on its 14 variants. They cover six categories, Terrestrial,
Host-plant and Host-animal in unit `taxa` and Alkaline, Hypersaline and Thermal in unit
`taxa (chemistry)`, with seven quantities each (n1, n2, a, b, and the standardised estimate
and its two bounds) and a tolerance per unit and quantity. **No step-4 run publishes unless
its grid reproduces all 630:** `publish!` in `scripts/habitat_sensitivity.jl` moves nothing
into the figure directory unless every gate in `PUBLISH_GATES` passes, its `anchor targets`
gate among them. The recorded run met all 630. Skipping the check has to be typed
(`--no-anchors`), the run records that it was skipped, and step 5 refuses to draw from it.

The targets are transcribed from independent plain-Python reimplementations of the estimates,
written by AI coding agents from a written specification without reading the analysis code.
Each specification was implemented three times: the standard library alone, pandas with
numpy, and numpy.

- The 42 primary cells come from `impl_4.csv`, the standard-library implementation of the
  registered specification; all three implementations of it agree on all 42 cells to 2.3e-13.
- The 588 variant cells come from `impl_w1.csv`, the standard-library implementation of the
  variants; all three agree on all 588 cells to 3.3e-12, and a variant row is written only
  where at least two of the three agree to 1e-9 on every one of the 588 cells.

Both source tables ship beside the targets file, which declares each one's sha256; a step-4
run checks both digests before publishing (its `target sources` gate). A run without those
sources is not refused: it records each source as not beside the targets file.

The analysis plan does not determine every rule the targets depend on; three come from this
repository's implementation. The thin-stratum thresholds 10 and 20, shared by every level,
which the plan names without stating the collapse rule that uses them or valuing both. The
rule definitions of four variant levels, `T:alphabetical`, `T:unresolved`, `C:share` and
`C:winner_fractional`. One edge case of a fifth, `I:iso_taxon`, affecting two taxa. For those
three the anchors check the wiring and the arithmetic, not the rule choice.
[`docs/DATA.md`](DATA.md) lists the same three.

*Detail:* [`docs/TOOLING.md`](TOOLING.md), *Sensitivity grid and the robustness flag* and
*Where the numbers come from, and where they do not*.

---
