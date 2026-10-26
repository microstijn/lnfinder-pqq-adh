# LnFinder

LnFinder finds PQQ-dependent alcohol dehydrogenases in large sequence sets and calls which
metal each one uses, a lanthanide or calcium, from the residues at one calibrated active-site
column. It also carries the analysis of where the lanthanide-type carriers live, and the
figure that reports it.

This is one enzyme family, not lanthanide-binding proteins in general. No LanM, LanD or LanP
homolog is searched anywhere here.

| residues at the calibrated active-site column | call |
|---|---|
| `D-x-D` — Asp, any residue, Asp | lanthanide |
| `D-x-[ATS]` — Asp, any residue, Ala/Thr/Ser | calcium |
| neither, or a gap | unknown |

The rule itself is published: the extra aspartate was proposed by Keltjens et al. 2014,
[doi:10.1007/s00253-014-5766-8](https://doi.org/10.1007/s00253-014-5766-8), resolved in the
crystal structures of Good et al. 2020 (cited in full under
[The cofactor call](#the-cofactor-call)), and applied at scale by Voutsinos et al. 2025,
[doi:10.1093/ismejo/wraf057](https://doi.org/10.1093/ismejo/wraf057). What this repository
adds is the machinery around it:

## What you need

Julia 1.12 or newer, which `Project.toml` pins. Python 3.12 or newer for two steps, standard
library only, no requirements file. HMMER 3 and network access only if you re-run the
acquisition steps; everything else is offline.

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
```

`Manifest.toml` is committed so `instantiate` reproduces the environment instead of resolving
a new one. Run the test file directly; `Pkg.test` is not wired up. The suite is offline and
needs no input the repository does not carry.

## Draw the figure

This needs nothing but a clone. The three tables it reads are tracked under `dat/figure/`.

```
julia --project=. scripts/fold_change_figure.jl notitle
```

It writes PDF, SVG and PNG to `$LNFINDER_FIGDIR`, or beside the tables if that is unset. Drop
`notitle` for the titled version; `pooled` and `twophylum` draw the two thermal variants. The
PDF is not byte-reproducible, because cairo stamps a timestamp inside a compressed stream, so
compare the SVG if you need to check a rebuild.

## Rebuild it from the data

Everything below needs the archived data set under `$LNFINDER_DATA`; see
[`docs/DATA.md`](docs/DATA.md), and check it with `julia scripts/verify_data.jl` before
running anything. The chain, in order:

```
scripts/prokaryote_background.jl        network; draws the reference-proteome pool, fetches habitats
src/scriptLn.jl                         HMMER; the census and the Asp-x-Asp call
scripts/habitat_reaggregate.jl          one habitat call per taxon
scripts/background_pqq_scan.jl          HMMER; carriers the census never saw
scripts/prevalence.py                   carrier prevalence per setting
scripts/prevalence_phylum_average.py    the same, each phylum weighted equally
scripts/fold_change_figure.jl           the figure
```

Two steps take a flag with no default and refuse without it: `habitat_enrichment.jl` needs
`--categories` pointed at `criteria/category_selection.csv` in the data set, and
`habitat_sensitivity.jl` needs `--anchor-targets dat/criteria/anchor_targets.csv`, which ships
here. `habitat_figure.jl` needs both that same `--categories` and a `--figdir` that exists.

The two Python steps stay in Python because their intervals come from a seeded
`random.Random` stream; reimplementing them elsewhere moves the bounds and so the figure.

## Reading further

[`docs/ANALYSIS.md`](docs/ANALYSIS.md) is the method. [`docs/DATA.md`](docs/DATA.md) is what
each step reads and where it expects to find it. [`docs/REPRODUCING.md`](docs/REPRODUCING.md)
says which steps run from what, and where the inputs came from.
[`docs/TOOLING.md`](docs/TOOLING.md) covers the external binaries and services.

## Licence and attribution

The software is **MIT licensed** — see [`LICENSE`](LICENSE) — with one exception, set out in
that file under the heading "EXCEPTIONS TO THE LICENCE ABOVE". In summary, and in this
README's own words rather than the licence's: `dat/models/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm`
is a third-party profile HMM that LnFinder did not build. It is redistributed unmodified and
attributed. Its upstream repository states no licence at all, so no permission was granted for
it and none is granted or implied here. Anyone wishing to reuse it should obtain terms from its
upstream authors. **The licence file is the operative text; read it there.**

Its source, upstream commit, checksums, citation and licence status are recorded in
[`dat/models/PROVENANCE.md`](dat/models/PROVENANCE.md), together with the checksum you can
verify against your own checkout.

`dat/control_set/known_cofactor_mdhs.tsv` is a compilation of published facts: every row
carries the citation the claim comes from, and a few rows are adopted from a third party's
supplementary table and cite it as their source. **That source is named in full in the file's
own header** — author list, title, preprint DOI, the DOI of the supplementary tables, and its
CC BY 4.0 licence — together with a statement of exactly which fields were taken from it and
which were not. The MIT licence covers the compilation; it does not replace those citations.
If you reuse the table, carry them with it.

### Data sources

Nothing in this repository redistributes these; the code reads them, and some test fixtures
are responses recorded from them (see "Run the tests"; one further fixture is a table recorded
from a local `hmmsearch` run, and HMMER is a tool rather than a data source, so it appears under
[External tools and services](#external-tools-and-services) instead of in the table below).
Citing LnFinder does not cite them.
If you use the data, cite each in its own right — each maintains a current citation on its
own site, which is why this list gives endpoints rather than a bibliography frozen at the
date this file was written.

| source | endpoint | what is read |
|---|---|---|
| UniProt | <https://rest.uniprot.org/> | sequences, TaxID resolution through UniProtKB and UniParc, reference proteomes, the InterPro-targeted prefilter |
| OntoBiotope / OmniCrobe | <https://omnicrobe.migale.inrae.fr/> | microbe–habitat relations and the habitat ontology the aggregation climbs |
| NCBI Taxonomy | <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/> | the taxonomy dump the census resolves TaxIDs against |
| NCBI E-utilities | <https://eutils.ncbi.nlm.nih.gov/> | one control-set sequence whose UniProtKB entry was deleted |
| InterPro | <https://www.ebi.ac.uk/interpro/> | the two entry accessions that define the step-2 prefilter, named in `docs/TOOLING.md` under *The frame, and the carriers hiding in it* |

## Citing this work

**The accompanying publication:** citation to be added at acceptance.

[`CITATION.cff`](CITATION.cff) carries the software metadata. **It is deliberately
incomplete**: it names the copyright holder and no one else, and it has no DOI and no release
date, because the accompanying publication is not out and no archived release exists. The
author list, the ORCIDs, the publication's own citation and the release DOI are added there
when the accompanying publication appears — and nothing is guessed in the meantime.

## How this code was written

Stijn Peeters wrote this software together with Claude Code (Anthropic), which drafted,
reviewed and tested a large part of it under his direction. Responsibility for what the code
does, and for every result obtained with it, rests with the named author. The commit history
records the model's contribution in `Co-Authored-By` trailers.
