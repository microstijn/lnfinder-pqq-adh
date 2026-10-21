# Third-party profile HMM

This directory holds one profile HMM that LnFinder did not build. It is redistributed here
unmodified, with attribution, and it is **not** covered by this repository's MIT licence.

## Source

| | |
|---|---|
| Repository | <https://github.com/banfieldlab/metabolic-hmms> |
| Commit | `92addd7d10` (2016-09-23) |
| File | `methanol_dehydrogenase_pqq_xoxF_mxaF.hmm` |
| Upstream git blob | `f50d4a542e156ba46d26ca53ae08b61eeceefc14` |
| SHA-256 | `50882c8390c33ad91a0447a81222cd72877c2f9b954bc41c9295a753c34e465f` |

The copy here is byte-identical to the upstream blob; both hashes above were re-checked
against the GitHub tree API on 2026-09-13.

That digest is also what a checkout of this repository gives you, on every platform.
`.gitattributes` marks `*.hmm` as non-text precisely so that line-ending normalisation
cannot rewrite these bytes on a Windows clone and falsify the checksum that stands behind
an unlicensed redistribution. Check it yourself:

```
sha256sum dat/models/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm
```

or, with no coreutils and nothing but the Julia this repository already needs:

```
julia -e 'using SHA; println(bytes2hex(open(sha256, "dat/models/methanol_dehydrogenase_pqq_xoxF_mxaF.hmm")))'
```

## What the model is

Built by the original authors with HMMER 3.0 `hmmbuild`. Its header records
`HMMER3/b [3.0 | March 2010]`, `NAME pqq_xoxF_mxaF_alignment`, `LENG 589`, `NSEQ 189`,
`DATE Tue Jul 14 13:16:10 2015`.

The model carries **no GA, TC or NC cutoffs**, so it defines no inclusion threshold of its
own: whatever threshold a search applies comes from the search command, not from this file.
The threshold used for a search is a property of that run and is recorded with the run's
own output, not here.

## Citation

Anantharaman K, et al. (2016) Thousands of microbial genomes shed light on interconnected
biogeochemical processes in an aquifer system. Volume 7, article 13219.
doi:[10.1038/ncomms13219](https://doi.org/10.1038/ncomms13219) — "Hmm databases used in this
study are available from <https://github.com/banfieldlab/metabolic-hmms>".

The same group's PQQ-dependent dehydrogenase HMM was reused by Diamond et al. (2019),
doi:[10.1038/s41564-019-0449-y](https://doi.org/10.1038/s41564-019-0449-y), and through
them by Voutsinos et al. (2025),
doi:[10.1093/ismejo/wraf057](https://doi.org/10.1093/ismejo/wraf057). That is reuse of a
model from the same source; it is not a claim that those studies used this identical file.

## Licence

**The upstream repository states no licence.** Re-checked on 2026-09-23: the repository
metadata reports its licence as null, the licence endpoint answers 404, the root holds no
`LICENSE` or equivalent, and the README sets out no terms. So no permission is granted
explicitly and none is implied here. This file is redistributed unmodified and attributed,
and it is excluded from the MIT licence that covers the rest of this repository — see the
"Exceptions" section of `LICENSE`. Anyone reusing it should look to the upstream authors
for terms.

## What is not here, and why

The copy of the upstream repository this one was built from carried 28 profile HMMs, of
which only `methanol_dehydrogenase_pqq_xoxF_mxaF.hmm` is used by LnFinder. (That count is
what was seen at the commit named above; check upstream for what it carries now.) The rest
are referenced by no file here, so they are not redistributed: they are not in this
repository — not in the working tree and not in its history.

The retained model is used by:

- `src/scriptLn.jl` — the census pipeline
- `src/validate_classifier.jl` — the control-set validation of the cofactor call
- `scripts/background_pqq_scan.jl` — the background carrier scan
- `docs/TOOLING.md` — the method documentation, which describes what it is searched with
