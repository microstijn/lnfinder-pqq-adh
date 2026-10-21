"""
prevalence_phylum_average.py -- the hit rate with every phylum given equal say.

THE PROBLEM IT FIXES. The pooled hit rate is dominated by Pseudomonadota: they are 37 % of the
sampled reference proteomes and carry a lanthanide-type PQQ-ADH at 47 % against 9 % elsewhere,
so they supply three quarters of all carriers. A pooled rate therefore mostly reports how many
Pseudomonadota happen to be in each setting, and a figure split by phylum group makes
Pseudomonadota the subject of the figure. Neither is the question.

WHAT THIS COMPUTES INSTEAD. Within each phylum separately, the hit rate for a setting; then the
unweighted mean across phyla. Each phylum gets one vote, so no phylum's abundance in the
reference set can drive the answer. This is direct standardisation to an equal-phylum reference.

WHY THE PHYLUM SET IS FIXED AND NOT PER-ROW. Averaging over "whichever phyla are big enough in
this row" compares different phylum sets between rows, and it is unstable: at a 20-genome floor
Aquatic averages 13.3 % over seven phyla, at a 30-genome floor 16.0 % over five. Both are
"the average" and they are not comparable to each other or to Host-plant's four. So the set is
FIXED -- the four phyla with at least MIN_CELL sampled genomes in every drawn row -- and a row
that cannot supply all four is DROPPED and named, rather than quietly averaged over fewer.

THE INTERVAL is a bootstrap: each phylum's count is redrawn from its own binomial at the
observed rate, the four rates are averaged, and the 2.5/97.5 percentiles are taken over 2000
replicates. It covers sampling within phyla. It does NOT cover the choice of the four phyla,
which is a structural choice, nor habitat misannotation, nor the reference-proteome frame.
"""
import csv, math, os, random, statistics


# ---------------------------------------------------------------------------
# Paths follow the repository's data contract (docs/DATA.md): every input resolves
# under the data root, and the step refuses before doing any work if one is missing.
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
ROOT = os.path.abspath(os.environ.get("LNFINDER_DATA") or os.path.join(_REPO, "dat"))
CALLS = os.path.join(ROOT, "habitat", "habitat_taxon_calls.csv")
SCAN = os.path.join(ROOT, "habitat", "background_pqq_scan.csv")
HERE = os.path.join(ROOT, "figure")
os.makedirs(HERE, exist_ok=True)
_NL = chr(10)

for _name, _p in (("habitat_taxon_calls.csv", CALLS), ("background_pqq_scan.csv", SCAN)):
    if not os.path.exists(_p):
        raise SystemExit(_NL.join([
            "missing input: " + _name,
            "  looked in : " + _p,
            "  the layout puts it at <data root>/habitat/" + _name,
            "  set $LNFINDER_DATA to the directory holding the archived data set."]))


ISOLATION = {"Food/substance", "Lab-medium"}
PHYLA = ["Pseudomonadota", "Actinomycetota", "Bacillota", "Bacteroidota"]
MIN_CELL = 20
NBOOT = 2000
SEED = 20260928


def truthy(s):
    return str(s).strip().lower() == "true"


def tally():
    scan_ln = {int(r["taxid"]) for r in csv.DictReader(open(SCAN, encoding="utf-8"))
               if "Lanthanide" in str(r.get("hit_classes", ""))}
    hab, chem, overall = {}, {}, {}
    for r in csv.DictReader(open(CALLS, encoding="utf-8")):
        if r["config"] != "primary" or not truthy(r["in_pool"]):
            continue
        ph = str(r["phylum"]).strip()
        ht = str(r["habitat_type"]).strip()
        if not (ph and ht and ht not in ISOLATION):
            continue
        ln = ((truthy(r["in_census"]) and int(r["n_ln_seqs"] or 0) > 0)
              or ((not truthy(r["in_census"])) and int(r["taxid"]) in scan_ln))
        a, b = overall.get(ph, (0, 0))
        overall[ph] = (a + int(ln), b + 1)
        # An unresolved tie ('Host-plant|Terrestrial') belongs to no single setting, so it is
        # left out of the per-phylum cells rather than counted in both.
        if "|" not in ht:
            a, b = hab.get((ht, ph), (0, 0))
            hab[(ht, ph)] = (a + int(ln), b + 1)
        for f in [x for x in str(r["chemistry_flags"]).split("+") if x]:
            a, b = chem.get((f, ph), (0, 0))
            chem[(f, ph)] = (a + int(ln), b + 1)
    return hab, chem, overall


def average_row(cells, rng):
    """Mean of the per-phylum rates, with a bootstrap interval over within-phylum sampling."""
    rates = [100 * k / n for _, k, n in cells]
    point = statistics.fmean(rates)
    draws = []
    for _ in range(NBOOT):
        draws.append(statistics.fmean(
            [100 * rng.binomialvariate(n, k / n) / n for _, k, n in cells]))
    draws.sort()
    return point, draws[int(0.025 * NBOOT)], draws[int(0.975 * NBOOT)]


def geomean(xs):
    return math.exp(statistics.fmean(math.log(max(x, 1e-9)) for x in xs))


def fold_row(cells, base_rate, rng):
    """
    Fold-change against each phylum's OWN overall rate, combined by geometric mean.

    WHY A RATIO AND NOT THE STANDARDISED PERCENTAGE. The equal-phylum percentage is the
    prevalence in a community where these four phyla are equally common, which is nowhere:
    Host-plant's 28.1 % is neither the observed 42.8 % nor any phylum's rate (75/16/12/9).
    A fold-change compares each phylum only against itself, so it needs no invented
    population, and the geometric mean is the right centre for ratios -- an arithmetic mean
    of 0.5x and 2.0x is 1.25x, which claims an increase where there is none.

    THE INTERVAL bootstraps the SETTING cell only, holding each phylum's overall rate fixed.
    The two are nested -- the setting's genomes are part of that phylum's total -- and the
    totals are 10-40x larger (2008, 1165, 923, 625), so their sampling error is the smaller
    term. It is an approximation, and it is the reason these intervals are narrower than a
    fully nested resampling would give.
    """
    # Haldane +0.5 on a ZERO CELL ONLY, the same correction crude() applies, and here it is
    # not cosmetic: Bacillota contributes about 4 carriers of 405 in Host-animal, so a
    # bootstrap draw empties it perfectly often, one zero sends the geometric mean to zero,
    # and the interval's lower bound prints as 0.00 for a row whose point estimate is 0.31.
    def rate(k, n):
        return 100 * (k + 0.5) / (n + 1) if k == 0 else 100 * k / n

    folds = [rate(k, n) / base_rate[p] for p, k, n in cells]
    point = geomean(folds)
    draws = []
    for _ in range(NBOOT):
        draws.append(geomean([rate(rng.binomialvariate(n, k / n), n) / base_rate[p]
                              for p, k, n in cells]))
    draws.sort()
    return point, draws[int(0.025 * NBOOT)], draws[int(0.975 * NBOOT)], folds


def main():
    rng = random.Random(SEED)
    hab, chem, overall = tally()

    out, dropped, rows_pending = [], [], []
    for axis, d in (("habitat", hab), ("chemistry", chem)):
        for setting in sorted({k[0] for k in d}):
            cells = [(p,) + d.get((setting, p), (0, 0)) for p in PHYLA]
            thin = [f"{p} ({n})" for p, _, n in cells if n < MIN_CELL]
            if thin:
                dropped.append((axis, setting, "; ".join(thin)))
                continue
            point, lo, hi = average_row(cells, rng)
            rows_pending.append((axis, setting, cells, point, lo, hi))

    base_cells = [(p,) + overall[p] for p in PHYLA]
    bpoint, blo, bhi = average_row(base_cells, rng)
    base_rate = {p: 100 * k / n for p, k, n in base_cells}

    for axis, setting, cells, point, lo, hi in rows_pending:
        f, flo, fhi, folds = fold_row(cells, base_rate, rng)
        out.append({"axis": axis, "setting": setting, "avg_pct": point,
                    "ci_lo_pct": lo, "ci_hi_pct": hi,
                    "fold": f, "fold_lo": flo, "fold_hi": fhi,
                    **{f"rate_{p}": 100 * k / n for p, k, n in cells},
                    **{f"fold_{p}": fd for (p, _, _), fd in zip(cells, folds)},
                    **{f"n_{p}": n for p, _, n in cells}})

    print(f"EQUAL-PHYLUM AVERAGE over {', '.join(PHYLA)}")
    print(f"(each phylum one vote; cells under {MIN_CELL} sampled genomes are not averaged)\n")
    print(f"overall {bpoint:.1f} % [{blo:.1f}, {bhi:.1f}]   "
          + "  ".join(f"{p[:6]} {100*k/n:.0f}%" for p, k, n in base_cells) + "\n")
    for axis in ("habitat", "chemistry"):
        rows = sorted([r for r in out if r["axis"] == axis],
                      key=lambda r: -r["avg_pct"])
        if not rows:
            continue
        print(axis.upper())
        for r in rows:
            per = "  ".join(f"{p[:6]} {r['rate_'+p]:.0f}%({r['n_'+p]})" for p in PHYLA)
            print(f"  {r['setting']:<15}{r['avg_pct']:>6.1f} % "
                  f"[{r['ci_lo_pct']:.1f}, {r['ci_hi_pct']:.1f}]"
                  f"   fold {r['fold']:.2f}x [{r['fold_lo']:.2f}, {r['fold_hi']:.2f}]"
                  f"   {per}")
        print()
    if dropped:
        print("NOT DRAWN — a row must supply all four phyla:")
        for axis, s, thin in dropped:
            print(f"  {s:<15} ({axis}) thin: {thin}")

    path = os.path.join(HERE, "prevalence_phylum_average.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        cols = (["axis", "setting", "avg_pct", "ci_lo_pct", "ci_hi_pct",
                 "fold", "fold_lo", "fold_hi"]
                + [f"rate_{p}" for p in PHYLA] + [f"fold_{p}" for p in PHYLA]
                + [f"n_{p}" for p in PHYLA])
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out:
            w.writerow(r)
        w.writerow({"axis": "baseline", "setting": "all sampled genomes",
                    "avg_pct": bpoint, "ci_lo_pct": blo, "ci_hi_pct": bhi,
                    "fold": 1.0, "fold_lo": 1.0, "fold_hi": 1.0,
                    **{f"rate_{p}": 100 * k / n for p, k, n in base_cells},
                    **{f"fold_{p}": 1.0 for p in PHYLA},
                    **{f"n_{p}": n for p, _, n in base_cells}})
    print(f"\nwrote {os.path.basename(path)}")


if __name__ == "__main__":
    main()
