"""
prevalence.py -- carrier PREVALENCE per habitat: if you sample here, what fraction of
genomes carry a lanthanide-type PQQ-ADH?

WHY THIS IS A DIFFERENT NUMBER FROM BOTH FIGURES SO FAR. The published figure and the
norm figure both answer "is this habitat over-represented among carriers". Neither answers
"what is the hit rate in this habitat", which is what decides where to prospect. Soil can
hold a quarter of all known carriers and still be a 3 % hit rate; an alkaline lake can hold
few and be 30 %.

THE FRAME, AND WHY IT IS THE ONLY HONEST ONE. Prevalence needs a numerator and a
denominator drawn from the SAME sample. The census is not a sample of anything -- it is a
UniProtKB-wide hmmsearch -- so carriers/(carriers+background) would be biased by however
completely the census caught carriers relative to the pool draw. The 6,000-taxon random
sample of reference proteomes IS a sample, drawn with MersenneTwister(20260810) from the
18,105 prokaryotic reference proteomes, so every rate below is computed INSIDE it:

    denominator  pool taxa with a usable habitat call
    numerator    those of them that carry a lanthanide-type enzyme

CARRIER STATUS COMES FROM TWO DIFFERENT MEASUREMENTS, and the split is reported rather
than blended away. A pool taxon the census registered is called from the census's own
hmmalign classification (n_ln_seqs > 0). A pool taxon the census never registered is called
from the targeted background scan's own model call (hit_classes contains Lanthanide). The
codebase is explicit that these are not the same measurement, so the output carries both
counts and a rate computed from each.

WHAT THIS IS STILL NOT. A rate among SEQUENCED REFERENCE PROTEOMES that have a habitat on
record -- not among organisms in that habitat in nature. Reference proteomes are roughly one
per species group and skew to cultured, studied organisms, and only ~88 % of the pool has a
habitat at all.
"""
import csv, os, math

Z = 1.959963984540054  # 95 %


def wilson(k, n):
    """
    Approximate 95 % Wilson interval on a weighted count.

    APPROXIMATE, AND HERE IS WHY. Taxa whose habitat evidence ties across types carry a
    fractional weight (the pipeline's own tie rule), so k and n are not integers and no exact
    binomial interval applies. Rounding them and using Wilson is a close approximation -- the
    winner-take-all counts, which ARE integers, put Host-plant at 43.5 % against the weighted
    42.8 %, so the choice moves the point estimate by well under the interval's width. What
    the interval does NOT cover is habitat misannotation or the sampling frame itself.
    """
    k, n = round(k), round(n)
    if n <= 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    h = Z / d * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n))
    return 100 * max(0.0, c - h), 100 * min(1.0, c + h)


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


# habitat_enrichment.jl: usable = habitat_type != "" and not an isolation source.
MIN_GROUP = 30   # a group-by-habitat cell below this is not a rate, it is noise
ISOLATION = {"Food/substance", "Lab-medium"}
PROK = {"Bacteria", "Archaea"}


def truthy(s):
    return str(s).strip().lower() == "true"


def type_weights(s):
    """'Aquatic:1.0;Terrestrial:0.5' -> {habitat: weight}, the pipeline's own format."""
    out = {}
    for part in str(s).split(";"):
        part = part.strip()
        if not part or ":" not in part:
            continue
        k, v = part.rsplit(":", 1)
        try:
            out[k] = float(v)
        except ValueError:
            pass
    return out


def flags(s):
    """'Hypersaline+Thermal' -> ['Hypersaline','Thermal']; non-exclusive, weight 1 each."""
    s = str(s).strip()
    return [f for f in s.split("+") if f] if s else []


def main():
    scan_ln = set()
    for r in csv.DictReader(open(SCAN, encoding="utf-8")):
        if "Lanthanide" in str(r.get("hit_classes", "")):
            scan_ln.add(int(r["taxid"]))

    rows = [r for r in csv.DictReader(open(CALLS, encoding="utf-8"))
            if r["config"] == "primary"]

    pool_den, pool_num = {}, {}          # habitat -> weight
    chem_den, chem_num = {}, {}
    # The same tallies split by phylum group. The overall rate is largely a composition
    # effect -- the pool is 37 % Pseudomonadota and they carry it at 47 % against 9 %
    # elsewhere -- so the split is reported rather than left for a reader to infer.
    grp_den, grp_num = {}, {}            # (group, axis, habitat) -> weight
    n_pool = n_pool_usable = 0
    n_ln_census_pool = n_ln_scan_pool = 0
    # the published Ln arm, recomputed here purely as a check on the masks
    ln_arm = 0
    ln_arm_w, ln_arm_chem = {}, {}

    for r in rows:
        taxid = int(r["taxid"])
        placed = str(r["phylum"]).strip() != ""
        ht = str(r["habitat_type"]).strip()
        usable = ht != "" and ht not in ISOLATION
        in_pool, in_cen = truthy(r["in_pool"]), truthy(r["in_census"])
        n_ln = int(r["n_ln_seqs"] or 0)

        if in_cen and n_ln > 0 and str(r["kingdom"]).strip() in PROK and placed and usable:
            ln_arm += 1
            for h, w in type_weights(r["type_weights"]).items():
                ln_arm_w[h] = ln_arm_w.get(h, 0.0) + w
            for f in flags(r["chemistry_flags"]):
                ln_arm_chem[f] = ln_arm_chem.get(f, 0.0) + 1.0

        if not in_pool:
            continue
        n_pool += 1
        if not (placed and usable):
            continue
        n_pool_usable += 1

        is_ln_census = in_cen and n_ln > 0
        is_ln_scan = (not in_cen) and taxid in scan_ln
        is_ln = is_ln_census or is_ln_scan
        n_ln_census_pool += int(is_ln_census)
        n_ln_scan_pool += int(is_ln_scan)

        grp = ("Pseudomonadota" if str(r["phylum"]).strip() == "Pseudomonadota"
               else "other phyla")
        for h, w in type_weights(r["type_weights"]).items():
            pool_den[h] = pool_den.get(h, 0.0) + w
            grp_den[(grp, "habitat", h)] = grp_den.get((grp, "habitat", h), 0.0) + w
            if is_ln:
                pool_num[h] = pool_num.get(h, 0.0) + w
                grp_num[(grp, "habitat", h)] = grp_num.get((grp, "habitat", h), 0.0) + w
        for f in flags(r["chemistry_flags"]):
            chem_den[f] = chem_den.get(f, 0.0) + 1.0
            grp_den[(grp, "chemistry", f)] = grp_den.get((grp, "chemistry", f), 0.0) + 1.0
            if is_ln:
                chem_num[f] = chem_num.get(f, 0.0) + 1.0
                grp_num[(grp, "chemistry", f)] = grp_num.get((grp, "chemistry", f), 0.0) + 1.0

    checks = [
        ("pool size == 6000", n_pool, 6000),
        ("Ln arm == 2397", ln_arm, 2397),
        ("Ln arm Terrestrial weight == 614.227", round(ln_arm_w.get("Terrestrial", 0), 3), 614.227),
        ("Ln arm Hypersaline flag == 621", round(ln_arm_chem.get("Hypersaline", 0)), 621),
        ("scan lanthanide-class taxa == 680", len(scan_ln), 680),
    ]
    print("CHECKS AGAINST KNOWN PIPELINE NUMBERS")
    ok = True
    for name, got, want in checks:
        good = got == want
        ok = ok and good
        print(f"  [{'OK ' if good else 'BAD'}] {name}: got {got}, want {want}")
    print()
    if not ok:
        print("A check failed: the masks do not match the pipeline, so the rates below")
        print("would be a different measurement wearing this one's name. Stopping.")
        return

    print(f"FRAME: {n_pool} pool taxa, {n_pool_usable} with a usable habitat and a phylum.")
    print(f"Carriers in the pool: {n_ln_census_pool} by census classification + "
          f"{n_ln_scan_pool} by the background scan = {n_ln_census_pool + n_ln_scan_pool}")
    print(f"Overall prevalence: "
          f"{100*(n_ln_census_pool+n_ln_scan_pool)/n_pool_usable:.1f} % of annotated pool taxa\n")

    out = []
    for label, num, den in (("habitat", pool_num, pool_den), ("chemistry", chem_num, chem_den)):
        print(f"{label.upper():<12}{'carriers':>10}{'sampled':>10}{'prevalence':>12}")
        for h in sorted(den, key=lambda k: -(num.get(k, 0.0) / den[k] if den[k] else 0)):
            d = den[h]
            if d < 20:
                continue
            n = num.get(h, 0.0)
            lo, hi = wilson(n, d)
            print(f"  {h:<26}{n:>8.1f}{d:>10.1f}{100*n/d:>11.1f} %"
                  f"  [{lo:>4.1f}, {hi:>4.1f}]")
            out.append((label, h, round(n, 3), round(d, 3), 100 * n / d, lo, hi))
        print()

    base_lo, base_hi = wilson(n_ln_census_pool + n_ln_scan_pool, n_pool_usable)
    print(f"baseline {100*(n_ln_census_pool+n_ln_scan_pool)/n_pool_usable:.1f} % "
          f"[{base_lo:.1f}, {base_hi:.1f}]")

    with open(os.path.join(HERE, "prevalence.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["axis", "habitat", "carriers_weighted", "pool_taxa_weighted",
                    "prevalence_pct", "ci_lo_pct", "ci_hi_pct"])
        w.writerows(out)
        w.writerow(["baseline", "all annotated pool taxa", n_ln_census_pool + n_ln_scan_pool,
                    n_pool_usable,
                    100 * (n_ln_census_pool + n_ln_scan_pool) / n_pool_usable,
                    base_lo, base_hi])

    # The split. Its own file: it answers a different question from the headline rate, and
    # mixing the two in one table is how a reader ends up quoting the wrong denominator.
    print("\nSPLIT BY PHYLUM GROUP")
    gb_num, gb_den = {}, {}
    for (g, ax, h), d in grp_den.items():
        gb_den[g] = gb_den.get(g, 0.0) + (d if ax == "habitat" else 0.0)
        gb_num[g] = gb_num.get(g, 0.0) + (grp_num.get((g, ax, h), 0.0) if ax == "habitat" else 0.0)
    for g in sorted(gb_den):
        lo, hi = wilson(gb_num[g], gb_den[g])
        print(f"  {g:<16}{gb_num[g]:>8.0f}/{gb_den[g]:<8.0f}"
              f"{100*gb_num[g]/gb_den[g]:>7.1f} %  [{lo:.1f}, {hi:.1f}]")
    with open(os.path.join(HERE, "prevalence_by_phylum.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["group", "axis", "habitat", "carriers_weighted", "pool_taxa_weighted",
                    "prevalence_pct", "ci_lo_pct", "ci_hi_pct"])
        for (g, ax, h) in sorted(grp_den):
            d = grp_den[(g, ax, h)]
            if d < MIN_GROUP:
                continue
            n = grp_num.get((g, ax, h), 0.0)
            lo, hi = wilson(n, d)
            w.writerow([g, ax, h, round(n, 3), round(d, 3), 100 * n / d, lo, hi])
        for g in sorted(gb_den):
            lo, hi = wilson(gb_num[g], gb_den[g])
            w.writerow([g, "baseline", "all annotated pool taxa", round(gb_num[g], 3),
                        round(gb_den[g], 3), 100 * gb_num[g] / gb_den[g], lo, hi])
    print("wrote prevalence_by_phylum.csv")
    print("wrote prevalence.csv")


if __name__ == "__main__":
    main()
