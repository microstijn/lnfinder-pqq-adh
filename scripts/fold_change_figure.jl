#=
    fold_change_figure.jl -- how much a setting raises the chance of carrying the gene,
    measured within each phylum against that phylum's own rate.

    WHY A FOLD-CHANGE AND NOT A PERCENTAGE. The equal-phylum percentage is the prevalence in
    a community where these four phyla are equally common, which is nowhere: Host-plant's
    28.1 % is neither the observed 42.8 % nor any phylum's rate (75/16/12/9). It reads as a
    fact about genomes and is a fact about an invented population. A fold-change compares
    each phylum only against itself, so nothing has to be invented, and the reference is 1x.

    WHAT IS DRAWN, one thing each:
      open circles  the four per-phylum fold-changes. They ARE the finding, not decoration:
                    two settings can share a mean and be completely different underneath.
      filled marker their geometric mean -- the right centre for ratios, since an arithmetic
                    mean of 0.5x and 2.0x is 1.25x, claiming a rise where there is none.
      whisker       its bootstrap interval (prevalence_phylum_average.py), which resamples
                    the setting cell with each phylum's overall rate held fixed.
      solid line    1x, the phylum's own average.
    The row label carries the raw OBSERVED percentage, unadjusted, as the concrete anchor a
    ratio cannot give.

    THE POINT OF THE FIGURE. Host-plant and Terrestrial have the same mean, 1.50x against
    1.49x, and are not the same result: Host-plant runs 0.81x to 3.56x and is BELOW 1x in
    Bacteroidota, while Terrestrial lifts all four (1.05, 1.59, 1.89, 1.54). Soil is the
    consistent signal; plants the large but erratic one. A bar chart of the means hides that
    exactly.
=#

using CSV, DataFrames, Printf, SHA
import CairoMakie as mk

# Paths follow the repository's data contract (docs/DATA.md): inputs under dat/, outputs
# where $LNFINDER_FIGDIR says, defaulting beside the inputs so a bare clone just works.
const REPO = normpath(joinpath(@__DIR__, ".."))
# No variable of its own: docs/DATA.md allows exactly two per-input overrides, and these
# tables are not among them. They resolve under the data root like everything else, and the
# root defaults to dat/ in the clone, which is where the tracked copies already sit.
const DATADIR = joinpath(get(ENV, "LNFINDER_DATA", joinpath(REPO, "dat")), "figure")
const OUTDIR = get(ENV, "LNFINDER_FIGDIR", DATADIR)
const HERE = OUTDIR
const SRC = joinpath(DATADIR, "prevalence_phylum_average.csv")
const THERM = joinpath(DATADIR, "thermal_options.csv")
# "" draws the figure as registered; "pooled" and "twophylum" append a Thermal row that is
# NOT on the same footing as the others, each marked on its face and named in the file name.
# "notitle" drops the in-figure title, for journals that set it from the caption.
const NOTITLE = "notitle" in ARGS
const MODE = (m = filter(a -> a != "notitle", ARGS); isempty(m) ? "" : m[1])
const RAW = joinpath(DATADIR, "prevalence.csv")
const PHYLA = ["Pseudomonadota", "Actinomycetota", "Bacillota", "Bacteroidota"]

# Display names. The pipeline's habitat_type values are joined on everywhere upstream, so the
# renaming lives here and the data keep their keys.
const DISPLAY = Dict("Host-plant"  => "Plant-associated",
                     "Host-animal" => "Animal-associated",
                     "Host-other"  => "Other host-associated")
show_name(s) = get(DISPLAY, s, s)

const COL_HI   = "#D55E00"
const COL_LO   = "#0072B2"
const COL_MID  = "#7F7F7F"
const COL_INK  = "#2E2E2E"
const COL_RULE = "#BFBEB6"
const PT_PER_MM = 72 / 25.4
const ROW_MM, FS = 8.5, 7
const LABEL_FS = FS - 1.5
const GAP, PNG_DPI = 1.0, 300

function main()
    mk.activate!(type = "pdf")
    d = CSV.read(SRC, DataFrame)
    d = d[string.(d.axis) .!= "baseline", :]
    raw = CSV.read(RAW, DataFrame)
    obs = Dict(string(r.habitat) => r.prevalence_pct
               for r in eachrow(raw) if string(r.axis) != "baseline")

    rows = NamedTuple[]
    for blk in ("habitat", "chemistry")
        sub = d[string.(d.axis) .== blk, :]
        sort!(sub, :fold, rev = true)
        for r in eachrow(sub)
            push!(rows, (block = blk, setting = string(r.setting), fold = r.fold,
                         lo = r.fold_lo, hi = r.fold_hi,
                         per = [r[Symbol("fold_", p)] for p in PHYLA],
                         obs = get(obs, string(r.setting), NaN)))
        end
    end
    # The Thermal row, when asked for. It fails the four-phylum rule -- Actinomycetota is 0 of
    # 13 there -- so neither variant is comparable with the rows above it, and both say so with
    # an open marker and a dagger rather than blending in.
    if MODE in ("pooled", "twophylum")
        th = CSV.read(THERM, DataFrame)
        want = MODE == "pooled" ? "pooled_fold" : "twophylum_fold"
        t = only(th[string.(th.option) .== want, :])
        push!(rows, (block = "chemistry", setting = "Thermal †", fold = t.value,
                     lo = t.lo, hi = t.hi, per = Float64[],
                     obs = MODE == "pooled" ?
                           only(th[string.(th.option) .== "pooled", :]).value : NaN))
    end
    nh = count(r -> r.block == "habitat", rows)
    ys = [-(i + (i > nh ? GAP : 0.0)) for i in 1:length(rows)]
    yrule = -(nh + GAP / 2)
    h = (length(rows) + 2.5) * ROW_MM * PT_PER_MM

    colour(r) = r.lo > 1 ? COL_HI : r.hi < 1 ? COL_LO : COL_MID
    label(r) = isnan(r.obs) ? show_name(r.setting) :
                @sprintf("%s (%.0f%%)", show_name(r.setting), r.obs)

    spread = [r for r in rows if !isempty(r.per)]
    lo_x = min(0.2, minimum(minimum(r.per) for r in spread) * 0.85,
               minimum(r.lo for r in rows) * 0.85)
    hi_x = max(5.0, maximum(maximum(r.per) for r in spread) * 1.15)
    ticks = [0.25, 0.5, 1, 2, 4]

    fig = mk.Figure(size = (150 * PT_PER_MM, h + 95), fontsize = FS)
    ax = mk.Axis(fig[1, 1];
                 title = NOTITLE ? "" : "Environment influences chances of finding Ln carriers.",
                 xlabel = "fold-change against the same phylum's own rate",
                 yticks = (ys, [label(r) for r in rows]),
                 xscale = log10, xticks = (ticks, ["0.25×", "0.5×", "1×", "2×", "4×"]),
                 ygridvisible = false, rightspinevisible = false, topspinevisible = false,
                 spinewidth = 0.6, titlesize = FS + 1, titlealign = :left, titlefont = :bold,
                 xlabelsize = FS, xticklabelsize = FS - 0.5, yticklabelsize = FS - 0.5,
                 xgridcolor = (COL_RULE, 0.5), xgridwidth = 0.4, xgridvisible = true,
                 height = h)
    mk.vlines!(ax, [1.0]; color = (COL_INK, 0.8), linewidth = 0.9)
    for (y, r) in zip(ys, rows)
        c = colour(r)
        if !isempty(r.per)
            mk.linesegments!(ax, [mk.Point2f(minimum(r.per), y), mk.Point2f(maximum(r.per), y)];
                             color = (COL_MID, 0.45), linewidth = 0.5)
            mk.scatter!(ax, r.per, fill(y, length(r.per)); markersize = 4.2,
                        color = (:white, 1.0), strokecolor = (COL_MID, 0.95), strokewidth = 0.6)
        end
        mk.errorbars!(ax, [r.fold], [y], [r.fold - r.lo], [r.hi - r.fold];
                      direction = :x, whiskerwidth = 3.0, linewidth = 1.0, color = c)
        mk.scatter!(ax, [r.fold], [y]; marker = :diamond, markersize = 6.0,
                    color = isempty(r.per) ? (:white, 1.0) : c,
                    strokecolor = isempty(r.per) ? c : COL_INK,
                    strokewidth = isempty(r.per) ? 1.1 : 0.4)
    end
    mk.hlines!(ax, [yrule]; color = (COL_RULE, 0.9), linewidth = 0.5, linestyle = :dot)
    for (txt, y) in (("HABITAT", 0.30), ("CHEMISTRY (non-exclusive)", yrule - 0.15))
        mk.text!(ax, lo_x * 1.04, y; text = txt, align = (:left, :top), fontsize = LABEL_FS,
                 color = (COL_MID, 0.95), font = :bold)
    end
    mk.xlims!(ax, lo_x, hi_x)
    mk.ylims!(ax, minimum(ys) - 0.8, 0.35)

    mk.Legend(fig[2, 1],
              [mk.MarkerElement(color = COL_HI, marker = :diamond, markersize = 6.0,
                                strokecolor = COL_INK, strokewidth = 0.4),
               mk.MarkerElement(color = COL_MID, marker = :diamond, markersize = 6.0,
                                strokecolor = COL_INK, strokewidth = 0.4),
               mk.MarkerElement(color = COL_LO, marker = :diamond, markersize = 6.0,
                                strokecolor = COL_INK, strokewidth = 0.4),
               mk.MarkerElement(color = (:white, 1.0), marker = :circle, markersize = 4.2,
                                strokecolor = (COL_MID, 0.95), strokewidth = 0.6),
               mk.LineElement(color = COL_INK, linewidth = 1.0)],
              ["raises the chance", "no clear effect", "lowers it",
               "one phylum", "95% interval on the mean"];
              orientation = :horizontal, framevisible = false, labelsize = FS - 0.5,
              patchlabelgap = 3, colgap = 10, nbanks = 2, patchsize = (7.0f0, 7.0f0),
              rowgap = 2)
    # No explanatory line under the figure. What the marks mean is caption material, and a
    # figure carrying its own methods note only competes with the caption that must repeat it.
    mk.rowgap!(fig.layout, 1, 6)
    mk.resize_to_layout!(fig)

    bp = joinpath(HERE, (isempty(MODE) ? "fig3_fold_change" :
                         "fig3_fold_change_thermal_" * MODE) * (NOTITLE ? "_notitle" : ""))
    written = String[]
    for (fmt, kw) in (("pdf", (pt_per_unit = 1,)), ("svg", (pt_per_unit = 1,)),
                      ("png", (px_per_unit = PNG_DPI / 72,)))
        p = bp * "." * fmt
        mk.save(p, fig; kw...)
        push!(written, p)
        println("wrote ", p)
    end
    open(joinpath(HERE, "FOLD_CHANGE_NOTE.txt"), "w") do io
        println(io, "Fold-change figure -- ",
                Base.Libc.strftime("%Y-%m-%dT%H:%M:%S", time()))
        println(io, "julia ", VERSION, "; CairoMakie ", pkgversion(mk))
        println(io, "rows: ", SRC, "  ", bytes2hex(open(sha256, SRC)))
        println(io, "raw : ", RAW, "  ", bytes2hex(open(sha256, RAW)))
        println(io, "\nreference 1x = that phylum's own overall rate\n")
        for r in rows
            @printf(io, "%-14s %5.2fx [%.2f, %.2f]  observed %4.1f %%  per-phylum %s\n",
                    r.setting, r.fold, r.lo, r.hi, r.obs,
                    join((@sprintf("%.2f", x) for x in r.per), "/"))
        end
        println(io)
        for p in written
            println(io, bytes2hex(open(sha256, p)), "  ", basename(p))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
