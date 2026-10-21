#=
    Offline self-test for scripts/habitat_figure.jl — the two-panel habitat enrichment figure,
    the gates between its five inputs, the values it exports and the staging that keeps a
    failed run out of the figure directory. No network, no repository data file and no
    renderer: every input is built here in a temporary directory, and the drawing step is
    replaced by a function that writes placeholder bytes.

    The script is included behind its own `PROGRAM_FILE` guard, so including it defines the
    functions and runs nothing; the suite checks that the include did not load the renderer.

    EVERY FIXTURE IS ABSTRACT. Categories are role names ("cat-a", "chem-x"), variants are
    "V:one".."V:three", digests are hashes of placeholder strings, and every estimate and
    interval is a number chosen so that one branch is decided. Nothing here is a quantity
    measured on the data the script analyses.

    THE FAILURES THIS SUITE EXISTS TO PREVENT.

    1. A COLOUR THAT DOES NOT COME FROM THE FLAG. A row whose primary interval excludes zero
       but whose flag is not a sign must be grey, and a flag that contradicts its own interval
       must stop the run rather than fall back to grey. The three inks are pinned by value and
       are pairwise different, and the ink of a drawn row of each class is checked, because
       the drawing step is deliberately not run here.
    2. A GATE THAT NEVER FIRES. Each named gate is driven to failure by one mutation of one
       input, and must be the only gate that fails under it. Each CONDITION inside a gate is
       driven separately, so a later edit cannot narrow one and keep the gate's name.
    2b. A FLAG ITS OWN GRID DOES NOT SUPPORT, including one written consistently into two
       files: the sign and the breakers are recomputed from the grid rows they summarise.
    2c. A PUBLISHED FILE THAT NAMES THE MACHINE IT RAN ON. The provenance records file names
       and digests; no value may look like a local path.
    3. A TICK FROM THE WRONG ROW. The comparator is keyed on unit AND category, so a chemistry
       flag never borrows a habitat category's tick of the same name.
    4. A MARK OUTSIDE ITS AXIS, OR UNDER ANOTHER MARK. Limits cover every drawn mark; the
       sub-row clearances are positive at the declared row pitch.
    5. A VALUE UNDER A KEY NOBODY DECLARED, OR A DECLARED KEY NOBODY BUILT.
    6. A HALF-PUBLISHED RUN. A failed gate leaves the figure directory byte-identical; a
       passing run publishes exactly the declared outputs and nothing else.

    Run:  julia --project=. test/test_habitat_figure.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Test, DataFrames, CSV, SHA

const CAIROMAKIE_ID = Base.PkgId(Base.UUID("13f3f980-e62b-5c42-98c6-ff1f3baf88f0"), "CairoMakie")
const RENDERER_LOADED_BEFORE = haskey(Base.loaded_modules, CAIROMAKIE_ID)

module HabitatFigure end
Base.include(HabitatFigure, joinpath(@__DIR__, "..", "scripts", "habitat_figure.jl"))
const HF = HabitatFigure

const LN   = HF.HE.LN_VS_BG
const NON  = HF.HE.NON_VS_BG
const ALL  = HF.HE.ALL_VS_BG
const CA   = HF.HE.CA_VS_BG
const HAB  = HF.HABITAT_UNIT
const CHEM = HF.CHEMISTRY_UNIT
# Two rule axes, so the declaration below reaches one unit with one axis and both with the
# other, and the per-unit variant counts are not the same number twice.
const VARIANTS = ["V:one", "V:two", "W:three"]
const N1, N2, NNON = 300, 500, 400
# THE SELECTION THE FIXTURE IS HELD TO. The script ships no category selection, so every
# run here is handed one; abstract role names, like everything else in this file. The key
# prefixes below are what `HE.key_slug` derives from these names, and the fixture's summary
# is keyed with them — so a run whose selection changed would export DIFFERENT keys, which is
# the property the assertions on them hold.
const SEL_TIER    = "cat-t"          # the depletion-tier rule's category
const SEL_GENUS   = "cat-g"          # the genus-clause rule's category
const SEL_SHARE   = "cat-c"          # the sharing rule's category: drawn, and read back
const SEL_GRID    = ("cat-b",)       # counted across every specification: drawn
const SEL_PHYLUM  = "phy-2"          # NOT the fixture's heaviest stratum ("phy-1"), so the
                                     # share identities are off until a test switches them on
const K_TIER   = HF.HE.key_slug(SEL_TIER)
const K_PHYLUM = HF.HE.key_slug(SEL_PHYLUM)
const GRID_CATS = SEL_GRID
const SHARING_CAT = SEL_SHARE

"""Write a category-selection criteria file. `over` replaces individual role values."""
function selection_csv(dir::AbstractString; grid = SEL_GRID, sharing = SEL_SHARE,
                       tier = SEL_TIER, genus = SEL_GENUS, phylum = SEL_PHYLUM,
                       name = "categories.csv")
    rows = vcat([("expect", "selections", string(4 + length(grid)), "fixture")],
                [("select", "depletion_tier", tier, "fixture"),
                 ("select", "genus_clause", genus, "fixture"),
                 ("select", "sharing", sharing, "fixture")],
                [("select", "grid_count", c, "fixture") for c in grid],
                [("select", "standardisation_phylum", phylum, "fixture")])
    path = joinpath(dir, name)
    CSV.write(path, DataFrame(kind = [r[1] for r in rows], role = [r[2] for r in rows],
                              value = [r[3] for r in rows], note = [r[4] for r in rows]))
    return path
end

const SEL_PATH = selection_csv(mktempdir())
const SELECTED = HF.HE.read_category_selection(SEL_PATH)

# ---------------------------------------------------------------- the fixture

mutable struct Fixture
    stats::DataFrame
    summary::Dict{String,String}
    rob::DataFrame
    grid::DataFrame
    prov::Dict{String,String}
    calls::Vector{UInt8}
    panels::Dict{String,Vector{UInt8}}
end

# unit, category, a, b, std, lo, hi, crude, crude_lo, crude_hi, tick, tick_lo, tick_hi,
# sign, pez, breakers, zero_cells, min_cell
const SPECS = [
    (HAB,  "cat-a", 40.5,  128.0, -1.10, -1.60, -0.60, -0.90, -1.30, -0.50, -0.70, -1.00, -0.40, "-",       true,  "",        0, 3.5),
    (HAB,  "cat-b", 90.25, 140.0, -0.10, -0.35,  0.15,  0.30,  0.10,  0.50,  0.21,  0.02,  0.40, "null",    false, "W:three", 0, 12.0),
    (HAB,  "cat-c", 70.0,  60.0,   0.80,  0.40,  1.20,  0.50,  0.20,  0.80,  0.40,  0.10,  0.70, "+",       true,  "",        0, 2.75),
    (HAB,  "cat-d", 25.75, 35.0,   0.35, -0.25,  0.95,  0.10, -0.30,  0.50,  0.36, -0.20,  0.92, "null",    false, "",        1, 0.0),
    (HAB,  "cat-e", 6.5,   9.0,    0.20, -1.50,  1.90,  0.05, -1.00,  1.10,  0.15, -1.20,  1.50, "null",    false, "",        0, 1.5),
    (CHEM, "chem-x", 60.0, 70.0,   0.45,  0.10,  0.80,  0.55,  0.25,  0.85,  0.62,  0.30,  0.94, "fragile", true,  "V:two",   0, 8.0),
    (CHEM, "chem-y", 33.0, 70.0,  -0.30, -0.90,  0.30, -0.50, -0.95, -0.05, -0.20, -0.70,  0.30, "null",    false, "W:three", 0, 4.0),
]

# THE GRID IS THE FLAGS' OWN EVIDENCE, so the fixture's variant rows are the rows that give
# the flag each SPECS line declares: a signed flag is reproduced by all three variants, a
# `null` is spanned by every variant that is not named as a breaker, and a named breaker
# carries an interval on the other side of the decision. `(unit, category) => three intervals`
# in VARIANTS order; `nothing` for a variant the grid did not run.
const NO_ROW = nothing      # the grid did not run this variant on this row
const SAME = :same          # the variant returned the registered specification, bit for bit
const VARIANT_CELLS = Dict(
    (HAB, "cat-a")       => [(-0.95, -1.40, -0.50), (-1.05, -1.55, -0.55), (-1.20, -1.73, -0.70)],
    (HAB, "cat-b")       => [(-0.15, -0.40,  0.10), (-0.05, -0.30,  0.20), ( 0.32,  0.05,  0.60)],
    (HAB, "cat-c")       => [( 0.75,  0.35,  1.15), ( 0.85,  0.45,  1.25), ( 0.70,  0.30,  1.10)],
    (HAB, "cat-d")       => [( 0.30, -0.30,  0.90), ( 0.40, -0.20,  1.00), ( 0.25, -0.35,  0.85)],
    (HAB, "cat-e")       => [( 0.18, -1.45,  1.81), ( 0.24, -1.38,  1.86), ( 0.16, -1.47,  1.84)],
    (HAB, "cat-z")       => [( 0.40,  0.10,  0.70), NO_ROW, NO_ROW],
    (CHEM, "chem-x")     => [( 0.52,  0.16,  0.88), ( 0.05, -0.20,  0.30), ( 0.48,  0.12,  0.84)],
    (CHEM, "chem-y")     => [(-0.28, -0.88,  0.32), (-0.34, -0.94,  0.26), (-0.55, -1.00, -0.10)],
    # no robustness row reads these: they are the per-sequence spread the values file counts,
    # and V:one returns the primary unchanged, so a count of specifications is not a count of
    # distinct answers.
    ("sequences", "cat-b") => [SAME, (-0.22, -0.50, 0.05), (-0.40, -0.65, -0.15)],
)

stat_row(contrast, unit, cat, a, n1, b, n2, crude, clo, chi, std, lo, hi; zero = 0, mincell = 5.0,
         tested = true, fisher = 0.002) =
    (contrast = contrast, unit = unit, habitat = cat, a = a, n1 = n1, b = b, n2 = n2,
     pct_fg = 100 * a / n1, pct_bg = 100 * b / n2, log2_rr = crude, ci_lo = clo, ci_hi = chi,
     log2_rr_std = std, ci_lo_std = lo, ci_hi_std = hi, n_strata = 4, fisher_p = fisher,
     q_bh = 2 * fisher, tested = tested, continuity_corrected = false, rounded_for_fisher = true,
     n_zero_stratum_cells = zero, min_bg_stratum_cell = mincell)

grid_row(role, variant, s) =
    (variant = variant, axis = role == "primary" ? "primary" : first(variant, 1),
     level = "lvl", role = role, contrast = s.contrast, unit = s.unit, category = s.habitat,
     a = s.a, b = s.b, n1 = s.n1, n2 = s.n2, crude = s.log2_rr, crude_lo = s.ci_lo,
     crude_hi = s.ci_hi, std = s.log2_rr_std, lo = s.ci_lo_std, hi = s.ci_hi_std, n_strata = 4,
     n_zero_stratum_cells = s.n_zero_stratum_cells, min_bg_stratum_cell = s.min_bg_stratum_cell,
     tested = s.tested, q_bh = s.q_bh, feeds_flag = s.contrast == LN)

function variant_row(variant, unit, cat, std, lo, hi)
    s = stat_row(LN, unit, cat, 10.0, N1, 10.0, N2, 0.1, 0.0, 0.2, std, lo, hi)
    return grid_row("variant", variant, s)
end

function rob_row(unit, cat, sign, pez, breakers, std, lo, hi, a; expected = 3, with_row = 3)
    return (unit = unit, category = cat, robust_sign = sign, primary_excludes_zero = pez,
            n_breaking_variants = length(HF.split_codes(breakers)), breaking_variants = breakers,
            breaking_variants_are = "", depletion_tier = "", depletion_tier_breakers = "",
            primary_std = std, primary_lo = lo, primary_hi = hi, primary_a = a,
            variants_expected = expected, variants_with_a_row = with_row,
            n_variants_that_moved_the_estimate = 2)
end

"Rebuild the provenance's flag keys and robust lists from the robustness rows."
function sync_provenance!(fx::Fixture)
    for k in collect(keys(fx.prov))
        (startswith(k, "flag_") || startswith(k, "breakers_")) && delete!(fx.prov, k)
    end
    for r in eachrow(fx.rob)
        fx.prov["flag_$(r.unit)_$(r.category)"] = r.robust_sign
        isempty(r.breaking_variants) || (fx.prov["breakers_$(r.unit)_$(r.category)"] = r.breaking_variants)
    end
    cats(sign) = join(sort(unique([r.category for r in eachrow(fx.rob) if r.robust_sign == sign])), ";")
    fx.prov["robust_enriched"] = cats("+")
    fx.prov["robust_depleted"] = cats("-")
    fx.prov["fragile_categories"] = cats("fragile")
    return fx
end

function base_fixture()
    stats = NamedTuple[]
    for s in SPECS
        unit, cat, a, b, std, lo, hi, crude, clo, chi, tick, tlo, thi, _, _, _, zero, mincell = s
        push!(stats, stat_row(LN, unit, cat, a, N1, b, N2, crude, clo, chi, std, lo, hi;
                              zero = zero, mincell = mincell))
        push!(stats, stat_row(NON, unit, cat, 1.25 * a, NNON, b, N2, tick - 0.1, tlo - 0.1,
                              thi - 0.1, tick, tlo, thi))
        push!(stats, stat_row(ALL, unit, cat, 2.25 * a, N1 + NNON, b, N2, crude, clo, chi,
                              (std + tick) / 2, lo, hi))
        unit == HAB &&
            push!(stats, stat_row(LN, "sequences", cat, 1.5 * a, 450, 1.5 * b, 700, crude, clo,
                                  chi, std - 0.35, lo - 0.35, hi - 0.35))
    end
    push!(stats, stat_row(CA, HAB, "cat-a", 3.0, 30, 130.0, N2, -0.2, -1.0, 0.6, -0.1, -1.1, 0.9))
    stats_df = DataFrame(stats)

    grid = NamedTuple[grid_row("primary", "primary", s) for s in stats]
    primary_cell(unit, cat) =
        (s = stats[findfirst(r -> r.contrast == LN && r.unit == unit && r.habitat == cat, stats)];
         (s.log2_rr_std, s.ci_lo_std, s.ci_hi_std))
    for ((unit, cat), cells) in sort(collect(VARIANT_CELLS), by = first)
        for (v, cell) in zip(VARIANTS, cells)
            cell === NO_ROW && continue
            push!(grid, variant_row(v, unit, cat, (cell === SAME ? primary_cell(unit, cat) : cell)...))
        end
    end

    rob = NamedTuple[rob_row(s[1], s[2], s[14], s[15], s[16], s[5], s[6], s[7], s[3]) for s in SPECS]
    push!(rob, rob_row(HAB, "cat-z", "undecidable", false, "V:one", NaN, NaN, NaN, NaN; with_row = 1))

    calls = Vector{UInt8}("placeholder call table\n")
    scan = bytes2hex(sha256("placeholder scan table"))
    summary = Dict{String,String}(
        "n1_taxa" => string(N1), "n2_taxa" => string(N2), "n_nonln_carrier_taxa" => string(NNON),
        "background_scan_sha256" => scan, "$(K_TIER)_tier_code" => "code-t",
        "n_sequences_fg" => "450", "n_bacteria_fg" => "288", "n_archaea_fg" => "12",
        "eukaryotic_ln_taxa_excluded" => "7", "eukaryotic_ln_sequences_excluded" => "11",
        "pct_$(K_PHYLUM)_fg" => "61.25", "pct_$(K_PHYLUM)_bg" => "12.5",
        # the sampled frame less what the exclusion rule removed is the analysed background:
        # 1000 - 360 = 640, an identity the gate checks
        "pct_$(K_PHYLUM)_bg_before_scan" => "30.125", "bg_removed_by_scan" => "360",
        "nonln_coverage_numerator" => "400", "nonln_coverage_denominator" => "460",
        "nonln_coverage_pct" => string(100 * 400 / 460),
        "fg_coverage_numerator" => "300", "fg_coverage_denominator" => "320",
        "fg_coverage_pct" => "93.75",
        "bg_analysed_coverage_numerator" => "500", "bg_analysed_coverage_denominator" => "640",
        "bg_analysed_coverage_pct" => "78.125", "coverage_fisher_p_after_scan" => "0.00064",
        "bg_frame_coverage_numerator" => "800", "bg_frame_coverage_denominator" => "1000",
        "bg_frame_coverage_pct" => "80.0", "coverage_fisher_p" => "0.0312",
        "fg_relations_text_mined_pct" => "65.5", "bg_relations_text_mined_pct" => "66.25",
        "fg_taxa_tied_on_distinct" => "40", "fg_taxa_tie_broken_by_mentions" => "25",
        "fg_taxa_still_tied_after_mentions" => "16",
        "bg_taxa_tied_on_distinct" => "50", "bg_taxa_tie_broken_by_mentions" => "30",
        "bg_taxa_still_tied_after_mentions" => "20",
        "fg_isolation_terms_total" => "1200", "fg_taxa_with_isolation_term" => "60",
        "fg_generic_terms_total" => "55", "fg_taxa_with_generic_term" => "44",
        "bg_isolation_terms_total" => "130", "bg_taxa_with_isolation_term" => "70",
        "bg_generic_terms_total" => "65", "bg_taxa_with_generic_term" => "51",
        "phylum_strata_used_fg_vs_bg" => "4", "standardisation_heaviest_stratum" => "phy-1",
        "standardisation_heaviest_stratum_weight" => "0.625",
        "standardisation_heaviest_stratum_bg_taxa" => "80",
        "n_genera_fg" => "150", "n_genera_bg" => "210",
        "background_licence_code" => "code-l", "omnicrobe_release" => "release-r",
        "obt_release" => "release-o", "uniprot_release" => "release-u",
        "sharing_ratio_nonln_over_ln" => string(0.40 / 0.80),
        "sharing_ratio_ci_lo" => "0.25", "sharing_ratio_ci_hi" => "0.75",
        "sharing_code" => "code-s", "sharing_rule_used" => "rule-r",
        "categories_legend_only" => "chem-x;cat-d")
    prov = Dict{String,String}(
        "calls_sha256" => bytes2hex(sha256(calls)), "scan_sha256" => scan,
        "anchor_checks_declined" => "false", "anchor_targets_failed" => "0",
        "publish_gates" => "gate-1;gate-2", "publish_gates_failed" => "",
        "panel_written" => "panel.pdf;panel.png",
        "variants_in_the_registered_set" => join(VARIANTS, ";"),
        "n_variants_in_the_registered_set" => "3",
        "axis_unit_declaration" => "V=$(HAB),$(CHEM);W=$(CHEM)",
        "variants_whose_axis_reaches_$(HAB)" => "2",
        "variants_whose_axis_reaches_$(CHEM)" => "3",
        "panel_categories" => "6", "panel_categories_not_drawn" => "cat-z")
    fx = Fixture(stats_df, summary, DataFrame(rob), DataFrame(grid), prov, calls,
                 Dict("panel.pdf" => Vector{UInt8}("pdf bytes"), "panel.png" => Vector{UInt8}("png bytes")))
    return sync_provenance!(fx)
end

write_kv(path, d) = CSV.write(path, DataFrame(key = collect(keys(d)), value = collect(values(d))))

function write_inputs(fx::Fixture, dir::AbstractString)
    src = joinpath(dir, "in"); mkpath(src)
    fig = joinpath(dir, "fig"); mkpath(fig)
    CSV.write(joinpath(src, "stats.csv"), fx.stats)
    write_kv(joinpath(src, "summary.csv"), fx.summary)
    CSV.write(joinpath(src, "rob.csv"), fx.rob)
    CSV.write(joinpath(src, "grid.csv"), fx.grid)
    write_kv(joinpath(src, "sprov.csv"), fx.prov)
    write(joinpath(src, "calls.csv"), fx.calls)
    for (n, b) in fx.panels
        write(joinpath(src, n), b)
    end
    args = ["--categories", SEL_PATH,
            "--figdir", fig, "--stats", joinpath(src, "stats.csv"),
            "--summary", joinpath(src, "summary.csv"), "--robustness", joinpath(src, "rob.csv"),
            "--sensitivity", joinpath(src, "grid.csv"),
            "--sensitivity-provenance", joinpath(src, "sprov.csv"),
            "--calls", joinpath(src, "calls.csv"), "--panel-dir", src]
    return (src = src, fig = fig, args = args)
end

const NO_ENV = Dict{String,String}()

function load(fx::Fixture)
    p = write_inputs(fx, mktempdir())
    opt = HF.parse_run_options(p.args; env = NO_ENV)
    return (paths = p, opt = opt, inp = HF.load_inputs(opt))
end

failed_gates(fx::Fixture) = [g.name for g in HF.input_gates(load(fx).inp) if !g.ok]

function set_stat!(fx, contrast, unit, cat, col, value)
    i = findfirst(r -> r.contrast == contrast && r.unit == unit && r.habitat == cat, eachrow(fx.stats))
    fx.stats[i, col] = value
    return fx
end
function set_grid!(fx, contrast, unit, cat, col, value; variant = "primary")
    i = findfirst(r -> r.variant == variant && r.contrast == contrast && r.unit == unit &&
                       r.category == cat, eachrow(fx.grid))
    fx.grid[i, col] = value
    return fx
end
function set_rob!(fx, unit, cat, col, value)
    i = findfirst(r -> r.unit == unit && r.category == cat, eachrow(fx.rob))
    fx.rob[i, col] = value
    return fx
end

"Every file under `dir` with its bytes, and every directory, as one comparable value."
function snapshot(dir::AbstractString)
    out = Dict{String,Vector{UInt8}}()
    for (root, dirs, files) in walkdir(dir)
        for d in dirs
            out[relpath(joinpath(root, d), dir) * "/"] = UInt8[]
        end
        for f in files
            out[relpath(joinpath(root, f), dir)] = read(joinpath(root, f))
        end
    end
    return out
end

function fake_render(plan, base)
    return [begin
                p = base * "." * f
                write(p, "placeholder $f bytes for $(length(plan.rows)) rows")
                p
            end for f in HF.IMAGE_FORMATS]
end

run_fixture(opt; render = fake_render) =
    HF.run_figure(opt; render = render, io = devnull, err = IOBuffer())

"""The loaded inputs with one role of the selection replaced, for the export checks that need
a selection the figure does not draw."""
with_selection(inp; kw...) =
    merge(inp, (selection = HF.HE.read_category_selection(
                                selection_csv(mktempdir(); kw...)),))

const GATE_MUTATIONS = [
    # the comparator arm size alone, with the summary left consistent with itself (its coverage
    # numerator is the same arm), so only the comparison against the grid's own rows can fire
    "arm sizes" => fx -> (fx.summary["n_nonln_carrier_taxa"] = string(NNON + 3);
                          fx.summary["nonln_coverage_numerator"] = string(NNON + 3);
                          fx.summary["nonln_coverage_pct"] = string(100 * (NNON + 3) / 460)),
    "summary identities" => fx -> (fx.summary["fg_coverage_pct"] = "90.0"),
    "robustness rows" => fx -> (deleteat!(fx.rob, findfirst(==("cat-d"), fx.rob.category));
                                sync_provenance!(fx)),
    "comparator rows" => fx -> deleteat!(fx.stats, findfirst(r -> r.contrast == NON && r.habitat == "cat-d",
                                                             eachrow(fx.stats))),
    "robustness values" => fx -> set_rob!(fx, HAB, "cat-d", :primary_std, 0.35 + 1e-6),
    "robustness flags" => fx -> set_rob!(fx, CHEM, "chem-x", :n_breaking_variants, 2),
    "flags against grid" => fx -> set_grid!(fx, LN, HAB, "cat-c", :lo, -0.05; variant = "V:one"),
    "sensitivity primary cells" => fx -> set_grid!(fx, LN, HAB, "cat-d", :std, 0.35 + 1e-6),
    "provenance flags" => fx -> (fx.prov["flag_$(HAB)_cat-d"] = "fragile"),
    "scan digest" => fx -> (fx.prov["scan_sha256"] = bytes2hex(sha256("another scan"))),
    "calls digest" => fx -> (fx.calls = Vector{UInt8}("another call table\n")),
    "sensitivity run state" => fx -> (fx.prov["anchor_checks_declined"] = "true"),
    "panel files" => fx -> delete!(fx.panels, "panel.png"),
    "drawn rows" => fx -> set_stat!(fx, LN, HAB, "cat-d", :tested, false),
]

@testset "habitat_figure.jl" begin

# ---------------------------------------------------------------- include
@testset "including the script defines the functions and loads no renderer" begin
    @test isdefined(HF, :draw_figure) && isdefined(HF, :cairomakie) && isdefined(HF, :main)
    if RENDERER_LOADED_BEFORE
        @warn "CairoMakie was already loaded before the include; the no-renderer check is vacuous here"
    end
    @test RENDERER_LOADED_BEFORE || !haskey(Base.loaded_modules, CAIROMAKIE_ID)
    @test all(c -> startswith(c, "#") && length(c) == 7,
              (HF.COL_ENRICHED, HF.COL_DEPLETED, HF.COL_NEUTRAL, HF.COL_INK, HF.COL_RULE))
end

# ---------------------------------------------------------------- colour
@testset "colour class comes from the robustness flag and nothing else" begin
    cc(sign, pez, lo, hi; std = NaN) = HF.colour_class(sign, pez, lo, hi; std = std, what = "u | c")
    @test cc("+", "true", 0.2, 0.9) === :enriched
    @test cc("-", "true", -0.9, -0.2) === :depleted
    # a primary interval that excludes zero does not colour a fragile row, on either side
    @test cc("fragile", "true", 0.2, 0.9) === :neutral
    @test cc("fragile", "true", -0.9, -0.2) === :neutral
    @test cc("null", "false", -0.2, 0.9) === :neutral
    @test cc("undecidable", "false", NaN, NaN) === :neutral
    @test cc("+", "true", 0.2, 0.9; std = 0.5) === :enriched

    inconsistent = [
        ("+", "false", 0.2, 0.9), ("+", "true", -0.2, 0.9), ("+", "true", -0.9, -0.2),
        ("+", "true", NaN, NaN),
        ("-", "false", -0.9, -0.2), ("-", "true", -0.9, 0.2), ("-", "true", 0.2, 0.9),
        ("fragile", "false", 0.2, 0.9), ("fragile", "true", -0.2, 0.9), ("fragile", "true", NaN, NaN),
        ("null", "true", -0.2, 0.9), ("null", "false", 0.2, 0.9), ("null", "false", NaN, NaN),
        ("undecidable", "true", NaN, NaN), ("undecidable", "false", -0.2, 0.9),
        ("maybe", "false", -0.2, 0.9), ("null", "yes", -0.2, 0.9), ("+", "", 0.2, 0.9),
    ]
    for (sign, pez, lo, hi) in inconsistent
        err = try
            cc(sign, pez, lo, hi); nothing
        catch e
            e
        end
        @test err isa ErrorException
        err isa ErrorException && @test occursin("u | c", err.msg)
    end
    # an estimate outside its own interval is not a flag to draw
    @test_throws ErrorException cc("+", "true", 0.2, 0.9; std = 1.5)
    @test_throws ErrorException cc("null", "false", -0.2, 0.9; std = -0.5)
end

# ---------------------------------------------------------------- gates
@testset "the fixture passes every gate, and the gates are the declared ones" begin
    gates = HF.input_gates(load(base_fixture()).inp)
    @test [g.name for g in gates] == collect(HF.INPUT_GATES)
    @test all(g -> g.ok, gates)
end


@testset "each gate fails alone under its own mutation" begin
    @test sort(first.(GATE_MUTATIONS)) == sort(collect(HF.INPUT_GATES))
    for (name, mutate!) in GATE_MUTATIONS
        fx = base_fixture()
        mutate!(fx)
        @test failed_gates(fx) == [name]
    end
end

@testset "further mutations reach the gate they belong to" begin
    # the sensitivity run failed one of its own publish gates, or missed an anchor target
    fx = base_fixture(); fx.prov["publish_gates_failed"] = "gate-2"
    @test failed_gates(fx) == ["sensitivity run state"]
    fx = base_fixture(); fx.prov["anchor_targets_failed"] = "1"
    @test failed_gates(fx) == ["sensitivity run state"]
    fx = base_fixture(); delete!(fx.prov, "publish_gates_failed")
    @test failed_gates(fx) == ["sensitivity run state"]
    # a sensitivity run that records no publish gate at all is not a gated run
    fx = base_fixture(); fx.prov["publish_gates"] = ""
    @test failed_gates(fx) == ["sensitivity run state"]
    # the foreground arm size on a Ln row, changed everywhere the other gates read it
    fx = base_fixture(); set_stat!(fx, LN, HAB, "cat-b", :n1, N1 + 1)
    set_grid!(fx, LN, HAB, "cat-b", :n1, N1 + 1)
    set_stat!(fx, LN, HAB, "cat-b", :pct_fg, 100 * 90.25 / (N1 + 1))
    @test failed_gates(fx) == ["arm sizes"]
    # a flag written consistently in both files but contradicting its own interval; the grid
    # it summarises still says what it says, so both gates name it
    fx = base_fixture(); set_rob!(fx, HAB, "cat-d", :robust_sign, "+")
    set_rob!(fx, HAB, "cat-d", :primary_excludes_zero, true); sync_provenance!(fx)
    @test failed_gates(fx) == ["robustness flags", "flags against grid"]
    # a robustness row with no enrichment row that is not undecidable
    fx = base_fixture(); set_rob!(fx, HAB, "cat-z", :robust_sign, "null")
    set_rob!(fx, HAB, "cat-z", :primary_lo, -0.5); set_rob!(fx, HAB, "cat-z", :primary_hi, 0.5)
    set_rob!(fx, HAB, "cat-z", :primary_std, 0.0); sync_provenance!(fx)
    @test failed_gates(fx) == ["robustness values", "flags against grid"]
    # a duplicated robustness row
    fx = base_fixture(); push!(fx.rob, fx.rob[findfirst(==("cat-a"), fx.rob.category), :])
    @test "robustness rows" in failed_gates(fx)
    # a category with an enrichment row on a primary unit but no flag, held back or not
    for cat in ("cat-d", "cat-e")
        fx = base_fixture(); deleteat!(fx.rob, findfirst(==(cat), fx.rob.category))
        sync_provenance!(fx)
        @test failed_gates(fx) == ["robustness rows"]
    end
    # the robust list in the provenance names a category the flags do not
    fx = base_fixture(); fx.prov["robust_enriched"] = "cat-c;cat-d"
    @test failed_gates(fx) == ["provenance flags"]
    # the fragile list, and the breaker list, are checked the same way
    fx = base_fixture(); fx.prov["fragile_categories"] = "chem-y"
    @test failed_gates(fx) == ["provenance flags"]
    fx = base_fixture(); fx.prov["breakers_$(CHEM)_chem-x"] = "V:one"
    @test failed_gates(fx) == ["provenance flags"]
    # a flag key for a row the robustness file does not hold
    fx = base_fixture(); fx.prov["flag_$(HAB)_cat-q"] = "null"
    @test failed_gates(fx) == ["provenance flags"]
    # a drawn comparator with no finite interval, consistent in the grid
    fx = base_fixture(); set_stat!(fx, NON, CHEM, "chem-y", :ci_hi_std, NaN)
    set_grid!(fx, NON, CHEM, "chem-y", :hi, NaN)
    @test failed_gates(fx) == ["comparator rows"]
    # an untested comparator row, and one that exists only on the other unit
    fx = base_fixture(); set_stat!(fx, NON, HAB, "cat-b", :tested, false)
    @test failed_gates(fx) == ["comparator rows"]
    fx = base_fixture(); set_stat!(fx, NON, CHEM, "chem-x", :unit, HAB)
    set_grid!(fx, NON, CHEM, "chem-x", :unit, HAB)
    @test failed_gates(fx) == ["comparator rows"]
    # a drawn crude estimate that is not a number, consistent in the grid
    fx = base_fixture(); set_stat!(fx, LN, CHEM, "chem-y", :log2_rr, NaN)
    set_grid!(fx, LN, CHEM, "chem-y", :crude, NaN)
    @test failed_gates(fx) == ["drawn rows"]
    # no silent crude fallback: a drawn standardised estimate that is not a number fails
    fx = base_fixture(); set_stat!(fx, LN, HAB, "cat-d", :log2_rr_std, NaN)
    @test "drawn rows" in failed_gates(fx)
    # a share that is not a number, and a share that is a number but not its own count over
    # its own arm size
    fx = base_fixture(); set_stat!(fx, LN, HAB, "cat-b", :pct_fg, NaN)
    @test failed_gates(fx) == ["drawn rows"]
    fx = base_fixture(); set_stat!(fx, LN, HAB, "cat-b", :pct_bg, 99.0)
    @test failed_gates(fx) == ["drawn rows"]
end

@testset "a value the robustness file states twice is checked against itself" begin
    # the count of breaking variants against the names beside it
    fx = base_fixture(); set_rob!(fx, CHEM, "chem-x", :n_breaking_variants, 2)
    @test failed_gates(fx) == ["robustness flags"]
    # a sign is a claim that no variant broke it: it may not also name one
    fx = base_fixture(); set_rob!(fx, HAB, "cat-c", :breaking_variants, "V:one")
    set_rob!(fx, HAB, "cat-c", :n_breaking_variants, 1); sync_provenance!(fx)
    @test failed_gates(fx) == ["robustness flags", "flags against grid"]
    # and 'fragile' is a claim that one did
    fx = base_fixture(); set_rob!(fx, CHEM, "chem-x", :breaking_variants, "")
    set_rob!(fx, CHEM, "chem-x", :n_breaking_variants, 0); sync_provenance!(fx)
    @test failed_gates(fx) == ["robustness flags", "flags against grid"]
    # completeness is checked on a depleted row too, not only on an enriched one: drop one
    # variant row and record it, so only the completeness condition is left to fire
    fx = base_fixture()
    deleteat!(fx.grid, findfirst(r -> r.variant == "W:three" && r.unit == HAB && r.category == "cat-a",
                                 eachrow(fx.grid)))
    set_rob!(fx, HAB, "cat-a", :variants_with_a_row, 2)
    @test failed_gates(fx) == ["robustness flags"]
end

@testset "a flag is recomputed from the grid rows it summarises" begin
    # the identity the fixture already satisfies, stated directly
    inp = load(base_fixture()).inp
    @test isempty(HF.check_flags_against_grid(inp))
    @test HF.grid_sign(0.2, 0.9) == "+" && HF.grid_sign(-0.9, -0.2) == "-"
    @test HF.grid_sign(-0.2, 0.9) == "0" && HF.grid_sign(NaN, NaN) == "undecidable"
    @test HF.flag_from_grid("+", [("a", "+"), ("b", "+")]) == ("+", String[])
    @test HF.flag_from_grid("+", [("b", "0"), ("a", "+")]) == ("fragile", ["b"])
    @test HF.flag_from_grid("0", [("a", "+"), ("b", "0")]) == ("null", ["a"])
    @test HF.flag_from_grid("undecidable", [("a", "0")]) == ("undecidable", ["a"])
    # A ONE-FILE EDIT OF THE GRID. The flags and the provenance still agree with each other.
    fx = base_fixture(); set_grid!(fx, LN, HAB, "cat-c", :lo, -0.05; variant = "V:two")
    @test failed_gates(fx) == ["flags against grid"]
    # A TWO-FILE EDIT THAT AGREES WITH ITSELF: the flag and its breakers moved from a fragile
    # row onto a signed one, consistently, in the robustness file and in the provenance.
    fx = base_fixture()
    set_rob!(fx, CHEM, "chem-x", :robust_sign, "+")
    set_rob!(fx, CHEM, "chem-x", :breaking_variants, ""); set_rob!(fx, CHEM, "chem-x", :n_breaking_variants, 0)
    set_rob!(fx, HAB, "cat-c", :robust_sign, "fragile")
    set_rob!(fx, HAB, "cat-c", :breaking_variants, "V:two"); set_rob!(fx, HAB, "cat-c", :n_breaking_variants, 1)
    sync_provenance!(fx)
    @test failed_gates(fx) == ["flags against grid"]
    # a variant row the grid does not hold is not a variant that ran
    fx = base_fixture()
    deleteat!(fx.grid, findfirst(r -> r.variant == "V:two" && r.unit == HAB && r.category == "cat-d",
                                 eachrow(fx.grid)))
    @test failed_gates(fx) == ["flags against grid"]
    # a registered set of a different size than the flag rows were decided on
    fx = base_fixture(); set_rob!(fx, HAB, "cat-d", :variants_expected, 4)
    @test "flags against grid" in failed_gates(fx)
    # with no registered set named, no flag can be checked and the gate says so
    fx = base_fixture(); fx.prov["variants_in_the_registered_set"] = ""
    @test "flags against grid" in failed_gates(fx)
end

@testset "the arm summary and the grid provenance are checked against themselves" begin
    for (key, value) in ("fg_coverage_denominator" => "321",
                         "bg_analysed_coverage_pct" => "78.2",
                         "bg_frame_coverage_numerator" => "801",
                         "n_bacteria_fg" => "287",
                         "variants_whose_axis_reaches_$(HAB)" => "3")
        fx = base_fixture()
        haskey(fx.summary, key) ? (fx.summary[key] = value) : (fx.prov[key] = value)
        @test failed_gates(fx) == ["summary identities"]
    end
    # what the exclusion rule removed, changed on its own: the sampled frame no longer accounts
    # for the background the figure draws, and the caption would print both halves
    for (key, value) in ("bg_removed_by_scan" => "361", "bg_frame_coverage_denominator" => "1001")
        fx = base_fixture()
        fx.summary[key] = value
        key == "bg_frame_coverage_denominator" &&
            (fx.summary["bg_frame_coverage_pct"] = string(100 * 800 / 1001))
        @test failed_gates(fx) == ["summary identities"]
    end
    # a coverage numerator that is not the arm the figure draws
    fx = base_fixture(); fx.summary["fg_coverage_numerator"] = "299"
    fx.summary["fg_coverage_pct"] = string(100 * 299 / 320)
    @test failed_gates(fx) == ["summary identities"]
    # The comparator arm's coverage is printed beside the figure as one fraction, so its two
    # halves are checked against each other and its numerator against the arm the ticks are
    # drawn from. Either half moving alone publishes a fraction of two different computations.
    fx = base_fixture(); fx.summary["nonln_coverage_denominator"] = "470"
    @test failed_gates(fx) == ["summary identities"]
    fx = base_fixture(); fx.summary["nonln_coverage_numerator"] = "390"
    fx.summary["nonln_coverage_pct"] = string(100 * 390 / 460)
    @test failed_gates(fx) == ["summary identities"]
    # an axis declaration that does not cover a registered variant, with the per-unit counts
    # left consistent, so the coverage condition is the only one able to fire
    fx = base_fixture()
    fx.prov["axis_unit_declaration"] = "V=$(HAB),$(CHEM);W2=$(CHEM)"
    fx.prov["variants_whose_axis_reaches_$(CHEM)"] = "2"
    bad = [g for g in HF.input_gates(load(fx).inp) if !g.ok]
    @test [g.name for g in bad] == ["summary identities"]
    @test occursin("W:three", bad[1].message) && !occursin("reaches", bad[1].message)
    # while the heaviest stratum is the phylum the two shares are about, they are its weight
    # and its background count; the fixture's heaviest stratum is another one, so the identity
    # is switched on here and then broken
    fx = base_fixture()
    fx.summary["standardisation_heaviest_stratum"] = SEL_PHYLUM
    fx.summary["pct_$(K_PHYLUM)_fg"] = "62.5"       # 100 * weight 0.625
    fx.summary["pct_$(K_PHYLUM)_bg"] = "16.0"       # 100 * 80 background taxa / n2 500
    @test isempty(failed_gates(fx))
    for (key, value) in ("pct_$(K_PHYLUM)_fg" => "61.25",
                         "pct_$(K_PHYLUM)_bg" => "12.5",
                         "standardisation_heaviest_stratum_bg_taxa" => "81")
        fy = base_fixture()
        fy.summary["standardisation_heaviest_stratum"] = SEL_PHYLUM
        fy.summary["pct_$(K_PHYLUM)_fg"] = "62.5"
        fy.summary["pct_$(K_PHYLUM)_bg"] = "16.0"
        fy.summary[key] = value
        @test failed_gates(fy) == ["summary identities"]
    end
end

@testset "cells of the grid's registered specification beyond the standardised estimate" begin
    for (col, value) in ((:b, 36.0), (:crude_lo, -0.31), (:crude_hi, 0.51), (:n2, N2 + 2),
                         (:n_zero_stratum_cells, 2.0), (:min_bg_stratum_cell, 0.5))
        fx = base_fixture(); set_grid!(fx, LN, HAB, "cat-d", col, value)
        @test failed_gates(fx) == ["sensitivity primary cells"]
    end
    # an enrichment row with no registered-specification row at all. The gate has to NAME the
    # row it could not find: a version that indexes the absent row instead fails too, with an
    # exception for a message, and a failure that is really a crash is not a gate.
    fx = base_fixture()
    deleteat!(fx.grid, findfirst(r -> r.variant == "primary" && r.contrast == NON &&
                                      r.unit == HAB && r.category == "cat-b", eachrow(fx.grid)))
    bad = [g for g in HF.input_gates(load(fx).inp) if !g.ok]
    @test [g.name for g in bad] == ["sensitivity primary cells"]
    @test occursin("0 registered-specification rows", bad[1].message)
    @test occursin("cat-b", bad[1].message) && !occursin("could not be evaluated", bad[1].message)
end

@testset "the robustness file's own copy of the primary interval is checked in full" begin
    for (col, value) in ((:primary_lo, -0.25 - 1e-6), (:primary_hi, 0.95 + 1e-6),
                         (:primary_a, 25.75 + 1e-6))
        fx = base_fixture(); set_rob!(fx, HAB, "cat-d", col, value)
        @test failed_gates(fx) == ["robustness values"]
    end
end

@testset "a digest that was never taken is not a digest" begin
    fx = base_fixture(); fx.summary["background_scan_sha256"] = ""; fx.prov["scan_sha256"] = ""
    @test failed_gates(fx) == ["scan digest"]
    fx = base_fixture(); fx.prov["calls_sha256"] = "not-a-digest"
    @test failed_gates(fx) == ["calls digest"]
end

@testset "a panel file that exists but holds nothing is not a panel file" begin
    fx = base_fixture(); fx.panels["panel.png"] = UInt8[]
    @test failed_gates(fx) == ["panel files"]
end

# ---------------------------------------------------------------- rows
@testset "the minimum count compares the weighted count, and the boundary is inclusive" begin
    fx = base_fixture()
    set_stat!(fx, LN, HAB, "cat-d", :a, 19.999)
    set_stat!(fx, LN, HAB, "cat-e", :a, 20.0)
    sel = HF.select_rows(load(fx).inp.stats)
    @test sort(sel.habitat.habitat) == ["cat-a", "cat-b", "cat-c", "cat-e"]
    @test sel.held_back.habitat == ["cat-d"]
    @test HF.MIN_N == 20
end

@testset "rows are sorted on the standardised estimate with a deterministic tie-break" begin
    fx = base_fixture()
    for c in ("cat-a", "cat-b", "cat-c", "cat-d")
        set_stat!(fx, LN, HAB, c, :log2_rr_std, 0.25)
    end
    st = load(fx).inp.stats
    orders = [HF.select_rows(st[perm, :]).habitat.habitat
              for perm in (collect(1:nrow(st)), reverse(collect(1:nrow(st))),
                           vcat(collect(2:2:nrow(st)), collect(1:2:nrow(st))))]
    @test all(==(["cat-a", "cat-b", "cat-c", "cat-d"]), orders)
    base = HF.select_rows(load(base_fixture()).inp.stats)
    @test base.habitat.habitat == ["cat-a", "cat-b", "cat-d", "cat-c"]
    @test base.chemistry.habitat == ["chem-y", "chem-x"]
end

@testset "a tick is keyed on unit and category, never linked across units" begin
    st = DataFrame([
        stat_row(LN, HAB, "cat-a", 30.0, N1, 40.0, N2, 0.1, -0.1, 0.3, 0.2, -0.2, 0.6),
        stat_row(LN, CHEM, "cat-a", 30.0, N1, 40.0, N2, 0.1, -0.1, 0.3, -0.4, -0.9, 0.1),
        stat_row(NON, HAB, "cat-a", 30.0, NNON, 40.0, N2, 0.1, -0.1, 0.3, 0.55, 0.1, 1.0),
        stat_row(NON, CHEM, "cat-a", 30.0, NNON, 40.0, N2, 0.1, -0.1, 0.3, -0.65, -1.1, -0.2)])
    rob = DataFrame([rob_row(HAB, "cat-a", "null", "false", "", 0.2, -0.2, 0.6, 30.0),
                     rob_row(CHEM, "cat-a", "null", "false", "", -0.4, -0.9, 0.1, 30.0)])
    rob.primary_excludes_zero = string.(rob.primary_excludes_zero)
    st.tested = string.(st.tested)
    plan = HF.figure_plan((stats = st, robustness = rob), HF.select_rows(st))
    tick(unit) = only(r.tick for r in plan.rows if r.unit == unit)
    @test tick(HAB) == 0.55
    @test tick(CHEM) == -0.65
    @test HF.family_tick(st, CHEM, "cat-a").log2_rr_std == -0.65
    @test_throws ErrorException HF.family_tick(st, "sequences", "cat-a")
    # the shared name also refuses to export colliding keys
    @test_throws ErrorException HF.require_drawn(plan, "cat-a")
end

# ---------------------------------------------------------------- limits and geometry
@testset "axis limits cover every drawn mark" begin
    fx = base_fixture()
    set_stat!(fx, LN, HAB, "cat-a", :b, 0.43 * N2)            # a background share above any old constant
    set_stat!(fx, LN, HAB, "cat-a", :pct_bg, 43.0)
    set_stat!(fx, LN, HAB, "cat-b", :log2_rr, -2.7)           # a crude ring beyond every interval
    set_stat!(fx, NON, CHEM, "chem-x", :log2_rr_std, 3.1)     # a tick beyond every interval
    inp = load(fx).inp
    plan = HF.figure_plan(inp, HF.select_rows(inp.stats))
    @test plan.share_max >= maximum(max(r.pct_fg, r.pct_bg) for r in plan.rows)
    @test plan.share_max == 50 && plan.share_step == 10
    @test plan.xb_lo < -2.7 && plan.xb_hi > 3.1
    for r in plan.rows
        @test plan.xb_lo < min(r.lo, r.crude, r.tick) && max(r.hi, r.crude, r.tick) < plan.xb_hi
        @test plan.ylo < r.y + HF.RING_DY && r.y + HF.tick_top() < plan.ymax
    end
    @test HF.share_axis(41.2) == (top = 50, step = 10)
    @test HF.share_axis(7.3) == (top = 8, step = 2)
    @test HF.share_axis(0.0).top > 0
    @test HF.enrichment_axis([0.3, 0.9]).lo < 0 < HF.enrichment_axis([-0.9, -0.3]).hi   # zero stays in view
    @test_throws ErrorException HF.enrichment_axis([0.3, NaN])
end

@testset "sub-row clearances are positive at the declared row pitch" begin
    @test HF.tick_dot_clearance_pt() > 0
    @test HF.dot_ring_clearance_pt() > 0
    @test HF.row_to_row_clearance_pt() > 0
    @test HF.bottom_clearance_pt() > 0
    @test HF.habitat_label_clearance_pt() > 0
    @test HF.chemistry_label_clearance_pt() > 0
    @test HF.TICK_DY - HF.TICK_HALF > HF.DOT_DY > HF.RING_DY       # tick above dot above ring
    inp = load(base_fixture()).inp
    plan = HF.figure_plan(inp, HF.select_rows(inp.stats))
    @test plan.axis_height_pt ≈ (plan.ymax - plan.ylo) * HF.ROW_MM * HF.PT_PER_MM
    @test plan.width_pt ≈ HF.WIDTH_MM * HF.PT_PER_MM
    # a tick at exactly the dot's estimate still clears the dot: the clearance does not depend on x
    r = first(plan.rows)
    tick_bottom_pt = (r.y + HF.tick_bottom()) * HF.row_pt()
    dot_top_pt = (r.y + HF.DOT_DY) * HF.row_pt() + HF.dot_half_pt()
    @test tick_bottom_pt - dot_top_pt > 0
    # block labels sit in data coordinates, under the axis top and under the rule
    labels = Dict(l.text => l for l in plan.block_labels)
    @test labels[HF.LABEL_HABITAT].y == plan.ymax - HF.LABEL_DY
    @test labels[HF.LABEL_CHEMISTRY].y == plan.yrule - HF.LABEL_DY
    @test all(l -> plan.xb_lo < l.x < plan.xb_hi, plan.block_labels)
end

@testset "colours and legend follow the classes the rows use" begin
    inp = load(base_fixture()).inp
    plan = HF.figure_plan(inp, HF.select_rows(inp.stats))
    cls = Dict(r.category => r.colour_class for r in plan.rows)
    @test cls == Dict("cat-a" => :depleted, "cat-b" => :neutral, "cat-c" => :enriched,
                      "cat-d" => :neutral, "chem-x" => :neutral, "chem-y" => :neutral)
    # THE THREE INKS, BY VALUE. The class-to-colour map is the last step the suite can reach:
    # draw_figure is deliberately not run here, so a swap inside the map or a hex collapsed
    # onto another would otherwise be visible for the first time on the page. They are the
    # Okabe-Ito pair and a neutral grey, and they are the sensitivity panel's three inks.
    @test HF.COL_ENRICHED == "#D55E00"
    @test HF.COL_DEPLETED == "#0072B2"
    @test HF.COL_NEUTRAL == "#7F7F7F"
    @test length(Set((HF.COL_ENRICHED, HF.COL_DEPLETED, HF.COL_NEUTRAL))) == 3
    @test HF.CLASS_COLOURS.enriched == HF.COL_ENRICHED
    @test HF.CLASS_COLOURS.depleted == HF.COL_DEPLETED
    @test HF.CLASS_COLOURS.neutral == HF.COL_NEUTRAL
    # and the ink each drawn row actually carries, one row of each class: enriched, depleted,
    # and a fragile row, whose primary interval excludes zero and which is grey all the same
    drawn_colour = Dict(r.category => r.colour for r in plan.rows)
    @test drawn_colour["cat-c"] == HF.COL_ENRICHED
    @test drawn_colour["cat-a"] == HF.COL_DEPLETED
    @test drawn_colour["chem-x"] == HF.COL_NEUTRAL
    @test drawn_colour["cat-b"] == HF.COL_NEUTRAL
    @test first.(plan.legend.marks) == [:enriched, :depleted, :neutral, :interval, :ring, :tick]
    @test first.(HF.legend_entries(Set([:neutral])).marks) == [:neutral, :interval, :ring, :tick]
    @test first.(HF.legend_entries(Set([:enriched, :neutral])).marks) == [:enriched, :neutral, :interval, :ring, :tick]
    labels = Dict(HF.legend_entries(Set([:enriched, :depleted, :neutral])).marks)
    @test labels[:enriched] == "robust, enriched" && labels[:depleted] == "robust, depleted" &&
          labels[:neutral] == "not robust"
    @test [r.label for r in HF.panel_order(plan)] ==
          ["cat-c (70)", "cat-d (25.8)", "cat-b (90.2)", "cat-a (40.5)", "chem-x (60)", "chem-y (33)"]
end

# ---------------------------------------------------------------- formatting
@testset "one formatter per kind of printed number" begin
    @test HF.fmt_count(12.0) == "12"
    @test HF.fmt_count(7.26) == "7.3"
    @test HF.fmt_count(3.04) == "3.0"
    @test_throws ErrorException HF.fmt_count(NaN)
    @test HF.fmt_int(1234567) == "1,234,567"
    @test HF.fmt_int(987) == "987"
    @test HF.fmt_int(-4321) == "-4,321"
    @test HF.fmt_pct(12.36) == "12.4"
    @test HF.fmt_log2(0.456) == "+0.46"
    @test HF.fmt_log2(-0.004) == "0.00"
    @test HF.fmt_interval(-0.123, 0.987) == "-0.12 to +0.99"
    @test HF.fmt_p(0.3172) == "0.32"
    @test HF.fmt_p(0.00064) == "6 × 10⁻⁴"
    @test HF.fmt_p(0.0000097) == "1 × 10⁻⁵"
    @test HF.fmt_p(2.4e-13) == "2 × 10⁻¹³"
    @test_throws ErrorException HF.fmt_p(1.5)
    @test HF.prose_list(String[]) == "none"
    @test HF.prose_list(["x"]) == "x"
    @test HF.prose_list(["x", "y"]) == "x and y"
    @test HF.prose_list(["x", "y", "z"]) == "x, y and z"
end

# ---------------------------------------------------------------- exported values
@testset "exported keys are exactly the declared keys, and every value is computed" begin
    inp = load(base_fixture()).inp
    d = HF.derive(inp)
    kv = Dict(d.values)
    order = [r.category for r in HF.panel_order(d.plan)]
    @test first.(d.values) == HF.caption_keys(order, SELECTED)
    # The keys that carry a selected name are DERIVED from the selection, not written down.
    @test haskey(kv, "$(K_TIER)_tier_code") && kv["$(K_TIER)_tier_code"] == "code-t"
    @test haskey(kv, "pct_$(K_PHYLUM)_fg")
    # A different selection declares keys under THAT selection's names, and neither
    # selection's names appear in the other's key set — so a name reaches a key only by
    # being supplied, whichever selection is supplied.
    alt = HF.HE.read_category_selection(selection_csv(mktempdir(); tier = "cat-w",
                                                      phylum = "phy-9", name = "alt.csv"))
    altk = HF.caption_keys(order, alt)
    @test "catw_tier_code" in altk && "pct_phy9_fg" in altk
    @test !("$(K_TIER)_tier_code" in altk) && !("pct_$(K_PHYLUM)_fg" in altk)
    @test !("catw_tier_code" in first.(d.values)) && !("pct_phy9_fg" in first.(d.values))
    @test length(unique(first.(d.values))) == length(d.values)
    @test all(k -> occursin(HF.KEY_PATTERN, k), first.(d.values))
    @test kv["n1"] == "300" && kv["n2"] == "500" && kv["min_n"] == "20"
    @test kv["drawn_habitat"] == "cat-c, cat-d, cat-b and cat-a"
    @test kv["drawn_chemistry"] == "chem-x and chem-y"
    @test kv["held_back"] == "cat-e" && kv["held_back_with_a"] == "cat-e (6.5)"
    @test kv["robust_enriched"] == "cat-c" && kv["robust_depleted"] == "cat-a"
    @test kv["any_robust"] == "true" && kv["robust_not_drawn"] == "none"
    @test kv["coloured_enriched"] == "cat-c" && kv["coloured_depleted"] == "cat-a"
    @test kv["fragile"] == "chem-x" && kv["fragile_with_breakers"] == "chem-x (V:two)"
    @test kv["legend_only"] == "chem-x and cat-d"          # read from the summary, in its order
    fx = base_fixture(); fx.summary["categories_legend_only"] = ""
    @test Dict(HF.derive(load(fx).inp).values)["legend_only"] == "none"
    @test kv["drawn_rows_with_a_zero_stratum_cell"] == "cat-d"
    @test kv["narrowest_ci_category"] == "cat-b" && kv["narrowest_ci_width"] == "0.50"
    @test kv["runner_up_ci_category"] == "chem-x" && kv["runner_up_ci_width"] == "0.70"
    @test kv["largest_fg_share_category"] == "cat-b"
    @test kv["largest_bg_share_category"] == "cat-b"
    # the drawn row whose ring is furthest from its dot, in either direction
    @test kv["largest_adjustment_shift_category"] == "cat-b"
    @test kv["largest_adjustment_shift"] == "0.40"
    @test kv["n_specifications"] == "4" && kv["registered_variants"] == "V:one;V:two;W:three"
    @test kv["variants_reaching_taxa"] == "2" && kv["variants_reaching_chemistry"] == "3"
    @test kv["fg_isolation_terms_total"] == "1,200"
    @test kv["coverage_fisher_p_after_scan"] == "6 × 10⁻⁴"
    @test kv["standardisation_heaviest_stratum_weight_pct"] == "62.5"
    @test kv["panel_categories_not_drawn"] == "cat-z"
    # the background is the frame minus what the scan removed, and the comparator arm has a
    # coverage of its own: both are printed beside the figure, so both are exported
    @test kv["bg_removed_by_scan"] == "360"
    @test kv["nonln_coverage_numerator"] == "400" && kv["nonln_coverage_denominator"] == "460"
    @test kv["cat-b.taxa_specs_excluding_zero"] == "1"
    @test kv["cat-b.taxa_specs_distinct"] == "4"
    @test kv["cat-b.sequences_specs_excluding_zero"] == "3"
    @test kv["cat-b.sequences_specs_excluding_zero_below"] == "3"
    @test kv["cat-b.sequences_specs_not_excluding_zero"] == "V:two"
    # V:one returns the primary bit for bit, so four specifications are three answers
    @test kv["cat-b.sequences_specs_distinct"] == "3"
    @test kv["cat-b.seq_ci"] == "-0.70 to -0.20"
    @test kv["cat-c.sharing_ratio"] == "0.50" && kv["cat-c.sharing_ratio_ci"] == "0.25 to 0.75"
    @test kv["cat-c.fam"] == "+0.40" && kv["cat-c.val"] == "+0.80" && kv["cat-c.colour"] == "enriched"
    # the exported colour of a depleted row and of a fragile one, not only of the coloured
    # enriched row: a fragile row's interval excludes zero and its colour key is still neutral
    @test kv["cat-a.colour"] == "depleted" && kv["cat-a.robust_sign"] == "-"
    @test kv["chem-x.colour"] == "neutral" && kv["chem-x.robust_sign"] == "fragile"
    @test kv["cat-c.rounded_for_fisher"] == "true"
    @test kv["cat-d.a"] == "25.8" && kv["cat-d.n_zero_stratum_cells"] == "1"
    # ranks within the block: the crude order and the adjusted order differ
    @test [kv["$c.rank_crude"] for c in ("cat-c", "cat-b", "cat-d", "cat-a")] == ["1", "2", "3", "4"]
    @test [kv["$c.rank_std"] for c in ("cat-c", "cat-d", "cat-b", "cat-a")] == ["1", "2", "3", "4"]
    @test kv["chem-x.rank_crude"] == "1" && kv["chem-y.rank_std"] == "2"
    @test kv["chem-x.breakers"] == "V:two"
    @test kv["chem-x.breaker_estimates"] == "V:two +0.05 (-0.20 to +0.30)"
    @test kv["chem-y.breaker_estimates"] == "W:three -0.55 (-1.00 to -0.10)"
    @test kv["cat-a.breaker_estimates"] == "none"
    @test !haskey(kv, "cat-e.val")                            # a held-back row exports no row keys

    # a category the named exports need but the figure does not draw is an error naming it
    err = try
        HF.caption_values(with_selection(inp; grid = ("cat-q",)), d.sel, d.plan)
    catch e
        e
    end
    @test err isa ErrorException && occursin("cat-q", err.msg)
    err = try
        HF.caption_values(with_selection(inp; sharing = "cat-e"), d.sel, d.plan)
    catch e
        e
    end
    @test err isa ErrorException && occursin("cat-e", err.msg)
    # a sharing ratio that is not the drawn tick over the drawn estimate is refused
    fx = base_fixture(); fx.summary["sharing_ratio_nonln_over_ln"] = "0.6"
    inp2 = load(fx).inp
    @test_throws ErrorException HF.derive(inp2)
    # the robust lists are checked against the provenance again at export, not only by a gate
    fx = base_fixture(); fx.prov["robust_enriched"] = "cat-c;cat-d"
    err = try
        HF.derive(load(fx).inp)
    catch e
        e
    end
    @test err isa ErrorException && occursin("robust_enriched", err.msg)
    # a missing summary key is named
    fx = base_fixture(); delete!(fx.summary, "obt_release")
    err = try
        HF.derive(load(fx).inp)
    catch e
        e
    end
    @test err isa ErrorException && occursin("obt_release", err.msg)
    # a grid category missing one registered variant is named. The per-sequence unit carries
    # no robustness row, so no gate reaches it: the export is the only thing that can.
    fx = base_fixture()
    deleteat!(fx.grid, findfirst(r -> r.variant == "V:two" && r.unit == "sequences", eachrow(fx.grid)))
    err = try
        HF.derive(load(fx).inp)
    catch e
        e
    end
    @test err isa ErrorException && occursin("V:two", err.msg)
    # keys outside the alphabet the text renderer reads are refused
    @test_throws ErrorException HF.caption_keys(["bad/name"], SELECTED)
end

@testset "the thinnest per-stratum background counts are computed from the drawn rows" begin
    vals(fx) = Dict(HF.derive(load(fx).inp).values)
    # The list is the ascending head of the DRAWN rows with their counts, and beside it the
    # smallest count among the drawn rows it leaves out, so a sentence built on it cannot name
    # one row while a thinner one goes unmentioned. cat-e's cell is thinner than three of the
    # listed ones and belongs in neither: it is not drawn.
    kv = vals(base_fixture())
    @test kv["thinnest_bg_stratum_rows"] == "cat-d (0), cat-c (2.8), cat-a (3.5) and chem-y (4)"
    @test kv["thinnest_bg_stratum_outside_the_list"] == "8"
    @test !occursin("cat-e", kv["thinnest_bg_stratum_rows"])

    # A rerun that makes ANOTHER drawn row the thin one changes the list and the bound; a
    # hand-typed sentence would have gone on naming the same rows.
    fx = base_fixture()
    set_stat!(fx, LN, HAB, "cat-b", :min_bg_stratum_cell, 0.5)
    set_grid!(fx, LN, HAB, "cat-b", :min_bg_stratum_cell, 0.5)
    kv = vals(fx)
    @test kv["thinnest_bg_stratum_rows"] == "cat-d (0), cat-b (0.5), cat-c (2.8) and cat-a (3.5)"
    @test kv["thinnest_bg_stratum_outside_the_list"] == "4"

    # A tie at the list's last place carries both rows into the list, so the bound is always
    # above every count the list names.
    fx = base_fixture()
    set_stat!(fx, LN, CHEM, "chem-x", :min_bg_stratum_cell, 4.0)
    set_grid!(fx, LN, CHEM, "chem-x", :min_bg_stratum_cell, 4.0)
    kv = vals(fx)
    @test kv["thinnest_bg_stratum_rows"] ==
          "cat-d (0), cat-c (2.8), cat-a (3.5), chem-x (4) and chem-y (4)"
    @test kv["thinnest_bg_stratum_outside_the_list"] == "12"

    # The bound is printed as "no other drawn row rests on fewer than this", so it is rounded
    # DOWN: a row resting on 5.875 must not be described by a bound of 5.9, which it is below.
    fx = base_fixture()
    set_stat!(fx, LN, CHEM, "chem-x", :min_bg_stratum_cell, 5.875)
    set_grid!(fx, LN, CHEM, "chem-x", :min_bg_stratum_cell, 5.875)
    kv = vals(fx)
    @test kv["thinnest_bg_stratum_rows"] == "cat-d (0), cat-c (2.8), cat-a (3.5) and chem-y (4)"
    @test kv["thinnest_bg_stratum_outside_the_list"] == "5.8"
    @test parse(Float64, kv["thinnest_bg_stratum_outside_the_list"]) <= 5.875

    # With every drawn row inside the list there is nothing outside it, and the bound says so
    # rather than naming a count that does not exist.
    fx = base_fixture()
    for (u, c) in ((HAB, "cat-b"), (CHEM, "chem-x"))
        set_stat!(fx, LN, u, c, :min_bg_stratum_cell, 4.0)
        set_grid!(fx, LN, u, c, :min_bg_stratum_cell, 4.0)
    end
    kv = vals(fx)
    @test kv["thinnest_bg_stratum_rows"] ==
          "cat-d (0), cat-c (2.8), cat-a (3.5), cat-b (4), chem-x (4) and chem-y (4)"
    @test kv["thinnest_bg_stratum_outside_the_list"] == "none"
end

@testset "with no robust row nothing is coloured and the values say so" begin
    fx = base_fixture()
    # a flag is fragile because a variant broke it, so the grid has to say so too
    for cat in ("cat-a", "cat-c")
        set_rob!(fx, HAB, cat, :robust_sign, "fragile")
        set_rob!(fx, HAB, cat, :breaking_variants, "V:one")
        set_rob!(fx, HAB, cat, :n_breaking_variants, 1)
        set_grid!(fx, LN, HAB, cat, :std, 0.15; variant = "V:one")
        set_grid!(fx, LN, HAB, cat, :lo, -0.10; variant = "V:one")
        set_grid!(fx, LN, HAB, cat, :hi, 0.40; variant = "V:one")
    end
    sync_provenance!(fx)
    inp = load(fx).inp
    @test all(g -> g.ok, HF.input_gates(inp))
    d = HF.derive(inp)
    kv = Dict(d.values)
    @test all(r -> r.colour_class === :neutral, d.plan.rows)
    @test kv["n_robust"] == "0" && kv["any_robust"] == "false"
    @test kv["robust_enriched"] == "none" && kv["robust_depleted"] == "none"
    @test kv["coloured_enriched"] == "none"
    @test kv["fragile"] == "cat-c, cat-a and chem-x"
    @test first.(d.plan.legend.marks) == [:neutral, :interval, :ring, :tick]
end

@testset "source data: every Ln, non-Ln and all-carrier row, with its role" begin
    inp = load(base_fixture()).inp
    d = HF.derive(inp)
    sd = d.source
    @test sort(unique(sd.contrast)) == sort([LN, NON, ALL])
    @test nrow(sd) == count(r -> r.contrast in (LN, NON, ALL), eachrow(inp.stats))
    role(contrast, unit, cat) = only(sd.role[(sd.contrast .== contrast) .& (sd.unit .== unit) .& (sd.habitat .== cat)])
    @test role(LN, HAB, "cat-b") == HF.ROLE_PLOTTED
    @test role(LN, HAB, "cat-e") == HF.ROLE_HELD_BACK
    @test role(LN, "sequences", "cat-b") == HF.ROLE_UNIT_CHECK
    @test role(NON, CHEM, "chem-x") == HF.ROLE_TICK
    @test role(NON, HAB, "cat-e") == HF.ROLE_TICK_HELD
    @test role(ALL, HAB, "cat-c") == HF.ROLE_ALL_CARRIER
    row(contrast, unit, cat) = sd[findfirst((sd.contrast .== contrast) .& (sd.unit .== unit) .& (sd.habitat .== cat)), :]
    @test row(NON, CHEM, "chem-x").family_tick_of == HF.rowkey(CHEM, "chem-x")
    @test row(LN, HAB, "cat-c").colour_class == "enriched" && row(LN, HAB, "cat-c").drawn == "true"
    # the class travels from the drawn row, so a depleted and a fragile row cannot be filed
    # under a colour the figure does not give them
    @test row(LN, HAB, "cat-a").colour_class == "depleted"
    @test row(LN, CHEM, "chem-x").colour_class == "neutral"
    @test row(LN, HAB, "cat-c").panel_row == "1" && row(LN, CHEM, "chem-y").panel_row == "6"
    @test row(LN, CHEM, "chem-x").robust_sign == "fragile" && row(LN, CHEM, "chem-x").breaking_variants == "V:two"
    @test row(LN, HAB, "cat-e").drawn == "false" && row(LN, HAB, "cat-e").colour_class == ""
    @test row(LN, HAB, "cat-d").a == "25.75"                  # cells travel as read
    @test all(==(""), sd.robust_sign[sd.contrast .!= LN])
    @test "min_bg_stratum_cell" in names(sd) && "n_zero_stratum_cells" in names(sd)
end

# ---------------------------------------------------------------- staging and publishing
@testset "a failed gate writes nothing; a passing run publishes exactly the declared outputs" begin
    # failing: the figure directory, holding an earlier output, is byte-identical afterwards
    fx = base_fixture(); fx.calls = Vector{UInt8}("another call table\n")
    l = load(fx)
    write(joinpath(l.paths.fig, "habitat_enrichment.png"), "earlier run")
    before = snapshot(l.paths.fig)
    @test run_fixture(l.opt) == 1
    @test snapshot(l.paths.fig) == before
    @test HF.main(l.paths.args; env = NO_ENV, render = fake_render, io = devnull, err = IOBuffer()) == 1
    @test snapshot(l.paths.fig) == before

    # check only: every gate passes, nothing is written
    l = load(base_fixture())
    before = snapshot(l.paths.fig)
    opt = HF.parse_run_options(vcat(l.paths.args, "--check-only"); env = NO_ENV)
    @test opt.check_only
    @test run_fixture(opt) == 0
    @test snapshot(l.paths.fig) == before

    # passing: exactly the declared outputs, the staging directory gone, earlier files replaced
    write(joinpath(l.paths.fig, "habitat_enrichment.png"), "earlier run")
    write(joinpath(l.paths.fig, "unrelated.txt"), "left alone")
    @test run_fixture(l.opt) == 0
    @test sort(readdir(l.paths.fig)) == sort(vcat(collect(HF.FIGURE_OUTPUTS), "unrelated.txt"))
    @test read(joinpath(l.paths.fig, "unrelated.txt"), String) == "left alone"
    @test startswith(read(joinpath(l.paths.fig, "habitat_enrichment.png"), String), "placeholder png")
    @test !isdir(joinpath(l.paths.fig, HF.STAGING_DIR))
    vals = CSV.read(joinpath(l.paths.fig, HF.VALUES_FILE), DataFrame; types = String)
    @test vals.key == HF.caption_keys([r.category for r in HF.panel_order(HF.derive(l.inp).plan)],
                                      SELECTED)
    prov = Dict(zip(CSV.read(joinpath(l.paths.fig, HF.PROVENANCE_FILE), DataFrame; types = String).key,
                    coalesce.(CSV.read(joinpath(l.paths.fig, HF.PROVENANCE_FILE), DataFrame; types = String).value, "")))
    @test prov["input_calls_sha256"] == bytes2hex(sha256(base_fixture().calls))
    @test prov["input_stats_given"] == "explicit" && prov["figdir_given"] == "explicit"
    @test prov["input_gates_failed"] == "" && prov["gate_calls digest"] == "pass"
    @test prov["julia_version"] == string(VERSION)
    @test prov["png_dpi"] == "300"
    # THE PUBLISHED PROVENANCE NAMES FILES, NOT LOCATIONS. The digest is what pins the bytes;
    # an absolute path pins nothing and would carry the run's directory layout into a file
    # that travels with the figure. Nothing here may name the directories this run used.
    @test prov["input_calls_file"] == "calls.csv" && prov["input_stats_file"] == "stats.csv"
    @test prov["input_sensitivity_provenance_file"] == "sprov.csv"
    @test !haskey(prov, "figdir") && !haskey(prov, "panel_dir")
    @test !any(v -> occursin(l.paths.src, v) || occursin(l.paths.fig, v), values(prov))
    @test !any(v -> HF.looks_like_a_local_path(v), values(prov))
    @test HF.looks_like_a_local_path(l.paths.fig) && HF.looks_like_a_local_path("/figs/x.csv")
    @test !HF.looks_like_a_local_path("scripts/habitat_figure.jl")
    @test !HF.looks_like_a_local_path("habitat_robustness.csv")
    # the guard itself, on rows it would have to refuse. It cannot be driven through a run:
    # on well-formed rows it is invisible, and the visible condition is the basenames above.
    @test HF.refuse_local_paths(["a" => "x.csv", "b" => "scripts/y.jl"]) ==
          ["a" => "x.csv", "b" => "scripts/y.jl"]
    for bad in (raw"C:\x\y.csv", "/figs/y.csv", raw"dir\y.csv")
        err = try
            HF.refuse_local_paths(["a" => "x.csv", "figdir" => bad]); nothing
        catch e
            e
        end
        @test err isa ErrorException && occursin("figdir", err.msg)
    end
end

@testset "a renderer that stages the wrong set publishes nothing" begin
    for render in ((plan, base) -> [write(base * ".pdf", "x"), write(base * ".png", "x")],      # no svg
                   (plan, base) -> (fake_render(plan, base); write(base * ".eps", "x")),     # an extra file
                   (plan, base) -> (fake_render(plan, base); write(base * ".svg", "")))      # an empty file
        l = load(base_fixture())
        write(joinpath(l.paths.fig, "habitat_enrichment.png"), "earlier run")
        err = try
            run_fixture(l.opt; render = render); nothing
        catch e
            e
        end
        @test err isa ErrorException && occursin("nothing was published", err.msg)
        @test read(joinpath(l.paths.fig, "habitat_enrichment.png"), String) == "earlier run"
        @test isdir(joinpath(l.paths.fig, HF.STAGING_DIR))                    # left for inspection
        @test sort(filter(f -> f != HF.STAGING_DIR, readdir(l.paths.fig))) == ["habitat_enrichment.png"]
    end
    # a renderer that fails outright publishes nothing either
    l = load(base_fixture())
    before = filter(p -> !startswith(first(p), HF.STAGING_DIR), snapshot(l.paths.fig))
    @test_throws ErrorException run_fixture(l.opt; render = (plan, base) -> error("renderer broke"))
    @test filter(p -> !startswith(first(p), HF.STAGING_DIR), snapshot(l.paths.fig)) == before
    # the output gate itself, on a staging directory built by hand
    stage = mktempdir()
    for f in HF.FIGURE_OUTPUTS
        write(joinpath(stage, f), "x")
    end
    @test HF.output_gate(stage).ok
    write(joinpath(stage, "stray.csv"), "x")
    g = HF.output_gate(stage)
    @test !g.ok && occursin("stray.csv", g.message) && g.name == HF.OUTPUT_GATE
    # staging is cleared when opened, and never creates a figure directory
    fig = mktempdir()
    s = HF.open_staging(fig); write(joinpath(s, "left over"), "x")
    @test isempty(readdir(HF.open_staging(fig)))
    absent = joinpath(mktempdir(), "absent")
    @test_throws ErrorException HF.open_staging(absent)
    @test !ispath(absent)
end

# ---------------------------------------------------------------- the command line
@testset "options are strict and the figure directory is never created" begin
    dir = mktempdir()
    cats = ["--categories", SEL_PATH]
    # THE SELECTION IS REQUIRED AND HAS NO DEFAULT: the sharing category, the grid-count
    # categories and the standardisation phylum were chosen from this data set's estimates, so
    # a run that names no criteria file is refused rather than given a set this file picked.
    err = try
        HF.parse_run_options(["--figdir", dir]; env = NO_ENV)
    catch e
        e
    end
    @test err isa ErrorException && occursin("--categories", err.msg)
    @test HF.parse_run_options([cats..., "--figdir", dir]; env = NO_ENV).categories == SEL_PATH
    @test HF.parse_run_options([cats..., "--figdir", dir]; env = NO_ENV).figdir_given == "explicit"
    o = HF.parse_run_options(cats; env = Dict("LNFINDER_FIGDIR" => dir))
    @test o.figdir == abspath(dir) && o.figdir_given == "environment"
    @test o.robustness.path == joinpath(abspath(dir), HF.ROBUSTNESS_FILE) && o.robustness.given == "default"
    @test o.panel_dir.path == abspath(dir) && !o.check_only
    @test HF.parse_run_options(["--help"]; env = NO_ENV).help
    @test_throws ErrorException HF.parse_run_options(String[]; env = NO_ENV)
    @test_throws ErrorException HF.parse_run_options(cats; env = NO_ENV)
    absent = joinpath(dir, "absent")
    @test_throws ErrorException HF.parse_run_options([cats..., "--figdir", absent]; env = NO_ENV)
    @test !ispath(absent)
    @test_throws ErrorException HF.parse_run_options([cats..., "--figdir", dir, "--figdir", dir];
                                                    env = NO_ENV)
    @test_throws ErrorException HF.parse_run_options([cats..., "--figdir", dir, "--stats"];
                                                    env = NO_ENV)
    @test_throws ErrorException HF.parse_run_options([cats..., "--figdir", dir, "--bogus"];
                                                    env = NO_ENV)
    err = IOBuffer()
    @test HF.main([cats..., "--figdir", absent]; env = NO_ENV, io = devnull, err = err) == 1
    @test occursin("does not exist", String(take!(err))) && !ispath(absent)
    out = IOBuffer()
    @test HF.main(["--help"]; env = NO_ENV, io = out) == 0
    @test occursin("--check-only", String(take!(out)))
end

@testset "from outside: --help loads nothing, a bad command line writes nothing" begin
    julia = Base.julia_cmd()
    clean_env = Dict(k => v for (k, v) in ENV if k != "LNFINDER_FIGDIR")
    function run_script(script, args)
        out, err = IOBuffer(), IOBuffer()
        cmd = setenv(`$julia --startup-file=no $script $args`, clean_env)
        p = run(pipeline(ignorestatus(cmd); stdout = out, stderr = err))
        return (code = p.exitcode, out = String(take!(out)), err = String(take!(err)))
    end
    # a copy in an empty tree: --help must exit before anything is loaded or written
    tree = mktempdir()
    mkpath(joinpath(tree, "scripts"))
    copy_ = joinpath(tree, "scripts", "habitat_figure.jl")
    cp(joinpath(@__DIR__, "..", "scripts", "habitat_figure.jl"), copy_)
    h = run_script(copy_, ["--help"])
    @test h.code == 0
    @test occursin("--figdir", h.out) && occursin("--sensitivity-provenance", h.out)
    @test !occursin("Activating", h.out * h.err)
    @test readdir(tree) == ["scripts"] && readdir(joinpath(tree, "scripts")) == ["habitat_figure.jl"]

    # the script itself: an unknown option, a missing criteria file and a missing figure
    # directory all exit 1, and none of them writes anything
    script = normpath(joinpath(@__DIR__, "..", "scripts", "habitat_figure.jl"))
    repo_figs = normpath(joinpath(@__DIR__, "..", "figs"))
    repo_figs_before = ispath(repo_figs)
    fig = mktempdir()
    cats = ["--categories", SEL_PATH]
    u = run_script(script, ["--figdir", fig, "--figdri", fig])
    @test u.code == 1 && occursin("unknown option '--figdri'", u.err)
    @test isempty(readdir(fig))
    # NO CATEGORIES, AND NO DEFAULT TO FALL BACK TO. The refusal names the option that
    # supplies them, and it comes before anything is read: a run held to no selection is not a
    # run with the selection this software would have picked, because it has none to pick.
    c = run_script(script, ["--figdir", fig])
    @test c.code == 1 && occursin("--categories", c.err)
    @test isempty(readdir(fig))
    absent = joinpath(fig, "absent")
    m = run_script(script, [cats..., "--figdir", absent])
    @test m.code == 1 && occursin("does not exist", m.err)
    n = run_script(script, cats)
    @test n.code == 1 && occursin("no figure directory", n.err)
    @test isempty(readdir(fig))
    @test ispath(repo_figs) == repo_figs_before
end

end
