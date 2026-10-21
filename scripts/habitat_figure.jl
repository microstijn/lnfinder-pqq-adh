#=
    habitat_figure.jl — where the carriers of lanthanide-type PQQ-dependent dehydrogenases
    live: the enrichment table drawn as two panels, coloured by the robustness flags.

    WHAT IS DRAWN. The two panels share one row per category.
      (a) The share of taxa in each category: the carrier arm as a solid bar inside the
          background arm's outline, on a common baseline, so an overshoot is read as a length
          against a fixed origin. One ink only: this panel runs no test, so it carries none of
          panel (b)'s colour.
      (b) log2 enrichment against the background, with three marks per row on three
          sub-rows. The crude estimate is an open ring below; the phylum-standardised estimate
          is a filled dot on its 95% interval; above both, a thin tick marks the same
          standardised estimate for the family comparator, the PQQ-ADH carriers with no
          lanthanide-type sequence. A connector joins ring and dot, so the distance between
          them is the contribution of phylum composition, drawn rather than asserted. It has
          no fixed direction: adjustment moves some rows toward zero and others away from it,
          and the connector is what tells them apart. Rows are sorted on the standardised
          estimate within each block, so the dots form a staircase and the rings are free not
          to.

    COLOUR IS THE ROBUSTNESS FLAG AND NOTHING ELSE. A dot is coloured only where the
    sensitivity grid's flag gives the row a sign that the primary specification and every
    one-at-a-time variant reproduce. An interval that excludes zero under the primary
    specification alone does not colour a row, so no interval test decides a colour here. A
    flag that contradicts its own interval stops the run; it is never drawn grey instead.

    THE TICK HAS A SUB-ROW OF ITS OWN. Where the comparator's estimate nearly coincides with
    the carriers' own, a tick drawn at the dot's height disappears under the dot and the
    figure shows one measurement where it made two. The offsets and marker sizes are named
    constants whose clearances the test suite checks as arithmetic, and the figure's height
    is derived from a row pitch, so the clearances that hold in the suite hold on the page.

    CHEMISTRY shares panel (b)'s axis as a second block below a rule: the same estimator in
    the same units. Its flags are not exclusive of each other, and the block label says so.

    A MINIMUM WEIGHTED COUNT decides which categories are drawn. It compares the fractional
    foreground count (a taxon whose terms tie is split between the tied categories), and it
    is blind to direction, so it cannot select on the outcome. The held-back categories are
    reported with their counts, and the sensitivity panel draws them beside the rest.

    NOTHING IS PUBLISHED UNLESS THE INPUTS DESCRIBE ONE ANALYSIS. The enrichment rows, their
    arm summary, the robustness flags, the sensitivity grid and the grid's provenance come
    from separate runs. Every identity between them is a named gate; every gate is evaluated
    and all failures are reported together, and a single failure stops the run before
    anything is staged. A flag is not taken on trust either: its sign and the variants it
    names as breaking it are recomputed from the grid rows it summarises, so a flag its own
    grid does not support stops the run whether one file was edited or two. The outputs are
    staged side by side and moved into the figure directory together, and only once the
    staged set is exactly the declared one.

    THE PUBLISHED PROVENANCE NAMES FILES, NOT LOCATIONS. It records each input's file name,
    whether the path was typed or defaulted, and its digest; the digests pin the bytes, so a
    local absolute path would add nothing and would travel with the file. A value that looks
    like one stops the run.

    WHICH CATEGORIES THIS FIGURE IS ABOUT IS AN INPUT, NOT A CONSTANT. The category whose
    sharing ratio is read back and the categories whose intervals are counted across every
    registered specification were chosen from this data set's own estimates, with the estimates
    in hand, after those estimates existed. The phylum whose two shares are quoted is not a
    choice of that kind: it is DETERMINED as the heaviest standardisation stratum, and the
    identity gate below asserts the weight identities only while the run's own heaviest stratum
    is still that phylum. It is supplied here for a different reason — the share keys are named
    after it and a reader has to be able to derive that name. So they are supplied per run in
    the criteria file --categories names — the same file step 3 is handed, read with step 3's
    own reader — this file holds no default for any of them, and a run without it is refused by
    name. The selection and the criteria file's digest go into the published provenance, and the
    keys that carry a selected name derive from it, so a run on other categories exports keys
    for THOSE categories rather than filing other numbers under names chosen here.

    Run:  julia --project=. scripts/habitat_figure.jl --figdir <dir> [options]
          (--help lists every option and exits before any package is loaded)
=#

const REPO = normpath(joinpath(@__DIR__, ".."))

const USAGE = """
Usage: julia --project=. scripts/habitat_figure.jl [options]

Options:
  --categories <path>              REQUIRED, and there is no default. The category-selection
                                   criteria file — the same one step 3 was given, so that the
                                   two steps cannot disagree about which categories the run is
                                   about. It names the category whose sharing ratio is read
                                   back, the categories whose intervals are counted across
                                   every specification, and the phylum whose shares are
                                   quoted. Those are selections made from this data set's own
                                   estimates, not constants of this software
  --figdir <path>                  the figure directory: every output goes there, and the
                                   inputs marked <figdir> are read from it by default
                                   (default: \$LNFINDER_FIGDIR). It must already exist; it is
                                   never created
  --stats <path>                   enrichment rows
                                   (default: <data>/habitat/habitat_enrichment.csv)
  --summary <path>                 arm summary, key,value
                                   (default: <data>/habitat/habitat_arm_summary.csv)
  --robustness <path>              robustness flags (default: <figdir>/habitat_robustness.csv)
  --sensitivity <path>             sensitivity grid (default: <figdir>/habitat_sensitivity.csv)
  --sensitivity-provenance <path>  provenance of that grid, key,value
                                   (default: <figdir>/habitat_sensitivity_provenance.csv)
  --calls <path>                   taxon-call table; only its digest is read
                                   (default: <data>/habitat/habitat_taxon_calls.csv)
  --panel-dir <path>               directory holding the sensitivity panel files the grid's
                                   provenance names (default: <figdir>)
  --check-only                     evaluate every input gate and derive every value, then
                                   stop: nothing is written, not even the staging directory
  -h, --help                       print this text and exit

Exit status: 0 on a published run or a passing check; 1 on a command-line error, a missing
input, a failed gate, or a value that could not be derived.
"""

# --help costs nothing: the enrichment module below takes seconds to load.
if abspath(PROGRAM_FILE) == (@__FILE__) && any(a -> a in ("-h", "--help"), ARGS)
    print(USAGE)
    exit(0)
end

using Pkg; Pkg.activate(REPO)
using DataFrames, CSV, Printf, Dates, SHA

# Where the step-3 tables live, by NAME; the layout is declared once, in
# src/ModuleDataInputs.jl, and documented once, in docs/DATA.md.
include(joinpath(REPO, "src", "ModuleDataInputs.jl"))

# The enrichment step, behind its own program-file guard: including it defines its functions
# and runs nothing. It supplies the strict option parser, the contrast names and the interval
# predicate, so none of them is spelled a second time here.
module Enrichment end
Base.include(Enrichment, joinpath(REPO, "scripts", "habitat_enrichment.jl"))
const HE = Enrichment
const OR = Enrichment.OutcomeRules

# ---------------------------------------------------------------------------
# Declared names. Nothing below is a quantity measured on any data.
# ---------------------------------------------------------------------------

const HABITAT_UNIT   = HE.PRIMARY_UNITS[1]
const CHEMISTRY_UNIT = HE.PRIMARY_UNITS[2]
endswith(CHEMISTRY_UNIT, "(chemistry)") ||
    error("the enrichment step's second primary unit is '$CHEMISTRY_UNIT', not a chemistry unit")

# The minimum weighted foreground count for a drawn row, compared on the fractional count.
const MIN_N = 20

# Two cells that are the same number written by two runs agree to this absolute tolerance.
const TOL = 1e-9

const FIGURE_BASE      = "habitat_enrichment"
const IMAGE_FORMATS    = ("pdf", "svg", "png")
const SOURCE_DATA_FILE = "habitat_enrichment_source_data.csv"
const VALUES_FILE      = "habitat_enrichment_caption_values.csv"
const PROVENANCE_FILE  = "habitat_figure_provenance.csv"
const FIGURE_OUTPUTS   = (Tuple(FIGURE_BASE * "." * f for f in IMAGE_FORMATS)...,
                          SOURCE_DATA_FILE, VALUES_FILE, PROVENANCE_FILE)

# A staging directory of this script's own, so a run of it never clears what another
# script's failed run left for inspection in the same figure directory.
const STAGING_DIR = "_staging_habitat_figure"

# Default file names of the inputs that live in the figure directory.
const ROBUSTNESS_FILE             = "habitat_robustness.csv"
const SENSITIVITY_FILE            = "habitat_sensitivity.csv"
const SENSITIVITY_PROVENANCE_FILE = "habitat_sensitivity_provenance.csv"

# The inputs, in the order they are read and recorded.
const INPUT_NAMES = (:stats, :summary, :robustness, :sensitivity, :sensitivity_provenance,
                     :calls)

# What each of those is and which option names it, so that a missing one is a sentence rather
# than a bare path. The three step-4 outputs live in the figure directory the reader chose, so
# the path in the message is theirs by construction; the other three are registry names and
# are refused by name in the program guard before this is reached.
const INPUT_WHAT = (stats = ("step 3's estimates table", "--stats"),
                    summary = ("step 3's arm summary", "--summary"),
                    robustness = ("step 4's robustness table", "--robustness"),
                    sensitivity = ("step 4's sensitivity grid", "--sensitivity"),
                    sensitivity_provenance = ("step 4's grid provenance",
                                              "--sensitivity-provenance"),
                    calls = ("step 1's habitat call table", "--calls"))

# The input gates, by name, in the order they are built. `input_gates` compares what it
# built against this list, so a gate cannot be dropped by one edit.
const INPUT_GATES = ("arm sizes", "summary identities", "robustness rows", "comparator rows",
                     "robustness values", "robustness flags", "flags against grid",
                     "sensitivity primary cells", "provenance flags", "scan digest",
                     "calls digest", "sensitivity run state", "panel files", "drawn rows")

# The name of the check made on the staged outputs just before they are published.
const OUTPUT_GATE = "staged outputs"

# THE SELECTED CATEGORIES ARE NOT CONSTANTS OF THIS FILE, AND IT SHIPS NONE OF THEM. Which
# category's sharing ratio is read back and which categories' intervals are counted across
# every specification are selections made from this data set's own estimates — chosen with the
# estimates in hand, not before. Which phylum's shares are quoted is not a selection of that
# kind: it is determined as the heaviest standardisation stratum, and the summary-identity gate
# asserts the weight identities only while the run's own heaviest stratum is still that phylum.
# It is supplied rather than written here because the share keys are named after it. They arrive
# in the criteria
# file --categories names (the same file step 3 is given, so the two steps cannot disagree
# about them), they are read by HE.read_category_selection, they travel with the inputs as
# `inp.selection`, and there is no default for any of them anywhere in this file. A run given
# no criteria file is refused by name.
#
# What that leaves here is one number that is NOT a category selection:

# How many drawn rows the thinnest-stratum list names. A listing convention, and outcome-
# independent BY CONSTRUCTION rather than by date: the list is extended through a tie at its
# last place, and the smallest count among the rows it leaves out is exported beside it, so the
# list plus that bound is a statement about every drawn row and not a selection. It was written
# during figure construction, after the results were known, which is why the defence is the
# construction and not a claim about when the number was fixed.
const N_THINNEST_LISTED = 4

# Roles of the rows in the source data.
const ROLE_PLOTTED     = "plotted estimate"
const ROLE_HELD_BACK   = "held back"
const ROLE_UNIT_CHECK  = "unit check"
const ROLE_TICK        = "family tick"
const ROLE_TICK_HELD   = "family comparator, row held back"
const ROLE_COMPARATOR  = "family comparator, unit not drawn"
const ROLE_ALL_CARRIER = "sensitivity only, circular"
const SOURCE_CONTRASTS = (HE.LN_VS_BG, HE.NON_VS_BG, HE.ALL_VS_BG)

# Okabe-Ito, as hex strings: a CVD-safe diverging pair with a neutral midpoint. Strings, not
# colour macros, so that including this file needs no renderer.
const COL_ENRICHED = "#D55E00"
const COL_DEPLETED = "#0072B2"
const COL_NEUTRAL  = "#7F7F7F"
const COL_INK      = "#2E2E2E"
const COL_RULE     = "#BFBEB6"
const CLASS_COLOURS = (enriched = COL_ENRICHED, depleted = COL_DEPLETED, neutral = COL_NEUTRAL)
const CLASS_ORDER   = (:enriched, :depleted, :neutral)

# Legend wording: the colour classes in the sensitivity panel's vocabulary, then the marks.
const LEGEND_LABELS = (enriched = "robust, enriched", depleted = "robust, depleted",
                       neutral = "not robust", interval = "95% interval",
                       ring = "crude, unadjusted", tick = "non-Ln PQQ-ADH carriers",
                       carriers = "Ln carriers", background = "background")
const LABEL_HABITAT   = "HABITAT"
const LABEL_CHEMISTRY = "CHEMISTRY (non-exclusive)"

# ---------------------------------------------------------------------------
# Geometry, in the house units: 180 mm wide, a fixed pitch per row. All offsets are in rows,
# all marker sizes in points. The clearances below are what the suite checks.
# ---------------------------------------------------------------------------

const WIDTH_MM    = 180
const PT_PER_MM   = 72 / 25.4
const ROW_MM      = 7.0
const FS          = 7
const LABEL_FS    = FS - 1.5
const DOT_DY      = 0.0      # the standardised estimate sits on the row's own line
const DOT_MS      = 4.6
const WHISKER_PT  = 3.0      # height of the interval's end caps
const RING_DY     = -0.30    # the crude estimate, one sub-row below
const RING_MS     = 3.6
const RING_STROKE = 0.6
const TICK_DY     = 0.33     # the comparator tick, one sub-row above
const TICK_HALF   = 0.12
const TICK_LW     = 0.9
const BAR_BG_WIDTH = 0.86
const BAR_FG_WIDTH = 0.46
const GAP         = 1.0      # blank rows between the chemistry and habitat blocks
const BOTTOM_PAD  = 0.7
const TOP_PAD     = 1.0
const LABEL_DY    = 0.08     # a block label's top edge sits this far below its anchor line
const LABEL_DX    = 0.01     # and this fraction of the axis in from its left edge
const B_PAD       = 0.12     # log2 units added on both sides of panel (b)'s extent
# The raster preview's resolution at the figure's own width. The vector files carry the figure;
# this is the number a raster submission is checked against, so it is stated rather than implied
# by a scale factor.
const PNG_DPI     = 300
const PT_PER_INCH = 72
const SHARE_STEPS = (1, 2, 5, 10, 20, 25, 50, 100)
const SHARE_MAX_TICK_INTERVALS = 5

row_pt() = ROW_MM * PT_PER_MM
figure_width_pt() = WIDTH_MM * PT_PER_MM
dot_half_pt() = max(DOT_MS, WHISKER_PT) / 2
ring_half_pt() = (RING_MS + RING_STROKE) / 2
tick_top() = TICK_DY + TICK_HALF
tick_bottom() = TICK_DY - TICK_HALF
# Bottom of the tick to the top of the dot (or of an interval cap) on the same row.
tick_dot_clearance_pt() = (tick_bottom() - DOT_DY) * row_pt() - dot_half_pt() - TICK_LW / 2
# Bottom of the dot to the top of the ring on the same row.
dot_ring_clearance_pt() = (DOT_DY - RING_DY) * row_pt() - dot_half_pt() - ring_half_pt()
# Top of a row's tick to the bottom of the ring of the row above it.
row_to_row_clearance_pt() = (1 + RING_DY - tick_top()) * row_pt() - ring_half_pt() - TICK_LW / 2
# Lowest ring to the bottom of the axis.
bottom_clearance_pt() = (BOTTOM_PAD + RING_DY) * row_pt() - ring_half_pt()
# A block label's lower edge to the highest mark under it, for the label under the top of the
# axis and for the label under the rule between the blocks.
habitat_label_clearance_pt() = (TOP_PAD - tick_top() - LABEL_DY) * row_pt() - LABEL_FS -
                               TICK_LW / 2
chemistry_label_clearance_pt() = ((GAP / 2 + 0.5) - tick_top() - LABEL_DY) * row_pt() -
                                 LABEL_FS - TICK_LW / 2

# ---------------------------------------------------------------------------
# Formatting: one function per kind of printed number, so a value is written the same way
# wherever it appears.
# ---------------------------------------------------------------------------

sv(x) = ismissing(x) ? "" : String(strip(string(x)))
fnum(x) = (s = sv(x); isempty(s) ? NaN : something(tryparse(Float64, s), NaN))
is_sha256(s::AbstractString) = length(s) == 64 && all(isxdigit, s)
rowkey(unit, category) = string(unit, " | ", category)

"A weighted count: an integer when it is one, otherwise one decimal."
function fmt_count(x::Real)
    isfinite(x) || error("not a finite count: $x")
    return isinteger(x) ? string(Int(x)) : @sprintf("%.1f", x)
end

"""
A weighted count printed as a LOWER bound: the same one decimal, rounded down rather than to
nearest, so a sentence of the form "no row rests on fewer than this" stays true of the row the
bound came from. Rounding to nearest would print 5.9 for a row resting on 5.875 and make the
sentence contradict the figure's own source data.
"""
function fmt_count_floor(x::Real)
    isfinite(x) || error("not a finite count: $x")
    return fmt_count(floor(x * 10) / 10)
end

"A number in a message: a count when finite, as read otherwise."
show_num(x::Real) = isfinite(x) ? fmt_count(x) : string(x)

"A whole number with thousands grouped by commas."
function fmt_int(n::Integer)
    s = string(abs(n))
    parts = String[]
    while length(s) > 3
        pushfirst!(parts, s[end-2:end])
        s = s[1:end-3]
    end
    pushfirst!(parts, s)
    return (n < 0 ? "-" : "") * join(parts, ",")
end

"A percentage, one decimal."
fmt_pct(x::Real) = isfinite(x) ? @sprintf("%.1f", x) : error("not a finite percentage: $x")

"A log2 estimate with its sign, two decimals; a value that rounds to zero prints unsigned."
function fmt_log2(x::Real)
    isfinite(x) || error("not a finite estimate: $x")
    s = @sprintf("%+.2f", x)
    return s in ("+0.00", "-0.00") ? "0.00" : s
end

fmt_interval(lo::Real, hi::Real) = fmt_log2(lo) * " to " * fmt_log2(hi)
fmt_ratio(x::Real) = isfinite(x) ? @sprintf("%.2f", x) : error("not a finite ratio: $x")

const SUPERSCRIPT = Dict(zip("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))

"A probability: two decimals from 0.01 up, otherwise one significant digit times a power of ten."
function fmt_p(x::Real)
    (isfinite(x) && 0 <= x <= 1) || error("not a probability: $x")
    x >= 0.01 && return @sprintf("%.2f", x)
    x == 0 && return "0"
    k = floor(Int, log10(x))
    m = round(Int, x / 10.0^k)
    if m >= 10
        m = 1
        k += 1
    end
    m < 1 && (m = 1)
    return string(m, " × 10", join(SUPERSCRIPT[c] for c in string(k)))
end

"Names as a list for reading: 'a', 'a and b', 'a, b and c', or 'none'."
function prose_list(xs)
    v = String[String(x) for x in xs]
    isempty(v) && return "none"
    length(v) == 1 && return v[1]
    return join(v[1:end-1], ", ") * " and " * v[end]
end

split_codes(s::AbstractString) = String[String(strip(x)) for x in split(s, ';') if !isempty(strip(x))]

# ---------------------------------------------------------------------------
# Options.
# ---------------------------------------------------------------------------

const OPTION_SPEC = Dict("--categories" => :value,
                         "--stats" => :value, "--summary" => :value, "--robustness" => :value,
                         "--sensitivity" => :value, "--sensitivity-provenance" => :value,
                         "--calls" => :value, "--panel-dir" => :value, "--figdir" => :value,
                         "--check-only" => :flag, "--help" => :flag, "-h" => :flag)

"""
    parse_run_options(args; env = ENV) -> NamedTuple

The command line as the settings of one run. Unknown, repeated or value-less options are
errors (the enrichment step's parser). The figure directory comes from `--figdir`, else from
LNFINDER_FIGDIR; with neither, or when it does not exist, the run is refused — a figure
directory is never created, so no default location can come into being by accident. Every
input records whether its path was typed or defaulted.
"""
function parse_run_options(args; env = ENV)
    cli = HE.parse_options(args, OPTION_SPEC)
    (haskey(cli, "--help") || haskey(cli, "-h")) && return (help = true,)
    # FIRST, and with no default: which category's sharing ratio this figure reads back, which
    # categories' intervals it counts across specifications and which phylum's shares it quotes
    # are selections made from this data set's estimates, not properties of this software. A
    # default would publish a choice as this file's own, and would be indistinguishable from a
    # selection somebody checked.
    categories = HE.option_string(cli, "--categories", "")
    isempty(categories) &&
        error("no --categories file was given, and this step has no default for what is in " *
              "it: the sharing category, the grid-count categories and the standardisation " *
              "phylum are SELECTIONS made from this data set's own estimates. Pass " *
              "--categories <path> — the same file step 3 was given, so that the two steps " *
              "cannot disagree about which categories the run is about")
    envdir = strip(get(env, "LNFINDER_FIGDIR", ""))
    figdir = HE.option_string(cli, "--figdir", String(envdir))
    isempty(figdir) &&
        error("no figure directory: pass --figdir <path> or set LNFINDER_FIGDIR")
    isdir(figdir) ||
        error("the figure directory '$figdir' does not exist, and this script never creates one")
    figdir = abspath(figdir)
    given(flag) = haskey(cli, flag) ? "explicit" : "default"
    input(flag, default) = (path = abspath(HE.option_string(cli, flag, default)),
                            given = given(flag))
    return (help = false,
            check_only = haskey(cli, "--check-only"),
            categories = abspath(categories),
            figdir = figdir,
            figdir_given = haskey(cli, "--figdir") ? "explicit" : "environment",
            stats = input("--stats", DataInputs.input_path("enrichment_table")),
            summary = input("--summary", DataInputs.input_path("arm_summary")),
            robustness = input("--robustness", joinpath(figdir, ROBUSTNESS_FILE)),
            sensitivity = input("--sensitivity", joinpath(figdir, SENSITIVITY_FILE)),
            sensitivity_provenance = input("--sensitivity-provenance",
                                           joinpath(figdir, SENSITIVITY_PROVENANCE_FILE)),
            calls = input("--calls", DataInputs.input_path("habitat_calls")),
            panel_dir = input("--panel-dir", figdir))
end

# ---------------------------------------------------------------------------
# Inputs. Each file is read once; its digest and its parse come from the same bytes.
# ---------------------------------------------------------------------------

function read_input(path::AbstractString; what::AbstractString = "", option::AbstractString = "")
    isfile(path) ||
        error("missing input" * (isempty(what) ? "" : " — " * what) * ": " * path *
              (isempty(option) ? "" : "\nPass " * option * " <path>, or put it where this " *
                                      "step looked. docs/DATA.md §7 says which step writes it."))
    bytes = read(path)
    return (path = String(path), bytes = bytes, sha256 = bytes2hex(sha256(bytes)))
end

parse_table(bytes::Vector{UInt8}) = CSV.read(bytes, DataFrame; types = String, pool = false)

function key_values(df::AbstractDataFrame, what::AbstractString)
    ("key" in names(df) && "value" in names(df)) || error("$what is not a key,value table")
    d = Dict{String,String}()
    for r in eachrow(df)
        k = sv(r.key)
        haskey(d, k) && error("$what carries the key '$k' twice")
        d[k] = sv(r.value)
    end
    return d
end

const STATS_TEXT    = (:contrast, :unit, :habitat, :tested, :rounded_for_fisher)
const STATS_NUMBERS = (:a, :n1, :b, :n2, :pct_fg, :pct_bg, :log2_rr, :ci_lo, :ci_hi,
                       :log2_rr_std, :ci_lo_std, :ci_hi_std, :fisher_p,
                       :n_zero_stratum_cells, :min_bg_stratum_cell)
const ROB_TEXT      = (:unit, :category, :robust_sign, :primary_excludes_zero,
                       :breaking_variants, :depletion_tier)
const ROB_NUMBERS   = (:n_breaking_variants, :primary_std, :primary_lo, :primary_hi, :primary_a,
                       :variants_expected, :variants_with_a_row,
                       :n_variants_that_moved_the_estimate)
const GRID_TEXT     = (:variant, :role, :contrast, :unit, :category, :tested)
const GRID_NUMBERS  = (:a, :b, :n1, :n2, :crude, :crude_lo, :crude_hi, :std, :lo, :hi,
                       :n_zero_stratum_cells, :min_bg_stratum_cell)

function typed(raw::AbstractDataFrame, text, numbers; what::AbstractString)
    absent = [String(c) for c in (text..., numbers...) if !(String(c) in names(raw))]
    isempty(absent) || error("$what lacks the column(s) " * join(absent, ", "))
    df = DataFrame()
    for c in text;    df[!, c] = String[sv(x) for x in raw[!, c]]; end
    for c in numbers; df[!, c] = Float64[fnum(x) for x in raw[!, c]]; end
    return df
end

"""
    load_inputs(opt) -> NamedTuple

Every input, read once: the raw bytes' digest, the parsed tables (text as read, numbers
typed), and the panel files the grid's provenance names, hashed where they exist.
"""
function load_inputs(opt)
    files = NamedTuple{INPUT_NAMES}(Tuple(
        read_input(getfield(opt, n).path;
                   what = getfield(INPUT_WHAT, n)[1], option = getfield(INPUT_WHAT, n)[2])
        for n in INPUT_NAMES))
    stats_raw = parse_table(files.stats.bytes)
    rob_raw = parse_table(files.robustness.bytes)
    grid_raw = parse_table(files.sensitivity.bytes)
    provenance = key_values(parse_table(files.sensitivity_provenance.bytes),
                            "the sensitivity provenance")
    panel_dir = opt.panel_dir.path
    panels = [(name = n, path = joinpath(panel_dir, n),
               sha256 = isfile(joinpath(panel_dir, n)) ?
                        bytes2hex(open(sha256, joinpath(panel_dir, n))) : "")
              for n in split_codes(get(provenance, "panel_written", ""))]
    # The category selection, read from the criteria file with step 3's own reader — one
    # reader, so a file step 3 accepted cannot be read differently here. Its digest comes from
    # the bytes it parsed, and reaches the provenance record.
    selection = HE.read_category_selection(opt.categories)
    return (files = files,
            selection = selection,
            stats_raw = stats_raw,
            stats = typed(stats_raw, STATS_TEXT, STATS_NUMBERS; what = "the enrichment rows"),
            summary = key_values(parse_table(files.summary.bytes), "the arm summary"),
            robustness_raw = rob_raw,
            robustness = typed(rob_raw, ROB_TEXT, ROB_NUMBERS; what = "the robustness flags"),
            grid = typed(grid_raw, GRID_TEXT, GRID_NUMBERS; what = "the sensitivity grid"),
            provenance = provenance,
            panel_dir = panel_dir,
            panels = panels)
end

# ---------------------------------------------------------------------------
# Rows and colour.
# ---------------------------------------------------------------------------

rows_of(stats, contrast, unit) = stats[(stats.contrast .== contrast) .& (stats.unit .== unit), :]
primary_rows(stats, contrast) =
    stats[(stats.contrast .== contrast) .& in.(stats.unit, Ref(HE.PRIMARY_UNITS)), :]

"""
    select_rows(stats; min_n = MIN_N) -> (habitat, chemistry, held_back)

The Ln-against-background rows on the two taxon units. A row is drawn when its weighted
foreground count reaches `min_n` (a count that is not a number does not). Each block is
sorted on the standardised estimate, ties broken by category name, so the order does not
depend on the order of the input rows.
"""
function select_rows(stats::AbstractDataFrame; min_n::Real = MIN_N)
    function block(unit)
        rows = rows_of(stats, HE.LN_VS_BG, unit)
        keep = Bool[r.a >= min_n for r in eachrow(rows)]
        return (drawn = sort(rows[keep, :], [:log2_rr_std, :habitat]),
                held = sort(rows[.!keep, :], :habitat))
    end
    h = block(HABITAT_UNIT)
    c = block(CHEMISTRY_UNIT)
    return (habitat = h.drawn, chemistry = c.drawn, held_back = vcat(h.held, c.held))
end

drawn_rows(sel) = vcat(sel.chemistry, sel.habitat)

"""
    colour_class(sign, primary_excludes_zero, lo, hi; std = NaN, what = "") -> Symbol

The colour class of one row, from its robustness flag alone: `:enriched` for `+` whose
primary interval excludes zero above it, `:depleted` for `-` whose interval excludes zero
below it, `:neutral` for `null`, `fragile` and `undecidable`. Every other combination of the
flag, the boolean beside it and the interval it rests on is an error naming the row: a flag
that disagrees with its own interval says the file is not what it claims to be, and drawing
the row grey would hide exactly that.
"""
function colour_class(sign::AbstractString, pez::AbstractString, lo::Real, hi::Real;
                      std::Real = NaN, what::AbstractString = "")
    bad(why) = error("robustness flag of $what: $why")
    p = pez == "true" ? true : pez == "false" ? false :
        bad("primary_excludes_zero is '$pez', not true or false")
    status = OR.interval_excludes_zero(lo, hi)
    if status !== :undecidable && isfinite(std) && !(lo <= std <= hi)
        bad("the estimate lies outside its own interval")
    end
    if sign == "+"
        p || bad("'+' with primary_excludes_zero false")
        (status === :excludes_zero && lo > 0) || bad("'+' on an interval that does not lie above zero")
        return :enriched
    elseif sign == "-"
        p || bad("'-' with primary_excludes_zero false")
        (status === :excludes_zero && hi < 0) || bad("'-' on an interval that does not lie below zero")
        return :depleted
    elseif sign == "fragile"
        p || bad("'fragile' with primary_excludes_zero false")
        status === :excludes_zero || bad("'fragile' on an interval that does not exclude zero")
        return :neutral
    elseif sign == "null"
        p && bad("'null' with primary_excludes_zero true")
        status === :includes_zero || bad("'null' on an interval that does not include zero")
        return :neutral
    elseif sign == "undecidable"
        p && bad("'undecidable' with primary_excludes_zero true")
        status === :undecidable || bad("'undecidable' on an interval that was computed")
        return :neutral
    end
    bad("unknown flag '$sign'")
end

colour_class(r; what = rowkey(r.unit, r.category)) =
    colour_class(r.robust_sign, r.primary_excludes_zero, r.primary_lo, r.primary_hi;
                 std = r.primary_std, what = what)

function one_row(df::AbstractDataFrame, what::AbstractString)
    nrow(df) == 1 || error("$what: $(nrow(df)) rows where exactly one is needed")
    return df[1, :]
end

robustness_row(rob, unit, category) =
    one_row(rob[(rob.unit .== unit) .& (rob.category .== category), :],
            "robustness flags for $(rowkey(unit, category))")

"""
    family_tick(stats, unit, category) -> DataFrameRow

The comparator row for one drawn row: the non-Ln carriers against the background, on the
SAME unit and category. Keyed on both, so a chemistry flag never borrows the tick of a
habitat category that happens to share its name.
"""
family_tick(stats, unit, category) =
    one_row(stats[(stats.contrast .== HE.NON_VS_BG) .& (stats.unit .== unit) .&
                  (stats.habitat .== category), :],
            "$(HE.NON_VS_BG) for $(rowkey(unit, category))")

# ---------------------------------------------------------------------------
# The input gates. Each returns the problems it found; an empty list passes.
# ---------------------------------------------------------------------------

function check_arm_sizes(inp)
    s = inp.summary
    problems = String[]
    function whole(k)
        v = tryparse(Int, get(s, k, ""))
        v === nothing && push!(problems, "the arm summary has no whole number under '$k'")
        return v
    end
    n1 = whole("n1_taxa"); n2 = whole("n2_taxa"); nn = whole("n_nonln_carrier_taxa")
    isempty(problems) || return problems
    ln = primary_rows(inp.stats, HE.LN_VS_BG)
    nrow(ln) == 0 && push!(problems, "the enrichment rows hold no $(HE.LN_VS_BG) row on a taxon unit")
    for r in eachrow(ln)
        (r.n1 == n1 && r.n2 == n2) ||
            push!(problems, "$(HE.LN_VS_BG) $(rowkey(r.unit, r.habitat)) has n1 " *
                            "$(show_num(r.n1)), n2 $(show_num(r.n2)) against the summary's $n1, $n2")
    end
    for r in eachrow(primary_rows(inp.stats, HE.NON_VS_BG))
        (r.n1 == nn && r.n2 == n2) ||
            push!(problems, "$(HE.NON_VS_BG) $(rowkey(r.unit, r.habitat)) has n1 " *
                            "$(show_num(r.n1)), n2 $(show_num(r.n2)) against the summary's $nn, $n2")
    end
    return problems
end

"""
    check_summary_identities(inp)

Statements the arm summary and the grid's provenance each make twice, checked against
themselves: a coverage percentage against its own numerator and denominator, those numerators
against the two arm sizes the figure draws, the sampled background against what the exclusion
rule removed from it and what it left, the foreground arm against its two domains, the
heaviest standardisation stratum against the phylum share printed beside it, and the per-unit
variant counts against the axis declaration they count. Every one of these is a number the
text beside the figure prints; a summary that contradicts itself would print both halves.
"""
function check_summary_identities(inp)
    s = inp.summary
    p = inp.provenance
    problems = String[]
    ARM = "the arm summary"
    PROV = "the sensitivity provenance"
    function num(d, k, what)
        v = tryparse(Float64, get(d, k, ""))
        v === nothing && push!(problems, "$what has no number under '$k'")
        return v
    end
    for arm in ("fg_coverage", "bg_analysed_coverage", "bg_frame_coverage", "nonln_coverage")
        n, d = num(s, "$(arm)_numerator", ARM), num(s, "$(arm)_denominator", ARM)
        pct = num(s, "$(arm)_pct", ARM)
        (n === nothing || d === nothing || pct === nothing) && continue
        if d == 0
            push!(problems, "$(arm)_denominator is zero")
        elseif !same_number(pct, 100 * n / d)
            push!(problems, "$(arm)_pct is $pct where its own numerator over its own " *
                            "denominator is $(100 * n / d)")
        end
    end
    for (ck, nk) in (("fg_coverage_numerator", "n1_taxa"),
                     ("bg_analysed_coverage_numerator", "n2_taxa"),
                     ("nonln_coverage_numerator", "n_nonln_carrier_taxa"))
        a, b = num(s, ck, ARM), num(s, nk, ARM)
        (a === nothing || b === nothing) && continue
        a == b || push!(problems, "$ck is $(show_num(a)) where $nk is $(show_num(b))")
    end
    # The background the figure compares against is what the sample had left after the
    # exclusion rule removed part of it. All three numbers are printed beside the figure, in
    # sentences a reader subtracts one from the other, so the subtraction is checked here.
    fr = num(s, "bg_frame_coverage_denominator", ARM)
    rm_ = num(s, "bg_removed_by_scan", ARM)
    an = num(s, "bg_analysed_coverage_denominator", ARM)
    if !any(x -> x === nothing, (fr, rm_, an)) && !same_number(fr - rm_, an)
        push!(problems, "bg_frame_coverage_denominator $(show_num(fr)) less bg_removed_by_scan " *
                        "$(show_num(rm_)) is $(show_num(fr - rm_)), not the " *
                        "bg_analysed_coverage_denominator $(show_num(an))")
    end
    nb, na, n1 = num(s, "n_bacteria_fg", ARM), num(s, "n_archaea_fg", ARM), num(s, "n1_taxa", ARM)
    if !(nb === nothing || na === nothing || n1 === nothing) && nb + na != n1
        push!(problems, "n_bacteria_fg $(show_num(nb)) and n_archaea_fg $(show_num(na)) do not " *
                        "sum to n1_taxa $(show_num(n1))")
    end
    # The two quoted phylum shares against the standardisation weights they must equal. WHICH
    # phylum is the run's own selection, and the identities hold only while that phylum is also
    # the heaviest stratum — so the phylum and its key prefix both come from the selection, and
    # the check makes no assertion when the run's heaviest stratum is some other phylum.
    phylum = inp.selection.standardisation_phylum
    kp = inp.selection.key.standardisation_phylum
    if get(s, "standardisation_heaviest_stratum", "") == phylum
        w = num(s, "standardisation_heaviest_stratum_weight", ARM)
        bt = num(s, "standardisation_heaviest_stratum_bg_taxa", ARM)
        n2 = num(s, "n2_taxa", ARM)
        pf = num(s, "pct_$(kp)_fg", ARM)
        pb = num(s, "pct_$(kp)_bg", ARM)
        if !any(x -> x === nothing, (w, bt, n2, pf, pb))
            same_number(pf, 100 * w) ||
                push!(problems, "pct_$(kp)_fg is $pf where the heaviest stratum's " *
                                "weight gives $(100 * w)")
            if n2 == 0
                push!(problems, "n2_taxa is zero")
            elseif !same_number(pb, 100 * bt / n2)
                push!(problems, "pct_$(kp)_bg is $pb where the heaviest stratum's " *
                                "background taxa give $(100 * bt / n2)")
            end
        end
    end
    declared = Dict{String,Vector{String}}()
    for part in split_codes(get(p, "axis_unit_declaration", ""))
        bits = split(part, '='; limit = 2)
        if length(bits) != 2
            push!(problems, "'$part' is not an axis=unit declaration")
            continue
        end
        declared[String(strip(bits[1]))] = String[String(strip(u)) for u in split(bits[2], ',')]
    end
    registered = split_codes(get(p, "variants_in_the_registered_set", ""))
    if isempty(declared) || isempty(registered)
        push!(problems, "$PROV records no axis declaration or no registered variant set")
        return problems
    end
    axis_of(v) = String(first(split(v, ':'; limit = 2)))
    for v in registered
        haskey(declared, axis_of(v)) ||
            push!(problems, "the axis declaration does not cover the registered variant '$v'")
    end
    for unit in HE.PRIMARY_UNITS
        k = "variants_whose_axis_reaches_$unit"
        want = count(v -> unit in get(declared, axis_of(v), String[]), registered)
        got = tryparse(Int, get(p, k, ""))
        if got === nothing
            push!(problems, "$PROV has no whole number under '$k'")
        elseif got != want
            push!(problems, "'$k' is $got where the axis declaration gives $want")
        end
    end
    return problems
end

"""
    check_robustness_rows(inp)

One robustness row per category on a primary unit, and one for every such category, drawn or
not. The reverse direction matters: a row the minimum count holds back today is a row the
figure may draw tomorrow, and a flag file that is short of it says the two runs disagree about
which categories exist, not that the figure needs less.
"""
function check_robustness_rows(inp)
    problems = String[]
    ks = [(r.unit, r.category) for r in eachrow(inp.robustness)]
    for k in unique(ks)
        n = count(==(k), ks)
        n > 1 && push!(problems, "$(rowkey(k...)) has $n robustness rows")
    end
    drawn = Set((r.unit, r.habitat) for r in eachrow(drawn_rows(select_rows(inp.stats))))
    for r in eachrow(primary_rows(inp.stats, HE.LN_VS_BG))
        k = (r.unit, r.habitat)
        k in ks && continue
        push!(problems, (k in drawn ? "drawn row " : "$(HE.LN_VS_BG) ") *
                        "$(rowkey(k...)) has no robustness row")
    end
    return problems
end

function check_comparator_rows(inp)
    problems = String[]
    for r in eachrow(drawn_rows(select_rows(inp.stats)))
        m = inp.stats[(inp.stats.contrast .== HE.NON_VS_BG) .& (inp.stats.unit .== r.unit) .&
                      (inp.stats.habitat .== r.habitat), :]
        if nrow(m) != 1
            push!(problems, "drawn row $(rowkey(r.unit, r.habitat)) has $(nrow(m)) " *
                            "$(HE.NON_VS_BG) rows where exactly one is needed")
        elseif !(all(isfinite, (m.log2_rr_std[1], m.ci_lo_std[1], m.ci_hi_std[1])) &&
                 m.tested[1] == "true")
            push!(problems, "the $(HE.NON_VS_BG) row for $(rowkey(r.unit, r.habitat)) has no " *
                            "tested, finite standardised estimate and interval")
        end
    end
    return problems
end

same_number(x::Real, y::Real) = (isnan(x) && isnan(y)) || abs(x - y) <= TOL

function check_robustness_values(inp)
    problems = String[]
    ln = primary_rows(inp.stats, HE.LN_VS_BG)
    byname = Dict{Tuple{String,String},Any}()
    for r in eachrow(ln)
        k = (r.unit, r.habitat)
        haskey(byname, k) && push!(problems, "the enrichment rows carry $(rowkey(k...)) twice")
        byname[k] = r
    end
    for r in eachrow(inp.robustness)
        k = (r.unit, r.category)
        if haskey(byname, k)
            e = byname[k]
            for (rc, ec) in ((:primary_std, :log2_rr_std), (:primary_lo, :ci_lo_std),
                             (:primary_hi, :ci_hi_std), (:primary_a, :a))
                same_number(r[rc], e[ec]) ||
                    push!(problems, "$(rowkey(k...)): robustness $rc $(r[rc]) against " *
                                    "enrichment $ec $(e[ec])")
            end
        elseif r.robust_sign != "undecidable"
            push!(problems, "$(rowkey(k...)) is flagged '$(r.robust_sign)' but has no " *
                            "enrichment row; only an undecidable flag may lack one")
        end
    end
    return problems
end

"""
    check_robustness_flags(inp)

Each flag against the rest of its own row: the interval it rests on (`colour_class`), the
count of variants it was decided on, and the breaker columns beside it. A sign means no
variant broke it and `fragile` means one did, so a row that carries both a sign and a list of
breakers is contradicting itself in writing, and that contradiction is exactly what would
otherwise colour a row the grid does not support.
"""
function check_robustness_flags(inp)
    problems = String[]
    for r in eachrow(inp.robustness)
        what = rowkey(r.unit, r.category)
        cls = try
            colour_class(r)
        catch e
            e isa ErrorException || rethrow()
            push!(problems, e.msg)
            continue
        end
        if cls !== :neutral && r.variants_with_a_row != r.variants_expected
            push!(problems, "$what is coloured on $(r.variants_with_a_row) of " *
                            "$(r.variants_expected) variants")
        end
        named = split_codes(r.breaking_variants)
        same_number(r.n_breaking_variants, length(named)) ||
            push!(problems, "$what counts $(show_num(r.n_breaking_variants)) breaking " *
                            "variant(s) and names $(length(named))")
        if r.robust_sign in ("+", "-") && !isempty(named)
            push!(problems, "$what is flagged '$(r.robust_sign)' and names the breaking " *
                            "variant(s) " * join(named, ", "))
        elseif r.robust_sign == "fragile" && isempty(named)
            push!(problems, "$what is flagged 'fragile' and names no breaking variant")
        end
    end
    return problems
end

"""
    grid_sign(lo, hi) -> String

The sign one specification's interval gives: `+` above zero, `-` below it, `0` across it, and
`undecidable` where the interval was not computed.
"""
function grid_sign(lo::Real, hi::Real)
    status = OR.interval_excludes_zero(lo, hi)
    status === :undecidable && return "undecidable"
    status === :includes_zero && return "0"
    return lo > 0 ? "+" : "-"
end

"""
    flag_from_grid(primary, variants) -> (sign, breakers)

The flag a set of specifications gives: a sign, and the sorted names of the registered
variants that do not reproduce it. A primary that lies on one side of zero and is
reproduced by every variant that ran is that sign; one a variant contradicts is `fragile`; one
that spans zero is `null`, broken by any variant that does not; one that was never computed is
`undecidable`, broken by any variant that was. `variants` holds only the variants with a row:
a variant the grid does not carry did not run, and a specification that did not run cannot
break anything.
"""
function flag_from_grid(primary::AbstractString, variants)
    if primary == "undecidable"
        return ("undecidable", sort(String[v for (v, s) in variants if s != "undecidable"]))
    elseif primary == "0"
        return ("null", sort(String[v for (v, s) in variants if s in ("+", "-")]))
    end
    breakers = sort(String[v for (v, s) in variants if s != primary])
    return (isempty(breakers) ? String(primary) : "fragile", breakers)
end

"""
    check_flags_against_grid(inp)

The robustness file is a summary of the sensitivity grid, and here it is recomputed from it:
the primary row and the registered variant rows of the same contrast, unit and category give a
sign and a breaker list, and both must be what the flag row says. This is the one identity
that survives an edit made consistently in two files, because the third file is the one that
holds the intervals. The declared variant counts are checked against the grid at the same
time, so `variants_with_a_row` cannot claim a specification the grid never ran.
"""
function check_flags_against_grid(inp)
    p = inp.provenance
    registered = split_codes(get(p, "variants_in_the_registered_set", ""))
    isempty(registered) &&
        return ["the sensitivity provenance names no registered variant set, so no flag can " *
                "be checked against the grid it summarises"]
    g = inp.grid[inp.grid.contrast .== HE.LN_VS_BG, :]
    index = Dict{Tuple{String,String,String},Vector{Int}}()
    for (i, r) in enumerate(eachrow(g))
        keep = (r.role == "primary" && r.variant == "primary") ||
               (r.role == "variant" && r.variant in registered)
        keep && push!(get!(index, (r.variant, r.unit, r.category), Int[]), i)
    end
    problems = String[]
    for r in eachrow(inp.robustness)
        what = rowkey(r.unit, r.category)
        hits(v) = get(index, (v, r.unit, r.category), Int[])
        for v in ("primary", registered...)
            length(hits(v)) > 1 &&
                push!(problems, "$what has $(length(hits(v))) '$v' rows in the sensitivity grid")
        end
        h = hits("primary")
        psign = length(h) == 1 ? grid_sign(g[h[1], :lo], g[h[1], :hi]) : "undecidable"
        if isempty(h) && r.robust_sign != "undecidable"
            push!(problems, "$what is flagged '$(r.robust_sign)' with no registered-specification " *
                            "row in the sensitivity grid")
        end
        seen = Tuple{String,String}[]
        for v in registered
            hv = hits(v)
            length(hv) == 1 && push!(seen, (v, grid_sign(g[hv[1], :lo], g[hv[1], :hi])))
        end
        want, breakers = flag_from_grid(psign, seen)
        want == r.robust_sign ||
            push!(problems, "$what is flagged '$(r.robust_sign)' where its grid rows give '$want'")
        have = sort(split_codes(r.breaking_variants))
        breakers == have ||
            push!(problems, "$what names the breaking variant(s) " * prose_list(have) *
                            " where its grid rows give " * prose_list(breakers))
        same_number(r.variants_expected, length(registered)) ||
            push!(problems, "$what expects $(show_num(r.variants_expected)) variants where the " *
                            "registered set holds $(length(registered))")
        same_number(r.variants_with_a_row, length(seen)) ||
            push!(problems, "$what records $(show_num(r.variants_with_a_row)) variant row(s) " *
                            "where the grid holds $(length(seen))")
    end
    return problems
end

# Every cell the two runs both hold, including the two thin-cell columns the text beside the
# figure discloses: a stratum count that differs between the grid and the enrichment rows is a
# different standardisation, whatever the estimates say.
const GRID_STATS_COLUMNS = ((:a, :a), (:b, :b), (:n1, :n1), (:n2, :n2), (:crude, :log2_rr),
                            (:crude_lo, :ci_lo), (:crude_hi, :ci_hi), (:std, :log2_rr_std),
                            (:lo, :ci_lo_std), (:hi, :ci_hi_std),
                            (:n_zero_stratum_cells, :n_zero_stratum_cells),
                            (:min_bg_stratum_cell, :min_bg_stratum_cell))

function check_sensitivity_primary_cells(inp)
    problems = String[]
    g = inp.grid[(inp.grid.role .== "primary") .& (inp.grid.variant .== "primary"), :]
    index = Dict{Tuple{String,String,String},Vector{Int}}()
    for (i, r) in enumerate(eachrow(g))
        push!(get!(index, (r.contrast, r.unit, r.category), Int[]), i)
    end
    for e in eachrow(inp.stats)
        k = (e.contrast, e.unit, e.habitat)
        hits = get(index, k, Int[])
        if length(hits) != 1
            push!(problems, "$(e.contrast) $(rowkey(e.unit, e.habitat)) has $(length(hits)) " *
                            "registered-specification rows in the grid")
            continue
        end
        r = g[hits[1], :]
        for (gc, ec) in GRID_STATS_COLUMNS
            same_number(r[gc], e[ec]) ||
                push!(problems, "$(e.contrast) $(rowkey(e.unit, e.habitat)): grid $gc " *
                                "$(r[gc]) against enrichment $ec $(e[ec])")
        end
    end
    return problems
end

function check_provenance_flags(inp)
    problems = String[]
    p = inp.provenance
    rob = inp.robustness
    for r in eachrow(rob)
        fk = "flag_$(r.unit)_$(r.category)"
        bk = "breakers_$(r.unit)_$(r.category)"
        if !haskey(p, fk)
            push!(problems, "the sensitivity provenance has no '$fk'")
        elseif p[fk] != r.robust_sign
            push!(problems, "'$fk' is '$(p[fk])' where the robustness flags say '$(r.robust_sign)'")
        end
        get(p, bk, "") == r.breaking_variants ||
            push!(problems, "'$bk' is '$(get(p, bk, ""))' where the robustness flags say " *
                            "'$(r.breaking_variants)'")
    end
    # and the reverse: a flag key that names no robustness row is a flag for a row that does
    # not exist, which is how a stale category survives a rerun.
    known = Set(String["$(pre)$(r.unit)_$(r.category)" for r in eachrow(rob)
                       for pre in ("flag_", "breakers_")])
    for k in sort(collect(keys(p)))
        (startswith(k, "flag_") || startswith(k, "breakers_")) || continue
        k in known || push!(problems, "'$k' names no robustness row")
    end
    for (k, sign) in (("robust_enriched", "+"), ("robust_depleted", "-"),
                      ("fragile_categories", "fragile"))
        want = sort(unique(String[r.category for r in eachrow(rob) if r.robust_sign == sign]))
        if !haskey(p, k)
            push!(problems, "the sensitivity provenance has no '$k'")
        elseif sort(split_codes(p[k])) != want
            push!(problems, "'$k' lists '$(p[k])' where the robustness flags give '" *
                            join(want, ";") * "'")
        end
    end
    return problems
end

function check_scan_digest(inp)
    a = get(inp.summary, "background_scan_sha256", "")
    b = get(inp.provenance, "scan_sha256", "")
    problems = String[]
    is_sha256(a) || push!(problems, "the arm summary's background_scan_sha256 is not a digest")
    is_sha256(b) || push!(problems, "the sensitivity provenance's scan_sha256 is not a digest")
    isempty(problems) && a != b &&
        push!(problems, "the arm summary was built on scan $a and the sensitivity grid on $b")
    return problems
end

function check_calls_digest(inp)
    d = get(inp.provenance, "calls_sha256", "")
    is_sha256(d) || return ["the sensitivity provenance's calls_sha256 is not a digest"]
    d == inp.files.calls.sha256 && return String[]
    return ["the sensitivity grid was built on call table $d and " *
            "$(inp.files.calls.path) is $(inp.files.calls.sha256)"]
end

function check_sensitivity_run_state(inp)
    p = inp.provenance
    problems = String[]
    get(p, "anchor_checks_declined", "") == "false" ||
        push!(problems, "anchor_checks_declined is '$(get(p, "anchor_checks_declined", "absent"))'")
    tryparse(Int, get(p, "anchor_targets_failed", "")) == 0 ||
        push!(problems, "anchor_targets_failed is '$(get(p, "anchor_targets_failed", "absent"))'")
    isempty(get(p, "publish_gates", "")) &&
        push!(problems, "the sensitivity run records no publish gates")
    if !haskey(p, "publish_gates_failed")
        push!(problems, "the sensitivity run does not record publish_gates_failed")
    elseif !isempty(p["publish_gates_failed"])
        push!(problems, "the sensitivity run failed its publish gate(s) $(p["publish_gates_failed"])")
    end
    return problems
end

function check_panel_files(inp)
    isempty(inp.panels) && return ["the sensitivity provenance names no panel file"]
    return String["panel file $(x.path) is missing or empty"
                  for x in inp.panels if !(isfile(x.path) && filesize(x.path) > 0)]
end

"""
    check_drawn_rows(inp)

Every cell a drawn row puts on the page: the marks are finite, the row was tested, and the two
shares panel (a) draws are the counts panel (b) tests divided by the arm sizes it was gated
on. The shares travel in their own columns, so without this they are the one pair of drawn
numbers no other gate reaches.
"""
function check_drawn_rows(inp)
    problems = String[]
    for r in eachrow(drawn_rows(select_rows(inp.stats)))
        what = rowkey(r.unit, r.habitat)
        bad = [String(c) for c in (:log2_rr, :log2_rr_std, :ci_lo_std, :ci_hi_std, :pct_fg, :pct_bg)
               if !isfinite(r[c])]
        isempty(bad) || push!(problems, "drawn row $what has no finite " * join(bad, ", "))
        r.tested == "true" || push!(problems, "drawn row $what is not tested")
        for (name, pct, num, den) in (("pct_fg", r.pct_fg, r.a, r.n1),
                                      ("pct_bg", r.pct_bg, r.b, r.n2))
            all(isfinite, (pct, num, den)) || continue
            if den == 0
                push!(problems, "drawn row $what has a zero arm size under $name")
            elseif !same_number(pct, 100 * num / den)
                push!(problems, "drawn row $what has $name $pct where its own count over its " *
                                "own arm size is $(100 * num / den)")
            end
        end
    end
    return problems
end

function gate(name::AbstractString, check, inp)
    problems = try
        check(inp)
    catch e
        ["could not be evaluated: " * sprint(showerror, e)]
    end
    shown = length(problems) > 8 ? vcat(problems[1:8], ["and $(length(problems) - 8) more"]) :
            problems
    return (name = String(name), ok = isempty(problems), message = join(shown, "; "))
end

"""
    input_gates(inp) -> Vector{NamedTuple}

Every statement that must hold before anything is staged, as `(name, ok, message)`. All are
evaluated and none short-circuits, so a run that fails several says so once. The built names
are compared with `INPUT_GATES`.
"""
function input_gates(inp)
    gates = [gate("arm sizes", check_arm_sizes, inp),
             gate("summary identities", check_summary_identities, inp),
             gate("robustness rows", check_robustness_rows, inp),
             gate("comparator rows", check_comparator_rows, inp),
             gate("robustness values", check_robustness_values, inp),
             gate("robustness flags", check_robustness_flags, inp),
             gate("flags against grid", check_flags_against_grid, inp),
             gate("sensitivity primary cells", check_sensitivity_primary_cells, inp),
             gate("provenance flags", check_provenance_flags, inp),
             gate("scan digest", check_scan_digest, inp),
             gate("calls digest", check_calls_digest, inp),
             gate("sensitivity run state", check_sensitivity_run_state, inp),
             gate("panel files", check_panel_files, inp),
             gate("drawn rows", check_drawn_rows, inp)]
    built = String[g.name for g in gates]
    built == String[n for n in INPUT_GATES] ||
        error("the input gates were built as " * join(built, ", ") * " where this file " *
              "declares " * join(INPUT_GATES, ", "))
    return gates
end

# ---------------------------------------------------------------------------
# The plan: every number the renderer needs, computed and checkable without it.
# ---------------------------------------------------------------------------

"Panel (a)'s upper limit on a tick of a round step, and that step."
function share_axis(maxval::Real)
    (isfinite(maxval) && maxval >= 0) || error("no finite share to scale panel (a) by: $maxval")
    for step in SHARE_STEPS
        top = step * max(1, ceil(Int, maxval / step))
        top / step <= SHARE_MAX_TICK_INTERVALS && return (top = top, step = step)
    end
    step = last(SHARE_STEPS)
    return (top = step * max(1, ceil(Int, maxval / step)), step = step)
end

"Panel (b)'s limits: every interval bound, crude value and tick, and zero, padded."
function enrichment_axis(values)
    v = Float64[x for x in values]
    all(isfinite, v) || error("a mark of panel (b) is not a finite number")
    push!(v, 0.0)
    return (lo = minimum(v) - B_PAD, hi = maximum(v) + B_PAD)
end

"""
    figure_plan(inp, sel) -> NamedTuple

The rows in drawing order with their positions, labels, colours and marks; the axis limits;
the block labels; the legend entries; the page geometry. Pure: the renderer maps it to marks
and computes nothing.
"""
function figure_plan(inp, sel)
    chem, hab = sel.chemistry, sel.habitat
    nchem, nhab = nrow(chem), nrow(hab)
    nchem + nhab == 0 && error("no category reaches the minimum count, so there is no row to draw")
    first_hab = nchem == 0 ? 1.0 : nchem + GAP + 1
    yrule = (nchem == 0 || nhab == 0) ? NaN : nchem + GAP / 2 + 0.5
    ylo = 1 - BOTTOM_PAD
    ymax = (nhab == 0 ? Float64(nchem) : first_hab + nhab - 1) + TOP_PAD
    rows = NamedTuple[]
    for (block, df, y0) in ((:chemistry, chem, 1.0), (:habitat, hab, first_hab))
        for (i, r) in enumerate(eachrow(df))
            rb = robustness_row(inp.robustness, r.unit, r.habitat)
            cls = colour_class(rb)
            t = family_tick(inp.stats, r.unit, r.habitat)
            push!(rows, (unit = r.unit, category = r.habitat, block = block, y = y0 + i - 1,
                         label = string(r.habitat, " (", fmt_count(r.a), ")"),
                         a = r.a, pct_fg = r.pct_fg, pct_bg = r.pct_bg,
                         crude = r.log2_rr, crude_lo = r.ci_lo, crude_hi = r.ci_hi,
                         std = r.log2_rr_std, lo = r.ci_lo_std, hi = r.ci_hi_std,
                         fisher_p = r.fisher_p, rounded_for_fisher = r.rounded_for_fisher,
                         n_zero_stratum_cells = r.n_zero_stratum_cells,
                         min_bg_stratum_cell = r.min_bg_stratum_cell,
                         tick = t.log2_rr_std, tick_lo = t.ci_lo_std, tick_hi = t.ci_hi_std,
                         tick_n_zero_stratum_cells = t.n_zero_stratum_cells,
                         tick_min_bg_stratum_cell = t.min_bg_stratum_cell,
                         colour_class = cls, colour = CLASS_COLOURS[cls],
                         robust_sign = rb.robust_sign, breakers = rb.breaking_variants))
        end
    end
    share = share_axis(maximum(max(r.pct_fg, r.pct_bg) for r in rows))
    xb = enrichment_axis(Iterators.flatten((r.lo, r.hi, r.crude, r.tick) for r in rows))
    label_x = xb.lo + LABEL_DX * (xb.hi - xb.lo)
    labels = NamedTuple[]
    nhab > 0 && push!(labels, (text = LABEL_HABITAT, x = label_x, y = ymax - LABEL_DY))
    nchem > 0 && push!(labels, (text = LABEL_CHEMISTRY, x = label_x,
                                y = (nhab == 0 ? ymax : yrule) - LABEL_DY))
    used = Set(r.colour_class for r in rows)
    return (rows = rows, yrule = yrule, ylo = ylo, ymax = ymax,
            share_max = share.top, share_step = share.step, xb_lo = xb.lo, xb_hi = xb.hi,
            block_labels = labels, legend = legend_entries(used),
            width_pt = figure_width_pt(), axis_height_pt = (ymax - ylo) * row_pt())
end

"""
    legend_entries(classes_used) -> (shares, marks)

The legend as `(key, label)` pairs: panel (a)'s two bars, then a swatch for each colour class
that some row uses, in a fixed order, then the interval, the ring and the tick. A class no row
uses gets no swatch, so the legend never names a colour the figure does not show.
"""
function legend_entries(classes_used)
    shares = [(:carriers, LEGEND_LABELS.carriers), (:background, LEGEND_LABELS.background)]
    marks = [(c, LEGEND_LABELS[c]) for c in CLASS_ORDER if c in classes_used]
    append!(marks, [(:interval, LEGEND_LABELS.interval), (:ring, LEGEND_LABELS.ring),
                    (:tick, LEGEND_LABELS.tick)])
    return (shares = shares, marks = marks)
end

"The drawn rows from top to bottom, as a reader meets them."
panel_order(plan) = sort(plan.rows, by = r -> -r.y)

# ---------------------------------------------------------------------------
# Source data.
# ---------------------------------------------------------------------------

"""
    source_data(inp, plan) -> DataFrame

Every enrichment row of the Ln, non-Ln and all-carrier contrasts on every unit, with its
cells exactly as read, and beside them: the row's role in the figure, its block, whether it
is drawn and where, its colour class, the robustness columns of the flag it was coloured by,
and for a comparator row the drawn row it ticks.
"""
function source_data(inp, plan)
    raw = inp.stats_raw
    st = inp.stats
    drawn = Dict((r.unit, r.category) => (i, r) for (i, r) in enumerate(panel_order(plan)))
    robraw = inp.robustness_raw
    robidx = Dict((sv(r.unit), sv(r.category)) => r for r in eachrow(robraw))
    rob_cols = (:robust_sign, :primary_excludes_zero, :n_breaking_variants,
                :breaking_variants, :depletion_tier, :variants_expected,
                :variants_with_a_row, :n_variants_that_moved_the_estimate)
    keep = Int[]
    extra = Dict{Symbol,Vector{String}}(c => String[] for c in
        (:role, :block, :drawn, :panel_row, :colour_class, :family_tick_of, rob_cols...))
    for i in 1:nrow(st)
        s = st[i, :]
        s.contrast in SOURCE_CONTRASTS || continue
        push!(keep, i)
        k = (s.unit, s.habitat)
        primary_unit = s.unit in HE.PRIMARY_UNITS
        is_drawn = s.contrast == HE.LN_VS_BG && haskey(drawn, k)
        role = s.contrast == HE.ALL_VS_BG ? ROLE_ALL_CARRIER :
               s.contrast == HE.LN_VS_BG ? (!primary_unit ? ROLE_UNIT_CHECK :
                                            is_drawn ? ROLE_PLOTTED : ROLE_HELD_BACK) :
               (!primary_unit ? ROLE_COMPARATOR : haskey(drawn, k) ? ROLE_TICK : ROLE_TICK_HELD)
        push!(extra[:role], role)
        push!(extra[:block], !primary_unit ? "" : s.unit == CHEMISTRY_UNIT ? "chemistry" : "habitat")
        push!(extra[:drawn], string(is_drawn))
        push!(extra[:panel_row], is_drawn ? string(drawn[k][1]) : "")
        push!(extra[:colour_class], is_drawn ? String(drawn[k][2].colour_class) : "")
        push!(extra[:family_tick_of], role == ROLE_TICK ? rowkey(k...) : "")
        rb = s.contrast == HE.LN_VS_BG && primary_unit ? get(robidx, k, nothing) : nothing
        for c in rob_cols
            push!(extra[c], rb === nothing ? "" : sv(rb[c]))
        end
    end
    out = raw[keep, :]
    for c in (:role, :block, :drawn, :panel_row, :colour_class, :family_tick_of, rob_cols...)
        out[!, c] = extra[c]
    end
    return out
end

# ---------------------------------------------------------------------------
# Values for the text beside the figure: key,value rows, each pre-formatted as printed.
# ---------------------------------------------------------------------------

"""
    global_keys(selected) -> Tuple{Vararg{String}}

The keys that do not depend on which categories are DRAWN, in export order.

Four of them are nevertheless named after the run's selection rather than written as literals:
the depletion-tier code carries the tier category's key prefix, and the three phylum shares
carry the phylum's. A run given other selections therefore exports keys for THOSE selections
instead of filing other numbers under names this file chose, and the order is still declared
in exactly one place.
"""
function global_keys(selected)
    t = selected.key.depletion_tier
    ph = selected.key.standardisation_phylum
    return (
    "n1", "n2", "min_n",
    "n_drawn_habitat", "drawn_habitat", "n_drawn_chemistry", "drawn_chemistry",
    "n_held_back", "held_back", "held_back_with_a",
    "n_robust", "any_robust", "n_robust_enriched", "robust_enriched",
    "n_robust_depleted", "robust_depleted", "robust_not_drawn",
    "coloured_enriched", "coloured_depleted",
    "n_fragile", "fragile", "fragile_with_breakers", "legend_only",
    "drawn_rows_with_a_zero_stratum_cell",
    "thinnest_bg_stratum_rows", "thinnest_bg_stratum_outside_the_list",
    "narrowest_ci_category", "narrowest_ci_width", "runner_up_ci_category", "runner_up_ci_width",
    "largest_fg_share_category", "largest_fg_share", "largest_bg_share_category",
    "largest_bg_share",
    "largest_adjustment_shift_category", "largest_adjustment_shift",
    "n_registered_variants", "n_specifications", "registered_variants",
    "variants_reaching_taxa", "variants_reaching_chemistry",
    "$(t)_tier_code",
    "n_nonln_carrier_taxa", "n_sequences_fg", "n_bacteria_fg", "n_archaea_fg",
    "eukaryotic_ln_taxa_excluded", "eukaryotic_ln_sequences_excluded",
    "pct_$(ph)_fg", "pct_$(ph)_bg", "pct_$(ph)_bg_before_scan",
    "bg_removed_by_scan",
    "nonln_coverage_numerator", "nonln_coverage_denominator", "nonln_coverage_pct",
    "fg_coverage_numerator", "fg_coverage_denominator", "fg_coverage_pct",
    "bg_analysed_coverage_numerator", "bg_analysed_coverage_denominator",
    "bg_analysed_coverage_pct", "coverage_fisher_p_after_scan",
    "bg_frame_coverage_numerator", "bg_frame_coverage_denominator", "bg_frame_coverage_pct",
    "coverage_fisher_p",
    "fg_relations_text_mined_pct", "bg_relations_text_mined_pct",
    "fg_taxa_tied_on_distinct", "fg_taxa_tie_broken_by_mentions",
    "fg_taxa_still_tied_after_mentions",
    "bg_taxa_tied_on_distinct", "bg_taxa_tie_broken_by_mentions",
    "bg_taxa_still_tied_after_mentions",
    "fg_isolation_terms_total", "fg_taxa_with_isolation_term",
    "fg_generic_terms_total", "fg_taxa_with_generic_term",
    "bg_isolation_terms_total", "bg_taxa_with_isolation_term",
    "bg_generic_terms_total", "bg_taxa_with_generic_term",
    "phylum_strata_used", "standardisation_heaviest_stratum",
    "standardisation_heaviest_stratum_weight_pct", "standardisation_heaviest_stratum_bg_taxa",
    "n_genera_fg", "n_genera_bg",
    "background_licence_code", "omnicrobe_release", "obt_release", "uniprot_release",
    "panel_files", "panel_categories", "panel_categories_not_drawn")
end

# Keys exported for every drawn row, as "<category>.<field>".
const ROW_FIELDS = (
    "block", "rank_crude", "rank_std", "a", "pct_fg", "pct_bg",
    "crude", "crude_lo", "crude_hi", "crude_ci", "fisher_p", "rounded_for_fisher",
    "val", "lo", "hi", "ci", "width", "adjustment_direction",
    "fam", "fam_lo", "fam_hi", "fam_ci",
    "colour", "robust_sign", "breakers", "breaker_estimates",
    "n_zero_stratum_cells", "min_bg_stratum_cell",
    "fam_n_zero_stratum_cells", "fam_min_bg_stratum_cell")

# Keys exported for each grid-count category the selection names, as "<category>.<field>".
const GRID_FIELDS = (
    "taxa_specs_excluding_zero", "taxa_specs_distinct",
    "sequences_specs_excluding_zero", "sequences_specs_excluding_zero_below",
    "sequences_specs_not_excluding_zero", "sequences_specs_distinct",
    "seq_val", "seq_lo", "seq_hi", "seq_ci")

# Keys exported for the selection's sharing category, as "<category>.<field>".
const SHARING_FIELDS = ("sharing_ratio", "sharing_ratio_lo", "sharing_ratio_hi",
                        "sharing_ratio_ci", "sharing_code", "sharing_rule_used")

const KEY_PATTERN = r"^[A-Za-z0-9_.\-]+$"

"""
    caption_keys(categories, selected) -> Vector{String}

The declared key list for a run that draws `categories` (in panel order) and was handed the
category selection `selected`. Both arguments are required: which categories are drawn is
measured from the data, and which categories fill the selection's roles is supplied per run —
neither is something this file may supply for itself.
"""
function caption_keys(categories, selected)
    keys_ = String[global_keys(selected)...]
    for c in categories, f in ROW_FIELDS
        push!(keys_, "$c.$f")
    end
    for c in selected.grid_count, f in GRID_FIELDS
        push!(keys_, "$c.$f")
    end
    for f in SHARING_FIELDS
        push!(keys_, "$(selected.sharing).$f")
    end
    bad = [k for k in keys_ if !occursin(KEY_PATTERN, k)]
    isempty(bad) || error("key(s) outside the key alphabet: " * join(bad, ", "))
    dup = [k for k in unique(keys_) if count(==(k), keys_) > 1]
    isempty(dup) || error("key(s) declared twice: " * join(dup, ", "))
    return keys_
end

function require_value(d::AbstractDict, k::AbstractString, what::AbstractString)
    haskey(d, k) || error("$what has no key '$k'")
    return d[k]
end

function require_int(d::AbstractDict, k::AbstractString, what::AbstractString)
    v = tryparse(Int, require_value(d, k, what))
    v === nothing && error("$what: '$k' is '$(d[k])', not a whole number")
    return v
end

function require_float(d::AbstractDict, k::AbstractString, what::AbstractString)
    v = tryparse(Float64, require_value(d, k, what))
    (v === nothing || !isfinite(v)) && error("$what: '$k' is '$(d[k])', not a finite number")
    return v
end

"""
    require_drawn(plan, category) -> NamedTuple

The drawn row of `category`, or an error naming it: a value about a row the figure does not
draw cannot be exported beside the figure.
"""
function require_drawn(plan, category::AbstractString)
    m = [r for r in plan.rows if r.category == category]
    isempty(m) && error("'$category' is not a drawn row, so its values cannot be exported")
    length(m) > 1 && error("'$category' is drawn in both blocks, so its keys would collide")
    return m[1]
end

"""
    across_specifications(grid, registered, unit, category) -> Vector{NamedTuple}

The standardised interval of `category` on `unit` under the registered specification and
under every registered one-at-a-time variant, one row each, or an error naming the missing
specification.
"""
function across_specifications(grid, registered, unit::AbstractString, category::AbstractString)
    out = NamedTuple[]
    for v in ("primary", registered...)
        role = v == "primary" ? "primary" : "variant"
        m = grid[(grid.variant .== v) .& (grid.role .== role) .& (grid.contrast .== HE.LN_VS_BG) .&
                 (grid.unit .== unit) .& (grid.category .== category), :]
        r = one_row(m, "the sensitivity grid's $v row for $(rowkey(unit, category))")
        push!(out, (variant = v, std = r.std, lo = r.lo, hi = r.hi))
    end
    return out
end

sequences_unit(unit) = unit == HABITAT_UNIT ? "sequences" :
                       unit == CHEMISTRY_UNIT ? "sequences (chemistry)" :
                       error("no per-sequence unit is defined for '$unit'")

function breaker_estimates(inp, row)
    names_ = split_codes(row.breakers)
    isempty(names_) && return "none"
    parts = String[]
    for b in names_
        m = inp.grid[(inp.grid.variant .== b) .& (inp.grid.role .== "variant") .&
                     (inp.grid.contrast .== HE.LN_VS_BG) .& (inp.grid.unit .== row.unit) .&
                     (inp.grid.category .== row.category), :]
        if nrow(m) == 1 && all(isfinite, (m.std[1], m.lo[1], m.hi[1]))
            push!(parts, string(b, " ", fmt_log2(m.std[1]), " (", fmt_interval(m.lo[1], m.hi[1]), ")"))
        else
            push!(parts, string(b, " no interval"))
        end
    end
    return join(parts, "; ")
end

"""
    caption_values(inp, sel, plan) -> Vector{Pair{String,String}}

Every value, computed from the inputs and formatted as printed, under exactly the keys
`caption_keys` declares for this run; a key built that is not declared, or declared and not
built, is an error. A required row or key that does not exist is an error naming it.

The category selection came in with the inputs (`inp.selection`) and is read from there, so
the keys built here and the keys declared there cannot be built from two different selections.
"""
function caption_values(inp, sel, plan)
    selected = inp.selection
    grid_categories = selected.grid_count
    sharing_category = selected.sharing
    s = inp.summary
    p = inp.provenance
    rob = inp.robustness
    ARM = "the arm summary"
    PROV = "the sensitivity provenance"
    out = Pair{String,String}[]
    put(k, v) = push!(out, String(k) => String(v))
    order = panel_order(plan)

    put("n1", fmt_int(require_int(s, "n1_taxa", ARM)))
    put("n2", fmt_int(require_int(s, "n2_taxa", ARM)))
    put("min_n", fmt_count(MIN_N))
    hab = [r.category for r in order if r.block === :habitat]
    chem = [r.category for r in order if r.block === :chemistry]
    put("n_drawn_habitat", string(length(hab)))
    put("drawn_habitat", prose_list(hab))
    put("n_drawn_chemistry", string(length(chem)))
    put("drawn_chemistry", prose_list(chem))
    held = sort(collect(eachrow(sel.held_back)), by = r -> (r.habitat, r.unit))
    put("n_held_back", string(length(held)))
    put("held_back", prose_list([r.habitat for r in held]))
    put("held_back_with_a", prose_list(["$(r.habitat) ($(fmt_count(r.a)))" for r in held]))

    # Robust rows: the drawn ones in panel order, then any others by unit and name.
    rank = Dict((r.unit, r.category) => i for (i, r) in enumerate(order))
    by_panel(rs) = sort(rs, by = r -> (get(rank, (r.unit, r.category), typemax(Int)), r.unit, r.category))
    robust(sign) = by_panel([r for r in eachrow(rob) if r.robust_sign == sign])
    enr, dep = robust("+"), robust("-")
    for (k, rs) in (("robust_enriched", enr), ("robust_depleted", dep))
        sort(split_codes(require_value(p, k, PROV))) == sort(unique([r.category for r in rs])) ||
            error("the robust rows of the robustness flags do not match $PROV's '$k'")
    end
    put("n_robust", string(length(enr) + length(dep)))
    put("any_robust", string(!isempty(enr) || !isempty(dep)))
    put("n_robust_enriched", string(length(enr)))
    put("robust_enriched", prose_list([r.category for r in enr]))
    put("n_robust_depleted", string(length(dep)))
    put("robust_depleted", prose_list([r.category for r in dep]))
    put("robust_not_drawn", prose_list([r.category for r in vcat(enr, dep)
                                        if !haskey(rank, (r.unit, r.category))]))
    put("coloured_enriched", prose_list([r.category for r in order if r.colour_class === :enriched]))
    put("coloured_depleted", prose_list([r.category for r in order if r.colour_class === :depleted]))
    frag = by_panel([r for r in eachrow(rob) if r.robust_sign == "fragile"])
    put("n_fragile", string(length(frag)))
    put("fragile", prose_list([r.category for r in frag]))
    put("fragile_with_breakers",
        prose_list(["$(r.category) ($(isempty(r.breaking_variants) ? "none" : r.breaking_variants))"
                    for r in frag]))
    put("legend_only", prose_list(split_codes(require_value(s, "categories_legend_only", ARM))))
    put("drawn_rows_with_a_zero_stratum_cell",
        prose_list([r.category for r in order if r.n_zero_stratum_cells > 0]))
    # The rows a standardised estimate rests on the fewest background taxa for: computed and
    # named here, so a sentence about them cannot name a row the data no longer supports. The
    # list is the ascending head of ALL drawn rows, carried past a tie at its last place, and
    # the bound is the smallest count it leaves out, rounded DOWN so the bound never claims
    # more than the row it came from - together they say something about every drawn row, so no
    # row can be thinner than the list without appearing in it.
    thin = sort([(r.min_bg_stratum_cell, r.unit, r.category) for r in plan.rows])
    cut = thin[min(N_THINNEST_LISTED, length(thin))][1]
    put("thinnest_bg_stratum_rows",
        prose_list(["$(t[3]) ($(fmt_count(t[1])))" for t in thin if t[1] <= cut]))
    rest = [t[1] for t in thin if t[1] > cut]
    put("thinnest_bg_stratum_outside_the_list",
        isempty(rest) ? "none" : fmt_count_floor(minimum(rest)))

    widths = sort([(r.hi - r.lo, r.unit, r.category) for r in plan.rows])
    put("narrowest_ci_category", widths[1][3])
    put("narrowest_ci_width", fmt_ratio(widths[1][1]))
    put("runner_up_ci_category", length(widths) > 1 ? widths[2][3] : "none")
    put("runner_up_ci_width", length(widths) > 1 ? fmt_ratio(widths[2][1]) : "none")
    fg = sort([(-r.pct_fg, r.unit, r.category, r.pct_fg) for r in plan.rows])[1]
    bg = sort([(-r.pct_bg, r.unit, r.category, r.pct_bg) for r in plan.rows])[1]
    put("largest_fg_share_category", fg[3])
    put("largest_fg_share", fmt_pct(fg[4]))
    put("largest_bg_share_category", bg[3])
    put("largest_bg_share", fmt_pct(bg[4]))
    # How far adjustment moved a row, in either direction: the drawn row whose ring is
    # furthest from its dot, computed rather than asserted.
    shift = sort([(-abs(r.std - r.crude), r.unit, r.category, abs(r.std - r.crude))
                  for r in plan.rows])[1]
    put("largest_adjustment_shift_category", shift[3])
    put("largest_adjustment_shift", fmt_ratio(shift[4]))

    registered = split_codes(require_value(p, "variants_in_the_registered_set", PROV))
    require_int(p, "n_variants_in_the_registered_set", PROV) == length(registered) ||
        error("$PROV's registered-set count disagrees with the names it lists")
    put("n_registered_variants", string(length(registered)))
    put("n_specifications", string(1 + length(registered)))
    put("registered_variants", join(registered, ";"))
    put("variants_reaching_taxa",
        string(require_int(p, "variants_whose_axis_reaches_$(HABITAT_UNIT)", PROV)))
    put("variants_reaching_chemistry",
        string(require_int(p, "variants_whose_axis_reaches_$(CHEMISTRY_UNIT)", PROV)))
    put("$(selected.key.depletion_tier)_tier_code",
        require_value(s, "$(selected.key.depletion_tier)_tier_code", ARM))

    for k in ("n_nonln_carrier_taxa", "n_sequences_fg", "n_bacteria_fg", "n_archaea_fg",
              "eukaryotic_ln_taxa_excluded", "eukaryotic_ln_sequences_excluded")
        put(k, fmt_int(require_int(s, k, ARM)))
    end
    for k in ("pct_$(selected.key.standardisation_phylum)_fg",
              "pct_$(selected.key.standardisation_phylum)_bg",
              "pct_$(selected.key.standardisation_phylum)_bg_before_scan")
        put(k, fmt_pct(require_float(s, k, ARM)))
    end
    put("bg_removed_by_scan", fmt_int(require_int(s, "bg_removed_by_scan", ARM)))
    put("nonln_coverage_numerator", fmt_int(require_int(s, "nonln_coverage_numerator", ARM)))
    put("nonln_coverage_denominator", fmt_int(require_int(s, "nonln_coverage_denominator", ARM)))
    put("nonln_coverage_pct", fmt_pct(require_float(s, "nonln_coverage_pct", ARM)))
    for arm in ("fg_coverage", "bg_analysed_coverage")
        put("$(arm)_numerator", fmt_int(require_int(s, "$(arm)_numerator", ARM)))
        put("$(arm)_denominator", fmt_int(require_int(s, "$(arm)_denominator", ARM)))
        put("$(arm)_pct", fmt_pct(require_float(s, "$(arm)_pct", ARM)))
    end
    put("coverage_fisher_p_after_scan", fmt_p(require_float(s, "coverage_fisher_p_after_scan", ARM)))
    put("bg_frame_coverage_numerator", fmt_int(require_int(s, "bg_frame_coverage_numerator", ARM)))
    put("bg_frame_coverage_denominator", fmt_int(require_int(s, "bg_frame_coverage_denominator", ARM)))
    put("bg_frame_coverage_pct", fmt_pct(require_float(s, "bg_frame_coverage_pct", ARM)))
    put("coverage_fisher_p", fmt_p(require_float(s, "coverage_fisher_p", ARM)))
    for k in ("fg_relations_text_mined_pct", "bg_relations_text_mined_pct")
        put(k, fmt_pct(require_float(s, k, ARM)))
    end
    for a in ("fg", "bg"), k in ("taxa_tied_on_distinct", "taxa_tie_broken_by_mentions",
                                 "taxa_still_tied_after_mentions")
        put("$(a)_$(k)", fmt_int(require_int(s, "$(a)_$(k)", ARM)))
    end
    for a in ("fg", "bg"), k in ("isolation_terms_total", "taxa_with_isolation_term",
                                 "generic_terms_total", "taxa_with_generic_term")
        put("$(a)_$(k)", fmt_int(require_int(s, "$(a)_$(k)", ARM)))
    end
    put("phylum_strata_used", string(require_int(s, "phylum_strata_used_fg_vs_bg", ARM)))
    put("standardisation_heaviest_stratum", require_value(s, "standardisation_heaviest_stratum", ARM))
    put("standardisation_heaviest_stratum_weight_pct",
        fmt_pct(100 * require_float(s, "standardisation_heaviest_stratum_weight", ARM)))
    put("standardisation_heaviest_stratum_bg_taxa",
        fmt_int(require_int(s, "standardisation_heaviest_stratum_bg_taxa", ARM)))
    put("n_genera_fg", fmt_int(require_int(s, "n_genera_fg", ARM)))
    put("n_genera_bg", fmt_int(require_int(s, "n_genera_bg", ARM)))
    for k in ("background_licence_code", "omnicrobe_release", "obt_release", "uniprot_release")
        put(k, require_value(s, k, ARM))
    end
    put("panel_files", require_value(p, "panel_written", PROV))
    put("panel_categories", string(require_int(p, "panel_categories", PROV)))
    nd = require_value(p, "panel_categories_not_drawn", PROV)
    put("panel_categories_not_drawn", isempty(nd) ? "none" : nd)

    # Rank within the row's own block, 1 for the highest estimate, ties broken by name.
    function ranks(field)
        out_ = Dict{Tuple{String,String},Int}()
        for blk in (:habitat, :chemistry)
            rs = sort([x for x in plan.rows if x.block === blk], by = x -> (-getfield(x, field), x.category))
            for (i, x) in enumerate(rs)
                out_[(x.unit, x.category)] = i
            end
        end
        return out_
    end
    rank_crude, rank_std = ranks(:crude), ranks(:std)
    for r in order
        c = r.category
        put("$c.block", String(r.block))
        put("$c.rank_crude", string(rank_crude[(r.unit, r.category)]))
        put("$c.rank_std", string(rank_std[(r.unit, r.category)]))
        put("$c.a", fmt_count(r.a))
        put("$c.pct_fg", fmt_pct(r.pct_fg))
        put("$c.pct_bg", fmt_pct(r.pct_bg))
        put("$c.crude", fmt_log2(r.crude))
        put("$c.crude_lo", fmt_log2(r.crude_lo))
        put("$c.crude_hi", fmt_log2(r.crude_hi))
        put("$c.crude_ci", fmt_interval(r.crude_lo, r.crude_hi))
        put("$c.fisher_p", fmt_p(r.fisher_p))
        put("$c.rounded_for_fisher", r.rounded_for_fisher)
        put("$c.val", fmt_log2(r.std))
        put("$c.lo", fmt_log2(r.lo))
        put("$c.hi", fmt_log2(r.hi))
        put("$c.ci", fmt_interval(r.lo, r.hi))
        put("$c.width", fmt_ratio(r.hi - r.lo))
        # Which way the connector runs on this row. Adjustment has no fixed direction, and a
        # sentence about it has to be able to say which one this row took.
        put("$c.adjustment_direction",
            abs(r.std) > abs(r.crude) ? "away from zero" :
            abs(r.std) < abs(r.crude) ? "toward zero" : "unchanged")
        put("$c.fam", fmt_log2(r.tick))
        put("$c.fam_lo", fmt_log2(r.tick_lo))
        put("$c.fam_hi", fmt_log2(r.tick_hi))
        put("$c.fam_ci", fmt_interval(r.tick_lo, r.tick_hi))
        put("$c.colour", String(r.colour_class))
        put("$c.robust_sign", r.robust_sign)
        put("$c.breakers", isempty(r.breakers) ? "none" : r.breakers)
        put("$c.breaker_estimates", breaker_estimates(inp, r))
        put("$c.n_zero_stratum_cells", fmt_count(r.n_zero_stratum_cells))
        put("$c.min_bg_stratum_cell", fmt_count(r.min_bg_stratum_cell))
        put("$c.fam_n_zero_stratum_cells", fmt_count(r.tick_n_zero_stratum_cells))
        put("$c.fam_min_bg_stratum_cell", fmt_count(r.tick_min_bg_stratum_cell))
    end

    for c in grid_categories
        r = require_drawn(plan, c)
        taxa = across_specifications(inp.grid, registered, r.unit, c)
        seqs = across_specifications(inp.grid, registered, sequences_unit(r.unit), c)
        excl(x) = OR.interval_excludes_zero(x.lo, x.hi) === :excludes_zero
        # How many of those specifications are distinct estimates: a variant whose rule cannot
        # touch this row returns the primary bit for bit, and a count of specifications is not
        # a count of independent answers.
        distinct(xs) = length(unique([(x.std, x.lo, x.hi) for x in xs]))
        put("$c.taxa_specs_excluding_zero", string(count(excl, taxa)))
        put("$c.taxa_specs_distinct", string(distinct(taxa)))
        put("$c.sequences_specs_excluding_zero", string(count(excl, seqs)))
        put("$c.sequences_specs_excluding_zero_below", string(count(x -> excl(x) && x.hi < 0, seqs)))
        put("$c.sequences_specs_not_excluding_zero",
            (v = [x.variant for x in seqs if !excl(x)]; isempty(v) ? "none" : join(v, ";")))
        put("$c.sequences_specs_distinct", string(distinct(seqs)))
        su = rows_of(inp.stats, HE.LN_VS_BG, sequences_unit(r.unit))
        e = one_row(su[su.habitat .== c, :], "$(HE.LN_VS_BG) $(rowkey(sequences_unit(r.unit), c))")
        put("$c.seq_val", fmt_log2(e.log2_rr_std))
        put("$c.seq_lo", fmt_log2(e.ci_lo_std))
        put("$c.seq_hi", fmt_log2(e.ci_hi_std))
        put("$c.seq_ci", fmt_interval(e.ci_lo_std, e.ci_hi_std))
    end

    sr = require_drawn(plan, sharing_category)
    ratio = require_float(s, "sharing_ratio_nonln_over_ln", ARM)
    abs(ratio - sr.tick / sr.std) <= TOL ||
        error("$ARM's sharing ratio is not the drawn tick over the drawn estimate of " *
              "'$sharing_category'")
    sc = sharing_category
    put("$sc.sharing_ratio", fmt_ratio(ratio))
    lo_ = require_float(s, "sharing_ratio_ci_lo", ARM)
    hi_ = require_float(s, "sharing_ratio_ci_hi", ARM)
    put("$sc.sharing_ratio_lo", fmt_ratio(lo_))
    put("$sc.sharing_ratio_hi", fmt_ratio(hi_))
    put("$sc.sharing_ratio_ci", fmt_ratio(lo_) * " to " * fmt_ratio(hi_))
    put("$sc.sharing_code", require_value(s, "sharing_code", ARM))
    put("$sc.sharing_rule_used", require_value(s, "sharing_rule_used", ARM))

    built = first.(out)
    declared = caption_keys([r.category for r in order], selected)
    built == declared ||
        error("the exported keys differ from the declared ones: built but not declared " *
              join(setdiff(built, declared), ", ") * "; declared but not built " *
              join(setdiff(declared, built), ", "))
    return out
end

# ---------------------------------------------------------------------------
# Provenance of this run.
# ---------------------------------------------------------------------------

"""
    looks_like_a_local_path(s) -> Bool

A value that names a place on the machine that ran rather than a thing that was read: a drive
or root prefix, or any backslash-separated path. The provenance this script publishes records
file names, digests and declared names, none of which can contain one.
"""
looks_like_a_local_path(s::AbstractString) =
    occursin(r"^[A-Za-z]:[\\/]", s) || startswith(s, "/") || occursin('\\', s)

"""
    refuse_local_paths(rows) -> rows

`rows` unchanged, or an error naming every key whose value names a place on this machine. It
guards the provenance against an edit that records a path again; on well-formed rows it is
invisible, which is the point of it.
"""
function refuse_local_paths(rows)
    leaked = String[String(k) for (k, v) in rows if looks_like_a_local_path(v)]
    isempty(leaked) ||
        error("the run provenance would record a local path under " * join(leaked, ", ") *
              "; it records file names and digests only")
    return rows
end

"""
    provenance_rows(opt, inp, gates, plan; run_at) -> Vector{Pair{String,String}}

What this run read, what it checked and what it drew. Inputs are recorded by FILE NAME beside
their digest and whether the path was typed or defaulted: the digest is what pins the bytes, an
absolute path pins nothing and carries the local directory layout into a published file. The
rows are refused if any value still looks like one.
"""
function provenance_rows(opt, inp, gates, plan; run_at::AbstractString)
    out = Pair{String,String}[]
    put(k, v) = push!(out, String(k) => string(v))
    put("script", "scripts/habitat_figure.jl")
    put("run_at", run_at)
    put("julia_version", VERSION)
    put("figdir_given", opt.figdir_given)
    # The criteria file this run was held to, and what it selected. A figure whose categories
    # came from a file records which file and which categories, so two runs can be told apart
    # by what they were asked about and not only by what came out.
    put("input_categories_file", basename(opt.categories))
    put("input_categories_sha256", inp.selection.digest)
    put("selection_roles", join([r.role for r in HE.SELECTION_ROLES], ";"))
    for r in HE.SELECTION_ROLES
        v = getfield(inp.selection, Symbol(r.role))
        put("selection_$(r.role)", v isa Tuple ? join(v, ";") : v)
    end
    for n in INPUT_NAMES
        put("input_$(n)_file", basename(getfield(opt, n).path))
        put("input_$(n)_given", getfield(opt, n).given)
        put("input_$(n)_sha256", getfield(inp.files, n).sha256)
    end
    put("panel_dir_given", opt.panel_dir.given)
    for x in inp.panels
        put("panel_file_$(x.name)_sha256", x.sha256)
    end
    put("input_gates", join(INPUT_GATES, ";"))
    put("input_gates_failed", join([g.name for g in gates if !g.ok], ";"))
    for g in gates
        put("gate_$(g.name)", g.ok ? "pass" : "fail")
    end
    put("output_gate", OUTPUT_GATE)
    put("outputs", join(FIGURE_OUTPUTS, ";"))
    put("staging_dir", STAGING_DIR)
    put("min_n", MIN_N)
    order = panel_order(plan)
    put("drawn_rows_top_to_bottom", join([rowkey(r.unit, r.category) for r in order], ";"))
    put("colour_classes_top_to_bottom", join([String(r.colour_class) for r in order], ";"))
    s = inp.summary
    # The first of these is named after the selected phylum, so it is derived rather than
    # written out; the rest are selection-independent.
    for k in ("pct_$(inp.selection.key.standardisation_phylum)_bg_before_scan",
              "frame_size", "frame_taxa_with_usable_term",
              "bg_frame_coverage_numerator", "bg_frame_coverage_denominator",
              "bg_frame_coverage_pct", "coverage_fisher_p", "bg_removed_by_scan")
        put("summary_$k", get(s, k, ""))
    end
    put("row_mm", ROW_MM)
    put("axis_height_pt", plan.axis_height_pt)
    put("tick_dot_clearance_pt", tick_dot_clearance_pt())
    put("dot_ring_clearance_pt", dot_ring_clearance_pt())
    put("row_to_row_clearance_pt", row_to_row_clearance_pt())
    put("png_dpi", PNG_DPI)
    return refuse_local_paths(out)
end

write_key_values(path, pairs) =
    CSV.write(path, DataFrame(key = first.(pairs), value = last.(pairs)))

# ---------------------------------------------------------------------------
# Staging and publishing.
# ---------------------------------------------------------------------------

"""
    open_staging(figdir) -> String

This script's staging directory inside an existing figure directory, emptied. The figure
directory itself is never created here.
"""
function open_staging(figdir::AbstractString)
    isdir(figdir) || error("the figure directory '$figdir' does not exist")
    stage = joinpath(figdir, STAGING_DIR)
    mkpath(stage)
    for f in readdir(stage)
        rm(joinpath(stage, f); force = true, recursive = true)
    end
    return stage
end

"""
    output_gate(stage) -> NamedTuple

The staged set is exactly `FIGURE_OUTPUTS`, and no staged file is empty.
"""
function output_gate(stage::AbstractString)
    staged = sort(readdir(stage))
    problems = String[]
    staged == sort(collect(FIGURE_OUTPUTS)) ||
        push!(problems, "staged " * join(staged, ", ") * " where the declared outputs are " *
                        join(sort(collect(FIGURE_OUTPUTS)), ", "))
    for f in staged
        filesize(joinpath(stage, f)) > 0 || push!(problems, "$f is empty")
    end
    return (name = OUTPUT_GATE, ok = isempty(problems), message = join(problems, "; "))
end

"""
    publish!(stage, figdir, gates) -> Vector{String}

Move every staged file into the figure directory, but only if every gate passed. On a
failure nothing moves and the staging directory stays, named in the error, for inspection.
"""
function publish!(stage::AbstractString, figdir::AbstractString, gates)
    failed = String["$(g.name): $(g.message)" for g in gates if !g.ok]
    isempty(failed) ||
        error(join(failed, "\n") * "\nnothing was published; the staged files are in " * stage)
    written = String[]
    for name in sort(readdir(stage))
        mv(joinpath(stage, name), joinpath(figdir, name); force = true)
        push!(written, joinpath(figdir, name))
    end
    # The empty staging directory goes too. The removal can be refused while this process
    # runs, so it is retried briefly; if it is still refused, the empty directory stays, the
    # run says so, and the next run clears it.
    last_error = ""
    for _ in 1:10
        try
            rm(stage; recursive = true)
            break
        catch e
            last_error = sprint(showerror, e)
            sleep(0.2)
        end
    end
    isdir(stage) && println(stderr, "note: the empty staging directory ", stage,
                            " was left in place (", last_error, "); the next run clears it")
    return written
end

# ---------------------------------------------------------------------------
# The renderer, loaded only when a figure is drawn.
# ---------------------------------------------------------------------------

"""
    cairomakie() -> Module

The renderer, loaded at the moment it is needed. Nothing but `draw_figure` uses it, so
including this file for its rules never pays for it.
"""
cairomakie() = Core.eval(@__MODULE__, :(import CairoMakie; CairoMakie))

"""
    draw_figure(mk, plan, path_base) -> Vector{String}

Map a plan to marks and write it as pdf, svg and png. Every position, limit, colour and label
comes from the plan; nothing is computed here.
"""
function draw_figure(mk, plan, path_base::AbstractString)
    mk.activate!(type = "pdf")
    rows = plan.rows
    ys = [r.y for r in rows]
    axis_style = (ygridvisible = false, rightspinevisible = false, topspinevisible = false,
                  spinewidth = 0.6, titlesize = FS + 1, titlealign = :left, titlefont = :bold,
                  xlabelsize = FS, xticklabelsize = FS - 0.5, yticklabelsize = FS - 0.5,
                  xgridcolor = (COL_RULE, 0.5), xgridwidth = 0.4, xgridvisible = true,
                  height = plan.axis_height_pt)
    fig = mk.Figure(size = (plan.width_pt, plan.axis_height_pt + 200), fontsize = FS)

    axA = mk.Axis(fig[1, 1]; title = "a   Share of taxa", xlabel = "% of taxa",
                  yticks = (ys, [r.label for r in rows]),
                  xticks = collect(0:plan.share_step:plan.share_max), axis_style...)
    for r in rows
        mk.barplot!(axA, [r.y], [r.pct_bg]; direction = :x, width = BAR_BG_WIDTH,
                    color = (:white, 0.0), strokecolor = (COL_INK, 0.55), strokewidth = 0.5)
        mk.barplot!(axA, [r.y], [r.pct_fg]; direction = :x, width = BAR_FG_WIDTH,
                    color = COL_INK)
    end
    isnan(plan.yrule) ||
        mk.hlines!(axA, [plan.yrule]; color = (COL_RULE, 0.9), linewidth = 0.5, linestyle = :dot)
    mk.xlims!(axA, 0, plan.share_max)
    mk.ylims!(axA, plan.ylo, plan.ymax)

    axB = mk.Axis(fig[1, 2]; title = "b   Enrichment vs background", xlabel = "log₂ enrichment",
                  yticks = (ys, fill("", length(rows))), yticksvisible = false, axis_style...)
    mk.vlines!(axB, [0.0]; color = (COL_INK, 0.5), linewidth = 0.8)
    isnan(plan.yrule) ||
        mk.hlines!(axB, [plan.yrule]; color = (COL_RULE, 0.9), linewidth = 0.5, linestyle = :dot)
    for r in rows
        y = r.y
        mk.lines!(axB, [r.crude, r.std], [y + RING_DY, y + DOT_DY];
                  color = (COL_NEUTRAL, 0.5), linewidth = 0.5)
        mk.linesegments!(axB, [mk.Point2f(r.tick, y + tick_bottom()),
                               mk.Point2f(r.tick, y + tick_top())];
                         color = (COL_INK, 0.8), linewidth = TICK_LW)
        mk.scatter!(axB, [r.crude], [y + RING_DY]; color = :white,
                    strokecolor = (COL_INK, 0.6), strokewidth = RING_STROKE, markersize = RING_MS)
        mk.errorbars!(axB, [r.std], [y + DOT_DY], [r.std - r.lo], [r.hi - r.std];
                      direction = :x, whiskerwidth = WHISKER_PT, linewidth = 1.1, color = r.colour)
        mk.scatter!(axB, [r.std], [y + DOT_DY]; color = r.colour, markersize = DOT_MS)
    end
    mk.xlims!(axB, plan.xb_lo, plan.xb_hi)
    mk.ylims!(axB, plan.ylo, plan.ymax)
    for l in plan.block_labels
        mk.text!(axB, l.x, l.y; text = l.text, align = (:left, :top), fontsize = LABEL_FS,
                 color = (COL_NEUTRAL, 0.95), font = :bold)
    end

    element(key) =
        key === :carriers   ? mk.PolyElement(color = COL_INK) :
        key === :background ? mk.PolyElement(color = (:white, 0.0), strokecolor = (COL_INK, 0.55),
                                              strokewidth = 0.5) :
        key === :interval   ? mk.LineElement(color = (COL_INK, 0.8), linewidth = 1.1) :
        key === :ring       ? mk.MarkerElement(color = :white, marker = :circle,
                                               markersize = RING_MS, strokecolor = (COL_INK, 0.6),
                                               strokewidth = RING_STROKE) :
        key === :tick       ? mk.LineElement(color = (COL_INK, 0.8), linewidth = TICK_LW,
                                             points = mk.Point2f[(0.5, 0.15), (0.5, 0.85)]) :
        mk.MarkerElement(color = CLASS_COLOURS[key], marker = :circle, markersize = DOT_MS)
    lg = plan.legend
    mk.Legend(fig[2, 1:2],
              [[element(k) for (k, _) in lg.shares], [element(k) for (k, _) in lg.marks]],
              [[l for (_, l) in lg.shares], [l for (_, l) in lg.marks]],
              ["a", "b"];
              orientation = :horizontal, framevisible = false, labelsize = FS - 0.5,
              titlesize = FS - 0.5, titlefont = :bold, titlegap = 2, groupgap = 16,
              patchlabelgap = 3, colgap = 8, titleposition = :left, nbanks = 2,
              patchsize = (7.0f0, 7.0f0), rowgap = 2, padding = (0.0f0, 0.0f0, 0.0f0, 0.0f0))
    mk.colsize!(fig.layout, 1, mk.Relative(0.34))
    mk.colgap!(fig.layout, 6)
    mk.rowgap!(fig.layout, 1, 6)
    mk.resize_to_layout!(fig)

    written = String[]
    for (fmt, kw) in (("pdf", (pt_per_unit = 1,)), ("svg", (pt_per_unit = 1,)),
                      ("png", (px_per_unit = PNG_DPI / PT_PER_INCH,)))
        path = path_base * "." * fmt
        mk.save(path, fig; kw...)
        push!(written, path)
    end
    return written
end

# ---------------------------------------------------------------------------
# The run.
# ---------------------------------------------------------------------------

"""
    derive(inp) -> NamedTuple

Everything computed from the inputs once the gates have passed: the selected rows, the plan,
the source data and the exported values. Pure; any missing row or key raises here, before
anything is written. The category selection is `inp.selection`, read with the other inputs.
"""
function derive(inp)
    sel = select_rows(inp.stats)
    plan = figure_plan(inp, sel)
    return (sel = sel, plan = plan, source = source_data(inp, plan),
            values = caption_values(inp, sel, plan))
end

"""
    run_figure(opt; render, io, err) -> Int

One run: read, gate, derive, and unless `opt.check_only`, stage, draw, check the staged set
and publish. Returns the exit status. `render(plan, path_base)` replaces the renderer when
given; otherwise the renderer is loaded here and called in the latest world.
"""
function run_figure(opt; render = nothing, io::IO = stdout, err::IO = stderr)
    println(io, "habitat_figure.jl")
    println(io, "  figdir  ", opt.figdir, " (", opt.figdir_given, ")")
    println(io, "  ", rpad("categories", 24), opt.categories, " (required, no default)")
    for n in INPUT_NAMES
        println(io, "  ", rpad(String(n), 24), getfield(opt, n).path, " (", getfield(opt, n).given, ")")
    end
    inp = load_inputs(opt)
    println(io, "  selected: rules speak for ", join(inp.selection.governed, ";"),
            ";  grid counts ", join(inp.selection.grid_count, ";"),
            ";  standardisation phylum ", inp.selection.standardisation_phylum)
    gates = input_gates(inp)
    println(io, "\nINPUT GATES")
    for g in gates
        println(io, "  ", g.ok ? "pass  " : "FAILED", "  ", g.name)
    end
    failed = [g for g in gates if !g.ok]
    if !isempty(failed)
        for g in failed
            println(err, "gate '", g.name, "' failed: ", g.message)
        end
        println(err, "nothing was written")
        return 1
    end
    d = derive(inp)
    order = panel_order(d.plan)
    println(io, "\nROWS (top to bottom)")
    for r in order
        println(io, "  ", rpad(rowkey(r.unit, r.category), 36), String(r.colour_class))
    end
    println(io, "  held back: ", prose_list([r.habitat for r in eachrow(d.sel.held_back)]))
    if opt.check_only
        println(io, "\ncheck only: every gate passed and every value was derived; nothing was written")
        return 0
    end
    stage = open_staging(opt.figdir)
    CSV.write(joinpath(stage, SOURCE_DATA_FILE), d.source)
    write_key_values(joinpath(stage, VALUES_FILE), d.values)
    write_key_values(joinpath(stage, PROVENANCE_FILE),
                     provenance_rows(opt, inp, gates, d.plan; run_at = string(now())))
    base = joinpath(stage, FIGURE_BASE)
    if render === nothing
        mk = cairomakie()
        Base.invokelatest(draw_figure, mk, d.plan, base)
    else
        render(d.plan, base)
    end
    for p in publish!(stage, opt.figdir, [output_gate(stage)])
        println(io, "wrote ", p)
    end
    return 0
end

"""
    main(args; env, render, io, err) -> Int

The command line: `--help`, strict options, then one run. A command-line error or a run
error prints its message and returns 1.
"""
function main(args; env = ENV, render = nothing, io::IO = stdout, err::IO = stderr)
    opt = try
        parse_run_options(args; env = env)
    catch e
        e isa ErrorException || rethrow()
        println(err, "habitat_figure.jl: ", e.msg, "\n")
        print(err, USAGE)
        return 1
    end
    if opt.help
        print(io, USAGE)
        return 0
    end
    try
        return run_figure(opt; render = render, io = io, err = err)
    catch e
        e isa ErrorException || rethrow()
        println(err, "habitat_figure.jl: ", e.msg)
        return 1
    end
end

# ---------------------------------------------------------------------------
# Everything above is definitions: including this file runs nothing.
# ---------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    # The inputs that come from the data root, by name, before `main` opens anything. Only
    # the ones still at their REGISTRY DEFAULT are checked: a path the reader typed is
    # theirs, and `read_input` reports it as given. Without this a missing data root reaches
    # the reader as `habitat_figure.jl: missing input: <absolute path>` — no input name, no
    # layout path, and a path off one machine.
    #
    # In the program guard and not in `main`, because `require_or_exit` ends the process and
    # `main` is a function the suite calls and reads the return code of. The option reader is
    # run here and its failure swallowed: a bad command line stays `main`'s to report, with
    # the usage, and stays FIRST — a reader who mistyped `--figdir` must be told that and not
    # told about a file they never asked for.
    # The three figure-directory inputs are step 4's outputs and are not under the data root
    # at all, so they are not registry names and stay where they are.
    let opt = try
                  parse_run_options(ARGS)
              catch
                  nothing
              end
        if opt !== nothing && !opt.help
            opt.stats.given   == "default" && DataInputs.require_or_exit("enrichment_table")
            opt.summary.given == "default" && DataInputs.require_or_exit("arm_summary")
            opt.calls.given   == "default" && DataInputs.require_or_exit("habitat_calls")
        end
    end
    exit(main(ARGS))
end
