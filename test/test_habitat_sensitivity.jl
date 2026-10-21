#=
    Offline self-test for scripts/habitat_sensitivity.jl — the one-at-a-time sensitivity
    grid, the robustness flag it decides, and the checks that keep both honest. No network
    and no repository data file: every input below is built here, at the size of the
    property under test.

    The script is included behind its own `PROGRAM_FILE` guard, so including it defines the
    functions and runs nothing.

    EVERY FIXTURE IS ABSTRACT. Categories are role names such as "cat-a", flags are
    "chem-x", phyla are "phy-1", taxids are small arbitrary integers, and every interval in
    the classifier cases is a pair of numbers chosen so one branch is decided. Nothing here
    is a quantity measured on the data this script analyses: a suite that encodes a run's own
    numbers passes only until the next run, and it publishes what the run found.

    THE FAILURES THIS SUITE EXISTS TO PREVENT.

    1. A VARIANT THAT AGREED BY BEING ABSENT. The whole point of the grid is that a result
       survives every one-at-a-time change of the rules. A variant that failed to compute, or
       a category that variant has no row for, must count as a variant that did NOT reproduce
       the primary sign — never as one that quietly agreed. The classifier cases below pin
       that, and the grid refuses to hand the flag a table that is missing a variant at all.

    2. A FLAG THAT IGNORES ITS BREAKERS. The registered rule has no majority tier: one
       variant crossing zero is enough to make a category fragile. A classifier that counted
       votes, or that looked only at the primary, would still colour most rows the same way,
       so the cases below name exactly one breaker and require it to decide the answer.

    3. A UNIT CHECK THAT FILTERED ONE ARM. The thin-evidence subset has to bite on the
       foreground AND the background. A subset applied to the foreground alone is not a unit
       check — it is a different foreground compared against the whole background, and it
       moves every estimate in one direction for a reason that has nothing to do with
       evidence.

    4. A ROW WRITTEN TWICE, OR TWO MEASUREMENTS UNDER ONE KEY. The long table carries several
       units and several contrasts per specification. Without the unit and the contrast in
       the key, a genus row and a taxon row of the same specification and category are one
       key with two values, and the table keeps whichever was written last while looking
       complete.

    5. AN ANCHOR THAT PASSED BY NOT BEING CHECKED. A target whose row this run never wrote
       has not been reproduced. Reporting it as untested would let a wiring error read as a
       missing input, which is the one failure an anchor exists to catch.

    6. A SEAM THAT CARRIES A DIFFERENT RULE FROM THE ONE IT WAS GIVEN. Both aggregation axes
       reach the estimators through a rewritten weight column, and both rewrites are silent
       when wrong: a type column rebuilt from the wrong tallies, or a chemistry weight handed
       to a parser that reads labels, produces a plausible table under a rule nobody asked
       for.

    7. A LEVEL THAT MOVED NOBODY. The move-in level is the only one that hands taxa to the
       foreground, so it is the only one that can break a result rather than move it — and a
       level that silently moved nobody reports the unmoved specification under its name,
       which no row of the output can show.

    Run:  julia --project=. test/test_habitat_sensitivity.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Test, DataFrames, CSV, Printf

# The script goes into a module of its own, the way it includes the step it calls: the
# aggregate runner includes every suite into one process, and a script's constants must not
# collide with another suite's.
module HabitatGrid end
Base.include(HabitatGrid, joinpath(@__DIR__, "..", "scripts", "habitat_sensitivity.jl"))
const SG = HabitatGrid
const HE = HabitatGrid.HE
const ET = HabitatGrid.ET

# ---------------------------------------------------------------------------
# Fixtures. Role names throughout; the isolation label is read off the library, so the
# fixture cannot name an isolation source the library does not recognise.
# ---------------------------------------------------------------------------

const ISO = first(sort(collect(ET.ISOLATION_SOURCE)))
const PHYLA = ("phy-1", "phy-2")

# The five tally shapes, each chosen so that a different pair of tie rules disagrees on it.
#   :clear    one outright winner — every rule agrees
#   :mention  a count tie the mention step settles — `fractional` splits it, the others do not
#   :dead     a tie the mention step cannot settle — `unresolved` declines to call it
#   :second   the winner is another category, so a rule cannot be right by always saying cat-a
#   :thin     five categories of one term each — no category reaches a quarter of the tally
const TALLIES = Dict(
    :clear   => (["cat-a" => 3, "cat-b" => 1], ["cat-a" => 4, "cat-b" => 1]),
    :mention => (["cat-a" => 2, "cat-b" => 2], ["cat-a" => 5, "cat-b" => 1]),
    :dead    => (["cat-a" => 2, "cat-b" => 2], ["cat-a" => 2, "cat-b" => 2]),
    :second  => (["cat-b" => 4, "cat-c" => 1], ["cat-b" => 6, "cat-c" => 1]),
    :thin    => (["cat-a" => 1, "cat-b" => 1, "cat-c" => 1, "cat-d" => 1, "cat-e" => 1],
                 ["cat-a" => 1, "cat-b" => 1, "cat-c" => 1, "cat-d" => 1, "cat-e" => 1]))

# The chemistry shapes, as the five columns the call table writes them in: the flag label,
# the exclusive winner, the share map over flag COMBINATIONS, the number of chemistry-bearing
# terms, and the per-flag counts. They are consistent by construction, which is the property
# `validate_combinations!` checks on a real table, and the label is what the flags rule
# `any_term` produces from the share map, which is what the route check compares.
#
# EACH SHAPE EXISTS TO SEPARATE A PAIR OF RULES, exactly as the tally shapes do for the type
# axis. `:two` gives one flag a fifth of the evidence, which is below the share threshold and
# above the any-term one, so the share rule and the any-term rule disagree on it; `:spread`
# ties two combinations, which is what the fractional winner rule settles differently; and
# `:only_y` names a different flag from `:two`, so a configuration that rewrites the chemistry
# text has something to rewrite it INTO. Without shapes like these a rule variant produces
# the registered specification's own rows, and a fixture in which a level cannot move
# anything cannot test the rule that refuses such a level.
const CHEMISTRY = Dict(
    :none   => ("", "", "", 0, ""),
    :one    => ("chem-x", "chem-x", "chem-x:1.0", 2, "chem-x:2"),
    :two    => ("chem-x+chem-y", "chem-x", "chem-x:0.8;chem-x+chem-y:0.2", 5,
                "chem-x:5;chem-y:1"),
    :spread => ("chem-x+chem-y", "chem-y", "chem-x:0.5;chem-y:0.5", 2, "chem-x:1;chem-y:1"),
    :only_y => ("chem-y", "chem-y", "chem-y:1.0", 2, "chem-y:2"))

const SHAPES = (:clear, :mention, :dead, :second, :thin)
const CHEMS  = (:one, :two, :spread, :none)

"""
    evidence() -> Vector{NamedTuple}

Twelve census carriers and sixteen frame taxa.

The last two carriers hold census sequences that are not lanthanide-type, so the comparator
arms are populated rather than empty. The first eight carriers and the first ten frame taxa
rest on several distinct terms and the rest do not, so the thin-evidence subset bites on BOTH
arms — which is the property that subset has to have.
"""
function evidence()
    out = NamedTuple[]
    for (k, t) in enumerate(11:22)
        push!(out, (taxid = t, census = true, n_ln = t >= 21 ? 0 : 2, n_cen = 2,
                    shape = SHAPES[mod1(k, 5)], n_distinct = t <= 18 ? 4 : 1,
                    chem = CHEMS[mod1(k, 4)], phylum = PHYLA[mod1(k, 2)],
                    genus = "gen-$(mod1(k, 5))"))
    end
    for (k, t) in enumerate(101:116)
        push!(out, (taxid = t, census = false, n_ln = 0, n_cen = 0,
                    shape = SHAPES[mod1(k + 2, 5)], n_distinct = t <= 110 ? 4 : 1,
                    chem = CHEMS[mod1(k + 1, 4)], phylum = PHYLA[mod1(k, 2)],
                    genus = "gen-$(mod1(k, 7))"))
    end
    return out
end

"""
    fixture_calls() -> DataFrame

The taxon-call table in the six configurations the grid reads, built so that each axis moves
what that axis is supposed to move and nothing else:

* the generic-term configuration changes the TYPE tally (one category gains terms) AND the
  chemistry evidence of the taxa that carry one flag, so it moves both units while the arms
  stay put — which is what keeping a generic term does: the term is read for its habitat and
  for its condition, so a filter on terms cannot reach one and not the other;
* the two chemistry-text configurations change only the CHEMISTRY columns, so they move the
  chemistry rows alone;
* the isolation configuration marks one taxon of each arm as an isolation source, so the
  usable mask — and therefore n1 and n2 — moves;
* the earlier configuration carries the exclusive chemistry winner the anchor reads.

`type_weights` is written by the library's own rule and formatter, so the column is what the
registered rule would have written. A fixture that wrote it by hand would be checking the
grid's reconstruction against a second guess at the rule.
"""
function fixture_calls()
    rows = NamedTuple[]
    for config in ("primary", "text_first", "text_object", "generic_none", "iso_taxon",
                   "legacy")
        for e in evidence()
            (cpairs, mpairs) = TALLIES[e.shape]
            counts = Dict{String,Float64}(k => float(v) for (k, v) in cpairs)
            mentions = Dict{String,Float64}(k => float(v) for (k, v) in mpairs)
            if config == "generic_none"
                counts["cat-c"] = get(counts, "cat-c", 0.0) + 2.0
                mentions["cat-c"] = get(mentions, "cat-c", 0.0) + 2.0
            end
            (label, weights, _, _) = ET.resolve(counts, mentions, :mentions_fractional)
            ht = (config == "iso_taxon" && e.taxid in (15, 111)) ? ISO : label
            shape = e.chem
            config == "text_first" && shape === :two && (shape = :one)
            # AND IT REWRITES THE TIED TAXA TOO, which the chemistry SHARE rule does not:
            # both a fifth of the evidence and a half of it are rewritten by a configuration
            # that reads different text, while the share rule only drops the fifth. Without
            # this the fixture's `X:text_first` and `C:share` are the same measurement under
            # two names — an accidental clone, which is precisely what the pairwise guard
            # exists to refuse, and a fixture that contains one cannot exercise it. Here they
            # stay two different measurements, for that same reason, so the guard has a margin
            # to measure rather than a clone to refuse.
            config == "text_first" && shape === :spread && (shape = :only_y)
            config == "text_object" && shape === :two && (shape = :only_y)
            # A generic term is read for its condition as well as for its habitat, so the
            # configuration that keeps them adds chemistry evidence too. Without this the
            # fixture's generic axis reaches one unit where the real one reaches both, and
            # the declared axis map could not be exercised on it.
            config == "generic_none" && shape === :one && (shape = :two)
            (flags, winner, profile, n_terms, flag_counts) = CHEMISTRY[shape]
            push!(rows, (config = config, taxid = e.taxid,
                         in_census = e.census, in_pool = !e.census,
                         n_census_seqs = e.n_cen, n_ln_seqs = e.n_ln,
                         habitat_type = ht,
                         type_weights = ET._weights_str(weights),
                         type_counts = ET._weights_str(counts),
                         type_mentions = ET._weights_str(mentions),
                         n_distinct = e.n_distinct,
                         kingdom = e.census ? "Bacteria" : "",
                         domain = "Bacteria",
                         phylum = e.phylum, genus = e.genus,
                         organism = "Org $(e.taxid)",
                         chemistry_flags = flags, chemistry = winner,
                         chemistry_profile = profile, n_chem_terms = n_terms,
                         chemistry_flag_counts = flag_counts))
        end
    end
    return DataFrame(rows)
end

"""
    fixture_scan(; taxa) -> DataFrame

The background carrier scan over the frame, with each background rule given a different,
non-nested set — so a level that read the wrong column could not accidentally agree with the
level it was supposed to differ from. The lanthanide-classed taxa are the pool the move-in
level moves, and none of them carries a hit the any-hit rule would already have removed, so
under that level the foreground GROWS and the background loses the same taxa.

The sets also differ in the HABITAT and the CHEMISTRY the taxa they hold carry. Two
background rules that swap one taxon for another of the same habitat shape, or of the same
chemistry profile, move one unit and leave the other exactly where it was — an accidental
orthogonality that no background rule of this analysis has, and one that would let a level
lose its wiring to a unit unnoticed.
"""
fixture_scan(; taxa = collect(101:116)) =
    DataFrame(taxid = taxa,
              n_candidates = [t in (101, 102, 103, 104, 105) ? 2 : 0 for t in taxa],
              n_hits = [t in (101, 102) ? 2 : 0 for t in taxa],
              hit_classes = [t in (103, 104, 105) ? "Lanthanide" : "" for t in taxa],
              in_ipr017512_ge500 = [t in (101, 108) for t in taxa],
              n_hits_ln_or_ca = [t in (101, 103) ? 1 : 0 for t in taxa])

fixture_context(; calls = fixture_calls(), scan = fixture_scan(),
                  movein = nothing) =
    SG.context(; calls = calls, scan = scan, ca_class_taxids = nothing,
               movein_taxids = movein === nothing ?
                               HE.exclusion_taxids(scan, "lntype") : movein)

# Interval shorthands for the classifier cases. Nothing here is measured: each pair is two
# numbers chosen so that one branch of the rule is decided by it.
const POS  = (0.4, 1.2)        # excludes zero, enriched
const NEG  = (-1.2, -0.4)      # excludes zero, depleted
const ZED  = (-0.3, 0.9)       # includes zero, point estimate above it
const ZEDN = (-0.9, 0.3)       # includes zero, point estimate below it
const GONE = (NaN, NaN)        # no interval at all

"""
    flag_rows(primary, variants; unit, category) -> Vector{NamedTuple}

The rows the flag reads, by hand: one category on one unit, with the intervals given. Used
where the property under test is what the flag does with a table, so the table is built
rather than measured.
"""
function flag_rows(primary, variants; unit = "taxa", category = "cat-a")
    # every quantity a grid row carries, because what makes two rows the same estimate is a
    # property of all of them and a fixture that carried only the interval would be answering
    # a narrower question than the run asks
    row(v, b) = (variant = v, unit = unit, category = category, contrast = HE.LN_VS_BG,
                 a = 1.0, b = 2.0, n1 = 1, n2 = 2,
                 crude = (b[1] + b[2]) / 2, crude_lo = b[1], crude_hi = b[2],
                 std = (b[1] + b[2]) / 2, lo = b[1], hi = b[2], feeds_flag = true)
    return vcat([row("primary", primary)], [row(v, b) for (v, b) in variants])
end

# THE REGISTERED SENSITIVITY SET, BY NAME. The script enumerates the set from its own tuples,
# so the set and the thing that checks the set are the same object: a rename of two members
# that keeps the count keeps every count-based check green, and a level that is the registered
# specification under another name then joins the set as a variant that cannot break anything.
# So the NAMES are pinned here, outside the script, sorted, and nothing below reads them from
# the script. Each is a rule name — an axis letter and a level of that axis — and none is a
# quantity measured on any data.
const REGISTERED_VARIANTS = ["B:hmmsearch_movein", "B:ipr017512", "B:none",
                             "C:share", "C:winner_fractional",
                             "G:generic_none", "I:iso_taxon",
                             "T:alphabetical", "T:any_type", "T:fractional", "T:share",
                             "T:unresolved",
                             "X:text_first", "X:text_object"]

"""
    respec(spec; field = ...) -> Spec

One specification with some of its fields replaced. Every attack below is a specification
that differs from a registered one in ONE place — a name moved to another axis, a level that
runs a rule from outside the set, a field the comparison is not looking at — so they are
written as edits of the registered set rather than as thirteen positional arguments each,
where a mistyped field would quietly make the attack a different attack.
"""
respec(s::SG.Spec; name = s.name, axis = s.axis, level = s.level, config = s.config,
       tie_rule = s.tie_rule, chem_col = s.chem_col, chem_rule = s.chem_rule,
       bg_mode = s.bg_mode, subset = s.subset, role = s.role, emit = s.emit,
       unit_suffix = s.unit_suffix, note = s.note) =
    SG.Spec(name, axis, level, config, tie_rule, chem_col, chem_rule, bg_mode, subset, role,
            emit, unit_suffix, note)

# What one of the guards said when it refused, or "" when it did not refuse at all. Every
# attack below asserts on the TEXT, because a guard that fires without naming the offender
# leaves a reader of the failure no way to tell which of fourteen levels is wrong.
refusal(f, args...; kw...) = try; f(args...; kw...); "" catch e; sprint(showerror, e) end

# A COMPLETE, PASSING SET OF INPUTS TO THE PUBLISH GATE, with any of them replaced. Every
# attack on the gate below moves exactly ONE of these off its good value, so the gate that
# answers for it is the only thing that differs between two runs — and a gate whose entry was
# deleted, or whose `ok` was replaced by `true`, then reports no failure where one input is
# bad. Nothing here is measured: the counts are small numbers chosen so one branch is decided.
const GOOD_GATE_INPUT = (target_failures = String[], unanchored_levels = String[],
                         source_checks = [(file = "an_implementation.csv",
                                           status = "matches", declared = "the declared bytes",
                                           found = "the declared bytes", ok = true)],
                         reference_ok = true, reference_rows_compared = 3,
                         reference_rows_expected = 3, stale_panels = String[],
                         figdir = "a figure directory")
gates_for(; kw...) = SG.publish_gates(; merge(GOOD_GATE_INPUT, NamedTuple(kw))...)

# The fixture grid, run once: every specification the script defines, on the fixture above.
const CTX   = fixture_context()
# The anchor levels are RETIRED (ANCHOR_LEVELS_RETIRED): they reproduce rows of a census other
# than the one a production run reads, so the production grid holds none of them. The suite
# builds its grid with `retired = false` on purpose: the anchor machinery must keep working,
# because a run against the census those levels describe produces them by flipping that one
# constant, and a retired level whose code has quietly rotted would be discovered only then.
# The production setting is pinned separately, in section 6.
const SPECS = SG.all_specs(retired = false)
const S     = SG.variant_names(SPECS)
const GRID  = SG.run_grid(CTX, SPECS)

"""
    md_cell(md, unit, category, order, variant) -> String

One cell of the written sensitivity table, read back out of the markdown. Used where the
property under test is what the TABLE says — a cell that repeats the primary's numbers and a
cell that says the axis never reached it are the same run and two different documents.
"""
function md_cell(md, unit, category, order, variant)
    lines = split(md, "\n")
    i = findfirst(==("## unit: " * unit), lines)
    i === nothing && return missing
    j = findnext(l -> startswith(l, "## unit: "), lines, i + 1)
    stop = j === nothing ? length(lines) : j - 1
    k = findnext(l -> startswith(l, "| " * category * " |"), lines, i)
    (k === nothing || k > stop) && return missing
    cells = strip.(split(lines[k], "|"))
    # "| category | c1 | ... |" splits to ["", category, c1, ..., ""]
    return String(cells[2 + findfirst(==(variant), order)])
end

# The source a fixture file declares its targets were transcribed from. The name is one no
# fixture directory holds, so the default case is the one a run meets when the named source is
# not distributed beside it — a source that is named and not present — and the digest is 64 hex
# characters that hash nothing, because what is under test is what the reader does with the
# declaration and not whether this suite can reproduce a checksum.
const FIXTURE_SOURCE = "an_implementation_that_never_read_this_script.csv"
const FIXTURE_DIGEST = repeat("ab", 32)

# The rows every targets file must carry beside its targets: the registered set by name, the
# three counts the file states about itself, and the file the targets were transcribed from
# with its digest. They are written here from the same pinned list the suite uses for
# everything else, so a fixture file is a file the run would accept. Each piece can be left
# out on its own, because each is refused for a different reason and a switch that removed all
# of them could not tell those reasons apart.
function declaration_rows(n_targets; names = REGISTERED_VARIANTS, expect_variants = nothing,
                          reference_rows = 0, source = FIXTURE_SOURCE,
                          digest = FIXTURE_DIGEST, with_expect = true, with_source = true,
                          with_variants = true)
    out = Any[]
    # `with_expect` is true, false, or the names of the counts to write — each `expect` row is
    # refused for its own reason, and a switch that could only remove all three could not tell
    # those reasons apart.
    wants(q) = with_expect isa Bool ? with_expect : (q in with_expect)
    wants(:targets) &&
        push!(out, ("expect", "", "", "", "targets", float(n_targets), missing,
                    "how many targets", ""))
    wants(:variants) &&
        push!(out, ("expect", "", "", "", "variants",
                    float(expect_variants === nothing ? length(names) : expect_variants),
                    missing, "how many levels the registered set holds", ""))
    wants(:reference_rows) &&
        push!(out, ("expect", "", "", "", "reference_rows", float(reference_rows), missing,
                    "how many rows the comparison against the earlier run covers", ""))
    with_source &&
        push!(out, ("source", source, "", "", "", missing, missing,
                    "an implementation written without reading this script", digest))
    with_variants &&
        append!(out, [("variant", v, "", "", "", missing, missing, "a registered level", "")
                      for v in names])
    return out
end

# A targets table in the shape the option reads. `value`/`tol` are arbitrary; what is under
# test is what the checker does with a target it meets, misses, or cannot find. `declare`
# appends the declaration the run requires — turned off, in whole or one piece at a time,
# where the property under test is what happens to a file that does not carry it. A row given
# with eight fields carries no digest, which is every kind of row but `source`.
function write_targets(dir; rows, declare = true, names = REGISTERED_VARIANTS,
                       expect = count(r -> r[1] == "target", rows), name = "targets.csv",
                       kw...)
    pad(r) = length(r) == 9 ? r : (r..., "")
    all = declare ?
          vcat(Any[pad(r) for r in rows],
               declaration_rows(expect; names = names, kw...)) :
          Any[pad(r) for r in rows]
    p = joinpath(dir, name)
    CSV.write(p, DataFrame(kind = [r[1] for r in all], anchor = [r[2] for r in all],
                           unit = [r[3] for r in all], category = [r[4] for r in all],
                           quantity = [r[5] for r in all], value = [r[6] for r in all],
                           tol = [r[7] for r in all], note = [r[8] for r in all],
                           digest = [r[9] for r in all]))
    return p
end

@testset "habitat_sensitivity grid" begin

# ---------------------------------------------------------------------------
# 1. The robustness rule.
# ---------------------------------------------------------------------------

@testset "a category is robust only when every variant reproduces its sign" begin
    c = SG.classify_category(POS..., ["v1" => POS, "v2" => POS, "v3" => POS])
    @test c.flag == "+"
    @test isempty(c.breakers)

    d = SG.classify_category(NEG..., ["v1" => NEG, "v2" => NEG])
    @test d.flag == "-"
    @test isempty(d.breakers)
end

@testset "exactly one variant crossing zero makes a category fragile" begin
    c = SG.classify_category(POS..., ["v1" => POS, "v2" => ZED, "v3" => POS])
    @test c.flag == "fragile"
    @test c.breakers == ["v2"]
    @test !isempty(c.meaning)

    # The same on the depleted side, so the rule cannot be right by only looking upward.
    d = SG.classify_category(NEG..., ["v1" => NEG, "v2" => ZED])
    @test d.flag == "fragile"
    @test d.breakers == ["v2"]
end

@testset "a variant that excludes zero on the OTHER side is a breaker" begin
    # A decided interval of the opposite sign has not reproduced the primary sign, and a
    # test that only asked "does it exclude zero" would call this robust.
    c = SG.classify_category(POS..., ["v1" => POS, "v2" => NEG])
    @test c.flag == "fragile"
    @test c.breakers == ["v2"]
end

@testset "a primary that includes zero is null, whatever the variants say" begin
    c = SG.classify_category(ZED..., ["v1" => ZED, "v2" => ZED])
    @test c.flag == "null"
    @test isempty(c.breakers)
    @test c.primary_excludes_zero == false

    # A variant that DOES distinguish the category does not make the row fragile. The
    # fragile shape is reserved for a claim the registered specification actually makes and
    # a variant does not reproduce; a primary that includes zero makes no claim, so the row
    # is neutral and the variants are reported beside it rather than in place of the word.
    d = SG.classify_category(ZED..., ["v1" => ZED, "v2" => POS])
    @test d.flag == "null"
    @test d.primary_excludes_zero == false
    @test d.breakers == ["v2"]            # still named: they are the informative part
    @test d.meaning != c.meaning          # the two senses of "breaker" are named apart

    # And the word cannot be reached from this side at all, in either direction.
    @test SG.classify_category(ZED..., ["v1" => NEG]).flag == "null"
    @test SG.classify_category(ZEDN..., ["v1" => POS, "v2" => NEG]).flag == "null"
end

@testset "primary_excludes_zero is the fact the word rests on" begin
    # A later step keys colour off this boolean rather than off the word, so it has to be
    # true exactly where the primary interval excludes zero and nowhere else.
    for (primary, want) in ((POS, true), (NEG, true), (ZED, false), (ZEDN, false),
                            (GONE, false))
        @test SG.classify_category(primary..., ["v1" => POS]).primary_excludes_zero == want
    end
    # fragile is reachable only where the boolean is true
    @test SG.classify_category(POS..., ["v1" => ZED]).flag == "fragile"
    @test SG.classify_category(POS..., ["v1" => ZED]).primary_excludes_zero == true
end

@testset "a variant with no interval is a breaker, never a variant that agreed" begin
    c = SG.classify_category(POS..., ["v1" => POS, "v2" => GONE])
    @test c.flag == "fragile"
    @test c.breakers == ["v2"]

    # ... and it is not counted as a variant that excludes zero either.
    @test SG.classify_category(ZED..., ["v1" => ZED, "v2" => GONE]).flag == "null"
end

@testset "no primary interval is undecidable, and names the variants that do have one" begin
    c = SG.classify_category(GONE..., ["v1" => POS, "v2" => GONE, "v3" => ZED])
    @test c.flag == "undecidable"
    # The variants that measured the category, not the ones that never saw it: naming all
    # three would report v2 as having disagreed about something it never computed.
    @test c.breakers == ["v1", "v3"]
end

@testset "the flag reads the table, and a missing row cannot pass as agreement" begin
    names = ["v1", "v2", "v3"]
    rows = flag_rows(POS, ["v1" => POS, "v2" => POS, "v3" => POS])
    f = SG.robustness(rows, names)[1]
    @test f.robust_sign == "+"
    @test f.variants_expected == 3
    @test f.variants_with_a_row == 3

    # Drop one variant's row entirely. It is now a variant with no interval, which is a
    # breaker — the failure the whole rule exists to prevent.
    thinned = [r for r in rows if r.variant != "v3"]
    g = SG.robustness(thinned, names)[1]
    @test g.robust_sign == "fragile"
    @test g.breaking_variants == "v3"
    @test g.variants_with_a_row == 2
end

@testset "the depletion tier is decided by the same grid" begin
    # The tier function has no grid of its own to consult; it is handed the variants
    # here, so a missing one leaves it undecidable rather than agreeing.
    t1 = SG.robustness(flag_rows(NEG, ["v1" => NEG, "v2" => NEG]), ["v1", "v2"])[1]
    @test t1.depletion_tier == "t1"

    t2 = SG.robustness(flag_rows(NEG, ["v1" => NEG, "v2" => ZED]), ["v1", "v2"])[1]
    @test t2.depletion_tier == "t2"
    @test t2.depletion_tier_breakers == "v2"

    t3 = SG.robustness(flag_rows(ZEDN, ["v1" => ZEDN, "v2" => ZEDN]), ["v1", "v2"])[1]
    @test t3.depletion_tier == "t3"

    gone = SG.robustness([r for r in flag_rows(NEG, ["v1" => NEG, "v2" => NEG])
                          if r.variant != "v2"], ["v1", "v2"])[1]
    @test gone.depletion_tier == "t2"      # an absent variant is a breaker, not an agreement
end

@testset "the depletion tier is left empty where there is no depletion to judge" begin
    # The tier's three levels say "state it", "state it tentatively" and "drop the claim"
    # about a category the registered specification puts BELOW the background. Printed on a
    # category it puts above, the worst of them reads as a verdict against a claim nobody
    # made — and it would sit in the column beside a row flagged robustly enriched.
    up = SG.robustness(flag_rows(POS, ["v1" => POS, "v2" => POS]), ["v1", "v2"])[1]
    @test up.robust_sign == "+"
    @test up.depletion_tier == ""
    @test up.depletion_tier_breakers == ""

    # A neutral row whose estimate is above zero is equally not a depletion.
    @test SG.robustness(flag_rows(ZED, ["v1" => ZED]), ["v1"])[1].depletion_tier == ""

    # Below zero the tier is emitted, because there is then something for it to be about.
    @test SG.robustness(flag_rows(ZEDN, ["v1" => ZEDN]), ["v1"])[1].depletion_tier == "t3"
    @test SG.robustness(flag_rows(NEG, ["v1" => NEG]), ["v1"])[1].depletion_tier == "t1"
end

# ---------------------------------------------------------------------------
# 2. The grid, on the fixture.
# ---------------------------------------------------------------------------

@testset "the registered set is exactly the fourteen levels, by name" begin
    # NOT the count. A count-preserving edit — two levels out, two re-runs of the registered
    # specification wearing variant names in — passes every count, meets every anchor, and
    # hands a category two variants that cannot break it. What has to be pinned is membership.
    @test sort(S) == REGISTERED_VARIANTS
    @test length(unique(S)) == length(S)
    @test sort(SG.variant_names(SG.all_specs())) == REGISTERED_VARIANTS

    # The registered levels of the two axes that carry a level by name are not among their
    # own variants — the script says so in a comment, and a comment is not a check.
    @test !(HE.PRIMARY_MODE in SG.BACKGROUND_LEVELS)
    @test !(SG.PRIMARY_CONFIG in SG.TEXT_CONFIGS)
    @test !(SG.PRIMARY_CONFIG in (SG.GENERIC_CONFIG, SG.ISOLATION_CONFIG))
    # and a disclosed background line is outside the set, not a fifteenth variant
    @test all(m -> !("B:$m" in S), SG.DISCLOSED_BACKGROUND_LEVELS)
end

@testset "no variant is the registered specification under another name" begin
    specs = SG.all_specs()
    @test SG.assert_variants_vary_the_primary(specs) == length(S)

    # Every variant differs from the primary in at least one of the four things a
    # specification is made of, and they are all different from each other.
    p = SG.defining_level(SG.primary_spec())
    levels = [SG.defining_level(x) for x in specs if x.role == SG.ROLE_VARIANT]
    @test all(!=(p), levels)
    @test length(unique(levels)) == length(levels)

    # THE COUNT-PRESERVING SWAP, built here: two levels leave and two re-runs of the
    # registered specification take their names. The count, the uniqueness and every anchor
    # survive it; this is the check that does not.
    doctored = [x for x in specs if !(x.name in ("B:none", "X:text_object"))]
    for name in ("B:none", "X:text_object")
        push!(doctored, SG.Spec(name, "B", "filler", SG.PRIMARY_CONFIG, nothing,
                                :chemistry_flags, nothing, HE.PRIMARY_MODE, :all,
                                SG.ROLE_VARIANT, :both, "", "the primary, renamed"))
    end
    @test length(SG.variant_names(doctored)) == length(S)      # the count survives
    @test sort(SG.variant_names(doctored)) == REGISTERED_VARIANTS   # so do the names
    err = try; SG.assert_variants_vary_the_primary(doctored); ""
          catch e; sprint(showerror, e); end
    @test occursin("B:none", err) && occursin("X:text_object", err)
end

@testset "a variant's name describes the rule it runs" begin
    # EVERY PUBLISHED STATEMENT ABOUT A LEVEL IS MADE UNDER ITS NAME — the breaker list
    # beside a category, the column heads of the written table, the levels named as excluding
    # zero where the primary does not. A name that has drifted from the computation is
    # therefore not a cosmetic defect: it publishes a correct number under the wrong rule.
    @test SG.assert_variant_names_describe_what_they_run(SPECS) == length(S)

    # the axes the set uses are exactly the axes the script declares — neither an undeclared
    # letter nor a declared axis nobody varies
    @test sort(unique([s.axis for s in SPECS if s.role == SG.ROLE_VARIANT])) ==
          sort([a.axis for a in SG.VARIANT_AXES])
    for s in SPECS
        s.role == SG.ROLE_VARIANT || continue
        a = SG.axis_declaration(s.axis)
        @test a !== nothing
        @test s.name == "$(s.axis):$(s.level)"
        @test s.level in a.levels
        @test SG.field_level(s, a.field) == s.level
    end

    # THE NAME/AXIS SWAP: two names are exchanged and the SET of names is untouched. The
    # count passes, the pinned name list passes, the comparison against the primary passes,
    # every number of the run is the same number — and the flag table names the other rule as
    # the one that broke a category.
    swap = Dict("B:none" => "T:alphabetical", "T:alphabetical" => "B:none")
    swapped = [haskey(swap, s.name) && s.role == SG.ROLE_VARIANT ?
               respec(s; name = swap[s.name]) : s for s in SPECS]
    @test sort(SG.variant_names(swapped)) == REGISTERED_VARIANTS
    @test SG.assert_variants_vary_the_primary(swapped) == length(S)
    err = refusal(SG.assert_variant_names_describe_what_they_run, swapped)
    @test occursin("B:none", err) && occursin("T:alphabetical", err)

    # THE IMPOSTOR LEVEL: a registered name running a background rule the recorded decision
    # deliberately keeps OUTSIDE the set. Nothing about the declaration is otherwise wrong —
    # it varies exactly one field, it differs from the primary, it is nobody's clone — and
    # every number it publishes travels under the name of the level it replaced.
    outside = first(SG.DISCLOSED_BACKGROUND_LEVELS)
    impostor = [s.name == "B:ipr017512" ? respec(s; bg_mode = outside) : s for s in SPECS]
    @test SG.assert_variants_vary_the_primary(impostor) == length(S)
    @test SG.assert_variants_differ_from_each_other(impostor) == length(S)
    e2 = refusal(SG.assert_variant_names_describe_what_they_run, impostor)
    @test occursin("B:ipr017512", e2) && occursin(outside, e2)

    # ... and the same rule renamed consistently onto its own level is refused too, because
    # the level itself is not one of that axis's registered levels
    renamed = [s.name == "B:ipr017512" ?
               respec(s; name = "B:$outside", level = outside, bg_mode = outside) : s
               for s in SPECS]
    @test occursin("not a registered level", refusal(
        SG.assert_variant_names_describe_what_they_run, renamed))

    # A FIELD THE COMPARISON DOES NOT LOOK AT. Both of these produce a full grid of plausible
    # rows under a registered name: one measured on a different population, one reading a
    # different column. Neither is the rule the name names, and the one-at-a-time property
    # is what refuses them.
    for (field, changed) in ((:subset, respec(SPECS[2]; subset = :min_distinct)),
                             (:chem_col, respec(SPECS[2]; chem_col = HE.ANCHOR_CHEM_COL)))
        @test SPECS[2].role == SG.ROLE_VARIANT
        e3 = refusal(SG.assert_variant_names_describe_what_they_run,
                     [changed; SPECS[3:end]])
        @test occursin(SPECS[2].name, e3)
        @test occursin(String(field), e3)
    end

    # an axis letter that is not a registered axis is refused by name, not ignored
    @test occursin("not one of the registered axes", refusal(
        SG.assert_variant_names_describe_what_they_run,
        [respec(SPECS[2]; axis = "Z", name = "Z:$(SPECS[2].level)")]))
end

@testset "two variants that are the same specification are refused" begin
    @test SG.assert_variants_differ_from_each_other(SPECS) == length(S)

    # CLONE-YOUR-SIBLING. Each of the two keeps its own name and runs a SIBLING's level.
    # Every guard that compares a variant with the primary still passes — each still differs
    # from the registered specification, each still moves estimates — the fourteen names are
    # present and unique, and every anchor is met. What the pair has done is delete two
    # levels from the set while leaving their columns standing, so a category whose only
    # breakers were those two levels is published as robust with no breakers at all.
    cloned = [s.name == "B:none" ? respec(s; bg_mode = "ipr017512") :
              s.name == "X:text_object" ? respec(s; config = "text_first") : s for s in SPECS]
    @test sort(SG.variant_names(cloned)) == REGISTERED_VARIANTS
    @test SG.assert_variants_vary_the_primary(cloned) == length(S)     # still silent
    err = refusal(SG.assert_variants_differ_from_each_other, cloned)
    @test occursin("B:none", err) && occursin("B:ipr017512", err)
    @test occursin("X:text_object", err) && occursin("X:text_first", err)

    # The count it returns is the number of DISTINCT specifications, so a clone shows up as a
    # set that is smaller than its own name list even where the error is not read.
    @test length(unique([SG.defining_level(s) for s in cloned
                         if s.role == SG.ROLE_VARIANT])) == length(S) - 2
end

@testset "the set the run enumerated is the set that was declared to it" begin
    # THE ONE CHECK THAT IS NOT MADE INSIDE THE SCRIPT. Everything else compares one part of
    # the script with another, and the axis table reads its levels off the same tuples the
    # specifications are built from — so a level deleted from a tuple is self-consistent: the
    # wiring and the vocabulary move together and nothing is left to disagree. What it costs
    # is a column, and a category whose breakers were the levels that left is published with a
    # shorter breaker list or none at all.
    @test SG.assert_the_set_is_the_declared_one(S, REGISTERED_VARIANTS) == length(S)
    @test SG.assert_the_set_is_the_declared_one(reverse(S), REGISTERED_VARIANTS) == length(S)

    # a level deleted from a tuple: the run enumerates thirteen and says so
    gone = refusal(SG.assert_the_set_is_the_declared_one,
                   filter(!=("B:none"), S), REGISTERED_VARIANTS)
    @test occursin("B:none", gone)
    @test occursin("missing from this run", gone)

    # and a rule from outside the set, substituted into one: both halves are named
    swapped = replace(copy(S), "B:ipr017512" => "B:lnca")
    both = refusal(SG.assert_the_set_is_the_declared_one, swapped, REGISTERED_VARIANTS)
    @test occursin("missing from this run", both) && occursin("B:ipr017512", both)
    @test occursin("absent from the declared set", both) && occursin("B:lnca", both)

    # A LEVEL ADDED is refused as well as one removed. A set that only refused shrinkage would
    # let a disclosed line — a rule the recorded decision deliberately keeps outside the set —
    # join it and decide a published breaker list.
    added = refusal(SG.assert_the_set_is_the_declared_one,
                    vcat(S, ["B:lntype"]), REGISTERED_VARIANTS)
    @test occursin("absent from the declared set", added) && occursin("B:lntype", added)
    @test !occursin("missing from this run", added)
end

@testset "the defining level carries every field a specification is made of" begin
    p = SG.primary_spec()
    @test length(SG.defining_level(p)) == length(SG.DEFINING_FIELDS)
    # The two fields a four-field tuple left out are exactly the two that let a variant be
    # measured on another population, or off another column, while comparing equal.
    @test :chem_col in SG.DEFINING_FIELDS
    @test :subset in SG.DEFINING_FIELDS
    @test SG.defining_level(respec(p; chem_col = HE.ANCHOR_CHEM_COL)) != SG.defining_level(p)
    @test SG.defining_level(respec(p; subset = :min_distinct)) != SG.defining_level(p)
    # "the level the call table already carries" is the empty string, which no named rule can
    # collide with, so a variant that varies a rule is never mistaken for one that does not
    @test SG.field_level(p, :tie_rule) == ""
    @test SG.field_level(SPECS[2], :tie_rule) == SPECS[2].level
end

@testset "the grid runs every registered variant and writes each row once" begin
    n = SG.assert_grid_complete(GRID.rows, S)
    @test n == length(GRID.rows)                  # no (variant, contrast, unit, category) twice
    for v in S
        @test any(r -> r.variant == v && r.feeds_flag, GRID.rows)
    end
end

@testset "a variant that moves nothing anywhere is refused, by name" begin
    mv = SG.estimate_movement(GRID.rows, S)
    @test SG.assert_every_variant_moves_something(mv, S) == length(S)
    @test all(v -> mv.per_variant[v] > 0, S)

    # An axis that cannot reach ONE unit is expected and must stay legal: the chemistry and
    # type rule families are orthogonal, so each leaves the other's unit exactly where it was.
    @test any(k -> mv.per_unit[k] == 0, keys(mv.per_unit))

    # A variant whose rows are the primary's on every unit is not. Built by overwriting one
    # variant's flag-feeding rows with the primary's numbers, which is what a re-run of the
    # registered specification under that name produces.
    victim = "B:none"
    primary = Dict((r.unit, r.category) => r for r in GRID.rows
                   if r.variant == "primary" && r.feeds_flag)
    doctored = [r.variant == victim && r.feeds_flag &&
                haskey(primary, (r.unit, r.category)) ?
                merge(primary[(r.unit, r.category)], (variant = victim,)) : r
                for r in GRID.rows]
    mv2 = SG.estimate_movement(doctored, S)
    @test mv2.per_variant[victim] == 0
    err = try; SG.assert_every_variant_moves_something(mv2, S); ""
          catch e; sprint(showerror, e); end
    @test occursin(victim, err)
end

@testset "the units a variant moved are the units its axis declares" begin
    mv = SG.estimate_movement(GRID.rows, S)
    @test SG.assert_axes_reach_the_units_they_declare(mv, SPECS) == length(S)

    # The declaration is not a restatement of the measurement: it says, per axis, which units
    # that family of rules can move at all, and the fixture reproduces the orthogonality the
    # real grid has — a tie rule moves no chemistry flag, a chemistry rule moves no habitat
    # call, and the three axes that change which taxa or which terms count move both.
    for s in SPECS
        s.role == SG.ROLE_VARIANT || continue
        a = SG.axis_declaration(s.axis)
        for u in mv.units
            @test (get(mv.per_unit, (s.name, u), 0) > 0) == (u in a.units)
        end
    end
    @test any(a -> length(a.units) == 1, SG.VARIANT_AXES)   # the orthogonality is real ...
    @test any(a -> length(a.units) > 1, SG.VARIANT_AXES)    # ... and is not total

    # A LEVEL THAT LOST ITS WIRING TO ONE UNIT. It still moves the other, so the "moves
    # nothing anywhere" guard stays silent — and the written table, which takes `n.r.` from
    # the same measurement, would print the loss as the sentence "that rule family does not
    # touch this unit at all". The declaration is what turns it into a contradiction.
    victim = "B:none"; unit = SG.chemistry_unit("taxa")
    lost_unit = copy(mv.per_unit); lost_variant = copy(mv.per_variant)
    lost_variant[victim] -= lost_unit[(victim, unit)]
    lost_unit[(victim, unit)] = 0
    lost = (per_cell = mv.per_cell, per_unit = lost_unit, per_variant = lost_variant,
            units = mv.units)
    @test lost.per_variant[victim] > 0
    @test SG.assert_every_variant_moves_something(lost, S) == length(S)     # silent
    e1 = refusal(SG.assert_axes_reach_the_units_they_declare, lost, SPECS)
    @test occursin(victim, e1) && occursin(unit, e1)
    @test occursin("declared to reach", e1)

    # ... and the other direction: an axis that moved a unit it cannot reach is varying
    # something other than the rule its name names.
    extra_unit = copy(mv.per_unit); extra_variant = copy(mv.per_variant)
    extra_unit[("T:alphabetical", unit)] = 1; extra_variant["T:alphabetical"] += 1
    extra = (per_cell = mv.per_cell, per_unit = extra_unit, per_variant = extra_variant,
             units = mv.units)
    e2 = refusal(SG.assert_axes_reach_the_units_they_declare, extra, SPECS)
    @test occursin("T:alphabetical", e2) && occursin("declared unable to reach", e2)
end

@testset "how far each variant moved is recorded, not only that it moved" begin
    mv = SG.estimate_movement(GRID.rows, S)
    byname = Dict(s.name => s for s in SPECS if s.role == SG.ROLE_VARIANT)
    for v in S
        @test mv.per_variant[v] == sum(get(mv.per_unit, (v, u), 0) for u in mv.units)
        # the denominator is a property of the axis and not of the numbers: how many cells
        # that level could have moved at all
        @test SG.declared_reachable_cells(byname[v], mv) >= mv.per_variant[v]
    end

    # THE ONE-CELL CLONE. One variant is reproduced on every flag cell but one — a
    # half-applied configuration, not a rule variant. It moves something, so
    # `assert_every_variant_moves_something` is silent, and its interval is the registered
    # specification's everywhere it matters. What makes it visible is the count beside its
    # sibling's, which is why the count is recorded per variant and not only per cell.
    victim = "B:none"; sibling = "B:hmmsearch_movein"
    prim = Dict((r.unit, r.category) => r for r in GRID.rows
                if r.variant == "primary" && r.feeds_flag)
    moved = sort([(String(r.unit), String(r.category)) for r in GRID.rows
                  if r.variant == victim && r.feeds_flag &&
                     haskey(prim, (r.unit, r.category)) &&
                     !SG.rows_agree(prim[(r.unit, r.category)], r)])
    @test length(moved) > 1
    kept = first(moved)
    doctored = [r.variant == victim && r.feeds_flag &&
                (String(r.unit), String(r.category)) != kept &&
                haskey(prim, (r.unit, r.category)) ?
                merge(prim[(r.unit, r.category)], (variant = victim,)) : r for r in GRID.rows]
    mv2 = SG.estimate_movement(doctored, S)
    @test mv2.per_variant[victim] == 1
    @test SG.assert_every_variant_moves_something(mv2, S) == length(S)      # silent
    @test mv2.per_variant[sibling] > mv2.per_variant[victim]                # and visible
    @test SG.declared_reachable_cells(byname[victim], mv2) > mv2.per_variant[victim]
    # A clone that loses a whole unit is caught outright by the declared reach; one that
    # keeps a cell on every unit is not, and for that one the recorded pair is the evidence.
    @test occursin(victim, refusal(SG.assert_axes_reach_the_units_they_declare, mv2, SPECS))
end

@testset "two levels that are one measurement are refused, and a near-clone is not" begin
    # THE DUPLICATE THAT NEEDS NO EDIT TO THE SCRIPT. The clone check above compares
    # DECLARATIONS, and a configuration of the call table written twice — a re-export that
    # reused a frame, a label assigned in the wrong loop iteration — leaves every declaration
    # distinct while making two levels the same measurement. The set has then lost a level
    # while keeping its column, and a category both of them broke keeps one breaker where a
    # reader counts two.
    @test SG.assert_no_two_variants_are_one_measurement(GRID.rows, S) == length(S)

    a, b = "X:text_first", "X:text_object"
    src = Dict((r.unit, r.category) => r for r in GRID.rows if r.variant == a && r.feeds_flag)
    dup = [r.variant == b && r.feeds_flag && haskey(src, (r.unit, r.category)) ?
           merge(src[(r.unit, r.category)], (variant = b,)) : r for r in GRID.rows]
    # every declaration still differs, and every other guard is silent
    @test SG.assert_variants_differ_from_each_other(SPECS) == length(S)
    dmv = SG.estimate_movement(dup, S)
    @test SG.assert_every_variant_moves_something(dmv, S) == length(S)
    @test SG.assert_axes_reach_the_units_they_declare(dmv, SPECS) == length(S)
    err = refusal(SG.assert_no_two_variants_are_one_measurement, dup, S)
    @test occursin(a, err) && occursin(b, err)          # both named, not just one
    @test occursin("one measurement", err)
    # and the error says what separates the two cases, so a near-clone is not read as this
    @test occursin("near-clone", err)

    # THE NEAR-CLONE STILL PASSES. One cell apart is a legitimate difference — it is reported
    # as fewer cells moved, not refused — so the guard is a test for NO difference and not for
    # a small one.
    cells = sort([(String(r.unit), String(r.category)) for r in GRID.rows
                  if r.variant == b && r.feeds_flag && haskey(src, (r.unit, r.category)) &&
                     !SG.rows_agree(src[(r.unit, r.category)], r)])
    @test !isempty(cells)
    keep = first(cells)
    near = [r.variant == b && r.feeds_flag &&
            (String(r.unit), String(r.category)) != keep &&
            haskey(src, (r.unit, r.category)) ?
            merge(src[(r.unit, r.category)], (variant = b,)) : r for r in GRID.rows]
    @test SG.assert_no_two_variants_are_one_measurement(near, S) == length(S)

    # A cell one of them wrote a row for and the other did not counts as a DIFFERENCE, for the
    # same reason it counts as movement: the two tables say different things about it. So one
    # row removed from the duplicate separates the pair again, and a guard that skipped the
    # cells where only one side has a row would still call these two one measurement.
    thin = [r for r in dup if !(r.variant == b && r.feeds_flag &&
                                (String(r.unit), String(r.category)) == keep)]
    @test SG.assert_no_two_variants_are_one_measurement(thin, S) == length(S)
end

@testset "how close the closest pair came is reported per unit, not only overall" begin
    # THE DENOMINATOR DOES LESS WORK THAN IT LOOKS. The guard fires at zero cells apart, so a
    # margin read against every cell the flag reads looks like far more headroom than it is.
    # Two levels of one rule family agree on every cell of the
    # other family's unit BY CONSTRUCTION, because their axis is declared unable to reach it.
    # The cells that could have separated them are the ones both axes reach.
    d = SG.variant_pair_differences(GRID.rows, S)
    @test length(d.pairs) == binomial(length(S), 2)
    @test length(d.units) > 1
    # the guard and the margin read ONE comparison, so the number the report prints is the
    # number the refusal was made against
    @test SG.assert_no_two_variants_are_one_measurement(GRID.rows, S, d) == length(S)

    m = SG.variant_pair_margin(d, SPECS)
    @test m.overall !== nothing
    @test m.overall.differing > 0                      # no pair is one measurement
    @test m.overall.compared > 0
    @test length(m.per_unit) == length(d.units)
    tighter = false
    for (u, x) in m.per_unit
        @test x.differing > 0
        @test x.compared <= m.overall.compared
        x.compared < m.overall.compared && (tighter = true)
        # the pair a unit's margin names is a pair BOTH of whose axes reach that unit — a pair
        # that agrees there because neither axis can touch it says nothing about either level
        for v in split(x.pair, " vs ")
            s = first([s for s in SPECS if s.name == String(v)])
            @test u in SG.axis_declaration(s.axis).units
        end
    end
    @test tighter                                      # the correction is not cosmetic

    # and the per-unit margin is a margin: a pair made identical on one unit it can reach is
    # reported there as zero cells apart, which is what the guard refuses on the whole grid
    a, b = "X:text_first", "X:text_object"
    src = Dict((r.unit, r.category) => r for r in GRID.rows if r.variant == a && r.feeds_flag)
    dup = [r.variant == b && r.feeds_flag && haskey(src, (r.unit, r.category)) ?
           merge(src[(r.unit, r.category)], (variant = b,)) : r for r in GRID.rows]
    md = SG.variant_pair_margin(SG.variant_pair_differences(dup, S), SPECS)
    @test md.overall.differing == 0
    @test any(x -> last(x).differing == 0, md.per_unit)
end

@testset "a category the library does not declare is refused" begin
    # A category is the name of an ontology branch or of a condition flag, not whatever the
    # column being read happened to hold. The fixture names its own vocabulary, because its
    # categories are role names; the run's default vocabulary is the library's own.
    stem_is_chem(u) = SG.is_chemistry_unit(SG.unit_stem(String(u)))
    types = Set(String[String(r.category) for r in GRID.rows if !stem_is_chem(r.unit)])
    chems = Set(String[String(r.category) for r in GRID.rows if stem_is_chem(r.unit)])
    @test !isempty(types) && !isempty(chems)
    @test isempty(intersect(types, chems))
    @test SG.assert_categories_are_registered(GRID.rows; type_categories = types,
                                              chemistry_categories = chems) ==
          length(GRID.rows)

    # A COMBINATION OF TWO FLAGS IS NOT A FLAG. The arm builder splits combinations before it
    # counts them, so a combination reaching a written row means the column being read is not
    # the column that was asked for — and the flag table would carry it as a category.
    combo = join(sort(collect(chems)), '+')
    @test !(combo in chems)
    spoiled = [stem_is_chem(r.unit) ? merge(r, (category = combo,)) : r for r in GRID.rows]
    err = refusal(SG.assert_categories_are_registered, spoiled;
                  type_categories = types, chemistry_categories = chems)
    @test occursin(combo, err)
    # and a label from the habitat vocabulary is not a condition flag either
    wrong = [stem_is_chem(r.unit) ? merge(r, (category = first(sort(collect(types))),)) : r
             for r in GRID.rows]
    @test occursin(first(sort(collect(types))),
                   refusal(SG.assert_categories_are_registered, wrong;
                           type_categories = types, chemistry_categories = chems))

    # The vocabulary the RUN uses is read off the library that assigns the categories, so it
    # cannot drift from it, and no combination is in it.
    @test !isempty(SG.REGISTERED_TYPE_CATEGORIES)
    @test !isempty(SG.REGISTERED_CHEMISTRY_CATEGORIES)
    @test all(c -> !occursin('+', c), SG.REGISTERED_CHEMISTRY_CATEGORIES)
    @test all(c -> !occursin('+', c), SG.REGISTERED_TYPE_CATEGORIES)
    @test !(join(sort(collect(SG.REGISTERED_CHEMISTRY_CATEGORIES))[1:2], '+') in
            SG.REGISTERED_CHEMISTRY_CATEGORIES)

    # The subset suffix must not make a chemistry unit look like a habitat one: the two
    # vocabularies are different, and a row sorted into the wrong one is checked against the
    # wrong list.
    @test stem_is_chem(SG.chemistry_unit("taxa") * SG.MIN_DISTINCT_SUFFIX)
    @test !stem_is_chem("taxa" * SG.MIN_DISTINCT_SUFFIX)
end

@testset "a table short of a contrast, a unit or a category is refused" begin
    # THE OTHER DIRECTION FROM EVERY GUARD ABOVE. Those ask whether a row that is here is
    # right; this asks whether a row that should be here is missing. A table short of a whole
    # contrast carries no wrong number anywhere, meets every anchor target, reproduces every
    # reference row and exits 0 — and what a downstream reader sees is a complete-looking
    # file, because "complete" is not a property any column of it states.
    stem_is_chem(u) = SG.is_chemistry_unit(SG.unit_stem(String(u)))
    types = Set(String[String(r.category) for r in GRID.rows if !stem_is_chem(r.unit)])
    chems = Set(String[String(r.category) for r in GRID.rows if stem_is_chem(r.unit)])
    shape(rows; kw...) = SG.assert_the_published_shape_is_the_declared_one(
        rows, SPECS; type_categories = types, chemistry_categories = chems, kw...)
    # The fixture context carries no Ca-class membership, which is exactly the state a run
    # whose census did not arrive is in, so the fixture's own declaration is the three
    # registered contrasts and the complete fixture grid satisfies it.
    fixture_contrasts = HE.REGISTERED_CONTRASTS
    @test shape(GRID.rows; contrasts = fixture_contrasts) == length(SPECS)

    # AND THE SAME GRID AGAINST THE DECLARATION THE RUN USES is a table short of the disclosed
    # comparator — one contrast gone from every specification, with every number that remains
    # correct. This is the shape a missing or stale census produces.
    @test HE.CA_VS_BG in SG.DECLARED_CONTRASTS
    @test all(c -> c in SG.DECLARED_CONTRASTS, HE.REGISTERED_CONTRASTS)
    short = refusal(shape, GRID.rows; contrasts = SG.DECLARED_CONTRASTS)
    @test occursin(HE.CA_VS_BG, short)
    @test occursin("no row at all", short)

    # A CONTRAST GONE FROM ONE SPECIFICATION is a specification that did not run, and it is
    # named rather than counted: the run-level check above cannot see it, because every other
    # specification still publishes that contrast.
    victim = "T:fractional"
    lost = [r for r in GRID.rows
            if !(r.variant == victim && r.contrast == HE.ALL_VS_BG &&
                 !endswith(String(r.unit), SG.MIN_DISTINCT_SUFFIX))]
    @test length(lost) < length(GRID.rows)
    one_short = refusal(shape, lost; contrasts = fixture_contrasts)
    @test occursin(victim, one_short) && occursin(HE.ALL_VS_BG, one_short)
    # and the subset twin of the same level, which is a different run under the same name, is
    # not what answered for it: it still carries the contrast and is not named
    @test !occursin(victim * SG.MIN_DISTINCT_SUFFIX, one_short)

    # A UNIT GONE takes a whole series out of the written table and out of the panel.
    no_genera = [r for r in GRID.rows if SG.unit_stem(String(r.unit)) != "genera"]
    @test length(no_genera) < length(GRID.rows)
    gone_unit = refusal(shape, no_genera; contrasts = fixture_contrasts)
    @test occursin("genera", gone_unit) && !occursin("genera (chemistry)", gone_unit)
    # and the chemistry twin is its own series, refused under its own name
    no_chem_genera = [r for r in GRID.rows
                      if SG.unit_stem(String(r.unit)) != SG.chemistry_unit("genera")]
    @test occursin(SG.chemistry_unit("genera"),
                   refusal(shape, no_chem_genera; contrasts = fixture_contrasts))

    # A CATEGORY GONE left no cell in the table and no line in the panel, and nothing in
    # either says the name was ever expected. This is the vocabulary check read from the end
    # `assert_categories_are_registered` does not cover.
    dropped = first(sort(collect(types)))
    no_cat = [r for r in GRID.rows
              if !(!stem_is_chem(r.unit) && String(r.category) == dropped)]
    @test length(no_cat) < length(GRID.rows)
    gone_cat = refusal(shape, no_cat; contrasts = fixture_contrasts)
    @test occursin(dropped, gone_cat) && occursin("published on no row", gone_cat)

    # The declaration the RUN uses is read off the library and off this file's own unit
    # helper, so a rename there moves it rather than leaving it describing something else.
    @test all(u -> u in SG.DECLARED_UNIT_SERIES, HE.PRIMARY_UNITS[1:1])
    @test length(SG.DECLARED_UNIT_SERIES) == length(unique(SG.DECLARED_UNIT_SERIES))
    @test !any(SG.is_chemistry_unit, SG.DECLARED_UNIT_SERIES)
end

@testset "the census is required, and an empty digest is not a digest" begin
    # THE ONE INPUT OF THE THREE WHOSE DROPPED PATH MUST ALSO STOP THE RUN: it is the input
    # that decides how much of the grid gets WRITTEN: the disclosed comparator arm is defined
    # by the census, and an arm with no members is left out of the arm builder rather than
    # emitted as a panel of zeros. That is the right decision for an arm and the wrong one
    # for a table.
    @test !SG.is_sha256("")                          # an absent file's empty digest field
    @test !SG.is_sha256(repeat("a", 63))
    @test !SG.is_sha256(repeat("z", 64))             # 64 characters, not 64 hex digits
    @test SG.is_sha256(repeat("ab", 32))
    mktempdir() do dir
        gone = joinpath(dir, "no_census_here.csv")
        # the library's reader is what makes this silent: it hands back no membership at all
        # for a file it cannot find, and the arm builder cannot tell that from a file that
        # named nobody
        @test HE.read_census_classes(gone).ca_class_taxids === nothing
        err = refusal(SG.require_census, gone)
        @test occursin("missing", err) && occursin("contrast", err)

        # A FILE THAT NAMES NOBODY is the same failure with bytes on disk, and is refused too:
        # a stale or half-written census is the ordinary way this happens.
        nobody = joinpath(dir, "nobody.csv")
        write(nobody, "classification,kingdom,TaxID\nLanthanide,Bacteria,11\n")
        @test isempty(HE.read_census_classes(nobody).ca_class_taxids)
        @test occursin("no members", refusal(SG.require_census, nobody))

        # and a census that names somebody is read, with a digest that is a digest
        good = joinpath(dir, "census.csv")
        write(good, "classification,kingdom,TaxID\nCalcium,Bacteria,11\n" *
                    "Lanthanide,Bacteria,12\n")
        c = SG.require_census(good)
        @test SG.is_sha256(c.sha256)
        @test c.census.ca_class_taxids == Set([11])   # calcium-class and never lanthanide
    end
end

@testset "the movement counts are what the flag table reports" begin
    mv = SG.estimate_movement(GRID.rows, S)
    flags = SG.robustness(GRID.rows, S; movement = mv)
    for f in flags
        @test f.n_variants_that_moved_the_estimate == mv.per_cell[(f.unit, f.category)]
        # it can never exceed the set, and a row nothing moved is possible but must be honest
        @test 0 <= f.n_variants_that_moved_the_estimate <= f.variants_expected
    end
    # The count is strictly below the size of the set wherever an axis cannot reach the unit,
    # which is the whole reason the column exists.
    @test any(f -> f.n_variants_that_moved_the_estimate < f.variants_expected, flags)
end

@testset "a missing variant is refused, and so is a duplicated row" begin
    @test_throws ErrorException SG.assert_grid_complete(GRID.rows, vcat(S, ["v:absent"]))
    @test_throws ErrorException SG.assert_grid_complete(vcat(GRID.rows, [GRID.rows[1]]), S)
end

@testset "a specification that computed nothing is an error" begin
    @test_throws ErrorException SG.spec_rows(SG.primary_spec(), (rows = NamedTuple[],),
                                             nothing)
end

@testset "only the registered rows feed the flag" begin
    for r in GRID.rows
        if r.feeds_flag
            @test r.contrast == HE.LN_VS_BG
            @test r.unit in HE.PRIMARY_UNITS
            @test r.role in (SG.ROLE_PRIMARY, SG.ROLE_VARIANT)
        end
    end
    # The unit checks, the anchors and the disclosed lines are all present AND all silent.
    for role in (SG.ROLE_UNIT, SG.ROLE_ANCHOR, SG.ROLE_AUDIT, SG.ROLE_DISCLOSED)
        rs = [r for r in GRID.rows if r.role == role]
        @test !isempty(rs)
        @test all(r -> !r.feeds_flag, rs)
    end
end

@testset "the genus series is complete over every graded specification" begin
    # A downstream decision rests on exactly this series, so what has to hold is a property of
    # every CELL — each graded specification measures each genus category the registered
    # specification measures — and not of the number of specifications with a genus row.
    graded = vcat(["primary"], S)
    g = SG.check_genus_series(GRID.rows, graded)
    @test g.specifications == length(graded)
    @test "genera" in g.units
    @test g.cells >= g.specifications        # at least one category per specification
    # every name it reports is a cell of the grid: a graded specification and a genus unit
    @test all(m -> any(v -> startswith(m, v * " | "), graded), g.missing_cells)
    @test all(m -> occursin(" | genera", m), g.missing_cells)

    # A cell CAN be legitimately absent — a tie rule that stops calling a category for any
    # taxon has nothing to measure there — which is why the check reports rather than
    # refuses. What it must not do is miss one.
    victim = last(sort(unique([r.category for r in GRID.rows
                               if r.variant == "primary" && r.unit == "genera" &&
                                  r.contrast == HE.LN_VS_BG])))
    @test !("B:none | genera | $victim" in g.missing_cells)      # present before
    thinned = [r for r in GRID.rows
               if !(r.variant == "B:none" && r.unit == "genera" && r.category == victim &&
                    r.contrast == HE.LN_VS_BG)]
    h = SG.check_genus_series(thinned, graded)
    @test length(h.missing_cells) == length(g.missing_cells) + 1
    @test "B:none | genera | $victim" in h.missing_cells

    # A count of variants carrying a genus row is unmoved by that deletion, which is why
    # this check counts cells instead: such a count reads the same whether the series is
    # whole or a cell short.
    countof(rs) = length(unique([r.variant for r in rs
                                 if r.unit == "genera" && r.contrast == HE.LN_VS_BG]))
    @test countof(thinned) == countof(GRID.rows)
end

@testset "a specification with an empty arm fails naming the arm" begin
    # A specification with nothing on one side writes undefined intervals everywhere, and the
    # complaint arrives one whole grid later as a thin-evidence subset that "left one of the
    # arms unchanged" — true, and about the wrong thing. A specification with no subset twin
    # (an anchor, a disclosed line) would not reach even that.
    calls = fixture_calls()
    calls.n_ln_seqs = zeros(Int, nrow(calls))          # nobody is a census carrier any more
    err = try; SG.build_arms(fixture_context(; calls = calls), SG.primary_spec()); ""
          catch e; sprint(showerror, e); end
    @test occursin("FOREGROUND", err)
    @test occursin("primary", err)

    # and the other arm, named as the other arm: every frame taxon excluded from it
    base = fixture_context()
    excl = Dict{String,Set{Int}}(k => v for (k, v) in base.exclusions)
    excl[HE.PRIMARY_MODE] = Set{Int}(fixture_calls().taxid)
    err2 = try; SG.build_arms(merge(base, (exclusions = excl,)), SG.primary_spec()); ""
           catch e; sprint(showerror, e); end
    @test occursin("BACKGROUND", err2)
    @test occursin("primary", err2)
end

# ---------------------------------------------------------------------------
# 3. The thin-evidence subset.
# ---------------------------------------------------------------------------

@testset "the thin-evidence subset filters BOTH arms" begin
    whole = GRID.runs[("primary", "")].type_arms
    sub = GRID.runs[("primary", SG.MIN_DISTINCT_SUFFIX)].type_arms
    SG.check_both_arms_filtered(sub, whole, "primary")     # the guard the run applies
    @test sub.n1 < whole.n1
    @test sub.n2 < whole.n2
    @test all(>=(SG.MIN_DISTINCT_TERMS), sub.fg.n_distinct)
    @test all(>=(SG.MIN_DISTINCT_TERMS), sub.bgm.n_distinct)

    # Its rows are marked by their unit, so a subset row cannot be lined up against a
    # whole-arm row of the same specification and category.
    subrows = [r for r in GRID.rows
               if r.variant == "primary" && endswith(r.unit, SG.MIN_DISTINCT_SUFFIX)]
    @test !isempty(subrows)
    @test all(r -> !r.feeds_flag, subrows)
    @test all(r -> r.role == SG.ROLE_UNIT, subrows)
end

@testset "a subset that filtered one arm only is refused" begin
    # The foreground filtered and the background whole: exactly the shape that moves every
    # estimate in one direction for a reason that is not evidence.
    calls = fixture_calls()
    calls.n_distinct = [r.in_census ? r.n_distinct : SG.MIN_DISTINCT_TERMS
                        for r in eachrow(calls)]
    ctx = fixture_context(; calls = calls)
    whole = SG.build_arms(ctx, SG.primary_spec())
    thin = SG.build_arms(ctx, SG.unit_check_specs([SG.primary_spec()])[1])
    @test thin.n1 < whole.n1
    @test thin.n2 == whole.n2
    @test_throws ErrorException SG.check_both_arms_filtered(thin, whole, "primary")
end

# ---------------------------------------------------------------------------
# 4. The two rule seams.
# ---------------------------------------------------------------------------

@testset "the type column is rebuilt from the table's own tallies" begin
    prim = fixture_calls()
    prim = prim[prim.config .== "primary", :]
    ok = SG.type_reconstruction_check(prim)
    @test ok.mismatches == 0
    @test ok.rows == nrow(prim)

    # The same check on a table whose weight column was written under ANOTHER rule: the
    # licence for every type variant has to fail here rather than pass quietly.
    wrong = copy(prim)
    wrong.type_weights = SG.type_weight_column(prim, :alphabetical)
    bad = SG.type_reconstruction_check(wrong)
    @test bad.mismatches > 0
    @test !isempty(bad.example)
end

@testset "both seams are checked over every configuration the grid reads" begin
    # THE SEAM CHECKS MUST NOT RUN ON THE PRIMARY ALONE. The type variants are recomputed
    # from the primary's tallies, so checking that one licenses the rewrite — but that is not
    # the whole of it. Four of the other configurations ARE registered levels: the text,
    # generic-term and isolation variants are their own committed weight columns read under a
    # level's name, and a primary-only check recomputes none of those from their own tallies.
    # A column re-exported from a sibling keeps that configuration's tallies, its counts and
    # its row count, so every other guard here stays silent.
    calls = fixture_calls()
    tconfigs = SG.registered_type_rule_configurations(SPECS)
    aconfigs = SG.grid_configurations(SPECS)
    @test length(tconfigs) > 1
    @test HE.ANCHOR_CONFIG in aconfigs                 # the grid does read it ...
    @test !(HE.ANCHOR_CONFIG in tconfigs)              # ... and its column is not the
                                                       #     registered rule's, on purpose
    t = SG.type_reconstruction_over(calls, tconfigs)
    @test t.mismatches == 0
    @test length(t.per_configuration) == length(tconfigs)
    # THE SCOPE IS RECORDED, not only the verdict: a check whose scope is not written down is
    # a check whose scope can shrink without the report changing.
    @test sum(last.(t.per_configuration)) == t.rows
    @test t.rows > nrow(calls[calls.config .== SG.PRIMARY_CONFIG, :])

    # A HALF-APPLIED RE-EXPORT OF ONE NON-PRIMARY CONFIGURATION. Its weight column becomes the
    # primary's for part of the table — the shape a regeneration that covered one arm and not
    # the other leaves behind — while its own tallies, counts and row count stay put.
    victim = "generic_none"
    @test victim in tconfigs
    prim = calls[calls.config .== SG.PRIMARY_CONFIG, :]
    spoiled = copy(calls)
    taxids = Dict(prim.taxid[i] => prim.type_weights[i] for i in 1:nrow(prim))
    half = [r.config == victim && haskey(taxids, r.taxid) && r.taxid < 105 ?
            taxids[r.taxid] : r.type_weights for r in eachrow(spoiled)]
    spoiled.type_weights = half
    @test any(spoiled.type_weights .!= calls.type_weights)
    # the primary is untouched, so the check as it was is silent ...
    @test SG.type_reconstruction_check(
        spoiled[spoiled.config .== SG.PRIMARY_CONFIG, :]).mismatches == 0
    # ... and the widened one names the configuration the substitution is in
    err = refusal(SG.type_reconstruction_over, spoiled, tconfigs)
    @test isempty(err)                                 # it returns a count, it does not throw
    @test SG.type_reconstruction_over(spoiled, tconfigs).mismatches > 0
    @test occursin(victim, SG.type_reconstruction_over(spoiled, tconfigs).example)

    # The combination tally is checked over every configuration INCLUDING the earlier
    # pipeline's, whose chemistry columns the audit rules are recomputed from, and a failure
    # names the configuration because a taxon id alone does not say which table holds it.
    c = SG.combinations_over(calls, aconfigs)
    @test length(c.per_configuration) == length(aconfigs)
    @test sum(last.(c.per_configuration)) == c.rows == nrow(calls)
    @test c.checked > 0
    broken = copy(calls)
    broken.n_chem_terms = [r.config == HE.ANCHOR_CONFIG && r.n_chem_terms > 0 ?
                           r.n_chem_terms + 1 : r.n_chem_terms for r in eachrow(broken)]
    e = refusal(SG.combinations_over, broken, aconfigs)
    @test occursin(HE.ANCHOR_CONFIG, e)

    # and a configuration the grid reads that the table does not hold is an error, not a check
    # made over nothing: "0 of 0 rows disagree" is the sentence a passing check writes too
    @test occursin("no rows for configuration",
                   refusal(SG.configuration_rows, calls, "a-configuration-nobody-exported"))
end

@testset "a type variant changes the weights and not the arms" begin
    base = GRID.runs[("primary", "")].type_arms
    for r in SG.TYPE_LEVELS
        v = GRID.runs[("T:$r", "")].type_arms
        @test v.n1 == base.n1
        @test v.n2 == base.n2
    end
    # `unresolved` declines to call a surviving tie, and an uncalled taxon stays one taxon of
    # its arm: it scores zero everywhere instead of leaving the denominator.
    prim = fixture_calls()
    prim = prim[prim.config .== "primary", :]
    @test any(isempty, SG.type_weight_column(prim, :unresolved))
end

@testset "the chemistry weights come from the counts, not from the rounded shares" begin
    prim = fixture_calls()
    prim = prim[prim.config .== "primary", :]
    v = SG.validate_combinations!(prim)
    @test v.rows == nrow(prim)
    @test v.checked > 0                              # some taxa do carry a chemistry term

    # A share map that does not account for the terms it was normalised by describes other
    # evidence than the counts beside it, and a rule computed from it is a rule applied to
    # the wrong tally.
    spoiled = copy(prim)
    spoiled.n_chem_terms = [x == 0 ? 0 : x + 1 for x in spoiled.n_chem_terms]
    @test_throws ErrorException SG.validate_combinations!(spoiled)

    # The other half of the same claim, and the one the first case cannot reach: a tally that
    # sums correctly while giving a flag a different mass from the one the table counted. A
    # rule keyed on a flag's share would then run on a flag the table does not have.
    drop_last(s) = (p = split(s, ';'); length(p) < 2 ? String(s) : join(p[1:end - 1], ";"))
    lost = copy(prim)
    lost.chemistry_flag_counts = drop_last.(lost.chemistry_flag_counts)
    @test any(lost.chemistry_flag_counts .!= prim.chemistry_flag_counts)
    @test_throws ErrorException SG.validate_combinations!(lost)
end

@testset "the registered chemistry rule agrees through both routes" begin
    # The weight column is the seam every chemistry variant crosses. Under the registered
    # rule it must reproduce the call table's own chemistry column cell for cell, or the
    # variants computed through it are different rules under the names that were asked for.
    base = GRID.runs[("primary", "")]
    check = SG.Spec("check", "check", "any_term", SG.PRIMARY_CONFIG, nothing,
                    :chemistry_flags, :any_term, HE.PRIMARY_MODE, :all, SG.ROLE_DISCLOSED,
                    :chemistry_only, "", "the registered chemistry rule via the weight column")
    derived = SG.spec_rows(check, base.type_arms, base.chem_arms)
    c = SG.chemistry_route_check([r for r in GRID.rows if r.variant == "primary"], derived)
    @test c.compared > 0
    @test c.worst <= SG.CHEMISTRY_ROUTE_TOL
end

@testset "a weight string survives the round trip at full precision" begin
    d = Dict("cat-a" => 1 / 3, "cat-b" => 2 / 3)
    back = SG.parse_tally(SG.weight_string(d))
    @test back["cat-a"] === 1 / 3
    @test back["cat-b"] === 2 / 3
    # The published formatter rounds for reading; this grid compares hit totals, so it must
    # not be the formatter the grid writes with.
    @test SG.weight_string(d) != ET._weights_str(d)
    @test isempty(SG.parse_tally(""))
    @test_throws ErrorException SG.parse_tally("cat-a")
end

# ---------------------------------------------------------------------------
# 5. The move-in level.
# ---------------------------------------------------------------------------

@testset "the move-in level grew the foreground by the taxa it moved" begin
    base = GRID.runs[("primary", "")].type_arms
    mv = GRID.runs[("B:$(HE.MOVEIN_MODE)", "")].type_arms
    g = SG.check_movein_growth(mv, base)
    @test g.moved > 0
    @test g.n1_after == g.n1_before + g.moved_annotated
    @test g.n2_after < g.n2_before                   # the same taxa left the background

    # A level that moved nobody is the unmoved specification wearing this level's name.
    @test_throws ErrorException SG.check_movein_growth(
        (n1 = base.n1, n2 = base.n2, n_moved = 0, n_moved_usable = 0), base)
    # ... and so is a level whose moved taxa never reached the arm.
    @test_throws ErrorException SG.check_movein_growth(
        (n1 = base.n1, n2 = base.n2, n_moved = 3, n_moved_usable = 3), base)
end

@testset "an empty move-in pool is refused before the arms are built" begin
    ctx = fixture_context(; movein = Set{Int}())
    spec = [s for s in SG.sensitivity_specs() if s.bg_mode == HE.MOVEIN_MODE][1]
    @test_throws ErrorException SG.build_arms(ctx, spec)
end

@testset "a background rule the scan cannot answer is an error, not an empty exclusion" begin
    scan = fixture_scan()
    @test_throws ErrorException fixture_context(; scan = select(scan, Not(:n_hits_ln_or_ca)))
end

# ---------------------------------------------------------------------------
# 6. The anchors.
# ---------------------------------------------------------------------------

@testset "the anchor levels are retired in production, and retirement is total" begin
    # anchor:all_legacy and the four audit:* levels ask a question about a census other than
    # the one a production run reads, so asked of that run they are a DIFFERENT test and a pass
    # could be misread as evidence that the two censuses give the same rows. What this pins is
    # that the retirement is COMPLETE - a half-retired grid that still emitted, say, the
    # audit:* rows would publish them under the registered set's name with nothing anchoring
    # them.
    @test SG.ANCHOR_LEVELS_RETIRED
    @test isempty(SG.anchor_specs())
    production = SG.all_specs()
    @test !any(s -> s.role in (SG.ROLE_ANCHOR, SG.ROLE_AUDIT), production)
    @test !any(s -> startswith(s.name, "anchor:") || startswith(s.name, "audit:"), production)
    # ... and that flipping the constant back is a real escape, not a dead branch.
    revived = SG.all_specs(retired = false)
    @test any(s -> s.name == "anchor:all_legacy", revived)
    @test count(s -> startswith(s.name, "audit:"), revived) == length(SG.AUDIT_CHEMISTRY_RULES)
    # The registered set is untouched by either setting: retiring an anchor may not quietly
    # shorten the grid that decides published colours.
    reg(v) = sort([s.name for s in v if s.role in (SG.ROLE_PRIMARY, SG.ROLE_VARIANT)])
    @test reg(production) == reg(revived)
end

@testset "an anchor target that is missed fails the check" begin
    row = [r for r in GRID.rows
           if r.variant == "anchor:all_legacy" && r.unit == "taxa" &&
              r.contrast == HE.LN_VS_BG][1]
    met = (anchor = row.variant, unit = row.unit, category = row.category,
           quantity = "a", value = float(row.a), tol = 0.001, note = "")

    @test SG.check_anchor_targets(GRID.rows, [met]).met == 1
    @test isempty(SG.check_anchor_targets(GRID.rows, [met]).failures)

    r = SG.check_anchor_targets(GRID.rows, [(met..., value = float(row.a) + 1.0)])
    @test r.met == 0
    @test length(r.failures) == 1

    # A target whose row this run never wrote has not been reproduced. Counting it as
    # untested would let a wiring error read as a missing input.
    a = SG.check_anchor_targets(GRID.rows,
                                [(met..., category = "cat-that-no-run-wrote")])
    @test a.met == 0
    @test length(a.failures) == 1
    @test occursin("no such row", a.failures[1])
end

@testset "every level this run publishes by name is pinned by a target" begin
    # A TARGET ON THE REGISTERED SPECIFICATION PINS ONE OF FIFTEEN INTERVALS. The word beside
    # a category is decided by that one and the fourteen levels together, and the levels that
    # reach print are the ones a breaker list names — that name is the published claim "this
    # rule, and not that one, is what the result does not survive". A wiring error confined to
    # such a level's own inputs rewrites the list while every count in the run agrees with
    # itself and every target on the primary is met.
    #
    # The requirement is read off what the run PUBLISHED rather than listed, so a level that
    # starts or stops breaking a category carries it along; a list written into the script
    # would go stale in exactly the run where it mattered.
    # The flag table is built here rather than measured, because the property under test is
    # what the requirement does with a breaker list and not which breakers this fixture
    # happens to produce.
    flags = [(unit = "taxa", category = "cat-a", robust_sign = "fragile",
              breaking_variants = "X:text_object;B:none"),
             (unit = SG.chemistry_unit("taxa"), category = "chem-x", robust_sign = "null",
              breaking_variants = "C:share"),
             (unit = "taxa", category = "cat-b", robust_sign = "+", breaking_variants = "")]
    named = SG.levels_named_in_a_breaker_list(flags)
    @test named == ["B:none", "C:share", "X:text_object"]   # sorted, de-duplicated, whole
    @test all(v -> v in S, named)
    # and every name a real run of this grid publishes is a level of the registered set, so
    # the requirement can never ask for a target on something that is not one
    @test all(v -> v in S,
              SG.levels_named_in_a_breaker_list(SG.robustness(GRID.rows, S)))

    target_on(v) = (anchor = v, unit = "taxa", category = "cat-a", quantity = "a",
                    value = 1.0, tol = 1.0, note = "")
    # with nothing anchored, every published name is unanchored ...
    @test SG.unanchored_published_levels(flags, NamedTuple[]) == named
    # ... a target on the registered specification anchors none of them ...
    @test SG.unanchored_published_levels(flags, [target_on("primary")]) == named
    # ... and a target on one level takes exactly that one off the list
    @test SG.unanchored_published_levels(flags, [target_on(first(named))]) ==
          filter(!=(first(named)), named)
    @test isempty(SG.unanchored_published_levels(flags, [target_on(v) for v in named]))
    # any cell of the level is enough: the failures this catches move a level's arms wholesale
    other = (anchor = first(named), unit = SG.chemistry_unit("taxa"),
             category = "chem-x", quantity = "n1", value = 1.0, tol = 1.0, note = "")
    @test !(first(named) in SG.unanchored_published_levels(flags, [other]))

    # a flag table with no breaker names nothing, and demands nothing
    quiet = [merge(f, (breaking_variants = "",)) for f in flags]
    @test isempty(SG.levels_named_in_a_breaker_list(quiet))
end

@testset "a targets file is read strictly, and an exception must carry its reason" begin
    mktempdir() do dir
        good = write_targets(dir; rows = [
            ("target", "anchor:all_legacy", "taxa", "cat-a", "std", 0.0, 1.0, ""),
            ("exception", "anchor:all_legacy", "taxa", "cat-b", "", 0.0, 0.0,
             "the unit was redefined on purpose")])
        t = SG.read_anchor_targets(good)
        @test length(t.targets) == 1
        @test length(t.exceptions) == 1

        for bad in (("exception", "a", "taxa", "cat-b", "", 0.0, 0.0, ""),       # no reason
                    ("something", "a", "taxa", "cat-b", "std", 0.0, 1.0, ""),    # neither kind
                    ("target", "a", "taxa", "cat-b", "pct_fg", 1.0, 1.0, ""),    # not a column
                    ("target", "a", "taxa", "cat-b", "std", 1.0, missing, ""))   # no tolerance
            @test_throws ErrorException SG.read_anchor_targets(write_targets(dir; rows = [bad]))
        end
    end
end

@testset "the targets file declares the set and its own size, or it is refused" begin
    one = [("target", "anchor:all_legacy", "taxa", "cat-a", "std", 0.0, 1.0, "")]
    mktempdir() do dir
        # THE REGISTERED SET, READ BY THE RUN. It is pinned in this suite too, and that pin
        # cannot stop a run: a suite makes its claims on a fixture and a run that publishes
        # never consults it. So the names also travel in the one input the run already refuses
        # to start without, and the two independent statements of the same fact is the point.
        t = SG.read_anchor_targets(write_targets(dir; rows = one))
        @test sort(t.variants) == REGISTERED_VARIANTS
        @test t.expected[:targets] == 1
        @test SG.assert_the_set_is_the_declared_one(S, t.variants) == length(S)

        # A FILE THAT DECLARES NO SET would let the run enumerate its own and compare it
        # against nothing, which is the circularity the declaration exists to break.
        @test_throws ErrorException SG.read_anchor_targets(
            write_targets(dir; rows = one, declare = false))
        @test occursin("no 'variant' rows", refusal(SG.read_anchor_targets,
            write_targets(dir; rows = one, names = String[], name = "noset.csv")))

        # A THINNED FILE IS A FAILURE, NOT A QUIETER SUCCESS. Without the count, a file that
        # lost rows is a run that checked fewer targets and reported it in the same words.
        @test occursin("expect", refusal(SG.read_anchor_targets,
            write_targets(dir; rows = one, with_expect = false, name = "nocount.csv")))
        @test occursin("thinned", refusal(SG.read_anchor_targets,
            write_targets(dir; rows = one, expect = 28, name = "short.csv")))

        # and the declaration itself is read strictly: no level twice, no count twice, and a
        # count that is a whole number
        @test_throws ErrorException SG.read_anchor_targets(
            write_targets(dir; rows = one, names = vcat(REGISTERED_VARIANTS, ["B:none"]),
                          name = "twice.csv"))
        for bad in (("expect", "", "", "", "targets", 2.0, missing, ""),        # counted twice
                    ("expect", "", "", "", "levels", 14.0, missing, ""),        # not a count
                    ("variant", "", "", "", "", missing, missing, ""))          # names nothing
            @test_throws ErrorException SG.read_anchor_targets(
                write_targets(dir; rows = vcat(one, [bad]), name = "bad.csv"))
        end
    end
end

@testset "the declared SIZE of the set is what a shortened pair of lists cannot keep" begin
    # THE ONE ATTACK THAT STILL CHANGED A PUBLISHED COLOUR. Everything inside the script that
    # could contradict a deleted level is enumerated from the tuples the level was deleted
    # from, and `assert_the_set_is_the_declared_one` compares two LISTS — which any two lists
    # that were shortened together satisfy. Delete the two levels that break the only fragile
    # row from the script's tuples, strike the same two names from the declaration, and the
    # run is self-consistent: it publishes a shorter grid under the registered set's name, with
    # that row coloured and its breaker list empty.
    #
    # A size is not a list. Striking a name does not change it, so the third statement
    # disagrees with the other two and the file is refused before the grid is computed.
    one = [("target", "anchor:all_legacy", "taxa", "cat-a", "std", 0.0, 1.0, "")]
    short = filter(v -> !(v in ("B:none", "X:text_object")), REGISTERED_VARIANTS)
    mktempdir() do dir
        # the two lists agree, so the name-for-name comparison is silent ...
        @test SG.assert_the_set_is_the_declared_one(short, short) == length(short)
        # ... and the file that declares them is refused all the same, because the size stayed
        err = refusal(SG.read_anchor_targets,
                      write_targets(dir; rows = one, names = short,
                                    expect_variants = length(REGISTERED_VARIANTS),
                                    name = "shortened.csv"))
        @test occursin("registered set of $(length(REGISTERED_VARIANTS))", err)
        @test occursin("names $(length(short))", err)
        # the size has to be lowered as well, which is a third edit and a bare number whose
        # only job is to be checked — the limit, stated rather than hidden
        t = SG.read_anchor_targets(write_targets(dir; rows = one, names = short,
                                                 name = "coherent.csv"))
        @test t.expected[:variants] == length(short)

        # and a file whose size is right while a name is missing is refused too, so the check
        # is on the pair and not on the direction that happens to be easier to make
        @test occursin("names", refusal(SG.read_anchor_targets,
            write_targets(dir; rows = one, names = REGISTERED_VARIANTS[1:13],
                          expect_variants = 14, name = "one_short.csv")))
    end
end

@testset "the guard's own message does not read as the edit that defeats it" begin
    # A failure that names the level the declaration holds and the run does not has written
    # half an edit: striking that name from the declaration clears the failure, makes every
    # count agree with itself again, and publishes the shorter grid. So the message has to
    # name the difference — a failure nobody can diagnose is a failure nobody fixes — and then
    # say outright which of the two repairs is the defect.
    err = refusal(SG.assert_the_set_is_the_declared_one,
                  filter(!=("B:none"), S), REGISTERED_VARIANTS)
    @test occursin("B:none", err)                       # still diagnosable
    @test occursin("not equivalent", err)               # and the two repairs are separated
    @test occursin("SMALLER grid", err)
    @test occursin("refuse", err)                       # named as the failure, not the fix
end

@testset "a target taken twice is one check counted twice, and is refused" begin
    # 28 rows that can fail, replaced by 28 more copies of rows that cannot: the file still
    # declares its number, still carries it, and still reports "N checked, N met". What moved
    # is the number of CELLS pinned, which no count in the output was taking.
    t = ("target", "anchor:all_legacy", "taxa", "cat-a", "std", 0.0, 1.0, "")
    u = ("target", "anchor:all_legacy", "taxa", "cat-b", "std", 0.0, 1.0, "")
    mktempdir() do dir
        @test length(SG.read_anchor_targets(write_targets(dir; rows = [t, u])).targets) == 2
        err = refusal(SG.read_anchor_targets,
                      write_targets(dir; rows = [t, t], name = "twice.csv"))
        @test occursin("taken twice", err)
        @test occursin("cat-a", err)
        # the count still adds up, which is exactly why counting rows cannot catch this
        @test occursin("one check counted twice", err)
        # a second target on the same row but a DIFFERENT quantity is not a duplicate: it is
        # another statement about the same cell and it can fail on its own
        @test length(SG.read_anchor_targets(
            write_targets(dir; rows = [t, (t[1:4]..., "lo", 0.0, 1.0, "")],
                          name = "other_quantity.csv")).targets) == 2
    end
end

@testset "the targets file says which bytes its numbers came from, and the run checks them" begin
    # A target transcribed from this grid's own output is byte-identical to an honestly
    # transcribed one, and the only record of which it is was a prose note that no program
    # reads and that a re-transcription does not update. The file therefore names the
    # implementation and carries its digest.
    one = [("target", "anchor:all_legacy", "taxa", "cat-a", "std", 0.0, 1.0, "")]
    mktempdir() do dir
        t = SG.read_anchor_targets(write_targets(dir; rows = one))
        @test length(t.sources) == 1
        @test t.sources[1].file == FIXTURE_SOURCE
        @test t.sources[1].sha256 == FIXTURE_DIGEST

        @test occursin("no 'source' row", refusal(SG.read_anchor_targets,
            write_targets(dir; rows = one, with_source = false, name = "nosource.csv")))
        # a digest that is not one, a source that is a path, a source declared twice, and a
        # source with no account of what it is — each refused on its own
        for (bad, why) in ((("source", "x.csv", "", "", "", missing, missing, "an impl", ""),
                            "sha256"),
                           (("source", "sub/x.csv", "", "", "", missing, missing, "an impl",
                             FIXTURE_DIGEST), "path"),
                           (("source", FIXTURE_SOURCE, "", "", "", missing, missing, "an impl",
                             FIXTURE_DIGEST), "twice"),
                           (("source", "x.csv", "", "", "", missing, missing, "",
                             FIXTURE_DIGEST), "no note"))
            e = refusal(SG.read_anchor_targets,
                        write_targets(dir; rows = vcat(one, [bad]), name = "badsource.csv"))
            @test occursin(why, e)
        end

        # A SOURCE THAT IS NOT THERE is reported, not refused: a targets file a reader writes
        # may name an implementation that does not travel with it, and a run of this script
        # must not need it, or its directory, present.
        p = write_targets(dir; rows = one, name = "absent.csv")
        s = SG.check_target_sources(p, SG.read_anchor_targets(p).sources)
        @test length(s) == 1 && s[1].ok
        @test occursin("not beside", s[1].status)

        # A SOURCE THAT IS THERE must be the bytes the file names. This is what stops the
        # attack that regenerates both halves of the wiring check from a broken run's own
        # output: the numbers still agree with each other, and the declared digest does not.
        body = "these are the bytes the targets were transcribed from\n"
        write(joinpath(dir, FIXTURE_SOURCE), body)
        real = bytes2hex(SG.sha256(Vector{UInt8}(body)))
        good = write_targets(dir; rows = one, digest = real, name = "match.csv")
        g = SG.check_target_sources(good, SG.read_anchor_targets(good).sources)
        @test g[1].ok && g[1].status == "matches" && g[1].found == real

        bad = SG.check_target_sources(p, SG.read_anchor_targets(p).sources)
        @test !bad[1].ok && bad[1].status == "differs"
        @test bad[1].declared == FIXTURE_DIGEST && bad[1].found == real
    end
end

@testset "the scope of the row-by-row comparison is declared, not discovered" begin
    # "N of N rows agree" is the same sentence at any N, so a reference file regenerated at
    # another scope reports a number that reads BETTER than the honest run's rather than
    # different from it. The count the comparison is supposed to cover is declared in the one
    # file the run already refuses to start without.
    one = [("target", "anchor:all_legacy", "taxa", "cat-a", "std", 0.0, 1.0, "")]
    mktempdir() do dir
        # an arbitrary count, like every other number in this suite: what is under test is
        # that the file states one and the run reads it back, not what the run's own scope is
        t = SG.read_anchor_targets(write_targets(dir; rows = one, reference_rows = 7))
        @test t.expected[:reference_rows] == 7
        @test :reference_rows in SG.EXPECTED_COUNTS
        # and a file that omits this count alone is refused on it by name, because a run
        # cannot discover from the reference on hand how much of the earlier run it was meant
        # to see
        @test occursin("reference_rows", refusal(SG.read_anchor_targets,
            write_targets(dir; rows = one, with_expect = (:targets, :variants),
                          name = "noref.csv")))
        @test occursin("variants", refusal(SG.read_anchor_targets,
            write_targets(dir; rows = one, with_expect = (:targets, :reference_rows),
                          name = "nosize.csv")))
    end
end

@testset "a named exception leaves the row-by-row comparison by name" begin
    ref = DataFrame(contrast = fill(HE.LN_VS_BG, 3), unit = ["taxa", "taxa", "genera"],
                    habitat = ["cat-a", "cat-b", "cat-a"],
                    ci_lo_std = [0.0, 0.0, 0.0], ci_hi_std = [1.0, 1.0, 1.0])
    ex = [(anchor = "anchor:all_legacy", unit = "taxa", category = "cat-b",
           quantity = "", value = NaN, tol = NaN, note = "one call changed on purpose")]
    s = SG.scoped_reference(ref, ex; units = ("taxa",))
    @test nrow(s.table) == 1                       # the genus row left with its unit ...
    @test s.table.habitat[1] == "cat-a"            # ... and cat-b left by name
    @test s.excluded == ["taxa | cat-b"]
end

@testset "the row-by-row comparison covers the scope it claims, and names what left" begin
    # TWO STATEMENTS, AND BOTH ARE THE RUN'S. How many rows the comparison was supposed to
    # cover is declared in the targets file and checked against what it did cover; which rows
    # left the comparison, and why, is declared here — and without the second the provenance
    # says "N of N rows agree" over a scope nobody can see, with the largest deliberate
    # difference between the two runs the one thing the output does not mention.
    @test !("genera" in SG.REFERENCE_SCOPE_UNITS)
    @test !(SG.chemistry_unit("genera") in SG.REFERENCE_SCOPE_UNITS)
    @test all(u -> u in SG.REFERENCE_SCOPE_UNITS, HE.PRIMARY_UNITS)
    # the reason travels with the list, so a later widening has to delete a stated reason
    # rather than edit a tuple
    @test occursin("genus", SG.REFERENCE_OUT_OF_SCOPE_REASON)
    @test occursin("symmetric", SG.REFERENCE_OUT_OF_SCOPE_REASON)
    @test occursin("not comparable", SG.REFERENCE_OUT_OF_SCOPE_REASON)

    ref = DataFrame(contrast = fill(HE.LN_VS_BG, 5),
                    unit = ["taxa", "sequences", SG.chemistry_unit("taxa"), "genera",
                            SG.chemistry_unit("genera")],
                    habitat = ["cat-a", "cat-a", "chem-x", "cat-a", "chem-x"],
                    ci_lo_std = zeros(5), ci_hi_std = ones(5))
    s = SG.scoped_reference(ref, NamedTuple[]; units = SG.REFERENCE_SCOPE_UNITS)
    @test nrow(s.table) == length(SG.REFERENCE_SCOPE_UNITS)
    @test isempty(s.excluded)
    # the rows that left come back BY NAME and counted, not absorbed into the scope
    @test s.out_of_scope == ["genera (chemistry) | chem-x", "genera | cat-a"]
    @test s.out_of_scope_units == ["genera", "genera (chemistry)"]
    # and the two declarations agree: nothing reported out of scope is in the declared scope
    @test all(u -> !(u in SG.REFERENCE_SCOPE_UNITS), s.out_of_scope_units)
    @test nrow(s.table) + length(s.out_of_scope) == nrow(ref)

    # THE DECLARED COUNT, which is the half a regenerated reference defeats: "N of N rows
    # agree" reads BETTER at the wrong N rather than different.
    @test SG.reference_scope_ok(nrow(s.table), nrow(s.table))
    @test !SG.reference_scope_ok(nrow(s.table), nrow(s.table) + 1)   # a thinned reference
    @test !SG.reference_scope_ok(nrow(s.table) + 1, nrow(s.table))   # a widened one
    # nothing was compared and nothing was declared: the refusal is the reference gate, not
    # this one, so a declined comparison does not fail on a scope it never had
    @test SG.reference_scope_ok(0, nothing)

    # AND IT IS THE GATE THAT READS IT. The predicate above is only a predicate until
    # something refuses to publish on it.
    scope_gate(compared, expected) =
        only([g for g in gates_for(; reference_rows_compared = compared,
                                   reference_rows_expected = expected)
              if g.name == "reference scope"])
    @test scope_gate(nrow(s.table), nrow(s.table)).ok
    @test !scope_gate(nrow(s.table), nrow(s.table) + 1).ok
end

# ---------------------------------------------------------------------------
# 7. The written tables.
# ---------------------------------------------------------------------------

@testset "the tables carry the unit, the contrast and both intervals" begin
    r = GRID.rows[1]
    for c in (:variant, :axis, :level, :role, :contrast, :unit, :category, :a, :b, :n1, :n2,
              :crude, :crude_lo, :crude_hi, :std, :lo, :hi, :n_strata,
              :n_zero_stratum_cells, :min_bg_stratum_cell, :tested, :q_bh, :feeds_flag)
        @test hasproperty(r, c)
    end
    @test nrow(DataFrame(GRID.rows)) == length(GRID.rows)

    flags = SG.robustness(GRID.rows, S)
    order = vcat(["primary"], S)
    md = SG.sensitivity_table(GRID.rows, order)
    @test occursin("| primary |", md)
    @test all(v -> occursin(v, md), S)

    unit_md = SG.unit_check_table(GRID.rows, order)
    @test occursin("genera", unit_md)
    @test occursin("NOT counted toward any flag", unit_md)
    @test occursin(SG.MIN_DISTINCT_SUFFIX, unit_md)
    @test occursin("unit", SG.columns_note()) && occursin("contrast", SG.columns_note())

    # Three different facts share the unnumbered state, and both tables name all three.
    for t in (md, unit_md)
        @test occursin("BLANK", t)
        @test occursin("n/a", t)
        @test occursin(SG.NOT_REACHED, t)
    end

    # A cell whose variant's axis cannot reach the unit is marked as that, not filled with
    # the primary's numbers under the variant's name. Otherwise a reader counts as many
    # independent perturbations as there are columns, and the effective number is smaller and
    # is stated nowhere.
    mv = SG.estimate_movement(GRID.rows, S)
    unreached = sort([k for k in collect(keys(mv.per_unit)) if mv.per_unit[k] == 0])
    @test !isempty(unreached)                 # the fixture reproduces the orthogonality
    for (v, unit) in unreached
        cats = sort(unique([r.category for r in GRID.rows
                            if r.variant == v && r.unit == unit && r.feeds_flag]))
        @test !isempty(cats)
        for c in cats
            p = first([r for r in GRID.rows if r.variant == "primary" && r.unit == unit &&
                                               r.category == c && r.feeds_flag])
            @test md_cell(md, unit, c, order, v) == SG.NOT_REACHED
            # the primary's own column still carries the number, so the marking is about the
            # variant's claim on the cell and not about the estimate being unavailable
            @test md_cell(md, unit, c, order, "primary") ==
                  SG.fmt_est(SG.fnum(p.std), SG.fnum(p.lo), SG.fnum(p.hi))
        end
    end
    # a variant that DOES reach the unit still prints its estimate
    reached = first(sort([k for k in collect(keys(mv.per_unit)) if mv.per_unit[k] > 0]))
    let (v, unit) = reached
        c = first(sort(unique([r.category for r in GRID.rows
                               if r.variant == v && r.unit == unit && r.feeds_flag])))
        @test md_cell(md, unit, c, order, v) != SG.NOT_REACHED
    end
    @test occursin("Variants whose axis reaches this unit", md)

    panel = SG.panel_rows(GRID.rows, flags)
    @test !isempty(panel)
    @test all(p -> p.unit in HE.PRIMARY_UNITS, panel)
    @test all(p -> p.flag in ("+", "-", "null", "fragile", "undecidable"), panel)
    # A category only some variant creates has no primary row, and is not drawn as one.
    @test all(p -> any(r -> r.variant == "primary" && r.unit == p.unit &&
                            r.category == p.category, GRID.rows), panel)

    # ... and the panel says so on the page, rather than leaving a reader holding the flag
    # table beside the picture to find rows in one that are missing from the other.
    note = SG.panel_note(flags, panel)
    @test occursin("weighted foreground count", note)
    drawn = Set((p.unit, p.category) for p in panel)
    off = [f for f in flags if !((f.unit, f.category) in drawn)]
    for f in off
        @test occursin(f.category, note)
    end
    if !isempty(off)
        @test occursin("Not drawn", note)
    end
    # a panel that drops nothing carries no such line, so the sentence cannot go stale
    @test !occursin("Not drawn", SG.panel_note(NamedTuple[], panel))

    # The renderer does not wrap, and an overrunning line is not clipped visibly: it ends
    # mid-sentence at the edge of the page. So the note arrives wrapped, on word boundaries.
    @test all(l -> length(l) <= SG.NOTE_WRAP, split(note, "\n"))
    long = join(fill("word", 80), " ")
    wrapped = SG.wrap_text(long, 20)
    @test all(l -> length(l) <= 20, split(wrapped, "\n"))
    @test replace(wrapped, "\n" => " ") == long      # nothing lost, nothing split
    @test SG.wrap_text("short enough", 40) == "short enough"
end

@testset "the two unit marks cannot print as one glyph" begin
    # Where the genus and the per-sequence estimate of a row agree closely, two marks at one
    # offset draw as a single glyph and the figure shows one measurement where it made two.
    # The offsets are far enough apart that the two markers' edges do not meet at this row
    # height — which is a property of four numbers, so it is checked rather than looked at.
    @test SG.GENUS_DY != SG.SEQUENCE_DY
    @test SG.unit_mark_separation_pt() > SG.unit_marks_touch_pt()
    # and both marks stay inside their own row, so the lower one cannot be read as the row
    # below's
    half_row = 0.5
    for (dy, ms) in ((SG.GENUS_DY, SG.GENUS_MS), (SG.SEQUENCE_DY, SG.SEQUENCE_MS))
        @test abs(dy) + (ms / 2) / (SG.ROW_MM * SG.PT_PER_MM) < half_row
    end
end

# ---------------------------------------------------------------------------
# 8. The gate between a finished run and the figure directory.
# ---------------------------------------------------------------------------

@testset "a run whose gate fails publishes nothing" begin
    # This is the check that stands between a failed run and the figure directory. Written
    # straight into that directory, a run whose anchors failed leaves a complete-looking set
    # of tables there — same names, same shapes, every number present — and is
    # distinguishable from a good run only by one provenance row and an exit status nobody
    # reads a month later.
    mktempdir() do dir
        figdir = joinpath(dir, "figs")
        mkpath(figdir)
        write(joinpath(figdir, "habitat_robustness.csv"), "the last good run\n")
        write(joinpath(figdir, "habitat_sensitivity.csv"), "the last good run\n")
        snapshot(d) = Dict(f => read(joinpath(d, f), String)
                           for f in readdir(d) if isfile(joinpath(d, f)))
        before = snapshot(figdir)

        stage = SG.open_staging(figdir)
        @test isdir(stage)
        @test isempty(readdir(stage))
        staged = ["habitat_robustness.csv", "habitat_sensitivity.csv"]
        for f in staged; write(joinpath(stage, f), "this run\n"); end

        err = refusal(SG.publish!, stage, figdir, staged,
                      [(ok = true, message = "the anchors were met"),
                       (ok = false, message = "an anchor target was not reproduced")])
        @test occursin("NOTHING WAS PUBLISHED", err)
        @test occursin("an anchor target was not reproduced", err)
        @test !occursin("the anchors were met", err)        # only the failures are reported
        @test snapshot(figdir) == before                    # the last good run is untouched
        @test all(f -> isfile(joinpath(stage, f)), staged)  # and the failure is inspectable

        # Every gate passed, and the same call publishes: the staged files replace the old
        # ones and the staging directory goes with them.
        written = SG.publish!(stage, figdir, staged, [(ok = true, message = "")])
        @test length(written) == length(staged)
        @test all(p -> read(p, String) == "this run\n", written)
        @test !isdir(stage) || isempty(readdir(stage))

        # A staging directory is cleared at the start of a run, so nothing an earlier attempt
        # left behind can be published by a later one.
        again = SG.open_staging(figdir)
        write(joinpath(again, "left_over.csv"), "from a failed run\n")
        @test isempty(readdir(SG.open_staging(figdir)))
        @test !isfile(joinpath(figdir, "left_over.csv"))
    end
end

@testset "every gate is built, and each one on its own stops the publish" begin
    # THE WIRING, NOT THE PREDICATES. Each predicate behind a gate is checked on its own
    # elsewhere in this file. What those checks do not cover is the six lines that HAND them
    # to `publish!`: written inline in the program body they would sit where no test can
    # reach them, so a gate entry deleted — or its `ok` quietly replaced by `true` — would
    # leave every predicate still tested, this whole suite still green, and the attack that
    # gate exists to refuse publishing at exit 0.
    gates = gates_for()
    @test [String(g.name) for g in gates] == String[String(n) for n in SG.PUBLISH_GATES]
    @test length(SG.PUBLISH_GATES) == length(unique(SG.PUBLISH_GATES))
    @test all(g -> g.ok, gates)
    @test all(g -> !isempty(String(g.message)), gates)   # a refusal a reader cannot diagnose

    # ONE BAD INPUT AT A TIME, and exactly one gate must answer for it. A gate that answers
    # for nothing is a gate that is no longer made; a gate that answers for someone else's
    # input is a gate wired to the wrong predicate.
    one_bad = Dict(
        "anchor targets"          => (target_failures = ["a target that was not reproduced"],),
        "published breaker lists" => (unanchored_levels = ["X:a_level_with_no_target"],),
        "target sources"          => (source_checks = [(file = "an_implementation.csv",
                                                        status = "differs",
                                                        declared = "the declared bytes",
                                                        found = "some other bytes",
                                                        ok = false)],),
        "reference rows"          => (reference_ok = false,),
        "reference scope"         => (reference_rows_compared = 2,),
        "stale panel"             => (stale_panels = ["a panel from the last run"],))
    @test sort(collect(keys(one_bad))) == sort(String[String(n) for n in SG.PUBLISH_GATES])

    for name in SG.PUBLISH_GATES
        g = gates_for(; one_bad[String(name)]...)
        @test String.([x.name for x in g if !x.ok]) == [String(name)]

        # AND THE GATE IS NOT A REPORT. `publish!` refuses on it, the figure directory keeps
        # the last good run byte for byte, and the staged tables stay where the failure is
        # inspectable. This is the half that has to hold for every gate, not only the one the
        # earlier case happened to build by hand.
        mktempdir() do dir
            write(joinpath(dir, "habitat_sensitivity.csv"), "the last good run\n")
            write(joinpath(dir, "habitat_robustness.csv"), "the last good run\n")
            snapshot() = Dict(f => read(joinpath(dir, f), String)
                              for f in readdir(dir) if isfile(joinpath(dir, f)))
            before = snapshot()
            staged = ["habitat_sensitivity.csv", "habitat_robustness.csv"]
            stage = SG.open_staging(dir)
            for f in staged; write(joinpath(stage, f), "this run\n"); end
            err = refusal(SG.publish!, stage, dir, staged, g)
            @test occursin("NOTHING WAS PUBLISHED", err)
            @test snapshot() == before
            @test all(f -> isfile(joinpath(stage, f)), staged)
            # the message reported is the one gate that failed, and no other
            failed = only([x for x in g if !x.ok])
            @test occursin(String(failed.message), err)
            for x in g; x.ok && @test !occursin(String(x.message), err); end
        end
    end

    # AND WITH EVERY GATE PASSING the same call publishes, so the refusals above are the gates
    # and not a `publish!` that never moves anything.
    mktempdir() do dir
        stage = SG.open_staging(dir)
        write(joinpath(stage, "habitat_sensitivity.csv"), "this run\n")
        written = SG.publish!(stage, dir, ["habitat_sensitivity.csv"], gates_for())
        @test length(written) == 1
        @test read(joinpath(dir, "habitat_sensitivity.csv"), String) == "this run\n"
    end
end

@testset "a run that draws no panel will not publish beside an old one" begin
    # `publish!` moves what was staged and leaves everything else where it is, so a
    # `--no-figure` run publishes its tables into a directory whose panel is the LAST run's.
    # The panel is the artefact a reader looks at: tables from run n beside the picture from
    # run n-1, under the names a reader and any downstream glob take as one result. One flag
    # reaches that state, and no file name in the directory distinguishes it from a complete
    # result.
    mktempdir() do dir
        @test SG.stale_panel_files(dir, false) == String[]      # nothing there to go stale
        pdf = SG.PANEL_BASE * ".pdf"; png = SG.PANEL_BASE * ".png"
        write(joinpath(dir, pdf), "last run\n")
        write(joinpath(dir, png), "last run\n")
        write(joinpath(dir, "habitat_sensitivity.csv"), "a table, not a panel\n")
        # A RUN THAT DRAWS ITS OWN PANEL replaces it, so there is nothing stale to refuse —
        # otherwise the check would refuse every ordinary re-run.
        @test SG.stale_panel_files(dir, true) == String[]
        @test SG.stale_panel_files(dir, false) == sort([pdf, png])
        # and it is the PANEL that is at stake, not everything in the directory: the tables
        # are the files this run replaces
        @test !any(f -> occursin("sensitivity.csv", f), SG.stale_panel_files(dir, false))

        # the gate is the publish gate, so nothing is published and the old files stand
        before = Dict(f => read(joinpath(dir, f), String)
                      for f in readdir(dir) if isfile(joinpath(dir, f)))
        stage = SG.open_staging(dir)
        write(joinpath(stage, "habitat_sensitivity.csv"), "this run\n")
        err = refusal(SG.publish!, stage, dir, ["habitat_sensitivity.csv"],
                      [(ok = isempty(SG.stale_panel_files(dir, false)),
                        message = "--no-figure was given and the directory already holds " *
                                  join(SG.stale_panel_files(dir, false), ";"))])
        @test occursin("NOTHING WAS PUBLISHED", err) && occursin(png, err)
        @test Dict(f => read(joinpath(dir, f), String)
                   for f in readdir(dir) if isfile(joinpath(dir, f))) == before
    end
end

@testset "the wiring check is refused by default and declined only on purpose" begin
    # A check that is skipped by default is a check nobody runs: given no targets this script
    # would write every table and exit 0, and its flag table would be byte-identical to a
    # fully checked run's.
    err = refusal(SG.parse_run_options, String[])
    @test occursin("--anchor-targets", err)
    @test occursin("--no-anchors", err)

    # WHILE THE ANCHOR LEVELS ARE RETIRED — the production setting — the targets are the one
    # half required. The grid holds no all-legacy rows, so a reference would be compared with
    # nothing whatever file were named, and requiring one would require a path nobody opens.
    @test SG.ANCHOR_LEVELS_RETIRED
    @test !occursin("--anchor-reference <path>", err)
    @test occursin("not required", err)
    alone = SG.parse_run_options(["--anchor-targets", "targets.csv"])
    @test !alone.no_anchors && alone.anchor_targets == "targets.csv"
    @test isempty(alone.anchor_reference)
    # A reference given anyway still arrives, so the run can record it by name.
    named = SG.parse_run_options(["--anchor-targets", "targets.csv",
                                  "--anchor-reference", "rows.csv"])
    @test named.anchor_reference == "rows.csv"
    # ... but it does not stand in for the half that IS required.
    other = refusal(SG.parse_run_options, ["--anchor-reference", "rows.csv"])
    @test occursin("--anchor-targets", other)
    # The production call and the explicit production setting are one rule.
    @test SG.parse_run_options(["--anchor-targets", "t.csv"];
                               retired = SG.ANCHOR_LEVELS_RETIRED).anchor_targets == "t.csv"

    # WITH THE LEVELS EVALUATED, BOTH HALVES, not only the first. The row-by-row comparison
    # against the earlier run is the half one dropped flag removes: the run records the
    # comparison as "not supplied" and publishes everything, and the only trace is a
    # provenance row nobody diffs.
    evaluated(args) = SG.parse_run_options(args; retired = false)
    both = refusal(evaluated, String[])
    @test occursin("--anchor-targets", both) && occursin("--anchor-reference", both)
    @test occursin("2 of the two wiring checks", both)
    half = refusal(evaluated, ["--anchor-targets", "targets.csv"])
    @test occursin("--anchor-reference <path>", half)
    @test occursin("'not supplied'", half)
    @test !occursin("--anchor-targets <path>", half)          # the half that WAS given
    @test occursin("--anchor-targets", refusal(evaluated, ["--anchor-reference", "rows.csv"]))
    t = evaluated(["--anchor-targets", "targets.csv", "--anchor-reference", "rows.csv"])
    @test !t.no_anchors && t.anchor_targets == "targets.csv"
    @test t.anchor_reference == "rows.csv"

    # --no-anchors declines the check in either setting, and nothing is required with it.
    for retired in (true, false)
        o = SG.parse_run_options(["--no-anchors"]; retired = retired)
        @test o.no_anchors && isempty(o.anchor_targets) && isempty(o.anchor_reference)
    end
    # the rest of the command line still arrives, so the refusal is about the check and not
    # about parsing
    @test SG.parse_run_options(["--no-anchors", "--no-figure"]).figure == false
    @test SG.parse_run_options(["--no-anchors", "--figdir", "somewhere"]).figdir == "somewhere"
end

@testset "a reference named while the levels are retired is recorded and not opened" begin
    recorded(path) = (out = Pair{String,Any}[];
                      ok = SG.record_reference_not_applicable!((k, v) -> push!(out, k => v),
                                                               path);
                      (ok = ok, out = Dict(out), keys = first.(out)))
    none = recorded("")
    @test none.ok
    @test none.out["anchor_reference_status"] == "not applicable: the anchor levels are retired"
    @test none.out["anchor_reference_rows_compared"] == 0
    @test !haskey(none.out, "anchor_reference_named_not_opened")
    mktempdir() do dir
        # NOT OPENED: a name that is no file at all is recorded exactly as one that is, and no
        # digest is taken of either.
        for p in (joinpath(dir, "absent_rows.csv"), joinpath(dir, "present_rows.csv"))
            endswith(p, "present_rows.csv") && write(p, "not,a,rows,table\n")
            r = recorded(p)
            @test r.ok
            @test r.out["anchor_reference_named_not_opened"] == basename(p)
            @test r.out["anchor_reference_status"] ==
                  "not applicable: the anchor levels are retired"
            @test !any(k -> occursin("sha256", k) || k == "anchor_reference_file", r.keys)
        end
    end
end

# THE RUN'S OWN BRANCH. The helper above decides what a run records only if the program calls
# it, and no include reaches the program block; so the branch is read as source: it records
# through the helper with the reference as typed, and nothing in it opens, stats or digests a
# file.
@testset "the run's retired branch records a named reference and opens nothing" begin
    src = read(joinpath(@__DIR__, "..", "scripts", "habitat_sensitivity.jl"), String)
    i = findfirst("if ANCHOR_LEVELS_RETIRED && isempty(anchor_rows)", src)
    j = i === nothing ? nothing : findnext("\nelseif !isempty(OPT.anchor_reference)", src, last(i))
    @test i !== nothing && j !== nothing
    branch = src[nextind(src, last(i)):prevind(src, first(j))]
    code = join(filter(l -> !startswith(lstrip(l), "#"), split(branch, '\n')), '\n')
    @test occursin("ref_ok = record_reference_not_applicable!(rec!, OPT.anchor_reference)", code)
    @test !occursin(r"\bread\(|\bisfile\(|\bopen\(|sha256|read_reference_rows", code)
end

end # testset
