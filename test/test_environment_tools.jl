#=
    Offline self-test for the habitat classification + aggregation. No network: the OBT
    ontology climb is dependency-injected (`path_of`), the same pattern contrastive_potts.jl
    uses for the DCA backend.

    The failures this suite exists to prevent, each of them reachable from one line of code:

    1. `json[1]`. OmniCrobe returns up to thousands of habitat relations per taxon, so
       taking the FIRST one assigns *Methylorubrum extorquens* — a model phyllosphere
       methylotroph — the habitat "blood", and everything downstream inherits that single
       arbitrary draw.

    2. An INVERTED ontology map. OBT:000013 is the natural-environment branch (soil, rock,
       forest, air), not "Aquatic", and OBT:000006 is the man-made branch (saltern, mine,
       sewage), not "Terrestrial". A map keyed at those two classifies `soil` and even
       `terrestrial habitat` as Aquatic, and leaves four of the six configured ids matching
       nothing at all. Nothing but a test on real paths catches it.

    3. A KEY THAT DOES NOT EXIST. A relation parser reading `obt_objects` never fires — the
       API returns the singular `obt_object` — and every habitat name then falls through to
       the fallback `obt_forms[1]`, the first surface form of the sentence rather than the
       ontology term ("spring, mesophilic saline spring" where the term is "spring"). Wrong
       names feed the term tallies and the representative-name column.

    4. `catch → empty vector`. Returning one value for both a failed fetch and a taxon with
       no habitat relation lets a network blip enter the analysis as the ecological claim
       "this organism has no known habitat". A failed fetch returns `nothing`.

    5. A DICT-ORDER WINNER. Picking the category with
       `first(sort(collect(tally), by = kv -> -kv[2]))` reads a `Dict`, which iterates in
       hash order, through a stable sort: among equal tallies the winner is whichever key
       the table happens to yield first — and the same stable sort on the term ranking
       decides which terms the 50-term cap keeps. Two runs over the same evidence can then
       disagree with no output column saying a tie existed at all. Ties are split
       explicitly instead.

    So the ANCHOR fixtures below use REAL ontology paths, copied verbatim from the live OBT
    API (https://omnicrobe.migale.inrae.fr/api/get/obt/<id>). A hand-invented path proves
    nothing —
    it would agree with an inverted table. The relation fixture is copied from a real
    /api/search/relations item for the same reason: a hand-invented payload would agree
    with an `obt_objects` parser too.

    Run:  julia --project=. test/test_environment_tools.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test, Logging, Random, DataFrames
include(joinpath(@__DIR__, "..", "src", "ModuleEnvironmentTools.jl"))
using .EnvironmentTools

# Expected warnings are silenced so a passing run stays readable; the one assertion that
# the warning is actually emitted uses @test_logs instead.
quiet(f) = with_logger(f, NullLogger())

# ---------------------------------------------------------------------------
# REAL paths, fetched from https://omnicrobe.migale.inrae.fr/api/get/obt/<id>
# ---------------------------------------------------------------------------
const REAL = Dict(
    "soil"                => "/OBT:000001/OBT:000013/OBT:000158/OBT:000427",
    "terrestrial habitat" => "/OBT:000001/OBT:000013/OBT:000158",
    "lake"                => "/OBT:000001/OBT:000013/OBT:000047/OBT:000326/OBT:000705/OBT:001171",
    "marine environment"  => "/OBT:000001/OBT:000013/OBT:000047/OBT:000338",
    "rock"                => "/OBT:000001/OBT:000013/OBT:000079/OBT:000350/OBT:000846",
    "solar saltern"       => "/OBT:000001/OBT:000006/OBT:000050/OBT:003164/OBT:003233 " *
                             "/OBT:000001/OBT:000009/OBT:000090/OBT:000302/OBT:000667/OBT:003164/OBT:003233 " *
                             "/OBT:000001/OBT:000013/OBT:000047/OBT:000326/OBT:000705/OBT:003164/OBT:003233",
    "mine"                => "/OBT:000001/OBT:000006/OBT:000100/OBT:000324/OBT:000624/OBT:001210",
    "sewage"              => "/OBT:000001/OBT:000006/OBT:000100/OBT:000325/OBT:000611/OBT:001307",
    "gut"                 => "/OBT:000001/OBT:000014/OBT:000039/OBT:000258/OBT:000642/OBT:001140",
    "rumen"               => "/OBT:000001/OBT:000014/OBT:000039/OBT:000465/OBT:000723/OBT:001289/OBT:001749",
    "cell"                => "/OBT:000001/OBT:000014/OBT:000061",
    "Arabidopsis"         => "/OBT:000001/OBT:000010/OBT:000089/OBT:000393/OBT:000912/OBT:001213/OBT:001390",
    "cattle"              => "/OBT:000001/OBT:000004/OBT:000035/OBT:002726/OBT:003295/OBT:003532",
    "salt"                => "/OBT:000001/OBT:000008/OBT:000094/OBT:000240/OBT:000477/OBT:001044/OBT:001753",
    # NOTE: these two carry ALL of the term's DAG lineages joined by a space, which is what
    # get_obt_path returns. Their first lineage is chemico-physical only; the aquatic one
    # comes second, and dropping it (path[1]) is what types `saline lake` as "Other".
    "hypersaline lake"    => "/OBT:000001/OBT:000009/OBT:002018/OBT:002592/OBT:003226/OBT:003261 " *
                             "/OBT:000001/OBT:000013/OBT:000047/OBT:000269/OBT:002018/OBT:002592/OBT:003226/OBT:003261 " *
                             "/OBT:000001/OBT:000013/OBT:000047/OBT:000326/OBT:000705/OBT:001171/OBT:003226/OBT:003261",
    "saline lake"         => "/OBT:000001/OBT:000009/OBT:002018/OBT:002592/OBT:003226 " *
                             "/OBT:000001/OBT:000013/OBT:000047/OBT:000269/OBT:002018/OBT:002592/OBT:003226 " *
                             "/OBT:000001/OBT:000013/OBT:000047/OBT:000326/OBT:000705/OBT:001171/OBT:003226",
    "soda lake"           => "/OBT:000001/OBT:000009/OBT:000090/OBT:000187/OBT:002940",
    "salt contaminated soil" => "/OBT:000001/OBT:000009/OBT:000141/OBT:000244/OBT:003043/OBT:003115 " *
                                "/OBT:000001/OBT:000009/OBT:002257/OBT:003043/OBT:003115 " *
                                "/OBT:000001/OBT:000013/OBT:000158/OBT:000427/OBT:002257/OBT:003043/OBT:003115",
    "anaerobic environment"  => "/OBT:000001/OBT:000009/OBT:000078/OBT:000191",
    "nutrient broth"      => "/OBT:000001/OBT:000007/OBT:000106/OBT:000360",
    "rhizosphere"         => "/OBT:000001/OBT:000014/OBT:000140/OBT:000410",
    # `living organism` (OBT:000010) in its OWN right: the lineage is the one every host
    # term above runs through, stopping at the generic node itself. Its type is Host-other,
    # which is exactly why leaving it in hands that pillar the taxa it says nothing about.
    "living organism"     => "/OBT:000001/OBT:000010",
    # `atmosphere part` (OBT:000052), the same lineage asserted inline above, kept here so
    # the label-neutrality sweep below reaches the Atmospheric type too.
    "atmosphere part"     => "/OBT:000001/OBT:000013/OBT:000052",
    # A term whose CONDITION is written down only in its name: the path is `lake`'s own
    # lineage, verbatim, so nothing thermal can enter through the ontology and the text
    # fallback is genuinely the thing under test (same shape as "lake, Salt Lake" above).
    "hot spring"          => "/OBT:000001/OBT:000013/OBT:000047/OBT:000326/OBT:000705/OBT:001171",
)
path_of(oid) = get(REAL, oid, "")     # fixtures key relations by term name for readability
rel(nm) = (obtid = nm, name = nm)

@testset "habitat_type — anchors against an inverted branch map" begin
    # These four are the exact cases an inverted map gets wrong (soil and terrestrial
    # habitat would read Aquatic; every engineered site would read Terrestrial).
    @test habitat_type(REAL["soil"])                == "Terrestrial"
    @test habitat_type(REAL["terrestrial habitat"]) == "Terrestrial"
    @test habitat_type(REAL["sewage"])              == "Anthropogenic"
    # A solar saltern is a hypersaline water body: OBT gives it an aquatic lineage
    # (/OBT:000013/OBT:000047/…) alongside its artificial one, so it types Aquatic and the
    # salinity rides the second axis. A map keyed at the man-made branch would call it
    # Terrestrial, which it is not.
    @test habitat_type(REAL["solar saltern"])       == "Aquatic"
    @test habitat_type(REAL["solar saltern"])       != "Terrestrial"

    @test habitat_type(REAL["lake"])               == "Aquatic"
    @test habitat_type(REAL["marine environment"]) == "Aquatic"
    @test habitat_type(REAL["rock"])               == "Geological"

    # OBT:000079 "environmental matter" must NOT be labelled wholesale. Only `rock` under
    # OBT:000350 "mineral matter" is geological; the branch also holds the bare term
    # `water` and aquatic sediments, which a parent-level "Geological" would have
    # classified as rock — the same failure shape as natural-environment→Aquatic.
    @testset "environmental-matter branch is split by child, not labelled at the parent" begin
        @test habitat_type("/OBT:000001/OBT:000013/OBT:000079/OBT:000468") == "Aquatic"      # water
        @test habitat_type("/OBT:000001/OBT:000013/OBT:000079/OBT:000350/OBT:000846") == "Geological"  # rock
        # Aquatic sediment routes through OBT:000428 "soil part" (Terrestrial) on its way to
        # OBT:001430 (Aquatic), so this also pins the precedence: Aquatic is tested first.
        @test habitat_type("/OBT:000001/OBT:000013/OBT:000079/OBT:000428/OBT:000885/OBT:001303/OBT:001430") == "Aquatic"
        # the parent alone resolves to nothing, which is the point
        @test habitat_type("/OBT:000001/OBT:000013/OBT:000079") == "Other"
    end

    @testset "atmosphere and whole-organism hosts" begin
        @test habitat_type("/OBT:000001/OBT:000013/OBT:000052") == "Atmospheric"          # air
        # Whole organisms, not just anatomy: `animal`/`plant` are separate branches from
        # `animal part`/`plant habitat`, so a table fed only by the part branches would
        # capture guts and roots while the organisms themselves fell to Host-other.
        @test habitat_type("/OBT:000001/OBT:000010/OBT:000089/OBT:000193") == "Host-animal"
        @test habitat_type("/OBT:000001/OBT:000010/OBT:000089/OBT:000393") == "Host-plant"
        # A saltern is aquatic BY THE ONTOLOGY (third lineage), not by our fiat.
        @test habitat_type("/OBT:000001/OBT:000006/OBT:000050/OBT:003164 " *
                           "/OBT:000001/OBT:000013/OBT:000047/OBT:000326/OBT:000705/OBT:003164") == "Aquatic"
    end
    @test habitat_type(REAL["mine"])               == "Anthropogenic"
    @test habitat_type(REAL["gut"])                == "Host-animal"
    @test habitat_type(REAL["rumen"])              == "Host-animal"
    @test habitat_type(REAL["rhizosphere"])        == "Host-plant"
    @test habitat_type(REAL["cell"])               == "Host-other"
    # Arabidopsis is a whole plant: its lineage runs through OBT:000393 "plant", so it must
    # type Host-plant. Without the whole-organism branches it falls to Host-other, i.e. a
    # table that captures roots and rhizospheres but not the plants they belong to.
    @test habitat_type(REAL["Arabidopsis"])        == "Host-plant"
    @test habitat_type(REAL["cattle"])             == "Agricultural"
    @test habitat_type(REAL["salt"])               == "Food/substance"
    @test habitat_type(REAL["nutrient broth"])     == "Lab-medium"

    @testset "isolation sources are flagged, not silently binned as habitat" begin
        @test is_isolation_source("Food/substance")
        @test is_isolation_source("Lab-medium")
        @test is_isolation_source(habitat_type(REAL["nutrient broth"]))
        @test is_isolation_source(habitat_type(REAL["salt"]))
        for t in ("Aquatic", "Terrestrial", "Host-plant", "Anthropogenic", "Geological")
            @test !is_isolation_source(t)
        end
        @test !is_isolation_source(missing)
    end
    @test habitat_type("/OBT:000001/OBT:999999")   == "Other"

    # Specific anatomy branches must be tested BEFORE the general living-organism one,
    # or gut/rhizosphere would collapse into Host-other.
    @test habitat_type(REAL["gut"])         != "Host-other"
    @test habitat_type(REAL["rhizosphere"]) != "Host-other"
end

@testset "coarse_pillar keeps the original five" begin
    @test coarse_pillar("Aquatic")        == "Aquatic"
    @test coarse_pillar("Geological")     == "Terrestrial"
    @test coarse_pillar("Host-plant")     == "Host-associated"
    @test coarse_pillar("Host-animal")    == "Host-associated"
    @test coarse_pillar("Agricultural")   == "Anthropogenic"
    @test coarse_pillar("Food/substance") == "Anthropogenic"
    @test coarse_pillar("Atmospheric")    == "Other"     # deliberate, not an unmapped gap
    @test coarse_pillar("Other")          == "Other"
    # every fine type must have an explicit coarse mapping — an unmapped one would fall to
    # "Other" silently, which is how a new pillar could quietly disappear from the χ².
    for (fine, _) in EnvironmentTools.OBT_TYPE
        @test haskey(EnvironmentTools.COARSE, fine)
    end
    @test category_from_path(REAL["soil"]) == "Terrestrial"     # end-to-end, through the table
end

@testset "habitat_chemistry — the orthogonal axis" begin
    # From the ontology's own chemico-physical branch (OBT:000009).
    @test habitat_chemistry(REAL["hypersaline lake"])       == "Hypersaline"
    @test habitat_chemistry(REAL["saline lake"])            == "Hypersaline"
    @test habitat_chemistry(REAL["salt contaminated soil"]) == "Hypersaline"
    @test habitat_chemistry(REAL["anaerobic environment"])  == "Anoxic"
    @test occursin("Alkaline", habitat_chemistry(REAL["soda lake"], "soda lake"))

    # OBT is a DAG. A saline lake hangs under BOTH the chemico-physical branch and the
    # aquatic branch; the chemico-physical lineage is returned first. Keeping only that one
    # (path[1] alone) leaves it typeless, so this is the assertion that both axes resolve.
    @testset "multi-parent terms resolve on BOTH axes" begin
        @test habitat_type(REAL["saline lake"])              == "Aquatic"
        @test habitat_chemistry(REAL["saline lake"])         == "Hypersaline"
        @test habitat_type(REAL["hypersaline lake"])         == "Aquatic"
        @test habitat_type(REAL["salt contaminated soil"])   == "Terrestrial"
        @test habitat_chemistry(REAL["salt contaminated soil"]) == "Hypersaline"
        # …and taking only the first lineage loses the second axis, so the both-axes walk
        # is load-bearing.
        @test habitat_type(first(split(REAL["saline lake"], ' '))) == "Other"
    end

    # Salinity is genuinely orthogonal to type: the SAME chemistry flag lands on different
    # habitat types, which one pillar list cannot express.
    @test habitat_type(REAL["hypersaline lake"])     == "Aquatic"
    @test habitat_type(REAL["salt contaminated soil"]) == "Terrestrial"
    @test habitat_type(REAL["hypersaline lake"]) != habitat_type(REAL["salt contaminated soil"])
    for (p, n) in [(REAL["hypersaline lake"], "hypersaline lake"),
                   (REAL["solar saltern"], "solar saltern"),
                   (REAL["salt contaminated soil"], "salt contaminated soil")]
        @test occursin("Hypersaline", habitat_chemistry(p, n))
    end

    # Text fallback: OmniCrobe gives one path, so a salt lake filed under Aquatic never
    # reaches the chemistry branch and must be caught by the term text.
    @test habitat_chemistry(REAL["lake"]) == ""                       # plain lake: no flag
    @test habitat_chemistry(REAL["lake"], "lake, Salt Lake") == "Hypersaline"
    @test habitat_chemistry(REAL["mine"], "mine, salt mine")  == "Hypersaline"

    # No false positives on ordinary habitats.
    for nm in ("soil", "lake", "marine environment", "gut", "rock")
        @test habitat_chemistry(REAL[nm], nm) == ""
    end

    @testset "bare `salt` is chemistry, and must not fire on basalt" begin
        # OmniCrobe can return habitat = "salt" for a halophile, including *salifodinae*
        # (salt-mine) strains. OBT files `salt` under "food", so without this a salt-mine
        # archaeon reads as an anthropogenic food habitat and nothing else.
        @test habitat_chemistry(REAL["salt"], "salt")       == "Hypersaline"
        @test habitat_chemistry(REAL["salt"], "salt, salt") == "Hypersaline"
        @test habitat_chemistry(REAL["salt"], "salt, salts") == "Hypersaline"
        # …but `basalt` is a rock and must not match.
        @test habitat_chemistry(REAL["rock"], "basalt")     == ""
        @test habitat_chemistry(REAL["rock"], "basaltic glass") == ""
    end

    # Anti-regression: four PARENT nodes look like the right ids but are far too broad,
    # and each was rejected during construction. If someone "simplifies" the table back to
    # a parent, these fail. Paths are real; only the leaf differs from the observed terms.
    @testset "over-broad parent nodes must not be used" begin
        # OBT:000141 "polluted environment" / OBT:000244 "contaminated site" — a contaminated
        # soil that is NOT salt-contaminated must not read as Hypersaline.
        @test habitat_chemistry("/OBT:000001/OBT:000009/OBT:000141/OBT:000244/OBT:003043",
                                "contaminated soil") == ""
        # OBT:000078 "environment wrt oxygen level" is the parent of aerobic AND anaerobic.
        @test habitat_chemistry("/OBT:000001/OBT:000009/OBT:000078", "aerobic environment") == ""
        # OBT:000090 "extreme environment" covers thermal/acidic too, not just alkaline.
        @test habitat_chemistry("/OBT:000001/OBT:000009/OBT:000090", "extreme environment") == ""
        # OBT:002018 "environmental water with chemical property" — any chemistry, not saline.
        @test habitat_chemistry("/OBT:000001/OBT:000009/OBT:002018", "acidic water") == ""
    end
end

@testset "aggregate_habitats" begin
    @testset "the json[1] failure mode is fixed" begin
        # First relation is a host term; the evidence is overwhelmingly soil/lake.
        rels = vcat([rel("gut")],
                    [rel("soil") for _ in 1:40],
                    [rel("lake") for _ in 1:10])
        @test habitat_type(path_of(rels[1].obtid)) == "Host-animal"   # what json[1] gave
        a = aggregate_habitats(rels, path_of; weighting = :mentions)
        @test a.habitat_type == "Terrestrial"
        @test a.name == "soil"
        @test a.n_relations == 51 && a.n_distinct == 3
    end

    @testset ":distinct vs :mentions weighting" begin
        rels = vcat([rel("lake") for _ in 1:100], [rel("soil"), rel("rock")])
        @test aggregate_habitats(rels, path_of; weighting = :mentions).category == "Aquatic"
        # Under :distinct each TERM counts once: 1 Aquatic vs 2 collapsing to Terrestrial.
        @test aggregate_habitats(rels, path_of; weighting = :distinct).category == "Terrestrial"
    end

    @testset "confidence reflects how split the evidence is" begin
        @test aggregate_habitats([rel("soil"), rel("soil")], path_of).confidence == 1.0
        @test aggregate_habitats([rel("soil"), rel("lake")], path_of).confidence == 0.5
    end

    @testset "chemistry is carried alongside type, not instead of it" begin
        # A haloarchaeon: saline habitats spanning aquatic, terrestrial and man-made.
        rels = [rel("hypersaline lake"), rel("saline lake"),
                rel("salt contaminated soil"), rel("solar saltern")]
        a = aggregate_habitats(rels, path_of)
        @test a.chemistry == "Hypersaline"
        @test a.chemistry_share == 1.0            # every term agrees on the chemistry …
        @test a.type_confidence < 1.0             # … while the TYPE is genuinely split
        @test a.chemistry_coverage == 1.0
    end

    @testset "chemistry_coverage separates 'not saline' from 'no evidence'" begin
        plain = aggregate_habitats([rel("soil"), rel("lake")], path_of)
        @test plain.chemistry == "" && plain.chemistry_coverage == 0.0
    end

    @testset "empty input and truncation are reported, not silent" begin
        e = aggregate_habitats(NamedTuple{(:obtid, :name), Tuple{String,String}}[], path_of)
        @test ismissing(e.category) && e.n_relations == 0
        t = aggregate_habitats([rel("soil"), rel("lake"), rel("rock")], path_of; max_distinct = 2)
        @test t.n_distinct == 3 && t.n_truncated == 1
    end

    @testset "profiles carry the full distributions" begin
        a = aggregate_habitats([rel("soil"), rel("lake")], path_of)
        @test occursin("Terrestrial:0.5", a.profile) && occursin("Aquatic:0.5", a.profile)
    end
end

@testset "the term filter and the tie rule, as units" begin
    # `is_generic_term` must read the term's OWN id, never its lineage. OBT:000010
    # "living organism" is the ANCESTOR of every plant and animal host term, so a substring
    # test against the path would delete the whole host-associated pillar while trying to
    # delete two uninformative nodes. The id is the LEAF of each lineage, hence both forms.
    @test is_generic_term("OBT:000010")
    @test is_generic_term("OBT:000061")
    @test is_generic_term("living organism", "/OBT:000001/OBT:000010")
    @test is_generic_term("cell", REAL["cell"])                    # lineage leaf OBT:000061
    @test !is_generic_term("Arabidopsis", REAL["Arabidopsis"])     # OBT:000010 is UPSTREAM
    @test !is_generic_term("gut", REAL["gut"])
    @test !is_generic_term("rhizosphere", REAL["rhizosphere"])
    @test !is_generic_term("soil", REAL["soil"])
    # …and the pillars those same paths carry are untouched by the filter.
    @test habitat_type(REAL["Arabidopsis"]) == "Host-plant"
    @test habitat_type(REAL["gut"]) == "Host-animal"

    # `resolve` stands where `first(sort(collect(t), by = kv -> -kv[2]))` would: that
    # expression returns whichever key the hash table yields first among equals — a coin toss
    # no output column records, and one that also decides which terms the max_distinct cap
    # drops.
    R = EnvironmentTools.resolve
    tie  = Dict("Aquatic" => 1, "Terrestrial" => 1)
    won  = Dict("Aquatic" => 1, "Terrestrial" => 2)
    @test R(won, Dict{String,Int}(), :mentions_fractional) ==
          ("Terrestrial", Dict("Terrestrial" => 1.0), 1, :none)
    # A residual tie is SPLIT 1/k, and labelled with a value no single category can take,
    # so downstream code cannot read a tie as a call.
    @test R(tie, tie, :mentions_fractional) ==
          ("Aquatic|Terrestrial", Dict("Aquatic" => 0.5, "Terrestrial" => 0.5), 2, :fractional)
    # Mentions break the tie before the split does.
    @test R(tie, Dict("Aquatic" => 3, "Terrestrial" => 1), :mentions_fractional) ==
          ("Aquatic", Dict("Aquatic" => 1.0), 1, :mentions)
    # :fractional ignores the mention count; :alphabetical picks one but still REPORTS that
    # there was a tie to pick from, which a bare single-winner pick cannot do.
    @test R(tie, Dict("Aquatic" => 3, "Terrestrial" => 1), :fractional)[1] ==
          "Aquatic|Terrestrial"
    @test R(tie, tie, :alphabetical) == ("Aquatic", Dict("Aquatic" => 1.0), 2, :alphabetical)
    @test R(Dict{String,Int}(), Dict{String,Int}(), :mentions_fractional)[1] == ""
    # Whatever the rule, the weights are a distribution: a tied taxon is still one taxon.
    for rule in (:mentions_fractional, :fractional, :alphabetical, :legacy_dict)
        @test sum(values(R(tie, tie, rule)[2])) ≈ 1.0
    end
end

@testset "agg_habitat keyword roll-up (unchanged by default)" begin
    @testset "default behaviour is unchanged" begin
        @test agg_habitat("rhizoplane")   == "Terrestrial"
        @test agg_habitat("phyllosphere") == "Other"
        @test agg_habitat("rhizosphere")  == "Host-associated"
        @test agg_habitat("human gut")    == "Host-associated"
        @test agg_habitat("marine water") == "Aquatic"
        @test agg_habitat("sewage")       == "Anthropogenic"
        @test agg_habitat(missing)        == "_Not Mapped_"
    end

    @testset "split_host separates plant from animal hosts" begin
        @test agg_habitat("rhizosphere";  split_host = true) == "Plant-associated"
        @test agg_habitat("root nodule";  split_host = true) == "Plant-associated"
        @test agg_habitat("human gut";    split_host = true) == "Animal-associated"
        @test agg_habitat("host";         split_host = true) == "Host-other"
        @test agg_habitat("rhizoplane";   split_host = true) == "Plant-associated"
        @test agg_habitat("phyllosphere"; split_host = true) == "Plant-associated"
        @test agg_habitat("marine water"; split_host = true) == "Aquatic"
        @test agg_habitat("soil";         split_host = true) == "Terrestrial"
    end
end

# ---------------------------------------------------------------------------
# Relation parsing. The fixture is one item of a real /api/search/relations
# response, keys and all; only the ids are the term names the REAL dict uses,
# so the injected `path_of` resolves them.
# ---------------------------------------------------------------------------
const REL_JSON = """
[
 {"docs": ["EU728658"],
  "obt_forms": ["spring, mesophilic saline spring"],
  "obt_object": {"name": "spring", "obtid": "OBT:001333"},
  "obtid": "OBT:001333",
  "obtroot": "ALL",
  "source": "GenBank",
  "taxid": "ncbi:560819",
  "taxon_forms": ["Tistlia consotensis USBA 355"],
  "taxon_object": {"name": "Tistlia consotensis USBA 355", "taxid": "ncbi:560819"},
  "taxroot": "ncbi:560819",
  "type": "habitat"}
]
"""

# Three distinct terms, two of which are Terrestrial, and `soil` asserted twice: both the
# winning type and the representative name are decided by a margin, so this fixture does
# not depend on how the aggregator breaks ties (that rule is a separate step).
const REL_MIX_JSON = """
[
 {"obtid": "soil", "obt_object": {"name": "soil"}, "obt_forms": ["soil", "soils"],
  "source": "PubMed", "type": "habitat"},
 {"obtid": "terrestrial habitat", "obt_object": {"name": "terrestrial habitat"},
  "obt_forms": ["terrestrial"], "source": "EPMC", "type": "habitat"},
 {"obtid": "lake", "obt_object": {"name": "lake"}, "obt_forms": ["lake"],
  "source": "DSMZ", "type": "habitat"},
 {"obtid": "soil", "obt_object": {"name": "soil"}, "obt_forms": ["soil"],
  "source": "EPMC", "type": "habitat"}
]
"""

@testset "parse_relations — the singular key the API really returns" begin
    r = EnvironmentTools.parse_relations(REL_JSON)
    @test length(r) == 1
    @test r[1].obtid  == "OBT:001333"
    @test r[1].name   == "spring"                                   # obt_object.name
    @test "spring, mesophilic saline spring" in r[1].forms
    @test r[1].source == "GenBank"
    # THE TRAP: with an `obt_objects` lookup the name comes from obt_forms[1], i.e. the
    # sentence's surface form and not the ontology term. If the two ever coincide the
    # regression is invisible, so the fixture is one where they differ.
    @test r[1].name != r[1].forms[1]
    @test eltype(r) == EnvironmentTools.Relation

    @testset "obt_objects survives as a legacy alias" begin
        obj  = """[{"obtid": "soil", "obt_objects": {"name": "soil"},
                    "obt_forms": ["soils"], "source": "PubMed"}]"""
        list = """[{"obtid": "lake", "obt_objects": [{"name": "lake"}],
                    "obt_forms": ["lakes"], "source": "PubMed"}]"""
        @test EnvironmentTools.parse_relations(obj)[1].name  == "soil"
        @test EnvironmentTools.parse_relations(list)[1].name == "lake"
    end

    @testset "forms are the fallback, not the source of truth" begin
        # No object at all: the first surface form is all there is.
        only_forms = """[{"obtid": "soil", "obt_forms": ["soil, paddy soil"]}]"""
        p = EnvironmentTools.parse_relations(only_forms)[1]
        @test p.name == "soil, paddy soil"
        @test p.source == ""                       # absent source, not a crash
        # EVERY form is kept: the chemistry axis reads the term text, and dropping the
        # forms throws away the only place where "saline" is written down.
        many = """[{"obtid": "lake", "obt_object": {"name": "lake"},
                    "obt_forms": ["lake", "saline lake", "Salt Lake"]}]"""
        @test EnvironmentTools.parse_relations(many)[1].forms ==
              ["lake", "saline lake", "Salt Lake"]
    end

    @testset "unusable items are skipped, not guessed at" begin
        # No obtid: nothing to climb the ontology with.
        @test isempty(EnvironmentTools.parse_relations(
            """[{"obt_object": {"name": "soil"}, "obt_forms": ["soil"]}]"""))
        @test isempty(EnvironmentTools.parse_relations("""[{"obtid": ""}]"""))
        # An obtid with no name anywhere: kept out rather than entered blank.
        @test isempty(EnvironmentTools.parse_relations("""[{"obtid": "soil"}]"""))
        # …while the usable neighbours of a skipped item still parse.
        mixed = """[{"obt_forms": ["soil"]}, {"obtid": "lake", "obt_object": {"name": "lake"}}]"""
        @test length(EnvironmentTools.parse_relations(mixed)) == 1
    end

    @testset "an empty response is an empty vector, a broken one raises" begin
        e = EnvironmentTools.parse_relations("[]")
        @test isempty(e) && e isa Vector{EnvironmentTools.Relation}
        # An error payload is an OBJECT. Parsing it as zero relations would read as "this
        # taxon has no habitat", which is the confusion this whole step removes.
        @test_throws ArgumentError EnvironmentTools.parse_relations("""{"detail": "nope"}""")
        @test_throws Exception EnvironmentTools.parse_relations("not json at all")
    end

    @testset "parsed rows still feed aggregate_habitats" begin
        a = aggregate_habitats(EnvironmentTools.parse_relations(REL_MIX_JSON), path_of)
        @test a.n_relations == 4 && a.n_distinct == 3
        @test a.habitat_type == "Terrestrial"
    end
end

@testset "fetch_relations — a failed fetch is not 'no habitat'" begin
    calls = Ref(0)
    function serve(status, body)
        return function (url; kwargs...)
            calls[] += 1
            return (status = status, body = body)
        end
    end
    function boom(url; kwargs...)
        calls[] += 1
        error("connection reset by peer")
    end

    @test length(EnvironmentTools.fetch_relations(560819; getter = serve(200, REL_JSON))) == 1
    # The two cases that must not collapse into one value.
    @test EnvironmentTools.fetch_relations(560819; getter = serve(200, "[]")) ==
          EnvironmentTools.Relation[]
    @test quiet(() -> EnvironmentTools.fetch_relations(560819; getter = serve(500, ""))) === nothing
    @test quiet(() -> EnvironmentTools.fetch_relations(560819; getter = boom)) === nothing

    # RETRY, AND ONLY WHERE RETRYING MEANS ANYTHING. In a fetch over thousands of taxa one
    # blip otherwise lands its taxon in fetch_failed and trips a tolerance-zero gate after the
    # whole network spend, so the long runs ask again -- but only on a transient failure.
    # retry_delay = 0 keeps the test instant; the default of 0 retries is what holds every
    # call site that does not ask for retries to a single request.
    calls[] = 0
    @test quiet(() -> EnvironmentTools.fetch_relations(560819; getter = serve(503, ""))) ===
          nothing
    @test calls[] == 1                      # no retries asked for, one request made

    calls[] = 0
    @test quiet(() -> EnvironmentTools.fetch_relations(
        560819; getter = serve(503, ""), retries = 2, retry_delay = 0)) === nothing
    @test calls[] == 3                      # the first attempt plus two retries

    # A request that raised before any status came back is transient too.
    calls[] = 0
    @test quiet(() -> EnvironmentTools.fetch_relations(
        560819; getter = boom, retries = 1, retry_delay = 0)) === nothing
    @test calls[] == 2

    # A 404 is an answer about this request; asking again only repeats it.
    calls[] = 0
    @test quiet(() -> EnvironmentTools.fetch_relations(
        560819; getter = serve(404, ""), retries = 2, retry_delay = 0)) === nothing
    @test calls[] == 1

    # A transient failure followed by a good response returns the good response.
    calls[] = 0
    flaky = function (url; kwargs...)
        calls[] += 1
        return calls[] == 1 ? (status = 502, body = "") : (status = 200, body = REL_JSON)
    end
    @test length(quiet(() -> EnvironmentTools.fetch_relations(
        560819; getter = flaky, retries = 3, retry_delay = 0))) == 1
    @test calls[] == 2

    # The request names its client in the User-Agent header and sets no other header of its
    # own, so the name is all it says about who is asking.
    sent = Ref{Any}(nothing)
    sniff = function (url; headers = Pair{String,String}[], kwargs...)
        sent[] = headers
        return (status = 200, body = REL_JSON)
    end
    EnvironmentTools.fetch_relations(560819; getter = sniff, user_agent = "probe/0")
    @test sent[] == ["User-Agent" => "probe/0"]

    # A 200 whose body is not a relation list is a failure too, not an empty habitat.
    @test quiet(() -> EnvironmentTools.fetch_relations(
        560819; getter = serve(200, """{"detail": "nope"}"""))) === nothing
    # …and the failure is LOUD. A silent `nothing` is as easy to ignore as a silent empty
    # vector: either way a lost fetch reads downstream as a taxon with no habitat.
    @test_logs (:warn,) match_mode = :any EnvironmentTools.fetch_relations(
        560819; getter = serve(500, ""))

    @testset "aggregate_habitats refuses a failed fetch instead of blanking it" begin
        @test_throws ArgumentError aggregate_habitats(nothing, path_of)
    end
end

@testset "fetch_relations — the on-disk cache is authoritative and raw" begin
    mktempdir() do dir
        calls = Ref(0)
        function counting(url; kwargs...)
            calls[] += 1
            return (status = 200, body = REL_JSON)
        end

        a = EnvironmentTools.fetch_relations(560819; cache_dir = dir, getter = counting)
        @test calls[] == 1
        @test isfile(joinpath(dir, "560819.json"))
        # The RAW body is stored, not the parsed rows: a later parser fix can be replayed
        # over the same bytes without refetching every taxon.
        @test read(joinpath(dir, "560819.json"), String) == REL_JSON
        @test !isfile(joinpath(dir, "560819.json.tmp"))     # write-then-rename, no litter

        b = EnvironmentTools.fetch_relations(560819; cache_dir = dir, getter = counting)
        @test calls[] == 1                                   # served from disk, not refetched
        @test a == b

        # A failed fetch must not poison the cache with an "empty" answer.
        function dead(url; kwargs...)
            return (status = 500, body = "")
        end
        @test quiet(() -> EnvironmentTools.fetch_relations(999; cache_dir = dir,
                                                           getter = dead)) === nothing
        @test !isfile(joinpath(dir, "999.json"))
        @test !isfile(joinpath(dir, "999.json.tmp"))

        # A corrupt cache entry (a half-written file from a killed run) is refetched,
        # rather than parsed as a taxon with no habitat.
        write(joinpath(dir, "7.json"), """[{"obtid": "soil", "obt_o""")
        r = quiet(() -> EnvironmentTools.fetch_relations(7; cache_dir = dir, getter = counting))
        @test calls[] == 2 && length(r) == 1
    end
end

@testset "OBT path cache: one shared loader/writer, and a strict offline mode" begin
    saved = copy(EnvironmentTools.OBT_CACHE)
    try
        empty!(EnvironmentTools.OBT_CACHE)
        EnvironmentTools.OBT_CACHE["soil"] = REAL["soil"]

        @test EnvironmentTools.get_obt_path("soil"; strict = true) == REAL["soil"]
        @test habitat_type(EnvironmentTools.get_obt_path("soil"; strict = true)) == "Terrestrial"
        # strict = true must RAISE on an uncached id. Returning "" (what the network path
        # does when the API is unreachable) types the term as "Other", which is
        # indistinguishable from a real answer in every downstream column.
        @test_throws ErrorException EnvironmentTools.get_obt_path("OBT:999999"; strict = true)
        @test habitat_type("") == "Other"          # why the silent "" is dangerous

        mktempdir() do dir
            csv = joinpath(dir, "obt_path_cache.csv")
            @test EnvironmentTools.load_obt_cache!(csv) == 0     # a cold start is not an error

            EnvironmentTools.OBT_CACHE["nutrient broth"] = REAL["nutrient broth"]
            EnvironmentTools.OBT_CACHE["OBT:000000"] = ""        # a term with no path at all
            @test EnvironmentTools.save_obt_cache(csv) == 3

            empty!(EnvironmentTools.OBT_CACHE)
            @test EnvironmentTools.load_obt_cache!(csv) == 3
            @test EnvironmentTools.OBT_CACHE["soil"] == REAL["soil"]
            @test EnvironmentTools.OBT_CACHE["nutrient broth"] == REAL["nutrient broth"]
            # An empty path round-trips as "" and not as `missing`, or the strict-mode
            # check above would throw on a term that IS cached.
            @test EnvironmentTools.OBT_CACHE["OBT:000000"] == ""
            @test EnvironmentTools.get_obt_path("OBT:000000"; strict = true) == ""

            # Sorted by id, so two runs that reached the same terms write the same file.
            ids = [first(split(l, ',')) for l in readlines(csv)[2:end]]
            @test ids == sort(ids)
        end
    finally
        empty!(EnvironmentTools.OBT_CACHE)
        merge!(EnvironmentTools.OBT_CACHE, saved)
    end
end

@testset "fetch_omnicrobe_env! runs offline from the relation cache" begin
    saved = copy(EnvironmentTools.OBT_CACHE)
    try
        for t in ("soil", "terrestrial habitat", "lake")
            EnvironmentTools.OBT_CACHE[t] = REAL[t]
        end
        mktempdir() do dir
            write(joinpath(dir, "42.json"), REL_MIX_JSON)
            write(joinpath(dir, "7.json"), "[]")          # genuinely unannotated taxon
            df = DataFrame(TaxID = [42, 7])
            fetch_omnicrobe_env!(df; cache_dir = dir)     # no network: every taxid is cached

            @test df.habitat_type[1] == "Terrestrial"
            @test df.habitat_specific[1] == "soil"
            @test df.habitat_n_relations[1] == 4 && df.habitat_n_distinct[1] == 3
            @test ismissing(df.habitat_category[2]) && df.habitat_n_relations[2] == 0
        end
    finally
        empty!(EnvironmentTools.OBT_CACHE)
        merge!(EnvironmentTools.OBT_CACHE, saved)
    end
end

# ---------------------------------------------------------------------------
# The tie rule and the term filters, driven end to end through aggregate_habitats.
#
# `resolve` is unit-tested above on bare tallies. These eight groups put the same rules
# through the aggregator, which is where each failure becomes reachable: a winner that
# depends on hash layout, a term cap whose victims depend on it too, isolation sources
# counted as habitats, two ontology nodes meaning "life in general" outvoting the one
# informative term, and a chemistry axis that can name only ONE condition per taxon while
# the flags are meant to be non-exclusive.
#
# All of them run on the REAL paths above. A hand-invented ontology path would agree with an
# inverted table, which is the standard the header sets.
# ---------------------------------------------------------------------------

# A relation carrying surface forms, which `rel` (obtid + name only) cannot express. The
# chemistry keyword fallback reads that text, so the forms are the difference between
# finding a condition and losing it.
relf(nm, forms...) = (obtid = nm, name = nm, forms = String[forms...])

@testset "1. TIES are split and labelled, never awarded" begin
    soil_lake = [rel("soil"), rel("lake")]

    f = aggregate_habitats(soil_lake, path_of; tie_rule = :fractional)
    @test f.type_weights == Dict("Terrestrial" => 0.5, "Aquatic" => 0.5)
    @test f.n_tied == 2
    @test f.resolved_by == :fractional
    @test f.habitat_type == "Aquatic|Terrestrial"
    # The label is a value no single category can take, so downstream code cannot read a
    # coin toss as a call.
    @test !(f.habitat_type in ("Aquatic", "Terrestrial"))

    # Mentions break the tie before the split does: the same two terms, one asserted twice.
    m = aggregate_habitats([rel("soil"), rel("soil"), rel("lake")], path_of;
                           tie_rule = :mentions_fractional)
    @test m.habitat_type == "Terrestrial"
    @test m.resolved_by == :mentions
    @test m.type_weights == Dict("Terrestrial" => 1.0)
    @test m.n_tied == 1                       # nothing was left for the split rule to do

    # and with nothing to break it, the primary rule falls back to the same 0.5/0.5.
    p = aggregate_habitats(soil_lake, path_of; tie_rule = :mentions_fractional)
    @test p.type_weights == Dict("Terrestrial" => 0.5, "Aquatic" => 0.5)
    @test p.resolved_by == :fractional
    @test p.habitat_type == "Aquatic|Terrestrial"

    a = aggregate_habitats(soil_lake, path_of; tie_rule = :alphabetical)
    @test a.habitat_type == "Aquatic"
    @test a.type_weights == Dict("Aquatic" => 1.0)
    @test a.n_tied == 2 && a.resolved_by == :alphabetical   # the tie is still REPORTED

    # `:legacy_dict` is the rule all of this is varied against: one label, full weight, and
    # no column saying that a coin was tossed. WHICH label it is depends on hash layout, so
    # only the shape is asserted here.
    l = aggregate_habitats(soil_lake, path_of; tie_rule = :legacy_dict)
    @test length(l.type_weights) == 1
    @test collect(values(l.type_weights)) == [1.0]
    @test l.habitat_type in ("Aquatic", "Terrestrial")
    @test l.n_tied == 2 && l.resolved_by == :none
end

#=
    LABEL NEUTRALITY: the test that catches a leftover first-max winner.

    A permutation test cannot do it. Julia's `Dict` iterates in HASH order, not insertion
    order, so shuffling the relations leaves the tally's iteration order untouched: a
    first-max winner answers every shuffle identically and looks deterministic. What it is
    not is NEUTRAL. The label it picks is a property of the two category strings, so the
    only question that exposes it has to be asked of every PAIR of labels: does a 1-vs-1
    tie come back as 0.5/0.5, whatever the two categories are called?

    Isolation and generic filtering are switched off here so that the pair under test is
    the whole evidence, and the question reaches the winner rule and nothing else.
=#
@testset "2. LABEL-NEUTRALITY across every pair of fine types with a real path" begin
    # One representative term per fine type, taken from the REAL paths above: whichever
    # term sorts first takes the slot, and any other term of the same type would serve.
    reps = Dict{String,String}()
    for nm in sort(collect(keys(REAL)))
        ft = habitat_type(REAL[nm])
        ft == "Other" && continue
        get!(reps, ft, nm)
    end
    types = sort(collect(keys(reps)))
    @test length(types) >= 11
    @test issubset(["Aquatic", "Terrestrial", "Geological", "Atmospheric", "Host-plant",
                    "Host-animal", "Host-other", "Agricultural", "Anthropogenic",
                    "Food/substance", "Lab-medium"], types)

    npairs = 0
    for i in 1:length(types), j in (i + 1):length(types)
        t1, t2 = types[i], types[j]
        npairs += 1
        want = Dict(t1 => 0.5, t2 => 0.5)
        lab  = join(sort([t1, t2]), "|")
        kw   = (isolation = :taxon, drop_generic = false)
        a = aggregate_habitats([rel(reps[t1]), rel(reps[t2])], path_of; kw...)
        b = aggregate_habitats([rel(reps[t2]), rel(reps[t1])], path_of; kw...)
        @test a.type_weights == want
        @test a.habitat_type == lab
        @test a.n_tied == 2
        @test a.resolved_by == :fractional
        @test isequal(a, b)
    end
    @test npairs >= 55
end

@testset "3. PERMUTATION: the call never depends on the input order" begin
    mix = [rel("soil"), rel("soil"), rel("lake"), rel("rock"), rel("gut"),
           rel("hypersaline lake"), rel("nutrient broth")]
    @test length(mix) == 7
    ref = aggregate_habitats(mix, path_of)
    # A DECIDED call, not a tie: two distinct aquatic terms against one of everything else,
    # so this group tests order-independence rather than the split rule a second time.
    @test ref.habitat_type == "Aquatic" && ref.n_tied == 1
    outs = [aggregate_habitats(shuffle(MersenneTwister(s), mix), path_of) for s in 1:50]
    @test length(outs) == 50
    @test all(o -> isequal(o, ref), outs)
    # Named separately so a failure says WHICH field moved.
    @test length(unique(o -> o.habitat_type, outs)) == 1
    @test length(unique(o -> o.name, outs)) == 1
    @test length(unique(o -> o.type_counts, outs)) == 1
    @test length(unique(o -> o.profile, outs)) == 1
    @test length(unique(o -> o.chemistry_flags, outs)) == 1
end

@testset "4. TRUNCATION drops the same terms whatever the order" begin
    # Equal mentions, so only the second sort key can decide: the obtids sort
    # lake < rock < soil, and the cap must therefore always spend its two slots on the
    # first two, whichever order they arrive in.
    terms  = ["lake", "rock", "soil"]
    orders = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
    for o in orders
        a = aggregate_habitats([rel(terms[i]) for i in o], path_of; max_distinct = 2)
        @test a.n_distinct == 3 && a.n_truncated == 1
        @test a.type_counts == "Aquatic:1;Geological:1"
        @test !occursin("Terrestrial", a.type_profile)      # `soil` is the one dropped
        @test a.habitat_type == "Aquatic|Geological"
    end
end

@testset "5. ISOLATION SOURCES leave the type axis term by term" begin
    a = aggregate_habitats([rel("salt"), rel("nutrient broth"), rel("soil")], path_of)
    @test a.habitat_type == "Terrestrial"
    @test a.n_isolation_terms == 2
    @test a.isolation_only == false
    @test a.type_counts == "Terrestrial:1"

    # Every term this taxon has records where a strain was OBTAINED. The call is blanked,
    b = aggregate_habitats([rel("salt"), rel("nutrient broth")], path_of)
    @test ismissing(b.habitat_type) && ismissing(b.category)
    @test b.isolation_only == true
    @test b.n_relations == 2      # and stays distinguishable from an unannotated taxon.

    # `salt` keeps speaking on the CHEMISTRY axis: OBT files it under food, but a habitat
    # statement of "salt" is a salinity statement.
    c = aggregate_habitats([rel("salt"), rel("soil")], path_of)
    @test occursin("Hypersaline", c.chemistry_flags)
    @test c.habitat_type == "Terrestrial"
    @test c.n_chem_terms == 1

    # A culture flask is neither a habitat nor a condition, so Lab-medium leaves BOTH axes
    # and its text cannot donate a flag.
    broth = [relf("nutrient broth", "alkaline minimal medium"), rel("soil")]
    d = aggregate_habitats(broth, path_of)
    @test d.chemistry_flags == ""
    @test d.n_chem_terms == 0
    @test d.habitat_type == "Terrestrial"
    # The filter is what removes it: the same text under the legacy taxon-level mode does
    # flag Alkaline, which is the medium's recipe and not the organism's environment.
    @test aggregate_habitats(broth, path_of; isolation = :taxon).chemistry_flags == "Alkaline"

    # Legacy taxon-level behaviour: nothing is filtered, so an isolation source can win.
    salty = [rel("salt"), rel("salt"), rel("soil")]
    g = aggregate_habitats(salty, path_of; isolation = :taxon)
    @test g.habitat_type == "Food/substance"
    @test g.resolved_by == :mentions
    @test g.n_isolation_terms == 1                     # counted under either mode
    @test aggregate_habitats(salty, path_of).habitat_type == "Terrestrial"
end

@testset "6. GENERIC TERMS are dropped, and the drop is counted" begin
    rels = [rel("soil"), rel("living organism"), rel("cell")]
    a = aggregate_habitats(rels, path_of)                       # drop_generic = true
    @test a.habitat_type == "Terrestrial"
    @test a.n_generic_terms == 2
    @test a.generic_only == false
    @test a.name == "soil"

    # Left in, the two uninformative nodes outvote the one informative term 2-1, which is
    # how the Host-other bucket fills up with taxa that say nothing about a host.
    b = aggregate_habitats(rels, path_of; drop_generic = false)
    @test b.habitat_type == "Host-other"
    @test b.type_counts == "Host-other:2;Terrestrial:1"
    @test b.n_generic_terms == 2                 # counted whether or not the filter acts

    c = aggregate_habitats([rel("cell")], path_of)
    @test ismissing(c.habitat_type) && ismissing(c.category)
    @test c.generic_only == true
    @test c.n_relations == 1                     # a blanked call, not an unannotated taxon
    @test aggregate_habitats([rel("cell")], path_of; drop_generic = false).habitat_type ==
          "Host-other"
    # The genuinely unannotated taxon must stay distinguishable from both blanked kinds.
    none = aggregate_habitats(EnvironmentTools.Relation[], path_of)
    @test none.generic_only == false && none.isolation_only == false
end

@testset "7. CHEMISTRY is non-exclusive" begin
    a = aggregate_habitats([rel("hypersaline lake"), rel("hot spring"), rel("soil")], path_of)
    @test a.chemistry_flags == "Hypersaline+Thermal"
    @test a.chemistry_flag_counts == "Hypersaline:1;Thermal:1"
    @test a.n_chem_terms == 2
    @test a.habitat_type == "Aquatic"            # the type axis is untouched by either flag
    # The back-compatible single-winner field CANNOT say "both": it reports a tie between
    # the two conditions, which is the expressiveness the union of flags exists to restore.
    @test a.chemistry == "Hypersaline|Thermal"
    @test a.chemistry_coverage ≈ 2 / 3           # two of the three kept terms said anything

    # ONE term, TWO flags: a combination is split into its parts before it is counted.
    soda = relf("soda lake", "soda lake", "hypersaline soda lake")
    b = aggregate_habitats([soda], path_of)
    @test b.chemistry_flags == "Alkaline+Hypersaline"
    @test b.chemistry_flag_counts == "Alkaline:1;Hypersaline:1"
    @test b.n_chem_terms == 1                    # terms, not flags

    c = aggregate_habitats([soda, rel("hypersaline lake"), rel("soil")], path_of)
    @test c.chemistry_flags == "Alkaline+Hypersaline"
    @test c.chemistry_flag_counts == "Hypersaline:2;Alkaline:1"
    @test c.n_chem_terms == 2
end

@testset "8. CHEM TEXT decides whether the condition is found at all" begin
    # The shape of the real API item in REL_JSON above: the ontology term is `spring` and
    # the salinity is written down only in a surface form.
    spring = (obtid = "OBT:001333", name = "spring",
              forms = ["spring", "mesophilic saline spring"], source = "GenBank")
    @test path_of(spring.obtid) == ""      # nothing can arrive from the ontology here

    u = aggregate_habitats([spring], path_of; chem_text = :forms_union)
    @test u.chemistry_flags == "Hypersaline"
    @test u.chemistry_coverage == 1.0

    o = aggregate_habitats([spring], path_of; chem_text = :object_name)
    @test o.chemistry_flags == ""
    @test o.chemistry_coverage == 0.0

    # `:first_form` is the surface-form text — forms[1] taken as the name — and it loses
    # the flag the moment the bare form is listed first.
    @test aggregate_habitats([spring], path_of; chem_text = :first_form).chemistry_flags == ""

    # The real payload behaves the same way through the parser.
    r = EnvironmentTools.parse_relations(REL_JSON)
    @test aggregate_habitats(r, path_of; chem_text = :forms_union).chemistry_flags ==
          "Hypersaline"
    @test aggregate_habitats(r, path_of; chem_text = :object_name).chemistry_flags == ""
end

# ---------------------------------------------------------------------------
# 9–13: the rules a one-at-a-time sensitivity run needs, which the aggregator could not
# express. Three of them decide the TYPE axis (leave a surviving tie uncalled; label every
# category holding a given share; label every category present at all) and two decide the
# CHEMISTRY axis (the same share test over flags; the exclusive winner with its tie split
# across the tied flag combinations).
#
# What each group has to establish, because each is a way the rules have failed elsewhere:
#
#   • the rule settles the case it was added for, on evidence small enough to check by hand,
#     INCLUDING the tie it exists to settle — a rule that only ever meets decided calls is a
#     rule that cannot be wrong;
#   • the threshold comparison sits exactly on the cut, tested from both sides with the same
#     two categories, so a `>` written for a `>=` moves an assertion;
#   • what one taxon's weights sum to under the rule, since the estimators take Real weights
#     per category and would read extra weight as extra hits with nothing to notice it;
#   • a type rule moves the type axis ONLY, or a variant meant to be one move off the
#     specification would quietly be two;
#   • an unknown rule name errors and names the valid rules, instead of falling through to
#     whichever branch an `else` happens to hold;
#   • and the rules that existed answer on these same fixtures exactly as they did before any
#     of this was added, which is the whole neutrality claim in one testset.
#
# All of it runs on the REAL ontology paths above, for the reason the header records.
# ---------------------------------------------------------------------------

# The sensitivity fixtures, each one relation set whose tally can be read off by eye.
#   S_TIE2     Terrestrial 1, Aquatic 1, equal mentions — the tie no rule can break on evidence
#   S_MENBREAK the same tally, with `soil` asserted twice — the tie the mention step DOES break
#   S_QUARTER  Aquatic 3, Terrestrial 1 — Terrestrial sits exactly on the default share
#   S_FIFTH    one Aquatic term more — the same Terrestrial term now sits just under it
#   S_SPREAD5  five categories of one term each — nothing reaches the share, nothing is tied off
#   S_HOTBRINE two conditions on two terms, so the chemistry WINNER is itself a tie
#   S_COMBO4   `Alkaline+Hypersaline` on 1 of 4 chemistry terms — the flag sits on the cut
#   S_COMBO5   one Hypersaline term more, and the same flag sits under it
const S_SODA     = relf("soda lake", "soda lake", "hypersaline soda lake")
const S_TIE2     = [rel("soil"), rel("lake")]
const S_MENBREAK = [rel("soil"), rel("soil"), rel("lake")]
const S_QUARTER  = [rel("lake"), rel("marine environment"), rel("hypersaline lake"), rel("soil")]
const S_FIFTH    = [rel("lake"), rel("marine environment"), rel("hypersaline lake"),
                    rel("saline lake"), rel("soil")]
const S_SPREAD5  = [rel("soil"), rel("lake"), rel("rock"), rel("gut"), rel("rhizosphere")]
const S_HOTBRINE = [rel("hypersaline lake"), rel("hot spring"), rel("soil")]
const S_COMBO4   = [S_SODA, rel("hypersaline lake"), rel("saline lake"), rel("solar saltern")]
const S_COMBO5   = [S_SODA, rel("hypersaline lake"), rel("saline lake"), rel("solar saltern"),
                    rel("salt")]
const S_FIXTURES = ["tie2" => S_TIE2, "menbreak" => S_MENBREAK, "quarter" => S_QUARTER,
                    "fifth" => S_FIFTH, "spread5" => S_SPREAD5, "hotbrine" => S_HOTBRINE,
                    "combo4" => S_COMBO4, "combo5" => S_COMBO5]

# "Cat:3;Cat:1" back into a tally, so the weight-sum checks below are driven by the raw count
# string rather than by the same numbers the rule under test worked from.
s_tally(s) = Dict(String(strip(p[1:prevind(p, findlast(==(':'), p))])) =>
                  parse(Int, strip(p[nextind(p, findlast(==(':'), p)):end]))
                  for p in split(s, ';') if !isempty(strip(p)))
s_flags(s) = isempty(s) ? String[] : sort(String.(split(s, '+')))

@testset "9. THE TYPE-AXIS SENSITIVITY RULES, on bare tallies" begin
    R = EnvironmentTools.resolve
    U = EnvironmentTools.UNRESOLVED
    tie     = Dict("Aquatic" => 1, "Terrestrial" => 1)
    men     = Dict("Aquatic" => 3, "Terrestrial" => 1)
    quarter = Dict("Aquatic" => 3, "Terrestrial" => 1)
    fifth   = Dict("Aquatic" => 4, "Terrestrial" => 1)
    spread  = Dict(c => 1 for c in ("Aquatic", "Terrestrial", "Geological",
                                    "Host-animal", "Host-plant"))
    one     = Dict("Aquatic" => 3)

    # The uncalled label must be a value no category can take, or a taxon the rule declined to
    # call would be scored as a habitat called "Unresolved".
    @test !(U in first.(EnvironmentTools.OBT_TYPE))
    @test !(U in values(EnvironmentTools.COARSE))
    # The registered cut. It is a parameter of the rule; a silent change here would move every
    # share row in the sensitivity run at once.
    @test EnvironmentTools.SHARE_THRESHOLD == 0.25

    @testset ":unresolved leaves a surviving tie uncalled, and keeps the taxon" begin
        # The tie it exists for: equal terms AND equal mentions, so the mention step cannot
        # settle it. No weight anywhere, and a label that is not a category.
        @test R(tie, tie, :unresolved) == (U, Dict{String,Float64}(), 2, :unresolved)
        # …but it is the PRIMARY rule minus the split, not minus the mention step as well: the
        # tie the mention count does settle must still be settled, or this would be two moves
        # off the primary rule rather than one.
        @test R(tie, men, :unresolved) == ("Aquatic", Dict("Aquatic" => 1.0), 1, :mentions)
        @test R(quarter, quarter, :unresolved) == ("Aquatic", Dict("Aquatic" => 1.0), 1, :none)
        @test R(Dict{String,Int}(), Dict{String,Int}(), :unresolved) ==
              ("", Dict{String,Float64}(), 0, :none)
    end

    @testset ":share labels every category AT the cut, and the cut is a parameter" begin
        @test R(quarter, quarter, :share) ==
              ("Aquatic+Terrestrial", Dict("Aquatic" => 1.0, "Terrestrial" => 1.0), 1, :share)
        # One term more on the other side and the same category drops out: 1 of 4 is in, 1 of
        # 5 is out. A `>` where the rule says `>=` fails the first of these.
        @test R(fifth, fifth, :share) == ("Aquatic", Dict("Aquatic" => 1.0), 1, :none)
        @test sort(collect(keys(R(fifth, fifth, :share; share_threshold = 0.2)[2]))) ==
              ["Aquatic", "Terrestrial"]
        @test collect(keys(R(quarter, quarter, :share; share_threshold = 0.5)[2])) == ["Aquatic"]
        @test isequal(R(quarter, quarter, :share),
                      R(quarter, quarter, :share;
                        share_threshold = EnvironmentTools.SHARE_THRESHOLD))
        # Five categories of one term each: nothing reaches a quarter, so the rule labels
        # NOTHING. The taxon is uncalled rather than dropped — it still carries a label, and it
        # still counts as one taxon of its arm, scoring zero everywhere.
        @test R(spread, spread, :share) == (U, Dict{String,Float64}(), 5, :share)
        # Multi-label is not "the tied set": nothing is tied in `quarter` at all.
        @test R(quarter, quarter, :share)[3] == 1
        # The mention count decides nothing here, which is what makes this one axis of change.
        @test isequal(R(quarter, Dict("Terrestrial" => 99), :share), R(quarter, quarter, :share))
    end

    @testset ":any_type labels every category present, tied or trailing" begin
        @test R(fifth, fifth, :any_type) ==
              ("Aquatic+Terrestrial", Dict("Aquatic" => 1.0, "Terrestrial" => 1.0), 1, :any_type)
        @test R(one, one, :any_type) == ("Aquatic", Dict("Aquatic" => 1.0), 1, :none)
        @test length(R(spread, spread, :any_type)[2]) == 5
        # The two multi-label rules differ exactly where the share bites: `:any_type` is the
        # share rule with no cut, and taking the tied set instead (a rule this module does
        # not offer) would answer "Aquatic" for both of these.
        @test issubset(keys(R(fifth, fifth, :share)[2]), keys(R(fifth, fifth, :any_type)[2]))
        @test keys(R(fifth, fifth, :share)[2]) != keys(R(fifth, fifth, :any_type)[2])
        @test R(spread, spread, :any_type)[1] != U
    end
end

@testset "10. THE TYPE-AXIS RULES through the aggregator, on real ontology paths" begin
    U = EnvironmentTools.UNRESOLVED

    u = aggregate_habitats(S_TIE2, path_of; tie_rule = :unresolved)
    @test u.habitat_type == U && u.category == U
    @test isempty(u.type_weights) && u.n_tied == 2 && u.resolved_by == :unresolved
    # Uncalled is not blank. A caller's usable-mask reads `habitat_type`, so a blank would
    # take the taxon out of the arm — and a rule that changes only the assignment would be
    # quietly changing the denominator too.
    @test !ismissing(u.habitat_type)
    @test !ismissing(aggregate_habitats(S_SPREAD5, path_of; tie_rule = :share).habitat_type)
    # No category was called, so no term may be put forward as representing one.
    @test ismissing(u.name)
    @test u.type_confidence == 0.0
    # The evidence itself is untouched: only the call changes.
    @test u.type_counts == aggregate_habitats(S_TIE2, path_of).type_counts

    m = aggregate_habitats(S_MENBREAK, path_of; tie_rule = :unresolved)
    @test m.habitat_type == "Terrestrial" && m.resolved_by == :mentions
    @test m.type_weights == Dict("Terrestrial" => 1.0)

    q = aggregate_habitats(S_QUARTER, path_of; tie_rule = :share)
    @test q.type_counts == "Aquatic:3;Terrestrial:1"          # 1 of 4 = the cut exactly
    @test q.habitat_type == "Aquatic+Terrestrial"
    @test q.type_weights == Dict("Aquatic" => 1.0, "Terrestrial" => 1.0)
    f = aggregate_habitats(S_FIFTH, path_of; tie_rule = :share)
    @test f.type_counts == "Aquatic:4;Terrestrial:1"          # …and 1 of 5 is under it
    @test f.habitat_type == "Aquatic" && f.type_weights == Dict("Aquatic" => 1.0)
    @test f.resolved_by == :none

    s5 = aggregate_habitats(S_SPREAD5, path_of; tie_rule = :share)
    @test s5.habitat_type == U && isempty(s5.type_weights) && s5.resolved_by == :share
    a5 = aggregate_habitats(S_SPREAD5, path_of; tie_rule = :any_type)
    @test length(a5.type_weights) == 5 && all(==(1.0), values(a5.type_weights))
    @test a5.habitat_type ==
          "Aquatic+Geological+Host-animal+Host-plant+Terrestrial"
    af = aggregate_habitats(S_FIFTH, path_of; tie_rule = :any_type)
    @test af.type_weights == Dict("Aquatic" => 1.0, "Terrestrial" => 1.0)   # cf. `f` above

    # A multi-label call is labelled with "+", the separator the chemistry axis uses for
    # conditions that hold together — never with "|", which downstream reads as a split tie.
    @test !occursin("|", af.habitat_type) && occursin("+", af.habitat_type)

    @testset "a type rule moves the type axis and nothing else" begin
        base = aggregate_habitats(S_HOTBRINE, path_of)
        # The chemistry winner on this fixture is itself a tie, so a type rule leaking into the
        # chemistry axis would show here instead of being hidden by a clear winner.
        @test base.chemistry == "Hypersaline|Thermal"
        for r in EnvironmentTools.TYPE_ONLY_RULES
            a = aggregate_habitats(S_HOTBRINE, path_of; tie_rule = r)
            @test a.chemistry == base.chemistry
            @test a.chemistry_flags == base.chemistry_flags
            @test a.chemistry_weights == base.chemistry_weights
            @test a.chemistry_profile == base.chemistry_profile
            @test a.chemistry_share == base.chemistry_share
            @test a.chemistry_coverage == base.chemistry_coverage
            @test a.n_chem_terms == base.n_chem_terms
            @test a.type_counts == base.type_counts
            @test a.type_mentions == base.type_mentions
            @test a.n_relations == base.n_relations && a.n_distinct == base.n_distinct
        end
    end
end

@testset "11. THE CHEMISTRY RULES: which conditions a taxon is scored for" begin
    C = EnvironmentTools.resolve_chemistry
    # The tally is over flag COMBINATIONS, keyed as `habitat_chemistry` writes them.
    tie  = Dict("Thermal" => 1, "Alkaline+Thermal" => 1)
    won  = Dict("Hypersaline" => 3, "Alkaline+Hypersaline" => 1)
    won5 = Dict("Hypersaline" => 4, "Alkaline+Hypersaline" => 1)
    five = Dict(f => 1 for f in ("Thermal", "Hypersaline", "Alkaline", "Anoxic", "Cold"))

    @testset ":any_term — one weight per flag any term carries" begin
        @test C(tie, :any_term) == Dict("Thermal" => 1.0, "Alkaline" => 1.0)
        @test C(won, :any_term) == Dict("Hypersaline" => 1.0, "Alkaline" => 1.0)
        @test C(Dict{String,Int}(), :any_term) == Dict{String,Float64}()
        # A term carrying no condition votes for nothing rather than for a blank flag.
        @test C(Dict("" => 3), :any_term) == Dict{String,Float64}()
    end

    @testset ":share — a combination counts toward EVERY flag it names" begin
        # Alkaline is named by 1 of the 4 chemistry-bearing terms, exactly the cut, and it is
        # named only inside a combination: a rule matching whole keys would score it 0.
        @test C(won, :share) == Dict("Hypersaline" => 1.0, "Alkaline" => 1.0)
        # One Hypersaline term more and the same flag is 1 of 5, under the cut.
        @test C(won5, :share) == Dict("Hypersaline" => 1.0)
        @test isempty(C(five, :share))                       # nothing reaches a quarter
        @test length(C(five, :share; share_threshold = 0.2)) == 5
        @test C(won5, :share; share_threshold = 0.2) == C(won5, :any_term)
    end

    @testset ":winner_fractional — split over the tied COMBINATIONS, not over flags" begin
        # Tied between `Thermal` and `Alkaline+Thermal`: Thermal is named by both tied
        # combinations and Alkaline by one, so the weights are 1 and 1/2 — not 1/2 each.
        @test C(tie, :winner_fractional) == Dict("Thermal" => 1.0, "Alkaline" => 0.5)
        # It is the EXCLUSIVE winner: 3 terms beat 1, and the trailing flag scores nothing,
        # which is what separates it from `:any_term` on the same tally.
        @test C(won, :winner_fractional) == Dict("Hypersaline" => 1.0)
        @test C(won, :any_term) != C(won, :winner_fractional)
        # A flag every tied combination names is a WHOLE hit, exactly: seven additions of 1/7
        # come to 0.9999999999999998, and a downstream parser that narrows an all-0/1 column to
        # a Bool vector would silently stop recognising it.
        seven = Dict("Thermal+F$i" => 1 for i in 1:7)
        w7 = C(seven, :winner_fractional)
        @test w7["Thermal"] === 1.0
        @test all(w7["F$i"] == 1 / 7 for i in 1:7)
        @test isempty(C(Dict{String,Int}(), :winner_fractional))
    end

    @testset "the chemistry rule decides the weights and touches nothing else" begin
        for (nm, fx) in S_FIXTURES
            base = aggregate_habitats(fx, path_of)
            for cr in EnvironmentTools.CHEMISTRY_RULES
                a = aggregate_habitats(fx, path_of; chemistry_rule = cr)
                @test isequal(Base.structdiff(a, (chemistry_weights = nothing,)),
                              Base.structdiff(base, (chemistry_weights = nothing,)))
            end
            # Under the primary rule the weights ARE the flag column the pipeline already
            # reads, so switching the pipeline onto the weights changes nothing by itself.
            @test sort(collect(keys(base.chemistry_weights))) == s_flags(base.chemistry_flags)
        end
        # Through the aggregator on real paths: one term names two conditions and one names
        # one, so the exclusive winner keeps only what the winning combination names.
        c4 = aggregate_habitats(S_COMBO4, path_of; chemistry_rule = :winner_fractional)
        @test c4.chemistry_profile == "Hypersaline:0.75;Alkaline+Hypersaline:0.25"
        @test c4.chemistry_weights == Dict("Hypersaline" => 1.0)
        @test aggregate_habitats(S_COMBO4, path_of; chemistry_rule = :share).chemistry_weights ==
              Dict("Hypersaline" => 1.0, "Alkaline" => 1.0)          # Alkaline is 1 of 4
        @test aggregate_habitats(S_COMBO5, path_of; chemistry_rule = :share).chemistry_weights ==
              Dict("Hypersaline" => 1.0)                             # …and 1 of 5 is out
        hb = aggregate_habitats(S_HOTBRINE, path_of; chemistry_rule = :winner_fractional)
        @test hb.chemistry_weights == Dict("Hypersaline" => 0.5, "Thermal" => 0.5)
    end
end

@testset "12. WEIGHTS: what one taxon may carry under each rule" begin
    U = EnvironmentTools.UNRESOLVED
    T = EnvironmentTools.SHARE_THRESHOLD
    # Downstream every category is estimated in a table of its own, over a denominator that
    # counts TAXA. So a weight above 1 in any category would make its hit total exceed the
    # number of taxa scored, and a rule that handed a taxon a second unit of one category
    # would read as an extra hit with nothing in the pipeline to notice it. The sums below are
    # checked against the raw count strings, not against the tallies the rules worked from.
    for (nm, fx) in S_FIXTURES
        base  = aggregate_habitats(fx, path_of)
        tc    = s_tally(base.type_counts)
        total = sum(values(tc))
        for r in EnvironmentTools.TIE_RULES
            a = aggregate_habitats(fx, path_of; tie_rule = r)
            w = collect(values(a.type_weights))
            @test all(x -> 0 < x <= 1, w)
            @test !ismissing(a.habitat_type)          # every fixture has usable evidence
            s = sum(w; init = 0.0)
            if r === :share
                @test s == count(v -> v / total >= T, values(tc))
                @test s <= floor(1 / T)
                @test (s == 0) == (a.habitat_type == U)
            elseif r === :any_type
                @test s == length(tc)
            elseif r === :unresolved
                @test s == (a.habitat_type == U ? 0.0 : 1.0)
            else
                @test s ≈ 1.0                          # a tied taxon is still one taxon
            end
        end

        nchem = base.n_chem_terms
        fc    = isempty(base.chemistry_flag_counts) ? Dict{String,Int}() :
                s_tally(base.chemistry_flag_counts)
        for cr in EnvironmentTools.CHEMISTRY_RULES
            a = aggregate_habitats(fx, path_of; chemistry_rule = cr)
            w = collect(values(a.chemistry_weights))
            @test all(x -> 0 < x <= 1, w)
            s = sum(w; init = 0.0)
            if nchem == 0
                @test s == 0.0                         # no condition stated, nothing scored
            elseif cr === :any_term
                @test s == length(fc)
            elseif cr === :share
                @test s == count(v -> v / nchem >= T, values(fc))
            else
                # Per combination the weights sum to 1; split into the flags a combination
                # names, one taxon can carry a whole unit of each.
                @test s >= 1.0
                @test s <= length(fc)
            end
        end
    end

    @testset "an unknown rule name errors, and the error names the rules" begin
        for (kw, rules) in ((:tie_rule, EnvironmentTools.TIE_RULES),
                            (:chemistry_rule, EnvironmentTools.CHEMISTRY_RULES))
            e = try
                aggregate_habitats(S_TIE2, path_of; NamedTuple{(kw,)}((:menfrac,))...)
                nothing
            catch err
                err
            end
            @test e isa ErrorException
            @test occursin(String(kw), e.msg)
            for r in rules
                @test occursin(":" * String(r), e.msg)       # the list cannot drift from the check
            end
        end
        # `resolve` refuses it too. Accepting anything lets an unrecognised symbol fall
        # through to the fractional split, so a misspelt rule produces a plausible table under
        # a rule nobody asked for.
        @test_throws ErrorException EnvironmentTools.resolve(
            Dict("Aquatic" => 1), Dict("Aquatic" => 1), :menfrac)
        @test_throws ErrorException EnvironmentTools.resolve_chemistry(
            Dict("Thermal" => 1), :anyflag)
    end

    @testset "a share is a fraction, and is checked where it is typed" begin
        # 25 for 25 % clears nothing and 0 clears everything; either would be a different rule
        # running under the name of this one. Refused even when no share rule is in play, so
        # the mistake is heard at the call that made it.
        for bad in (0, 25, 1.5, -0.1)
            @test_throws ErrorException aggregate_habitats(S_QUARTER, path_of;
                                                           tie_rule = :share,
                                                           share_threshold = bad)
            @test_throws ErrorException aggregate_habitats(S_QUARTER, path_of;
                                                           share_threshold = bad)
        end
        @test aggregate_habitats(S_QUARTER, path_of; tie_rule = :share,
                                 share_threshold = 0.75).habitat_type == "Aquatic"
    end
end

@testset "13. NEUTRALITY: the rules that existed answer exactly as they did" begin
    # Transcribed from a run of the aggregator as it stood BEFORE any of the above was added,
    # on the very fixtures the new rules are tested on. If a new rule reached into a shared
    # path — the mention step, the tally, the label, the representative name — this is where it
    # shows, and the replication configuration that every downstream anchor rests on is exactly
    # these four rules.
    prior = [
        # fixture      rule                  habitat_type            type_weights                       n_tied by             chemistry
        ("tie2",     :mentions_fractional, "Aquatic|Terrestrial", ["Aquatic" => 0.5, "Terrestrial" => 0.5], 2, :fractional,   ""),
        ("tie2",     :fractional,          "Aquatic|Terrestrial", ["Aquatic" => 0.5, "Terrestrial" => 0.5], 2, :fractional,   ""),
        ("tie2",     :alphabetical,        "Aquatic",             ["Aquatic" => 1.0],                       2, :alphabetical, ""),
        ("menbreak", :mentions_fractional, "Terrestrial",         ["Terrestrial" => 1.0],                   1, :mentions,     ""),
        ("menbreak", :fractional,          "Aquatic|Terrestrial", ["Aquatic" => 0.5, "Terrestrial" => 0.5], 2, :fractional,   ""),
        ("menbreak", :alphabetical,        "Terrestrial",         ["Terrestrial" => 1.0],                   1, :mentions,     ""),
        ("quarter",  :mentions_fractional, "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline"),
        ("quarter",  :fractional,          "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline"),
        ("quarter",  :alphabetical,        "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline"),
        ("fifth",    :mentions_fractional, "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline"),
        ("fifth",    :fractional,          "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline"),
        ("fifth",    :alphabetical,        "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline"),
        ("spread5",  :mentions_fractional, "Aquatic|Geological|Host-animal|Host-plant|Terrestrial",
                                                                  ["Aquatic" => 0.2, "Geological" => 0.2, "Host-animal" => 0.2,
                                                                   "Host-plant" => 0.2, "Terrestrial" => 0.2],
                                                                                                            5, :fractional,   ""),
        ("spread5",  :fractional,          "Aquatic|Geological|Host-animal|Host-plant|Terrestrial",
                                                                  ["Aquatic" => 0.2, "Geological" => 0.2, "Host-animal" => 0.2,
                                                                   "Host-plant" => 0.2, "Terrestrial" => 0.2],
                                                                                                            5, :fractional,   ""),
        ("spread5",  :alphabetical,        "Aquatic",             ["Aquatic" => 1.0],                       5, :alphabetical, ""),
        ("hotbrine", :mentions_fractional, "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline|Thermal"),
        ("hotbrine", :fractional,          "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline|Thermal"),
        ("hotbrine", :alphabetical,        "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline"),
        ("combo4",   :mentions_fractional, "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline"),
        ("combo4",   :alphabetical,        "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline"),
        ("combo5",   :mentions_fractional, "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline"),
        ("combo5",   :alphabetical,        "Aquatic",             ["Aquatic" => 1.0],                       1, :none,         "Hypersaline"),
    ]
    fixture = Dict(S_FIXTURES)
    for (nm, rule, ht, tw, nt, by, chem) in prior
        a = aggregate_habitats(fixture[nm], path_of; tie_rule = rule)
        @test a.habitat_type == ht
        @test a.type_weights == Dict(tw)
        @test a.n_tied == nt
        @test a.resolved_by == by
        @test a.chemistry == chem
    end

    # The raw evidence is the same under every rule, new or old: it is what the whole design
    # rests on — a variant rule is recomputable from these strings without refetching.
    counted = [("tie2", "Aquatic:1;Terrestrial:1", "", 0),
               ("menbreak", "Aquatic:1;Terrestrial:1", "", 0),
               ("quarter", "Aquatic:3;Terrestrial:1", "Hypersaline", 1),
               ("fifth", "Aquatic:4;Terrestrial:1", "Hypersaline", 2),
               ("spread5", "Aquatic:1;Geological:1;Host-animal:1;Host-plant:1;Terrestrial:1", "", 0),
               ("hotbrine", "Aquatic:2;Terrestrial:1", "Hypersaline+Thermal", 2),
               ("combo4", "Aquatic:3;Other:1", "Alkaline+Hypersaline", 4),
               ("combo5", "Aquatic:3;Other:1", "Alkaline+Hypersaline", 5)]
    for (nm, tc, flags, nchem) in counted, r in EnvironmentTools.TIE_RULES
        a = aggregate_habitats(fixture[nm], path_of; tie_rule = r)
        @test a.type_counts == tc
        @test a.chemistry_flags == flags
        @test a.n_chem_terms == nchem
    end

    # `:legacy_dict` is asserted on SHAPE only, as it is above: its winner is whichever key the
    # hash table yields first, which is the property it exists to reproduce.
    for (nm, fx) in S_FIXTURES
        l = aggregate_habitats(fx, path_of; tie_rule = :legacy_dict)
        @test length(l.type_weights) == 1
        @test collect(values(l.type_weights)) == [1.0]
        @test l.resolved_by == :none
        @test l.habitat_type in keys(s_tally(l.type_counts))
    end
end
