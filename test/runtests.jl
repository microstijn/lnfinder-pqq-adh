#=
    Aggregate offline test runner for the Ca-vs-Ln tooling. No network, no HMMER, no
    external binary and no input this repository does not carry: every suite makes its
    claim on a synthetic fixture it builds itself.

    Every suite in test/ is included here. A suite that ships but is not run is a suite
    nobody notices breaking, so the two sets are kept equal.

    Run:  julia --project=. test/runtests.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Test

@testset "LnFinder Ca/Ln tooling" begin
    include("test_source_layout.jl")         # every included file exists, spelled exactly
    include("test_data_inputs.jl")           # one data root: every input resolves under it, a missing one refuses by name
    include("test_validate_classifier.jl")   # classifier validation + calibration guard
    include("test_environment_tools.jl")     # OmniCrobe habitat aggregation + pillar roll-up
    include("test_taxonomy_tools.jl")        # NCBI dump vintages, retired TaxIDs, centroid filter
    include("test_accession_mapper.jl")      # UniProt/UniParc rescue, offline against live shapes
    include("test_habitat_stats.jl")         # enrichment estimators, fractional hits, symmetric genus unit
    include("test_outcome_rules.jl")         # result-to-decision predicates: three-valued intervals, anchor-relative categories
    include("test_background_scan.jl")       # background PQQ-ADH prefilter: tblout parsing, frame and exclusion joins
    include("test_habitat_enrichment.jl")    # enrichment consumer: scan evidence, licence wiring, anchor ladder
    include("test_habitat_sensitivity.jl")   # one-at-a-time rule grid: variant set, robustness flags, anchors
    include("test_habitat_reaggregate.jl")   # habitat call table: gate before write, declared gates and configs
    include("test_habitat_figure.jl")        # figure: input gates across five files, colour from the flag, staged publish
    include("test_package_surface.jl")       # the names the pipeline calls must be exported
end
