#=
    The names `using LnFinder` must actually bring into scope.

    The failure this suite exists to prevent is one no other suite here can see. Every
    other suite loads the module file it tests directly:

        include(joinpath(@__DIR__, "..", "src", "ModuleTaxonomyTools.jl"))
        using .TaxonomyTools

    That is deliberate — it keeps each suite standalone and offline — but it bypasses the
    package boundary completely. A function can be written, documented, unit-tested and
    green while never being re-exported from `module LnFinder`, because `using .SomeModule`
    inside LnFinder brings a name into LnFinder's namespace without re-exporting it. The
    census driver does plain `using LnFinder` and then calls the name unqualified, so the
    gap surfaces only at runtime, after the expensive part of the pipeline has already run.

    So this suite goes the other way round: it loads the package and asserts that the names
    the driver calls are exported. It is cheap, and it fails the moment an export is
    dropped rather than an hour into a rebuild.

    Run:  julia --project=. test/test_package_surface.jl
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
using LnFinder

@testset "LnFinder package surface" begin
    exported = Set(names(LnFinder))

    # Called unqualified by src/scriptLn.jl under `using LnFinder`. Each of these was, or
    # could have been, defined in a submodule and left unexported.
    @testset "the census driver's calls resolve" begin
        for f in (:parse_fasta_headers, :calibrate_active_site, :parse_hmmalign_results,
                  :filter_centroids, :resolve_taxids, :load_taxonomy_db, :append_taxonomy!,
                  :taxdump_provenance, :fetch_omnicrobe_env!, :agg_habitat)
            @test f in exported
        end
    end

    # A name being exported is not the same as it being callable: an export list can name a
    # symbol the package never defines, and `names()` will still list it.
    @testset "every exported name is defined" begin
        for f in names(LnFinder)
            @test isdefined(LnFinder, f)
        end
    end

    # The back-compat wrapper src/scriptLn.jl's older call sites still rely on.
    @testset "fetch_taxids_async survives alongside resolve_taxids" begin
        @test :fetch_taxids_async in exported
        @test :resolve_taxids in exported
    end
end
