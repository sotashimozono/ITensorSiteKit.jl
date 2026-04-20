using ITensors, ITensorMPS

N = 3
χ = 10

@testset "Index constructors" begin
    @testset "PhysInds" begin
        sites = PhysInds(N, χ; SiteType="S=1/2")
        @test length(sites) == N
        @test all(i -> hastags(i, PhysSite), sites)
        @test all(i -> dim(i) == 2, sites)
    end

    @testset "EnvInds" begin
        left, right = EnvInds(N, χ; SiteType="Environment")
        @test dim(left) == χ
        @test dim(right) == χ
        @test hastags(left, LeftEnvSite)
        @test hastags(right, RightEnvSite)

        left, right = EnvInds(N, χ; SiteType="S=1/2")
        @test dim(left) == 2
        @test dim(right) == 2
        @test hastags(left, LeftEnvSite)
        @test hastags(right, RightEnvSite)

        # Standard site types should still accept "Sz" operators after tagging.
        @test op("Sz", left) isa ITensor
        @test op("Sz", right) isa ITensor
    end

    @testset "AuxInds" begin
        left, right = AuxInds(N, χ)
        @test dim(left) == χ
        @test dim(right) == χ
        @test hastags(left, LeftAuxSite)
        @test hastags(right, RightAuxSite)
    end

    @testset "LinkInds" begin
        sites = PhysInds(N, χ; SiteType="S=1/2")
        links = LinkInds(χ, sites)
        @test length(links) == N - 1
        @test all(i -> hastags(i, "Link"), links)
        @test all(i -> dim(i) == χ, links)

        envs = EnvInds(N, χ; SiteType="Environment")
        full = [envs[1], sites..., envs[2]]
        links = LinkInds(χ, full)
        @test length(links) == length(full) - 1

        auxs = AuxInds(N, χ)
        full2 = [auxs[1], full..., auxs[2]]
        links = LinkInds(χ, full2)
        @test length(links) == length(full2) - 1
    end
end
