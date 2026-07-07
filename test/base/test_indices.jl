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

    # Dense is still the default when a plain `Int` space is passed.
    @testset "dense default (backward compat)" begin
        @test !hasqns(EnvInds(χ)[1])
        @test !hasqns(AuxInds(χ)[1])
        @test !hasqns(first(LinkInds(χ, PhysInds(N; SiteType="S=1/2"))))
        @test !hasqns(PhysInds(N; SiteType="S=1/2")[1])
    end
end

@testset "QN (symmetric) index constructors" begin
    # U(1): total-Sz sectors.
    u1 = [QN("Sz", 0) => 3, QN("Sz", 2) => 2]
    # Z(2): parity carrier (modulus 2).
    z2 = [QN("P", 0, 2) => 2, QN("P", 1, 2) => 2]

    @testset "PhysInds conserve_qns" begin
        sites = PhysInds(N; SiteType="S=1/2", conserve_qns=true)
        @test length(sites) == N
        @test all(hasqns, sites)                       # QN turned on...
        @test all(i -> dim(i) == 2, sites)             # ...dim unchanged
        @test all(i -> hastags(i, PhysSite), sites)    # ...tag preserved
        # kwargs also flow through the 3-arg form.
        @test all(hasqns, PhysInds(N, χ; SiteType="S=1/2", conserve_qns=true))
    end

    @testset "EnvInds QN carrier" begin
        left, right = EnvInds(u1)
        @test hasqns(left) && hasqns(right)
        @test dim(left) == 5 && dim(right) == 5        # dim = Σ degeneracies
        @test hastags(left, LeftEnvSite) && hastags(right, RightEnvSite)

        # Z(2) parity carrier.
        lp, rp = EnvInds(z2)
        @test hasqns(lp) && dim(lp) == 4 && hastags(rp, RightEnvSite)

        # 3-arg form (N ignored) forwards the QN space too.
        @test hasqns(EnvInds(N, u1)[1])

        # Physical SiteType branch forwards `conserve_qns` to `siteind`.
        ls, rs = EnvInds(χ; SiteType="S=1/2", conserve_qns=true)
        @test hasqns(ls) && dim(ls) == 2 && hastags(ls, LeftEnvSite)
        @test op("Sz", ls) isa ITensor
    end

    @testset "AuxInds QN carrier" begin
        left, right = AuxInds(u1)
        @test hasqns(left) && hasqns(right)
        @test dim(left) == 5
        @test hastags(left, LeftAuxSite) && hastags(right, RightAuxSite)
        @test hasqns(AuxInds(N, u1)[1])
    end

    @testset "LinkInds QN carrier" begin
        sites = PhysInds(N; SiteType="S=1/2", conserve_qns=true)
        links = LinkInds(u1, sites)
        @test length(links) == N - 1
        @test all(hasqns, links)
        @test all(i -> dim(i) == 5, links)
        @test all(i -> hastags(i, "Link"), links)
    end
end
