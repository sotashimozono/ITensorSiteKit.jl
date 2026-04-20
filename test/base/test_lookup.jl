using ITensors, ITensorMPS
using Random

N = 20
χ = 10

phys = PhysInds(N, χ; SiteType="S=1/2")
envs = EnvInds(N, χ; SiteType="S=1/2")
auxs = AuxInds(N, χ)

sites = [envs[1], phys..., envs[2]]
sites_aux = [auxs[1], sites..., auxs[2]]

# Build an MPS on the env+phys chain so we can exercise the MPS-dispatched
# lookup methods.
rng = MersenneTwister(0)
ψ = random_mps(rng, sites; linkdims=χ)

@testset "extract_siteinds" begin
    @test extract_siteinds(sites, PhysSite) == phys
    @test extract_siteinds(sites, LeftEnvSite) == [envs[1]]
    @test extract_siteinds(sites, RightEnvSite) == [envs[2]]
    @test extract_siteinds(ψ, PhysSite) == phys

    # ITensor dispatch: tensor with tagged indices
    T = ITensor(phys[1], phys[2])
    @test extract_siteinds(T, PhysSite) == [phys[1], phys[2]]
end

@testset "phys_site_index / phys_site_position" begin
    @test phys_site_index(sites, 3) == phys[3]
    @test phys_site_index(ψ, 3) == phys[3]

    # In `sites`, phys site n=1 sits at position 2 (after left env).
    @test phys_site_position(sites, 1) == 2
    @test phys_site_position(sites, N) == N + 1
    @test phys_site_position(ψ, 5) == 6
end

@testset "from_right" begin
    @test from_right(sites, 1) == N
    @test from_right(sites, N) == 1
    @test from_right(ψ, 1) == N

    for _ in 1:10
        l = rand(1:N)
        @test from_right(sites, l) == N - l + 1
        pos = phys_site_position(sites, from_right(sites, l))
        @test sites[pos] == phys[end - l + 1]
    end
end

@testset "env_site_position" begin
    @test env_site_position(sites, :left) == 1
    @test env_site_position(sites, :right) == N + 2
    @test env_site_position(ψ, :left) == 1
    @test_throws ErrorException env_site_position(sites, :top)
end

@testset "aux_site_position" begin
    @test aux_site_position(sites_aux, :left) == 1
    @test aux_site_position(sites_aux, :right) == length(sites_aux)
    @test_throws ErrorException aux_site_position(sites_aux, :top)
    # On an env-only chain there are no aux tags, so lookup should return nothing.
    @test aux_site_position(sites, :left) === nothing
end
