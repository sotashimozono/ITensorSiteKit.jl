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

@testset "site_position (index → position)" begin
    # every site index resolves to its own position; round-trips against the
    # tag-based phys_site_position (independent method).
    for (pos, i) in enumerate(sites)
        @test site_position(sites, i) == pos
    end
    @test site_position(sites, phys[5]) == 6          # [envL, phys_1..N, envR]
    @test site_position(ψ, envs[2]) == N + 2
    @test site_position(sites, envs[1]) == 1
    # phys site n's position matches the tag-based lookup
    @test site_position(sites, phys[7]) == phys_site_position(sites, 7)
    # absent index → nothing
    @test site_position(sites, auxs[1]) === nothing
end

@testset "classify_site_indices" begin
    c = classify_site_indices(sites)
    @test c[:physical] == phys
    @test c[:left_env] == [envs[1]]
    @test c[:right_env] == [envs[2]]
    @test isempty(c[:left_aux]) && isempty(c[:right_aux]) && isempty(c[:other])
    # partition: buckets are disjoint and cover every input index
    total = sum(length, values(c))
    @test total == length(sites)

    ca = classify_site_indices(sites_aux)
    @test ca[:left_aux] == [auxs[1]]
    @test ca[:right_aux] == [auxs[2]]
    @test ca[:physical] == phys
    @test sum(length, values(ca)) == length(sites_aux)

    # MPS dispatch agrees with the site-vector form
    @test classify_site_indices(ψ)[:physical] == phys
end
