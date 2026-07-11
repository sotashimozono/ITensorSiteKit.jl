# High-spin SiteType registration (S=3/2, S=2). Every check is against an
# INDEPENDENT closed-form expectation — the su(2) algebra the operators MUST
# satisfy regardless of how they were built — not the construction itself:
#   [Sz, S±] = ±S± ,  [S+, S-] = 2Sz ,  S² = Sx²+Sy²+Sz² = S(S+1)·I ,
#   Sz spectrum = {S, S-1, …, -S},  Sx=(S₊+S₋)/2,  Sy=(S₊−S₋)/2i.

using ITensorSiteKit, Test
using ITensors, ITensorMPS
using LinearAlgebra: I, diag, eigvals

mat(name, s) = matrix(op(name, s), prime(s), s)   # ⟨s'|O|s⟩ as a dense matrix

@testset "high-spin SiteType $st (S=$S)" for (st, S) in (("S=3/2", 1.5), ("S=2", 2.0))
    s = siteinds(st, 1)[1]
    d = Int(2S + 1)
    @test dim(s) == d

    Sz = mat("Sz", s)
    Sp = mat("S+", s)
    Sm = mat("S-", s)
    Sx = mat("Sx", s)
    Sy = mat("Sy", s)

    comm(A, B) = A * B - B * A
    Id = Matrix{ComplexF64}(I, d, d)

    # Sz spectrum = S, S-1, …, -S (basis ordered highest-m first)
    @test sort(real(diag(Sz)); rev=true) ≈ collect(S:-1:(-S))

    # su(2) ladder algebra — the falsifiable core
    @test comm(Sz, Sp) ≈ Sp atol = 1e-12
    @test comm(Sz, Sm) ≈ -Sm atol = 1e-12
    @test comm(Sp, Sm) ≈ 2Sz atol = 1e-12

    # Cartesian ↔ ladder
    @test Sx ≈ (Sp + Sm) / 2 atol = 1e-12
    @test Sy ≈ (Sp - Sm) / (2im) atol = 1e-12

    # Casimir S² = S(S+1)·I (independent of basis / construction)
    @test Sx^2 + Sy^2 + Sz^2 ≈ S * (S + 1) * Id atol = 1e-12

    # Sx, Sy, Sz Hermitian; S+ = (S-)†
    @test Sx ≈ Sx' atol = 1e-12
    @test Sy ≈ Sy' atol = 1e-12
    @test Sz ≈ Sz' atol = 1e-12
    @test Sp ≈ Sm' atol = 1e-12

    # Sz2 = Sz²
    @test mat("Sz2", s) ≈ Sz^2 atol = 1e-12

    # PhysInds routes the SiteType through siteinds
    ph = PhysInds(3; SiteType=st)
    @test length(ph) == 3 && all(dim.(ph) .== d)
end

# Integration: the new SiteType must work through the full
# OpSum → MPO → DMRG pipeline. Two spin-S Heisenberg dimer H = S₁·S₂ has the
# total-spin-0 singlet ground state with the closed-form energy
# E₀ = ½[S_tot(S_tot+1) − 2S(S+1)]|_{S_tot=0} = −S(S+1) — independent of the
# operator construction.
@testset "high-spin Heisenberg dimer GS = -S(S+1) (DMRG) $st" for (st, S) in (
    ("S=3/2", 1.5), ("S=2", 2.0)
)
    sites = siteinds(st, 2)
    os = OpSum()
    os += "Sz", 1, "Sz", 2
    os += 0.5, "S+", 1, "S-", 2
    os += 0.5, "S-", 1, "S+", 2
    H = MPO(os, sites)
    ψ0 = random_mps(sites; linkdims=4)
    E, _ = dmrg(H, ψ0; nsweeps=10, maxdim=[10, 20, 40, 40], cutoff=1e-13, outputlevel=0)
    @test E ≈ -S * (S + 1) atol = 1e-6
end

@testset "high-spin conserve_qns blocks (QN(\"Sz\", 2m))" begin
    using ITensors: SiteType, val
    # `space` returns the QN => degeneracy blocks; sectors = 2m for m=3/2..-3/2
    blocks = space(SiteType("S=3/2"); conserve_qns=true)
    @test all(deg == 1 for (_, deg) in blocks)
    @test sort([val(qn, "Sz") for (qn, _) in blocks]) == [-3, -1, 1, 3]
    # a real QN S=3/2 site index builds, and its densified ops obey the algebra
    s = siteinds("S=3/2", 1; conserve_qns=true)[1]
    @test hasqns(s) && dim(s) == 4
    Sz = mat("Sz", s)
    Sp = mat("S+", s)
    Sm = mat("S-", s)
    @test Sp * Sm - Sm * Sp ≈ 2Sz atol = 1e-12
end
