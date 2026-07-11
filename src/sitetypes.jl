# High-spin site types. ITensors ships only `SiteType"S=1/2"` and `SiteType"S=1"`;
# this registers `S=3/2` and `S=2` (dim 4 and 5) so `siteinds("S=3/2", N)` — and
# therefore [`PhysInds`](@ref)`(N; SiteType="S=3/2")` — work, unblocking high-spin
# lattice models (ITensorModels S1/S3-2 chains, QAtlas high-spin references,
# ThermalMPS S≥1 thermal benchmarks).
#
# Conventions match the ITensors built-in spin sites exactly:
# - basis ordered by decreasing magnetic number `m = S, S-1, …, -S` (index k ↔ m
#   = S-(k-1)), so `Sz = diag(S, S-1, …, -S)`;
# - `conserve_qns` blocks carry `QN("Sz", 2m)` (integer half-units, so S=1/2 is
#   `QN("Sz", ±1)` — the same units ITensors uses for S=1/2 and S=1).

using ITensors: ITensors, SiteType, OpName, StateName, QN

# Spin-S operator matrices in the `m = S, S-1, …, -S` basis. Closed-form ladder
# elements ⟨m±1|S±|m⟩ = √(S(S+1) − m(m±1)); Sx=(S₊+S₋)/2, Sy=(S₊−S₋)/2i.
function _spin_matrices(S::Rational{Int})
    d = Int(2S + 1)
    ms = [S - (k - 1) for k in 1:d]
    Sz = zeros(ComplexF64, d, d)
    for k in 1:d
        Sz[k, k] = ms[k]
    end
    Sp = zeros(ComplexF64, d, d)          # S₊: |m⟩ → |m+1⟩ (raises m, lowers index)
    for k in 2:d
        m = ms[k]
        Sp[k - 1, k] = sqrt(float(S * (S + 1) - m * (m + 1)))
    end
    Sm = collect(Sp')                     # S₋ = S₊†
    Sx = (Sp + Sm) / 2
    Sy = (Sp - Sm) / (2im)
    return (; Sz, Sp, Sm, Sx, Sy, Sz2=Sz^2)
end

# QN("Sz", 2m) blocks, one per basis state, highest m first.
function _sz_qn_blocks(S::Rational{Int}, qnname)
    d = Int(2S + 1)
    return [QN(qnname, Int(2 * (S - (k - 1)))) => 1 for k in 1:d]
end

const _OPS_S32 = _spin_matrices(3//2)
const _OPS_S2 = _spin_matrices(2//1)

# --- S=3/2 (dim 4) ---
function ITensors.space(
    ::SiteType"S=3/2"; conserve_qns=false, conserve_sz=conserve_qns, qnname_sz="Sz"
)
    conserve_sz && return _sz_qn_blocks(3//2, qnname_sz)
    return 4
end
ITensors.op(::OpName"Sz", ::SiteType"S=3/2") = copy(_OPS_S32.Sz)
ITensors.op(::OpName"S+", ::SiteType"S=3/2") = copy(_OPS_S32.Sp)
ITensors.op(::OpName"S-", ::SiteType"S=3/2") = copy(_OPS_S32.Sm)
ITensors.op(::OpName"Sx", ::SiteType"S=3/2") = copy(_OPS_S32.Sx)
ITensors.op(::OpName"Sy", ::SiteType"S=3/2") = copy(_OPS_S32.Sy)
ITensors.op(::OpName"Sz2", ::SiteType"S=3/2") = copy(_OPS_S32.Sz2)

# --- S=2 (dim 5) ---
function ITensors.space(
    ::SiteType"S=2"; conserve_qns=false, conserve_sz=conserve_qns, qnname_sz="Sz"
)
    conserve_sz && return _sz_qn_blocks(2//1, qnname_sz)
    return 5
end
ITensors.op(::OpName"Sz", ::SiteType"S=2") = copy(_OPS_S2.Sz)
ITensors.op(::OpName"S+", ::SiteType"S=2") = copy(_OPS_S2.Sp)
ITensors.op(::OpName"S-", ::SiteType"S=2") = copy(_OPS_S2.Sm)
ITensors.op(::OpName"Sx", ::SiteType"S=2") = copy(_OPS_S2.Sx)
ITensors.op(::OpName"Sy", ::SiteType"S=2") = copy(_OPS_S2.Sy)
ITensors.op(::OpName"Sz2", ::SiteType"S=2") = copy(_OPS_S2.Sz2)

# Named states for the extreme sublevels (the universal `Up`/`Dn` aliases).
# Arbitrary sublevels are addressable by 1-based basis index (m = S, S-1, …, -S
# order): `MPS(siteinds("S=3/2", N; conserve_qns=true), [1, 2, 4, …])` builds a
# definite-flux product state with no StateName needed (integer basis-index path).
ITensors.state(::StateName"Up", ::SiteType"S=3/2") = [1.0, 0.0, 0.0, 0.0]
ITensors.state(::StateName"Dn", ::SiteType"S=3/2") = [0.0, 0.0, 0.0, 1.0]
ITensors.state(::StateName"Up", ::SiteType"S=2") = [1.0, 0.0, 0.0, 0.0, 0.0]
ITensors.state(::StateName"Dn", ::SiteType"S=2") = [0.0, 0.0, 0.0, 0.0, 1.0]
