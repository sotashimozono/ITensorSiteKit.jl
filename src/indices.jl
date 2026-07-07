"""
    const IndexSpace = Union{Integer,AbstractVector{<:Pair{QN,<:Integer}}}

Accepted `space` argument for the index constructors in this file: either

- an `Integer` dimension → a **dense** `Index`, or
- a vector of `QN => degeneracy` blocks → an **abelian-symmetric**
  (block-sparse) `Index`,

mirroring `ITensors.Index`. This lets every constructor below build both
dense and quantum-number-conserving (U(1) / Z(N)) indices through the same
call. Non-abelian symmetries (e.g. SU(2)) are not representable in the
ITensors backend and are out of scope.
"""
const IndexSpace = Union{Integer,AbstractVector{<:Pair{QN,<:Integer}}}
export IndexSpace

"""
    PhysInds(N::Int; SiteType="S=1/2", kwargs...)

Construct `N` physical site indices of the given `SiteType`, tagged with
[`PhysSite`](@ref). `SiteType` follows ITensors.jl site-type conventions.

Any `kwargs` are forwarded to `ITensors.siteinds`, so quantum-number
conservation is opt-in in the usual way:

```julia
PhysInds(N; SiteType="S=1/2", conserve_qns=true)   # U(1) (total Sz) sites
```

Tagging with `addtags` leaves the QN block structure untouched, so the
returned indices stay block-sparse.

    PhysInds(N::Int, χ::Int; kwargs...)

Three-argument form accepted for symmetry with [`EnvInds`](@ref) /
[`AuxInds`](@ref); `χ` is ignored (physical site dim is fixed by
`SiteType`). `kwargs` are still forwarded.
"""
function PhysInds(N::Int; SiteType="S=1/2", kwargs...)
    sites_physical = siteinds(SiteType, N; kwargs...)
    sites_physical = addtags(sites_physical, PhysSite)
    return sites_physical
end
PhysInds(N::Int, χ::Int; kwargs...) = PhysInds(N; kwargs...)
export PhysInds

"""
    EnvInds(space::IndexSpace; SiteType="Environment", kwargs...)

Construct a pair `[left, right]` of environment site indices, tagged with
[`LeftEnvSite`](@ref) / [`RightEnvSite`](@ref). `space` is an
[`IndexSpace`](@ref): an `Integer` (dense carrier of that dimension) or a
vector of `QN => degeneracy` blocks (symmetric carrier).

If `SiteType == "Environment"` (default), the indices are raw
`"Environment,Site"` carriers built directly from `space` — pass QN blocks
here to give the env leg the sector structure `χ = Σ_q χ_q` of the cut
Schmidt basis it represents. Otherwise the indices are built via
`siteind(SiteType; kwargs...)` (a standard physical site type, with `kwargs`
such as `conserve_qns` forwarded) and then tagged — useful when the env
block shares the bulk local Hilbert space; in that case `space` is ignored
(the dim is fixed by `SiteType`).

```julia
EnvInds(χ)                                        # dense, dim χ
EnvInds([QN("Sz",0)=>3, QN("Sz",2)=>2])           # U(1), dim 5
EnvInds(χ; SiteType="S=1/2", conserve_qns=true)   # physical S=1/2 QN site
```

    EnvInds(N::Int, space::IndexSpace; kwargs...)

Three-argument form accepted for API symmetry; `N` is ignored.
"""
function EnvInds(space::IndexSpace; SiteType="Environment", kwargs...)
    if SiteType == "Environment"
        sites_left_environment = Index(space, LeftEnvSite * ",Site")
        sites_right_environment = Index(space, RightEnvSite * ",Site")
    else
        sites_left_environment = siteind(SiteType; kwargs...)
        sites_right_environment = siteind(SiteType; kwargs...)
        sites_left_environment = addtags(sites_left_environment, LeftEnvSite)
        sites_right_environment = addtags(sites_right_environment, RightEnvSite)
    end
    return [sites_left_environment, sites_right_environment]
end
EnvInds(N::Int, space::IndexSpace; kwargs...) = EnvInds(space; kwargs...)
export EnvInds

"""
    AuxInds(space::IndexSpace)

Construct a pair `[left, right]` of auxiliary site indices, tagged with
[`LeftAuxSite`](@ref) / [`RightAuxSite`](@ref). `space` is an
[`IndexSpace`](@ref): an `Integer` for a dense carrier of that dimension, or
a vector of `QN => degeneracy` blocks for a symmetric one.

    AuxInds(N::Int, space::IndexSpace)

Three-argument form accepted for API symmetry; `N` is ignored.
"""
function AuxInds(space::IndexSpace)
    sites_left_aux = Index(space, LeftAuxSite * ",Site")
    sites_right_aux = Index(space, RightAuxSite * ",Site")
    return [sites_left_aux, sites_right_aux]
end
AuxInds(N::Int, space::IndexSpace) = AuxInds(space)
export AuxInds

"""
    LinkInds(space::IndexSpace, sites)

Return a vector of `length(sites) - 1` freshly-generated link indices built
from `space` (an [`IndexSpace`](@ref)), each tagged `"Link, l=<i>"`.

Passing an `Integer` gives dense links of that dimension. Passing QN blocks
gives symmetric links that all carry the *same* sectors — this is bare
scaffolding: a physically valid quantum-number-conserving MPS needs its link
sectors chosen so the per-tensor flux is consistent, which
`ITensorMPS.random_mps(sites; ...)` handles automatically. Prefer that for
real symmetric initial states; use `LinkInds` with QN blocks only when you
are supplying the sector layout yourself.
"""
function LinkInds(space::IndexSpace, sites)
    return [Index(space, "Link, l=$l") for l in 1:(length(sites) - 1)]
end
export LinkInds
