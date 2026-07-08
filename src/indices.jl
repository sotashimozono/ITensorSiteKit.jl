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
    _space_index(space::IndexSpace, tag; dir=nothing)

Internal helper shared by the `space`-based constructors below: build a raw
`Index` from `space`/`tag`, optionally forwarding an explicit `dir` to
`Index`.

`dir=nothing` (the default) reproduces the plain `Index(space, tag)` call
with no `dir` keyword at all — so omitting `dir` is exactly backward
compatible. A non-`nothing` `dir` is meaningful for a QN (symmetric) `space`
(where `ITensors.Index` accepts a `dir` keyword) and errors for a dense
(`Integer`) `space`, mirroring `Index` itself: a dense index carries no
arrow, so `Index` has no `dir`-accepting method for it.
"""
function _space_index(space::IndexSpace, tag; dir=nothing)
    return dir === nothing ? Index(space, tag) : Index(space, tag; dir=dir)
end

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
    EnvInds(space::IndexSpace; SiteType="Environment", dir=nothing, kwargs...)

Construct a pair `[left, right]` of environment site indices, tagged with
[`LeftEnvSite`](@ref) / [`RightEnvSite`](@ref). `space` is an
[`IndexSpace`](@ref): an `Integer` (dense carrier of that dimension) or a
vector of `QN => degeneracy` blocks (symmetric carrier).

If `SiteType == "Environment"` (default), the indices are raw
`"Environment,Site"` carriers built directly from `space` — pass QN blocks
here to give the env leg the sector structure `χ = Σ_q χ_q` of the cut
Schmidt basis it represents. `dir`, if given, is forwarded to `Index` on
this branch (see [`IndexSpace`](@ref)'s `_space_index` helper); the default
`dir=nothing` reproduces the previous fixed behavior (`Index`'s own default,
`Out`, for a QN `space`). Otherwise the indices are built via
`siteind(SiteType; kwargs...)` (a standard physical site type, with `kwargs`
such as `conserve_qns` forwarded) and then tagged — useful when the env
block shares the bulk local Hilbert space; in that case `space` and `dir`
are both ignored (the dim is fixed by `SiteType`, and `siteind` has no `dir`
keyword).

```julia
EnvInds(χ)                                        # dense, dim χ
EnvInds([QN("Sz",0)=>3, QN("Sz",2)=>2])           # U(1), dim 5, dir=Out
EnvInds([QN("Sz",0)=>3, QN("Sz",2)=>2]; dir=ITensors.In)  # same, dir=In
EnvInds(χ; SiteType="S=1/2", conserve_qns=true)   # physical S=1/2 QN site
```

    EnvInds(N::Int, space::IndexSpace; kwargs...)

Three-argument form accepted for API symmetry; `N` is ignored. `kwargs`
(including `dir`) are forwarded.

    EnvInds(left_bond::Index, right_bond::Index)

Construct `[left, right]` environment indices *derived from an existing MPS
cut bond* via `sim`, so each output index has exactly the same `space` (QN
sectors, or dense dimension) **and** the same `dir` (arrow) as the
corresponding bond — merely re-tagged to `LeftEnvSite`/`RightEnvSite` (`sim`
replaces the tag set outright, so the bond's own `"Link"` tag is dropped,
not just supplemented).

This is the robust way to build an environment carrier that must
dag-contract cleanly against a real cut bond of a QN-MPS: unlike the
`space`-based constructor above — whose QN-blocks branch always returns a
fixed `dir=Out` handle unless you pass a matching `dir` yourself — this
preserves whatever `dir` the bond actually has (`ITensorMPS.linkinds`
bonds are typically `dir=In`), with no risk of a block-order/arrow mismatch.
Works identically for dense and QN bonds.

```julia
sites = siteinds("S=1/2", N; conserve_qns=true)
ψ = random_mps(sites, states)
bonds = linkinds(ψ)
left, right = EnvInds(bonds[i], bonds[j])   # space & dir both match the cut
```
"""
function EnvInds(space::IndexSpace; SiteType="Environment", dir=nothing, kwargs...)
    if SiteType == "Environment"
        sites_left_environment = _space_index(space, LeftEnvSite * ",Site"; dir=dir)
        sites_right_environment = _space_index(space, RightEnvSite * ",Site"; dir=dir)
    else
        sites_left_environment = siteind(SiteType; kwargs...)
        sites_right_environment = siteind(SiteType; kwargs...)
        sites_left_environment = addtags(sites_left_environment, LeftEnvSite)
        sites_right_environment = addtags(sites_right_environment, RightEnvSite)
    end
    return [sites_left_environment, sites_right_environment]
end
EnvInds(N::Int, space::IndexSpace; kwargs...) = EnvInds(space; kwargs...)
function EnvInds(left_bond::Index, right_bond::Index)
    sites_left_environment = sim(left_bond; tags=LeftEnvSite * ",Site")
    sites_right_environment = sim(right_bond; tags=RightEnvSite * ",Site")
    return [sites_left_environment, sites_right_environment]
end
export EnvInds

"""
    AuxInds(space::IndexSpace; dir=nothing)

Construct a pair `[left, right]` of auxiliary site indices, tagged with
[`LeftAuxSite`](@ref) / [`RightAuxSite`](@ref). `space` is an
[`IndexSpace`](@ref): an `Integer` for a dense carrier of that dimension, or
a vector of `QN => degeneracy` blocks for a symmetric one. `dir`, if given,
is forwarded to `Index` (see `_space_index`); the default `dir=nothing`
reproduces the previous fixed behavior (`Index`'s own default, `Out`, for a
QN `space`).

    AuxInds(N::Int, space::IndexSpace; kwargs...)

Three-argument form accepted for API symmetry; `N` is ignored. `kwargs`
(including `dir`) are forwarded.

    AuxInds(left_bond::Index, right_bond::Index)

Construct `[left, right]` auxiliary indices *derived from an existing MPS
cut bond* via `sim`: same mechanism as the [`EnvInds`](@ref) bond-derived
method — each output has the same `space` and `dir` as the corresponding
bond, re-tagged to `LeftAuxSite`/`RightAuxSite` (the bond's `"Link"` tag is
dropped since `sim` replaces the tag set). Works for dense and QN bonds.
"""
function AuxInds(space::IndexSpace; dir=nothing)
    sites_left_aux = _space_index(space, LeftAuxSite * ",Site"; dir=dir)
    sites_right_aux = _space_index(space, RightAuxSite * ",Site"; dir=dir)
    return [sites_left_aux, sites_right_aux]
end
AuxInds(N::Int, space::IndexSpace; kwargs...) = AuxInds(space; kwargs...)
function AuxInds(left_bond::Index, right_bond::Index)
    sites_left_aux = sim(left_bond; tags=LeftAuxSite * ",Site")
    sites_right_aux = sim(right_bond; tags=RightAuxSite * ",Site")
    return [sites_left_aux, sites_right_aux]
end
export AuxInds

"""
    LinkInds(χ::Int, sites)

Return a vector of `length(sites) - 1` freshly-generated **dense** link
indices of dimension `χ`, each tagged `"Link, l=<i>"`.

There is no QN (symmetric) form of this constructor: a quantum-number-
conserving MPS assigns a *different* accumulated flux sector to every bond
(built up tensor-by-tensor so the total flux is consistent), so handing
every link the same fixed `space` would describe a physically invalid MPS.
Build real symmetric link indices from an actual QN-conserving MPS instead —
`ITensorMPS.linkinds(ITensorMPS.random_mps(sites; ...))` — which handles the
per-bond sector layout automatically. If you then need an env/aux carrier
matching one of those real bonds exactly (same sectors, same arrow), pass
the bond straight to the `sim`-based [`EnvInds`](@ref)/[`AuxInds`](@ref)
methods rather than reconstructing it from blocks here.
"""
function LinkInds(χ::Int, sites)
    return [Index(χ, "Link, l=$l") for l in 1:(length(sites) - 1)]
end
export LinkInds
