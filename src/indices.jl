"""
    PhysInds(N::Int; SiteType="S=1/2")

Construct `N` physical site indices of the given `SiteType`, tagged with
[`PhysSite`](@ref). `SiteType` follows ITensors.jl site-type conventions.

    PhysInds(N::Int, χ::Int; kwargs...)

Three-argument form accepted for symmetry with [`EnvInds`](@ref) /
[`AuxInds`](@ref); `χ` is ignored (physical site dim is fixed by
`SiteType`).
"""
function PhysInds(N::Int; SiteType="S=1/2")
    sites_physical = siteinds(SiteType, N)
    sites_physical = addtags(sites_physical, PhysSite)
    return sites_physical
end
PhysInds(N::Int, χ::Int; kwargs...) = PhysInds(N; kwargs...)
export PhysInds

"""
    EnvInds(χ::Int; SiteType="Environment")

Construct a pair `[left, right]` of environment site indices of dimension
`χ`, tagged with [`LeftEnvSite`](@ref) / [`RightEnvSite`](@ref).

If `SiteType == "Environment"` (default), the indices carry a raw
`"Environment,Site"` tag set (generic env carrier of dim `χ`). Otherwise
the indices are built via `siteind(SiteType)` (a standard physical site
type) and then tagged — useful when the env block is constructed from a
full physical DMRG ground state and should share its local Hilbert space.

    EnvInds(N::Int, χ::Int; kwargs...)

Three-argument form accepted for API symmetry; `N` is ignored.
"""
function EnvInds(χ::Int; SiteType="Environment")
    if SiteType == "Environment"
        sites_left_environment = Index(χ, LeftEnvSite * ",Site")
        sites_right_environment = Index(χ, RightEnvSite * ",Site")
    else
        sites_left_environment = siteind(SiteType)
        sites_right_environment = siteind(SiteType)
        sites_left_environment = addtags(sites_left_environment, LeftEnvSite)
        sites_right_environment = addtags(sites_right_environment, RightEnvSite)
    end
    return [sites_left_environment, sites_right_environment]
end
EnvInds(N::Int, χ::Int; kwargs...) = EnvInds(χ; kwargs...)
export EnvInds

"""
    AuxInds(χ_aux::Int)

Construct a pair `[left, right]` of auxiliary site indices of dimension
`χ_aux`, tagged with [`LeftAuxSite`](@ref) / [`RightAuxSite`](@ref).

    AuxInds(N::Int, χ_aux::Int)

Three-argument form accepted for API symmetry; `N` is ignored.
"""
function AuxInds(χ_aux::Int)
    sites_left_aux = Index(χ_aux, LeftAuxSite * ",Site")
    sites_right_aux = Index(χ_aux, RightAuxSite * ",Site")
    return [sites_left_aux, sites_right_aux]
end
AuxInds(N::Int, χ_aux::Int) = AuxInds(χ_aux)
export AuxInds

"""
    LinkInds(χ::Int, sites)

Return a vector of `length(sites) - 1` freshly-generated link indices of
dimension `χ`, each tagged `"Link, l=<i>"`.
"""
function LinkInds(χ::Int, sites)
    return [Index(χ, "Link, l=$l") for l in 1:(length(sites) - 1)]
end
export LinkInds
