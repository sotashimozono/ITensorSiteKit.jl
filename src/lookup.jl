"""
    extract_siteinds(T, tag::AbstractString)

Return the indices from `T` that carry `tag`. `T` can be a
`Vector{<:Index}`, an `MPS`, or an `ITensor`.
"""
function extract_siteinds end

function extract_siteinds(si::Vector{<:Index}, tag::AbstractString)
    return filter(i -> hastags(i, tag), si)
end
function extract_siteinds(ψ::MPS, tag::AbstractString)
    return extract_siteinds(siteinds(ψ), tag)
end
function extract_siteinds(tensor::ITensor, tag::AbstractString)
    return filter(i -> hastags(i, tag), inds(tensor))
end
export extract_siteinds

"""
    phys_site_index(T, n::Int)

Return the physical site index labelled `n=<n>` in `T` (a
`Vector{<:Index}` or `MPS`).
"""
function phys_site_index(si::Vector{<:Index}, n::Int)
    si_phys = extract_siteinds(si, PhysSite)
    return filter(i -> hastags(i, "n=$n"), si_phys)[1]
end
phys_site_index(ψ::MPS, n::Int) = phys_site_index(siteinds(ψ), n)
export phys_site_index

"""
    phys_site_position(T, n::Int)

Return the position (1-based integer) of the physical site labelled
`n=<n>` within the full site ordering of `T`.
"""
function phys_site_position(si::Vector{<:Index}, n::Int)
    return findfirst(i -> hastags(i, "n=$n"), si)
end
phys_site_position(ψ::MPS, n::Int) = phys_site_position(siteinds(ψ), n)
export phys_site_position

"""
    from_right(T, n::Int)

Given a physical site number `n` counting from the left, return the
corresponding number counted from the right end of the physical sites.
"""
function from_right(si::Vector{<:Index}, n::Int)
    L = length(extract_siteinds(si, PhysSite))
    return L - n + 1
end
from_right(ψ::MPS, n::Int) = from_right(siteinds(ψ), n)
export from_right

"""
    env_site_position(T, direction::Symbol)

Position of the environment site (`:left` or `:right`) within `T`.
"""
function env_site_position(si::Vector{<:Index}, direction::Symbol)
    if direction === :left
        return findfirst(i -> hastags(i, LeftEnvSite), si)
    elseif direction === :right
        return findfirst(i -> hastags(i, RightEnvSite), si)
    else
        error("Unsupported direction: $direction. Use :left or :right.")
    end
end
env_site_position(ψ::MPS, direction::Symbol) = env_site_position(siteinds(ψ), direction)
export env_site_position

"""
    aux_site_position(T, direction::Symbol)

Position of the auxiliary site (`:left` or `:right`) within `T`.
"""
function aux_site_position(si::Vector{<:Index}, direction::Symbol)
    if direction === :left
        return findfirst(i -> hastags(i, LeftAuxSite), si)
    elseif direction === :right
        return findfirst(i -> hastags(i, RightAuxSite), si)
    else
        error("Unsupported direction: $direction. Use :left or :right.")
    end
end
aux_site_position(ψ::MPS, direction::Symbol) = aux_site_position(siteinds(ψ), direction)
export aux_site_position
