𝐁(x,α,β) = pdf(Beta(α,β),x)

struct BetaCopula{N} <:AbstractBivariateCopula 
    U::AbstractVector{<:Real}
    V::AbstractVector{<:Real}
    h::Real

    function BetaCopula(U::AbstractVector{<:Real}, V::AbstractVector{<:Real}; h = 0.05)
        @assert length(U) == length(V)
        @assert all(0 .≤ U .≤ 1)
        @assert all(0 .≤ V .≤ 1)

        N = length(U)

        new{N}(U,V,h)
    end
end

function pdf(C::BetaCopula{N}, x::AbstractVector{<:Real}) where {N}
    @assert length(x) == 2
    @assert all(0 .≤ x .≤ 1)

    u,v = x
    U,V = (C.U,C.V)
    h = C.h

    c = 0.
    @fastmath @simd for i in 1:N
        @inbounds c+= 𝐁(U[i], u/h + 1, (1-u)/h +1)*𝐁(V[i], v/h + 1, (1-v)/h +1)
    end

    return c/N
end

