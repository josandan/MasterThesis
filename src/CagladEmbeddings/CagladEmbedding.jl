struct LeftOpenInterval{R}
    a::R
    b::R
end

∈(x::Q, I::LeftOpenInterval{R}) where {Q <: Real, R <: Real} = begin
    return I.a < x ≤ I.b
end

endpoint(I::LeftOpenInterval) = I.b

struct CagladEmbedding{R,Q} <: Function
    p::Vector{R}
    q::Vector{LeftOpenInterval{Q}}
end

function CagladEmbedding(p::Vector{R}, q::Vector{Q}; lower = zero(Q)) where {R, Q <: Real} 
    @assert length(p) == length(q) 

    𝐪 = LeftOpenInterval{Q}[]
    for i in eachindex(q)
        if i == 1 
            push!(𝐪, LeftOpenInterval(lower, q[1]))
        else
            push!(𝐪, LeftOpenInterval(q[i-1],q[i]))
        end
    end

    CagladEmbedding{R,Q}(p, 𝐪)
end

function (C::CagladEmbedding{R,Q})(t::S) where {R <: Real,Q <: Real, S <: Real}
    q = C.q
    p = C.p
     
    for i in eachindex(q)
        if t ∈ q[i]
            return p[i] 
        end
    end

    return zero(R)
end

endpoints(C::CagladEmbedding) = begin
    map(endpoint, C.q)
end
