struct Curve{S,R,Q} <: Function 
    f::CagladEmbedding{R,Q}
end

domain(C::Curve{S,R,Q}) where {S,R,Q <: Real} = [0., maximum(C.f.q[end].b)]

function (C::Curve{S,R,Q})(q::T) where {R,Q,S, T <: Real}
    q = Q(q)
    𝐐 = length(C.f.q) > 0 ? C.f.q[end].b : 0 
    @assert 0 ≤ q ≤ 𝐐

    C.f(q) 
end

function Curve(dataframe::DataFrame) 
    side = dataframe[1,:Curve] |> Symbol
    I = side == :Supply ? sortperm(dataframe[:,:Price]) : sortperm(dataframe[:,:Price], rev = true )
    
    p = dataframe[I,:Price] 
    q = dataframe[I,:Quantity] |> cumsum

    
    C = CagladEmbedding(p,q)
    Curve{side, eltype(p), eltype(q)}(C)
end

function Curve(B::Vector{Tuple{R,Q}}, side::Union{String,Symbol}) where {R,Q} 
    side = Symbol(side)
    p = [b[1] for b in B]
    I = side == :Supply ? sortperm(p) : sortperm(p, rev = true)

    p = [b[1] for b in B[I]]
    q = [b[2] for b in B[I]] |> cumsum
    C = CagladEmbedding(p,q)

    Curve{side, eltype(p), eltype(q)}(C)
end

function BiddingCollection(dataframe::DataFrame)
    p = dataframe[:,:Price]
    q = dataframe[:,:Quantity]
    (collect ∘ zip)(p,q)
end

endpoints(C::Curve) = endpoints(C.f)

bids(C::Curve{S,R,Q}) where {S,R,Q} = begin
    q = [zero(Q), endpoints(C)...] |> diff 
    p = C.f.p 

    (collect ∘ zip)(p,q)
end