function ϕ(C::Curve{T,R,D}) where {T,R <:Real,D<:Real}
    Q = domain(C)[2]
    t::Real -> Q*t
end

function ξ(C::Curve{T,R,D}) where {T,R,D}
    t::Real -> [(C ∘ ϕ(C))(t), ϕ(C)(t)]
end

