#= Abstract Types =# 
const _DOMERR = "x ∉ [0,1]²"
abstract type AbstractBivariateCopula <: ContinuousMultivariateDistribution end

length(::AbstractBivariateCopula) = 2
eltype(C::AbstractBivariateCopula) = Real
@memoize function known(A::Type{C}) where {C <: AbstractBivariateCopula} 
    methodswith(A) |> x -> map(string,x) |> x-> map(u -> match(r"^[^\(]+", u).match, x)
end

@memoize function known(A::C) where {C <: AbstractBivariateCopula} 
    known(typeof(A).name.wrapper)
end

(C::AbstractBivariateCopula)(x::Vector{R}) where {R} = cdf(C,x)
(C::AbstractBivariateCopula)(x::Vararg{R}) where {R} = cdf(C,collect(x))

pdf(C::AbstractBivariateCopula, x::AbstractVector{R}) where {R} = "cdf" ∈ known(C) ? ForwardDiff.hessian(u -> cdf(C,u), x)[1,2] : throw("CDF is not known!")

_logpdf(C::AbstractBivariateCopula, x::AbstractArray) = log(pdf(C,x))

cdf(C::AbstractBivariateCopula, x::AbstractVector{R}) where {R} = begin
    any(x .≈ zero(R)) ? (return 0.) : nothing
    prod(x) ≈ one(R) ? (return 1.) : nothing

    "pdf" ∈ known(C) ? (size(x) == (2,)) & all(zero(R) .≤ x .≤ one(R)) ? hcubature(u -> pdf(C,u),zeros(R,2),x)[1] : throw(_DOMERR) : throw("PDF is not known!")
end

function gradient(C::AbstractBivariateCopula, x::AbstractVector{R}) where {R}  
    initialized = known(C)
    if "cdf" ∈ initialized 
        ∇C = ForwardDiff.gradient(u -> cdf(C,u),x)
    elseif "pdf" ∈ initialized
        ∂ᵤC = hquadrature(y -> pdf(C,[x[1],y]), 0., x[2])[1]
        ∂ᵥC = hquadrature(y -> pdf(C,[y,x[2]]), 0., x[1])[1]
        ∇C = [∂ᵤC, ∂ᵥC]
    end

    ∇C
end

function partial(C::AbstractBivariateCopula, x::AbstractVector{R}, i::Int) where {R}
    gradient(C,x)[i]
end

function h(C::AbstractBivariateCopula, x::AbstractVector{R}) where {R}
    gradient(C,x)[2]
end

function h(C::AbstractBivariateCopula, u::R,v::R) where {R}
    h(C,[u,v])
end

function h⁻¹(C::AbstractBivariateCopula, w::R, v::R; ε = 1e-8) where {R <: Real}
    (w ≈ zero(R)) | (w ≈ one(R)) ? (return w) : nothing

    a = 0.
    cₐ = 0.
    b = 1.
    cᵦ = 1.

    while b - a > ε

        m = (a+b)/2
        c = h(C,m,v)

        if w > c
            a = m
            cₐ = c
        else
            b = m
            cᵦ = c
        end
    end

    (a+b)/2
end


function _h_newton(C::AbstractBivariateCopula, x::R, p::R ; i::Int = 2, ε = 1e-8) where {R}
    q′ = Inf
    q = 0.5 #should be mode. 
    count = 0

    (p ≈ zero(R)) | (p ≈ one(R)) ? (return p) : nothing

    while (count < 100) & (abs(q - q′) > ε)
        count += 1
        q′ = q
        q = q + (p - partial(C,[q,x],i))/pdf(C,[q,x])
    end

    q
end

function _h_bisection(C::AbstractBivariateCopula, x::R, p::R ; i::Int = 2, ε = 1e-8) where {R}
    (p ≈ zero(R)) | (p ≈ one(R)) ? (return p) : nothing


    a = 0.
    cₐ = 0.
    b = 1.
    cᵦ = 1.


    while b - a > ε

        m = (a+b)/2
        c = partial(C,[m,x], i)

        if p > c
            a = m
            cₐ = c
        else
            b = m
            cᵦ = c
        end
    end

    (a+b)/2
end

Distributions._rand!(rng::AbstractRNG, C::AbstractBivariateCopula, x::AbstractVector{R}) where {R <: Real} = begin 
    v, p = rand(rng, 2)
    x[1] = _h_bisection(C,v,p)
    x[2] = v
    return x
end

conditional_sampler(rng::AbstractRNG, C::AbstractBivariateCopula, v::Real) = begin
    p = rand(rng)
    x = [_h_bisection(C,v,p), v]
    return x
end

conditional_sampler(rng::AbstractRNG, C::AbstractBivariateCopula, V::Vector{<:Real}) = begin
    U = map(v -> conditional_sampler(rng,C,v),V)
    return (U,V)
end

conditional_sampler(C::AbstractBivariateCopula, v::Real) = begin
    p = rand()
    _h_bisection(C,v,p)
end

conditional_sampler(C::AbstractBivariateCopula, V::Vector{<:Real}) = begin
    map(v -> conditional_sampler(C,v),V)
end

