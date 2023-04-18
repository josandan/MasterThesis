#= 

DiscretizedDistributions: 

    Discretizes ContinuousUnivariateDistributions to allow for fast evaluations, sampling, etc. using interpolations. 

=# 

# Imports
import Distributions: pdf, cdf, quantile, rand
import ForwardDiff: derivative, gradient
import Dates: month, year
import BivariateCopulas: h


struct DiscretizedDistribution <: ContinuousUnivariateDistribution
    F::LinearEmbedding
end

# function Discretize(D::ContinuousUnivariateDistribution, N::Int; ε = 1e-6)
#     L, U = quantile(D, [ε, 1-ε])
    
#     I = LinRange(L,U,N)

#     F = map(i -> cdf(D,i), I)
#     F = (F .- minimum(F))/(maximum(F) - minimum(F))

#     itp = interpolate(F, BSpline(Quadratic(Line(OnGrid()))))
#     itp = extrapolate(itp, Interpolations.Flat())
#     itp = Interpolations.scale(itp,I)

#     DiscretizedDistribution(itp)
# end

pdf(D::DiscretizedDistribution,x::Real) = derivative(D.F,x)
logpdf(D::DiscretizedDistribution,x) = log(pdf(D,x))
cdf(D::DiscretizedDistribution,x::Real) = D.F(x)
quantile(D::DiscretizedDistribution,p::Real; ε = 1e-6) = begin
    F = D.F
    a, b = extrema(D.F.x)
    Fᵃ, Fᵇ = extrema(D.F.f)

    c = 0.5*(a+b)
    Fᶜ = F(c)

    while abs(Fᶜ - p) > ε

        if Fᶜ - p > 0.
            b = c 
            Fᵇ = Fᶜ
        else
            a = c
            Fᵃ = Fᶜ
        end

        c = 0.5*(a+b)
        Fᶜ = F(c)
    end

    return c
end


rand(D::DiscretizedDistribution,n::Int) = map(p -> quantile(D,p), rand(n))


struct CopulaWrapper <: AbstractBivariateCopula 
    C::ScaledInterpolation
end

cdf(DC::CopulaWrapper, x::AbstractVector{<:Real}) = begin
    @assert length(x) == 2
    @assert all(0 .≤ x .≤ 1)

    DC.C(x...)
end

gradient(DC::CopulaWrapper, x::AbstractVector{<:Real}) = begin
    @assert length(x) == 2
    @assert all(0 .≤ x .≤ 1)
    u,v = x

    [min(max(derivative(u -> DC.C(u,v), u),0.),1), min(max(derivative(v -> DC.C(u,v),v),0.),1.)]
end


pdf(DC::CopulaWrapper, x::AbstractVector{<:Real}; ε = nothing) = begin
    @assert length(x) == 2
    @assert all(0 .≤ x .≤ 1)

    u,v = x
    if isnothing(ε)
        derivative(u -> derivative(v -> DC.C(u,v),v),u)
    else
        Dᵤ = v -> (DC.C(u + ε, v) - DC.C(u - ε, v))/(2*ε)

        return (Dᵤ(v + ε) - Dᵤ(v - ε))/(2*ε)
    end
end

h(DC::CopulaWrapper,x::AbstractVector{<:Real}) = begin
    @assert length(x) == 2
    @assert all(0 .≤ x .≤ 1)

    u,v = x

    return min(max(derivative(v -> DC.C(u,v),v),0.),1.)
end
