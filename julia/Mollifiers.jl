using ForwardDiff

struct Mollifier{R} 
    ε::R
end

function (ϕᵋ::Mollifier)(x::R) where {R <: Real}
    ε = ϕᵋ.ε
    1/ε*(exp(-1/(1-abs(x/ε)^2))/(0.443994)*(abs(x/ε) < 1))
end

function (ϕᵋ::Mollifier)(f::Function, Π; m = 100)
    # Get mollification epsilon
    ε = ϕᵋ.ε

    # Set mollification range 
    S = LinRange(-ε, ε, m)
    Δs = diff(S)[1]

    # Set up extension 
    a,b = Π[1], Π[end]
    n = length(Π)
    f̄(t) = a ≤ t ≤ b ? f(t) : t < a ? f(a) : f(b)

    ξ = fill(0., n)

    for (i, t) in enumerate(Π)
        χ = 0.
        for s in S
            χ += f̄(t - s)*ϕᵋ(s)*Δs
        end
        ξ[i] = χ
    end

    return Embeddings.LinearEmbedding(ξ,Π)
end

function ∂(ϕᵋ::Mollifier) 
    (f::Function, Π; m = 100) -> begin
       
        # Get mollification epsilon
        ε = ϕᵋ.ε

        # Set Mollifer Derivative 
        Dϕᵋ = t -> ForwardDiff.derivative(s -> ϕᵋ(s),t)

        # Set mollification range 
        S = LinRange(-ε, ε, m)
        Δs = diff(S)[1]

        # Set up extension 
        a,b = Π[1], Π[end]
        n = length(Π)
        f̄(t) = a ≤ t ≤ b ? f(t) : t < a ? f(a) : f(b)

        ξ = fill(0., n)

        for (i, t) in enumerate(Π)
            χ = 0.
            for s in S
                χ += f̄(t - s)*Dϕᵋ(s)*Δs
            end
            ξ[i] = χ
        end

        return Embeddings.LinearEmbedding(ξ,Π)  
    end
end
