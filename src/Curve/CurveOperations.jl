# Price and Quantity Functionals for the Supply Curve
function 𝐏(S::Curve{:Supply,R,Q},D::Curve{:Demand,R,Q}) where {R,Q}
    𝐪 = [endpoints(S)..., endpoints(D)...] |> sort

    if length(𝐪) == 0
        return 0.
    else
        i = 1

        while (D(𝐪[i]) - S(𝐪[i]) ≥ 0 ) & (i < length(𝐪))
            i += 1
        end

        return S(𝐪[i-1])
    end
end

function 𝐐(S::Curve{:Supply,R,K},D::Curve{:Demand,R,K}) where {R,K}
    𝐪 = [endpoints(S)..., endpoints(D)...] |> sort

    if length(𝐪) == 0
        return 0.
    else
        i = 1

        while (D(𝐪[i]) - S(𝐪[i]) ≥ 0 ) & (i < length(𝐪))
            i += 1
        end

        return 𝐪[i-1]
    end
end

# Sorting Operation for the Supply Curve
function ≤ˢ(b::Tuple{R,Q}, c::Tuple{R,Q}) where {R <: Real, Q <: Real}
    if b[1] < c[1]
        return true
    elseif (b[1] == c[1]) & (b[2] < c[2])
        return true
    else
        return false
    end
end

# Add a supply bid
function ⊕ˢ(B::Vector{Tuple{R,Q}}, b::Tuple{R,Q}) where {R,Q}
    i = searchsortedlast(B,b, lt = ≤ˢ)
    i == 0 ? (return [b, B...]) : i == length(B) ? [B..., b] : [B[1:i]..., b, B[(i+1):end]...]
end

function ⊕ˢ(S::Curve{:Supply, R, Q}, b::Tuple{R,Q}) where {R <: Real,Q <:Real}
    B = bids(S)
    B′ = B ⊕ˢ b
    return Curve(B′, "Supply")
end

# Obtain the Impact from adding a new bid to the supply curve 
function impact(S::Curve{:Supply,R,Q},D::Curve{:Demand,R,Q},b::Tuple{R,Q}) where {R <: Real, Q <: Real}
    S′ = S ⊕ˢ b    
    return 𝐏(S′,D) - 𝐏(S,D)
end

# Sorting Operations for the Demand Curve 
function ≤ᵈ(b::Tuple{R,Q}, c::Tuple{R,Q}) where {R <: Real, Q <: Real}
    if b[1] > c[1]
        return true
    elseif (b[1] == c[1]) & (b[2] < c[2])
        return true
    else
        return false
    end
end
