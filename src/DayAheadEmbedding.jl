module DayAheadEmbedding

# Packages 
using Plots 
using DataFrames

# Includes 
include("CagladEmbeddings/CagladEmbedding.jl")
include("Curve/Curve.jl")
include("Curve/CurveOperations.jl")
include("Curve/CurveParameterization.jl")
include("Curve/CurvePlot.jl")

# Exports 
export
    CagladEmbedding,
    Curve,
    BiddingCollection,
    domain,
    𝐏,
    𝐐,
    ≤ˢ,
    ≤ᵈ, 
    ⊕ˢ,
    impact,
    ϕ,
    ξ,
    endpoints,
    bids

end # module
