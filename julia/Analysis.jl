using Pkg
Pkg.activate(".")

# Packages and functions: 
using Distributions, Embeddings, StatsBase, LaTeXStrings
using EmpiricalCopulas, Chain, DataFramesMeta, ForwardDiff, Interpolations, BivariateCopulas

include("SampleScript\\Includes.jl")
include("Mollifiers.jl")
include("functions.jl")
include("DiscretizedDistributions.jl")

# Nonparametric estimation of the intensity

from_date = Date(2022,1,1)
to_date = Date(2022,1,1)
dates = from_date:Day(1):to_date
hours = 12
side = "Sell"
mollifier_tolerance = 10

comb_prices = CombinePriceDF(dates, hours, side)
price_pp = comb_prices.DF.Price
CumuIntensity = PWLinearEstimator((pp = price_pp, k = comb_prices.k, n = comb_prices.n))
EstIntensity = DifferentiatePWLinear(price_pp, CumuIntensity)
Intensity = GetIntensities(price_pp, comb_prices.k, comb_prices.n, mollifier_tolerance)

Λₚ = Intensity.Λ
λₚ = Intensity.λ
Λₚ⁻¹ = Intensity.Λ⁻¹

# Plot intensity

pₙ = price_pp[price_pp .!== -0.0]

p8 = plot(pₙ, CumuIntensity, label=L"\hat\Lambda(p)", legend=:right)
plot!(Λₚ, label=L"\hat\Lambda(p)*\varphi_\varepsilon")
scatter!(pₙ, zero, color=1, label = L"p_1,\ldots,p_n", marker = (:circle, 3), markerstrokewidth=0.5)
p9 = plot(EstIntensity, xlims=extrema(price_pp).*1.05, label=L"\hat\lambda(p)")
plot!(λₚ, label=L"\partial(\hat\lambda(p)*\varphi_\varepsilon)")
xlabel!("Price")
# xlims!(-5,500)

plot9 = plot(p8,p9, layout=(2,1), size=(600,600), ylabel="Intensity")

# savefig(plot9, "Figures/price_intensity.pdf")

# Test intensity

s = Λₚ.(pₙ)
Δs = [s[1], diff(s)...]
h20 = histogram(Δs, normalize = true, label=L"\Delta \tau")
plot!(t -> (t ≥ 0)*exp(-t), label=L"PDF\,\, of\,\, Exp(1)")
plot!(title=L"Testing\,\, \hat\Lambda(p)*\varphi_\epsilon",title_align=:left)
q20 = qqplot(Δs, Exponential(1), xlabel="Δτ", ylabel="Exp(1)", markerstrokewidth=0.5)
plot20 = plot(h20,q20, layout=(1,2), size=(700,300), bottom_margin=3mm)

# savefig(plot20, "Figures/price_residual_analysis.pdf")

# Resample and validate intensity

ogata = OgataThinning(pₙ, λₚ)
inv = SimulateByInversion(Λₚ⁻¹, pₙ[end])

plot21 = plot(yaxis=false, ylims=(-0.3,2.5), size=(600,200), xlabel="Price", bottom_margin=3mm)
scatter!(pₙ, one, color=1, marker=(:circle,3), label="True price process")
scatter!(ogata.sim, half, color=2, marker = (:circle,3), label = "Ogata simulation")
scatter!(inv, zero, color=3, marker = (:circle,3), label = "Simulation by inversion")

bins = pₙ[begin]:((pₙ[end]+1-pₙ[begin])/100):(pₙ[end]+1)
h21 = histogram(pₙ, bins = bins, alpha = 0.7, label = "Prices", ylabel="Frequency", bottom_margin=-1mm)
histogram!(ogata.sim, bins = bins, alpha = 0.7, label = "Ogata simulation")

h22 = histogram(pₙ, bins = bins, alpha = 0.7, label = "Prices", xlabel="Price", ylabel="Frequency", top_margin=-1mm)
histogram!(inv, color=3, bins = bins, alpha = 0.7, label = "Simulation by inversion")

# plot22 = plot(h21,h22, layout=(1,2), size=(700,300), bottom_margin=3mm, left_margin=3mm)
plot22 = plot(h21,h22, layout=(2,1), size=(600,500))

# savefig(plot21, "Figures/price_simulations.pdf")
# savefig(plot22, "Figures/price_simulations_hist.pdf")

ogata_length = Int64[]
inv_length = Int64[]
for _ in 1:1000
    ogata = OgataThinning(pₙ, λₚ)
    inv = SimulateByInversion(Λₚ⁻¹, pₙ[end])

    push!(ogata_length, length(ogata.sim))
    push!(inv_length, length(inv))
end
mean(ogata_length)
mean(inv_length)
1.96*std(ogata_length)
1.96*std(inv_length)
length(pₙ)/(length(dates)*length(hours))
