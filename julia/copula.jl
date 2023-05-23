using Pkg
Pkg.activate(".")

# Packages and functions: 
using Distributions, Embeddings, StatsBase, LaTeXStrings, Latexify, Base.Threads, ProgressMeter
using EmpiricalCopulas, Chain, DataFramesMeta, ForwardDiff, Interpolations, BivariateCopulas, DiscretizedCopulas

include("SampleScript\\Includes.jl")
include("Mollifiers.jl")
include("functions.jl")
include("DiscretizedDistributions.jl")

from_date = Date(2022,1,1)
to_date = Date(2022,1,1)
dates = from_date:Day(1):to_date
hours = 12

comb_prices = CombinePriceDF(dates, hours, "Sell")
price_pp = comb_prices.DF.Price
quantity_pp = comb_prices.DF.Quantity
bid_df = comb_prices.DF
fixed_bid = comb_prices.fixed_bid
S = Curve(push!(copy(bid_df), fixed_bid[1,:]))
D = GetOtherCurve(dates, hours, "Buy")

Supply = Curve(bid_df)
𝐛 = bids(Supply)
n = convert(Int, round(comb_prices.n/comb_prices.k))

# Plot Bids 
scatter(𝐛, label = "") 
xlabel!("Price") 
ylabel!("Quantity")

mollifier_tolerance = 10.

F_p = GetDensity(price_pp, comb_prices.k, comb_prices.n, mollifier_tolerance)
G_q = ecdf(quantity_pp)
# F_p⁻¹ = InverseDensity(cdf(F_p, price_pp), price_pp)
G_q⁻¹ = x -> quantile(quantity_pp, x)

p10 = plot(price_pp, x -> cdf(F_p,x), label=L"\hat F(p)", ylabel="Cumulative probability", xlabel="Price")
p11 = plot(sort(quantity_pp), x -> G_q(x), label=L"\hat G(q)", xlabel="Quantity")
plot_CDFs = plot(p10,p11, layout=(1,2), size=(700,300), left_margin=3mm, bottom_margin=3mm)
# savefig(plot_CDFs, "Figures/estimated_cdfs.pdf")

# p12 = plot([0.01:0.01:1...], F_p⁻¹, label=L"\hat F^{-1}(p)")
# p13 = plot([0.01:0.01:1...], G_q⁻¹, label=L"\hat G^{-1}(q)")
# plot!(xlabel="Cumulative probability", ylabel="Price and quantity")

U = cdf(F_p, price_pp)
V = G_q(quantity_pp)

# Fit copula
C = BetaCopula(U, V)
@time C = DiscretizedCopula{:PDF}(C, 300)

Z = [pdf(C,[u,v]) for u in LinRange(0.01,0.99,101), v in LinRange(0.01,0.99,101)]
plot_hmap = heatmap(Z)
plot!(xticks=(1:20:101, 0:0.2:1), yticks=(1:20:101, 0:0.2:1))
plot!(xlabel="u", ylabel="v")
surface(Z)
# savefig(plot_hmap, "Figures/heatmap_copula_density.pdf")

Intensities = GetIntensities(price_pp, comb_prices.k, comb_prices.n, mollifier_tolerance)
λₚ = Intensities.λ
Λₚ⁻¹ = Intensities.Λ⁻¹
ogata_p = OgataThinning(price_pp, λₚ).sim
inv_p = SimulateByInversion(Λₚ⁻¹, price_pp[end])

ogata_q = SimulateQuantity(ogata_p, F_p, G_q⁻¹)
inv_q = SimulateQuantity(inv_p, F_p, G_q⁻¹)

p20 = scatter(𝐛, label = "True supply bids", xlabel="Price", ylabel="Quantity") 
scatter!(ogata_p,ogata_q, color=2, label = "Ogata simulation")
p21 = scatter(𝐛, label = "True supply bids", xlabel="Price") 
scatter!(inv_p,inv_q, color=3, label = "Simulation by inversion")
plot_sim_bids = plot(p20,p21, layout=(1,2), size=(600,300))
# savefig(plot_sim_bids, "Figures/sim_supply_bids.pdf")

Supply_ogata = DataFrame(:Price => ogata_p, :Quantity => ogata_q, :Curve => :Supply) |> Curve
Supply_inv = DataFrame(:Price => inv_p, :Quantity => inv_q, :Curve => :Supply) |> Curve

# plot(Supply, color = 1, label = "True supply curve")
# plot!(Supply_ogata, color = 2, label = "Ogata simulation")
# plot!(Supply_inv, color = 3, label = "Simulation by inversion")
# plot!(xlabel="Quantity", ylabel="Price")

ogata_sum_q = Float64[]
inv_sum_q = Float64[]
for _ in 1:1000
    ogata_p = OgataThinning(price_pp, λₚ).sim
    inv_p = SimulateByInversion(Λₚ⁻¹, price_pp[end])
    ogata_q = SimulateQuantity(ogata_p, F_p, G_q⁻¹)
    inv_q = SimulateQuantity(inv_p, F_p, G_q⁻¹)
    push!(ogata_sum_q, sum(ogata_q))
    push!(inv_sum_q, sum(inv_q))
end
mean(ogata_sum_q)
mean(inv_sum_q)
sum(quantity_pp)


plot_ogata_curves = plot(ylabel="Price", bottom_margin=-1mm)
plot!(Supply_ogata, color = 2, alpha = 0.1, label = "Ogata simulation")
map(1:100) do _
    ogata_p = OgataThinning(price_pp, λₚ).sim
    ogata_q = SimulateQuantity(ogata_p, F_p, G_q⁻¹)
    Supply_ogata = DataFrame(:Price => ogata_p, :Quantity => ogata_q, :Curve => :Supply) |> Curve
    plot!(Supply_ogata, color = 2, alpha = 0.1, label = "")
end
plot!(Supply, color = 1, label = "True supply curve", markerstrokewidth=0.5)

plot_inv_curves = plot(xlabel="Quantity", ylabel="Price", top_margin=-1mm)
plot!(Supply_inv, color = 3, alpha = 0.1, label = "Simulation by inversion")
map(1:100) do _
    inv_p = SimulateByInversion(Λₚ⁻¹, price_pp[end])
    inv_q = SimulateQuantity(inv_p, F_p, G_q⁻¹)
    Supply_inv = DataFrame(:Price => inv_p, :Quantity => inv_q, :Curve => :Supply) |> Curve
    plot!(Supply_inv, color = 3, alpha = 0.1, label = "")
end
plot!(Supply, color = 1, label = "True supply curve", markerstrokewidth=0.5)

plot_sim_curves = plot(plot_ogata_curves,plot_inv_curves, layout=(2,1), size=(600,600))
# savefig(plot_sim_curves, "Figures/100_simulated_curves.pdf")

# ======================================
# Monte carlo

simulated_curves = map(1:1000) do _
    # ogata = OgataThinning(price_pp, λₚ).sim
    inversion = SimulateByInversion(Λₚ⁻¹, price_pp[end])

    X̂ = inversion
    Û = cdf(F_p, X̂)
    get_v = u -> h⁻¹(C, rand(), u)
    V̂ = get_v.(Û)
    Ŷ = G_q⁻¹(V̂)

    push!(X̂, fixed_bid.Price[1])
    push!(Ŷ, fixed_bid.Quantity[1])
    DataFrame(:Price => X̂, :Quantity => Ŷ, :Curve => :Supply) |> Curve
end


p_range = -50.:5:250 |> collect
q_range = 100.:10:1000 |> collect

MeanMatrix = zeros(Float64, length(p_range), length(q_range))
# StdMatrix = zeros(Float64, length(q_range), length(p_range))
ExeMatrix = zeros(Float64, length(q_range), length(p_range))

P = Progress(length(p_range) * length(q_range))
@threads for p_index in eachindex(p_range)
    for q_index in eachindex(q_range)
        Price_Change = Float64[]
        Bid_executed = Bool[]
        for i = 1:1000
            pc = 𝐏(simulated_curves[i] ⊕ˢ (p_range[p_index],q_range[q_index]),D) - 𝐏(simulated_curves[i],D)
            push!(Price_Change, pc)
            be = p_range[p_index] ≤ 𝐏(simulated_curves[i] ⊕ˢ (p_range[p_index],q_range[q_index]),D)
            push!(Bid_executed, be)
        end
        MeanMatrix[p_index,q_index] = mean(Price_Change)
        # StdMatrix[p_index,q_index] = std(Price_Change)
        ExeMatrix[p_index,q_index] = mean(Bid_executed)
        next!(P)
    end
end
MeanMatrix
# StdMatrix
ExeMatrix
SaveMeanMatrix
SaveExeMatrix

plot_impact = heatmap(MeanMatrix', xlabel="Quantity", ylabel="Price", right_margin=3mm)
plot!(xticks=(LinRange(1,length(q_range),5), convert.(Int, round.(LinRange(q_range[begin],q_range[end],5)))))
plot!(yticks=(LinRange(1,length(p_range),5), convert.(Int, round.(LinRange(p_range[begin],p_range[end],5)))))
# savefig(plot_impact, "Figures/impact_heatmap.pdf")

function get_surf(q,p)
    p_index = findfirst(isequal(p), p_range)
    q_index = findfirst(isequal(q), q_range)
    MeanMatrix'[p_index,q_index]
end
surf_impact = surface(q_range, p_range, get_surf, xlabel="Quantity", ylabel="Price", zlabel="Impact", size=(500,400), right_margin=3mm)
# savefig(surf_impact, "Figures/impact_surface.pdf")


heatmap(ExeMatrix)

plot(p_range, ExeMatrix[:,1], label=L"\mathbb{P}(Bid \,\, executed)", xlabel="Price", ylabel="Probability")




p_range = -50.:50:250 |> collect
q_range = 100.:100:1000 |> collect

MeanMatrix = zeros(Float64, length(p_range), length(q_range))
# StdMatrix = zeros(Float64, length(q_range), length(p_range))
ExeMatrix = zeros(Float64, length(p_range), length(q_range))

P = Progress(length(p_range) * length(q_range))
@threads for p_index in eachindex(p_range)
    for q_index in eachindex(q_range)
        Price_Change = Float64[]
        Bid_executed = Bool[]
        for i = 1:1000
            pc = 𝐏(simulated_curves[i] ⊕ˢ (p_range[p_index],q_range[q_index]),D) - 𝐏(simulated_curves[i],D)
            push!(Price_Change, pc)
            be = p_range[p_index] ≤ 𝐏(simulated_curves[i] ⊕ˢ (p_range[p_index],q_range[q_index]),D)
            push!(Bid_executed, be)
        end
        MeanMatrix[p_index,q_index] = mean(Price_Change)
        # StdMatrix[p_index,q_index] = std(Price_Change)
        ExeMatrix[p_index,q_index] = mean(Bid_executed)
        next!(P)
    end
end
MeanMatrix
# StdMatrix
ExeMatrix


1000 - 𝐐(simulated_curves[1] ⊕ˢ (0.,1000.),D) + 𝐐(simulated_curves[1],D)



plot_impact = heatmap(MeanMatrix, xlabel="Quantity", ylabel="Price", right_margin=3mm)
plot!(xticks=(LinRange(1,length(q_range),5), convert.(Int, round.(LinRange(q_range[begin],q_range[end],5)))))
plot!(yticks=(LinRange(1,length(p_range),5), convert.(Int, round.(LinRange(p_range[begin],p_range[end],5)))))


latexify(round.(MeanMatrix, digits=2), env=:table)
latexify(round.(StdMatrix, digits=2), env=:table)
latexify(round.(ExeMatrix, digits=4), env=:table)






# Intensity of marked point processes

F = x -> cdf(F_p, x)

ϕᵋ = Mollifier(10.)
Moll_G_q = ϕᵋ(x -> G_q(x), price_pp)
g_q = DifferentiatePWLinear(quantity_pp, Moll_G_q)
plot(quantity_pp, g_q)

λ_pq = (p,q) -> λₚ(p)*pdf(C,[F(p),Moll_G_q(q)])*g_q(q)
Z_intensity = [λ_pq(p,q) for p in -500:5:3000, q in 0:1:801]
heatmap(Z_intensity')
surface(-400:10:1000, 1:10:701, λ_pq)