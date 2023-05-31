using Pkg
Pkg.activate(".")

# Packages and functions: 
include("Includes.jl")
include("Mollifiers.jl")
include("DiscretizedDistributions.jl")
include("BetaCopula.jl")
include("functions.jl")

from_date = Date(2022,1,1)
to_date = Date(2022,1,1)
dates = from_date:Day(1):to_date
hours = 12
side = "Sell"

comb_prices = CombinePriceDF(dates, hours, side)
price_pp = comb_prices.DF.Price
quantity_pp = comb_prices.DF.Quantity
bid_df = comb_prices.DF
fixed_bid = comb_prices.fixed_bid
S = Curve(push!(copy(bid_df), fixed_bid[1,:]))
D = GetOtherCurve(dates, hours, side == "Sell" ? "Buy" : "Sell")

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
plot_hmap = heatmap(Z')
scatter!(100*U, 100*V, label="", xlims=(0.3,100), ylims=(0.3,99.8))
plot!(xticks=(LinRange(0.5,100,6), 0:0.2:1), yticks=(LinRange(0.5,99,6), 0:0.2:1))
plot!(xlabel="u", ylabel="v")
# surface!(Z')
plot_contour = contour(Z', fill=true)
scatter!(100*U, 100*V, label="", xlims=(0.7,100), ylims=(0.7,99.8))
plot!(xticks=(LinRange(0.7,100,6), 0.0:0.2:1), yticks=(LinRange(0.7,99,6), 0:0.2:1))
plot!(xlabel="u", ylabel="v")

# savefig(plot_hmap, "Figures/heatmap_copula_density.pdf")
# savefig(plot_contour, "Figures/contour_copula_density.pdf")

Intensities = GetIntensities(price_pp, comb_prices.k, comb_prices.n, mollifier_tolerance)
λₚ = Intensities.λ
Λₚ⁻¹ = Intensities.Λ⁻¹
ogata_p = OgataThinning(price_pp, λₚ).sim
inv_p = SimulateByInversion(Λₚ⁻¹, price_pp[end])

ogata_q = SimulateQuantity(ogata_p, F_p, G_q⁻¹, C)
inv_q = SimulateQuantity(inv_p, F_p, G_q⁻¹, C)

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
    ogata_q = SimulateQuantity(ogata_p, F_p, G_q⁻¹, C)
    inv_q = SimulateQuantity(inv_p, F_p, G_q⁻¹, C)
    push!(ogata_sum_q, sum(ogata_q))
    push!(inv_sum_q, sum(inv_q))
end
mean(ogata_sum_q)
mean(inv_sum_q)
quantile(ogata_sum_q, 0.025)
quantile(ogata_sum_q, 0.975)
quantile(inv_sum_q, 0.025)
quantile(inv_sum_q, 0.975)
sum(quantity_pp)


plot_ogata_curves = plot(ylabel="Price", bottom_margin=-1mm)
plot!(Supply_ogata, color = 2, alpha = 0.1, label = "Ogata simulation")
map(1:100) do _
    ogata_p = OgataThinning(price_pp, λₚ).sim
    ogata_q = SimulateQuantity(ogata_p, F_p, G_q⁻¹, C)
    Supply_ogata = DataFrame(:Price => ogata_p, :Quantity => ogata_q, :Curve => :Supply) |> Curve
    plot!(Supply_ogata, color = 2, alpha = 0.1, label = "")
end
plot!(Supply, color = 1, label = "True supply curve", markerstrokewidth=0.5)

plot_inv_curves = plot(xlabel="Quantity", ylabel="Price", top_margin=-1mm)
plot!(Supply_inv, color = 3, alpha = 0.1, label = "Simulation by inversion")
map(1:100) do _
    inv_p = SimulateByInversion(Λₚ⁻¹, price_pp[end])
    inv_q = SimulateQuantity(inv_p, F_p, G_q⁻¹, C)
    Supply_inv = DataFrame(:Price => inv_p, :Quantity => inv_q, :Curve => :Supply) |> Curve
    plot!(Supply_inv, color = 3, alpha = 0.1, label = "")
end
plot!(Supply, color = 1, label = "True supply curve", markerstrokewidth=0.5)

plot_sim_curves = plot(plot_ogata_curves,plot_inv_curves, layout=(2,1), size=(600,600))
savefig(plot_sim_curves, "Figures/100_simulated_curves.pdf")


# ======================================
# Monte carlo study on market impact

simulated_curves = map(1:1000) do _
    # ogata = OgataThinning(price_pp, λₚ).sim
    X̂ = SimulateByInversion(Λₚ⁻¹, price_pp[end])
    Ŷ = SimulateQuantity(X̂, F_p, G_q⁻¹)

    push!(X̂, fixed_bid.Price[1])
    push!(Ŷ, fixed_bid.Quantity[1])
    DataFrame(:Price => X̂, :Quantity => Ŷ, :Curve => :Supply) |> Curve
end

p_range = -50.:5:250 |> collect
q_range = 10.:10:1000 |> collect

MeanMatrix = zeros(Float64, length(p_range), length(q_range))
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
        ExeMatrix[p_index,q_index] = mean(Bid_executed)
        next!(P)
    end
end
MeanMatrix
ExeMatrix
SaveMeanMatrix 
SaveExeMatrix

plot_impact = heatmap(MeanMatrix, xlabel="Quantity", ylabel="Price", right_margin=3mm)
plot!(xticks=(LinRange(1,length(q_range),5), convert.(Int, round.(LinRange(q_range[begin],q_range[end],5)))))
plot!(yticks=(LinRange(1,length(p_range),5), convert.(Int, round.(LinRange(p_range[begin],p_range[end],5)))))
# savefig(plot_impact, "Figures/impact_heatmap.pdf")

surf_impact = surface(MeanMatrix, size=(500,400), right_margin=3mm)
plot!(xlabel="Quantity", ylabel="Price", zlabel="Impact")
# plot!(xticks=(LinRange(1,length(q_range),5), convert.(Int, round.(LinRange(q_range[begin],q_range[end],5)))))
plot!(xticks=(LinRange(1,length(q_range),5), convert.(Int, [q_range[begin],250:250:q_range[end]...])))
plot!(yticks=(LinRange(1,length(p_range),5), convert.(Int, round.(LinRange(p_range[begin],p_range[end],5)))))
# savefig(surf_impact, "Figures/impact_surface.pdf")

plot_bid_exe = plot(p_range, ExeMatrix[:,1], label=L"\mathbb{P}(Bid \,\, executed)")
plot!(xlabel="Price", ylabel="Probability")
# savefig(plot_bid_exe, "Figures/bids_executed.pdf")


p_range = -50.:50:250 |> collect
q_range = 100.:100:1000 |> collect

MeanMatrix = zeros(Float64, length(p_range), length(q_range))
StdMatrix = zeros(Float64, length(p_range), length(q_range))
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
        StdMatrix[p_index,q_index] = std(Price_Change)
        ExeMatrix[p_index,q_index] = mean(Bid_executed)
        next!(P)
    end
end
MeanMatrix
StdMatrix
ExeMatrix

latexify(round.(MeanMatrix', digits=2), env=:table)
latexify(round.(StdMatrix', digits=2), env=:table)
latexify(round.(ExeMatrix', digits=5), env=:table)



1000 - 𝐐(simulated_curves[1] ⊕ˢ (0.,1000.),D) + 𝐐(simulated_curves[1],D)





# Intensity of marked point processes

mollifier_tolerance = 100.

F_p = GetDensity(price_pp, comb_prices.k, comb_prices.n, mollifier_tolerance)
G_q = ecdf(quantity_pp)
plot(G_q)
G_q_moll = ϕᵋ(x -> G_q(x), quantity_pp)
plot!(G_q_moll)

U = cdf(F_p, price_pp)
V = G_q_moll.(quantity_pp)

C = BetaCopula(U, V)
@time C = DiscretizedCopula{:PDF}(C, 300)

F = x -> cdf(F_p, x)
G = x -> G_q(x)
plot(price_pp, F)
plot(sort(quantity_pp), G)

ϕᵋ = Mollifier(mollifier_tolerance)
Moll_G = ϕᵋ(G, sort(quantity_pp))
plot!(Moll_G)

g1 = ∂(ϕᵋ)(G, sort(quantity_pp))
g2 = DifferentiatePWLinear(quantity_pp, Moll_G_q)

plot(sort(quantity_pp), g1)
plot!(g2)

λₚ = GetIntensities(price_pp, comb_prices.k, comb_prices.n, mollifier_tolerance).λ
plot(price_pp, λₚ)

λ_pq = (p,q) -> λₚ(p)*pdf(C,[F(p),G_q(q)])*g1(q)
Z_intensity = [λ_pq(p,q) for p in -500:5:3000, q in 0:1:801]
heatmap(Z_intensity')
surface(-400:10:3000, 1:10:701, λ_pq)
contour(Z_intensity')

heatmap(price_pp, quantity_pp, λ_pq)
surface(price_pp, quantity_pp, λ_pq)
contour(-500:3000, 1:800, λ_pq)
