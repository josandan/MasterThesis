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

mollifier_tolerance = 10.

F_p = GetDensity(price_pp, comb_prices.k, comb_prices.n, mollifier_tolerance)
# F_p = ecdf(price_pp)
F_q = ecdf(quantity_pp)
F_p⁻¹ = InverseDensity(cdf(F_p, price_pp), price_pp)
# F_p⁻¹ = x -> quantile(price_pp, x)
F_q⁻¹ = x -> quantile(quantity_pp, x)

p10 = plot(price_pp, x -> cdf(F_p,x), label=L"\hat F(p)")
plot!(sort(quantity_pp), x -> F_q(x), label=L"\hat F(q)")
plot!(ylabel="Cumulative probability", xlabel="Price and quantity")
p11 = plot([0.01:0.01:1...], F_p⁻¹, label=L"\hat F^{-1}(p)")
plot!(F_q⁻¹, label=L"\hat F^{-1}(q)")
plot!(xlabel="Cumulative probability", ylabel="Price and quantity")
plot_CDFs = plot(p10,p11, layout=(1,2), size=(600,300))
# savefig(plot_CDFs, "Figures/estimated_cdfs_w_inv.pdf")

U = cdf(F_p, price_pp)
V = F_q(quantity_pp)

# Fit copula
C = BetaCopula(U, V)
@time C = DiscretizedCopula{:PDF}(C, 300)

# gradient(C, [0.5,0.5])

# Plot
Z = [pdf(C,[u,v]) for u in LinRange(0.01,0.99,101), v in LinRange(0.01,0.99,101)]
plot_hmap = heatmap(Z)
# plot!(xticks=(1:100/7:101, xticks), yticks=(1:100/7:101, yticks))
plot!(xticks=(1:20:101, 0:0.2:1), yticks=(1:20:101, 0:0.2:1))
# plot!(xlabel="Price", ylabel="Quantity")
plot!(xlabel="u", ylabel="v")
# savefig(plot_hmap, "Figures/heatmap_copula_density.pdf")

W = rand(C, n)' 
Û = W[:,1] 
V̂ = W[:,2] 

X̂ = convert.(Float64, F_p⁻¹.(Û))
Ŷ = F_q⁻¹.(V̂) 

plot_sim_bids = scatter(𝐛, label = "True supply bids") 
xlabel!("Price") 
ylabel!("Quantity") 
scatter!(X̂,Ŷ, label = "Inverse transform simulation")
# savefig(plot_sim_bids, "Figures/sim_supply_bids_alt.pdf")

Supply₀ = DataFrame(:Price => X̂, :Quantity => Ŷ, :Curve => :Supply) |> Curve

plot(Supply, color = 1, label = "True supply curve")
plot!(Supply₀, color = 2, label = "Inverse transform simulation")
xlabel!("Quantity")
ylabel!("Price")

# Plot 100 realizations

plot_its_curves = plot(xlabel="Quantity", ylabel="Price")
plot!(Supply₀, color = 2, alpha = 0.1, label = "Inverse transform simulation")
map(2:100) do _
    W = rand(C, n)' 
    Û = W[:,1] 
    V̂ = W[:,2] 
    X̂ = F_p⁻¹.(Û) 
    Ŷ = F_q⁻¹.(V̂) 
    Supply₀ = DataFrame(:Price => X̂, :Quantity => Ŷ, :Curve => :Supply) |> Curve
    plot!(Supply₀, color = 2, alpha = 0.1, label = "")
end
plot!(Supply, color = 1, label = "True supply curve", markerstrokewidth=0.5)
# savefig(plot_its_curves, "Figures/100_simulated_curves_alt.pdf")


sum_q = Float64[]
for _ in 1:1000
    W = rand(C, n)' 
    Û = W[:,1] 
    V̂ = W[:,2] 
    X̂ = F_p⁻¹.(Û) 
    Ŷ = F_q⁻¹.(V̂) 
    push!(sum_q, sum(Ŷ))
end
m = mean(sum_q)
quantile(sum_q, 0.025)
quantile(sum_q, 0.975)
sum(quantity_pp)


# Calculate error between true curve and simulated curves
# DF = transform(bid_df, :Quantity => cumsum => :Quantity)
# error = Float64[]
# Supply_Curves = Curve[]
# for _ in 1:1000
#     W = rand(C, n)' 
#     Û = W[:,1] 
#     V̂ = W[:,2] 
#     X̂ = F_p⁻¹.(Û) 
#     Ŷ = F_q⁻¹.(V̂) 

#     Supply₀ = DataFrame(:Price => X̂, :Quantity => Ŷ, :Curve => :Supply) |> Curve
#     push!(Supply_Curves, Supply₀)

#     DF_sim = @chain DataFrame(:Price => X̂, :Quantity => Ŷ) begin
#         sort!(:Price)
#         transform!(:Quantity => cumsum => :Quantity)
#     end

#     e = sqrt.((DF_sim.Price - DF.Price).^2 + (DF_sim.Quantity - DF.Quantity).^2) |> mean
#     push!(error, e)
# end

# error
# histogram(error)


# ======================================
# Monte carlo

Supply_Curves = Curve[]
for _ in 1:1000
    W = rand(C, n)' 
    Û = W[:,1] 
    V̂ = W[:,2] 
    X̂ = convert.(Float64, F_p⁻¹.(Û))
    Ŷ = F_q⁻¹.(V̂) 

    push!(X̂, fixed_bid.Price[1])
    push!(Ŷ, fixed_bid.Quantity[1])
    Supply₀ = DataFrame(:Price => X̂, :Quantity => Ŷ, :Curve => :Supply) |> Curve
    push!(Supply_Curves, Supply₀)
end
Supply_Curves

p,q = (0., 100.)
Price_Change = Float64[]
for i = 1:1000
    pc = 𝐏(Supply_Curves[i] ⊕ˢ (p,q),D) - 𝐏(Supply_Curves[i],D)
    push!(Price_Change, pc)
end
Price_Change
Price_Change |> mean
Price_Change |> std

histogram(Price_Change, label = "Change in MCP")
Price_Change[Price_Change .!= 0] |> length
Price_Change[Price_Change .<= -1] |> length

plot(D, color = 1, label = "Demand curve")
plot!(S, color = 2, label = "Supply curve")
plot!(Supply₀, color = 3, label = "Simulated supply curve")


MeanMatrix = zeros(Float64, 10, 7)
StdMatrix = zeros(Float64, 10, 7)
ExeMatrix = zeros(Float64, 10, 7)
@time for (p_index,p) in pairs([-50.:50:250...])
    for (q_index,q) in pairs([100.:100:1000...])
        Price_Change = Float64[]
        Bid_executed = Bool[]
        for i = 1:1000
            pc = 𝐏(Supply_Curves[i] ⊕ˢ (p,q),D) - 𝐏(Supply_Curves[i],D)
            push!(Price_Change, pc)
            be = p ≤ 𝐏(Supply_Curves[i] ⊕ˢ (p,q),D)
            push!(Bid_executed, be)
        end
        MeanMatrix[q_index,p_index] = mean(Price_Change)
        StdMatrix[q_index,p_index] = std(Price_Change)
        ExeMatrix[q_index,p_index] = mean(Bid_executed)
    end
end
MeanMatrix
StdMatrix
ExeMatrix

# (p,q) = (100., 1000.)
# i = 50
# 𝐏(Supply_Curves[i] ⊕ˢ (p,q),D) - 𝐏(Supply_Curves[i],D)
# 𝐏(Supply_Curves[i] ⊕ˢ (p,q),D)

# plot(Supply_Curves[i], ylims=(0,200), color=1)
# plot!(D, color=2)
# plot!(Supply_Curves[i] ⊕ˢ (p,q), ylims=(0,200), color=3)

latexify(round.(MeanMatrix, digits=2), env=:table)
latexify(round.(StdMatrix, digits=2), env=:table)
latexify(round.(ExeMatrix, digits=4), env=:table)

TrueMatrix = zeros(Float64, 10, 7)
TrueExeMatrix = zeros(Float64, 10, 7)
for (p_index,p) in pairs([-50.:50:250...])
    for (q_index,q) in pairs([100.:100:1000...])
        pc = 𝐏(S ⊕ˢ (p,q),D) - 𝐏(S,D)
        be = p ≤ 𝐏(S ⊕ˢ (p,q),D)
        TrueMatrix[q_index,p_index] = pc
        TrueExeMatrix[q_index,p_index] = be
    end
end
TrueMatrix
TrueExeMatrix
