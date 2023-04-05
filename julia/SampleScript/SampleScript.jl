using Pkg
Pkg.activate(".")

# Packages: 
include("Includes.jl")

# Get Random Date 
date = rand(Date(2022,1,1):Day(1):Date(2022,11,30))
hour = rand(1:24)

# Get the curves: 
S = Curve(get_data(date,hour, "Sell"))
D = Curve(get_data(date,hour, "Buy"))

# Plot the curves 
plot_1 = plot(S, color = 1, label = "Supply", legend = :outerright)
plot!(D, color = 2, label = "Demand")
scatter!([𝐐(S,D)],[𝐏(S,D)], label = "Clearing", color = :red)
xlabel!("Quantity")
ylabel!("Price")

savefig(plot_1, "Figures/plot_1.png")

# Plot curves after our bid 
p,q = (50., 1000.)
plot_2 = plot(S, color = 1, label = "Supply Before Bid", legend = :outerright, alpha = 0.25)
plot!(S ⊕ˢ (p,q), color = 3, label = "Supply After Bid", alpha = 0.25)
plot!(D, color = 2, label = "Demand", alpha = 0.25)
scatter!([𝐐(S,D)],[𝐏(S,D)], label = "Clearing Before Bid", color = :red)
scatter!([𝐐(S ⊕ˢ (p,q),D)],[𝐏(S ⊕ˢ (p,q),D)], label = "Clearing After Bid", color = :yellow)
title!("Impact of Bid Price 50 and Quantity 100")
xlabel!("Quantity")
ylabel!("Price")

savefig(plot_2, "Figures/plot_2.png")

# Plot the impact
𝐩 = LinRange(-1000, 𝐏(S,D)*1.3, 101)
𝐪 = LinRange(0, 10000, 101)
Z = [impact(S,D, (p,q)) for p ∈ 𝐩, q ∈ 𝐪]
plot_3 = heatmap(𝐪,𝐩, Z)
ylabel!("Price of Bid")
xlabel!("Quantity of Bid")
title!("Impact on Day-Ahead Price (Before EUPHEMIA)")
hline!([𝐏(S,D)], label = "Unimpacted Price", ls = :dash, color = :black)

savefig(plot_3, "Figures/plot_3.png")