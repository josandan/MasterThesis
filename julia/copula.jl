using Pkg
Pkg.activate(".")

# Packages and functions: 
using Distributions, Embeddings, StatsBase
using EmpiricalCopulas, Chain, DataFramesMeta, ForwardDiff, Interpolations, BivariateCopulas, DiscretizedCopulas

include("SampleScript\\Includes.jl")
include("Mollifiers.jl")
include("functions.jl")
include("DiscretizedDistributions.jl")



from_date = Date(2022,1,1)
to_date = Date(2022,1,5)
dates = from_date:Day(1):to_date
hours = 12
side = "Sell"

comb_prices = CombinePriceDF(dates, hours, side)
price_pp = comb_prices.DF.Price
quantity_pp = comb_prices.DF.Quantity
bid_df = comb_prices.DF

Supply = Curve(bid_df)
𝐛 = bids(Supply)
n = comb_prices.n/comb_prices.k |> round
n = convert(Int, n)

# Plot Bids 
scatter(𝐛, label = "") 
xlabel!("Price") 
ylabel!("Quantity")

mollifier_tolerance = 10.

F_p = GetDensity(price_pp, comb_prices.k, comb_prices.n, mollifier_tolerance)
# F_p = ecdf(price_pp)
F_q = ecdf(quantity_pp)

F_p⁻¹ = InverseDensity(cdf(F_p, price_pp), price_pp)
# F_p⁻¹ = x -> quantile(price_pp, x)
F_q⁻¹ = x -> quantile(quantity_pp, x)

plot(price_pp, x -> cdf(F_p,x))
# plot(price_pp, x -> F_p(x))
plot!(sort(quantity_pp), x -> F_q(x))

# rand(F_q,1000)

# Fit copula
C = BetaCopula(cdf(F_p, price_pp), F_q(quantity_pp))
# C = BetaCopula(F_p(price_pp), F_q(quantity_pp))
C = DiscretizedCopula{:PDF}(C, 300)

gradient(C, [0.5,0.5])

# Plot
Z = [pdf(C,[u,v]) for u in LinRange(0.01,0.99,101), v in LinRange(0.01,0.99,101)]
heatmap(Z)

W = rand(C, 739)' 
Û = W[:,1] 
V̂ = W[:,2] 

X̂ = F_p⁻¹.(Û) 
Ŷ = F_q⁻¹.(V̂) 

scatter(𝐛, label = "") 
xlabel!("Price") 
ylabel!("Quantity") 
scatter!(X̂,Ŷ, label = "")

Supply₀ = DataFrame(:Price => X̂, :Quantity => Ŷ, :Curve => :Supply) |> Curve

plot(Supply, color = 1)
plot!(Supply₀, color = 2)
xlabel!("Quantity")
ylabel!("Price")

