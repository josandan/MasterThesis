using Pkg
Pkg.activate(".")

# Packages and functions: 
using Distributions, Embeddings, StatsBase, LaTeXStrings
using EmpiricalCopulas, Chain, DataFramesMeta, ForwardDiff, Interpolations, BivariateCopulas, DiscretizedCopulas

include("SampleScript\\Includes.jl")
include("Mollifiers.jl")
include("functions.jl")
include("DiscretizedDistributions.jl")



from_date = Date(2022,1,1)
to_date = Date(2022,1,1)
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

p10 = plot(price_pp, x -> cdf(F_p,x), label=L"\hat F(p)", ylabel="Cumulative probability", xlabel="Price")
p11 = plot(sort(quantity_pp), x -> F_q(x), label=L"\hat F(q)", xlabel="Quantity")
plot(p10,p11, layout=(1,2), size=(600,300))

U = cdf(F_p, price_pp)
V = F_q(quantity_pp)

# Fit copula
C = BetaCopula(U, V)
# C = BetaCopula(F_p(price_pp), F_q(quantity_pp))
@time C = DiscretizedCopula{:PDF}(C, 300)

gradient(C, [0.5,0.5])

# Plot
Z = [pdf(C,[u,v]) for u in LinRange(0.01,0.99,101), v in LinRange(0.01,0.99,101)]
heatmap(Z)
# xticks = convert.(Int, round.([LinRange(price_pp[begin],price_pp[end],8)...]))
# yticks = convert.(Int, round.([LinRange(minimum(quantity_pp),maximum(quantity_pp),8)...]))
# plot!(xticks=(1:100/7:101, xticks), yticks=(1:100/7:101, yticks))
plot!(xticks=(1:20:101, 0:0.2:1), yticks=(1:20:101, 0:0.2:1))
plot!(xlabel="Price", ylabel="Quantity")


W = rand(C, n)' 
Û = W[:,1] 
V̂ = W[:,2] 

X̂ = F_p⁻¹.(Û) 
Ŷ = F_q⁻¹.(V̂) 

scatter(𝐛, label = "True supply bids") 
xlabel!("Price") 
ylabel!("Quantity") 
scatter!(X̂,Ŷ, label = "Simulated supply bids")

Supply₀ = DataFrame(:Price => X̂, :Quantity => Ŷ, :Curve => :Supply) |> Curve

plot(Supply, color = 1, label = "True supply curve")
plot!(Supply₀, color = 2, label = "Simulated supply curve")
xlabel!("Quantity")
ylabel!("Price")

# Plot 10 realizations
plot(xlabel="Quantity", ylabel="Price")
for i in 2:11
    W = rand(C, n)' 
    Û = W[:,1] 
    V̂ = W[:,2] 
    X̂ = F_p⁻¹.(Û) 
    Ŷ = F_q⁻¹.(V̂) 
    Supply₀ = DataFrame(:Price => X̂, :Quantity => Ŷ, :Curve => :Supply) |> Curve
    plot!(Supply₀, color = i, alpha = 0.3, label = "Simulated supply curve")
end
plot!(Supply, color = 1, label = "True supply curve")


# Calculate error between true curve and simulated curves
DF = transform(bid_df, :Quantity => cumsum => :Quantity)
error = Float64[]
Supply_Curves = Curve[]
for _ in 1:1000
    W = rand(C, n)' 
    Û = W[:,1] 
    V̂ = W[:,2] 
    X̂ = F_p⁻¹.(Û) 
    Ŷ = F_q⁻¹.(V̂) 

    Supply₀ = DataFrame(:Price => X̂, :Quantity => Ŷ, :Curve => :Supply) |> Curve
    push!(Supply_Curves, Supply₀)

    DF_sim = @chain DataFrame(:Price => X̂, :Quantity => Ŷ) begin
        sort!(:Price)
        transform!(:Quantity => cumsum => :Quantity)
    end

    e = sqrt.((DF_sim.Price - DF.Price).^2 + (DF_sim.Quantity - DF.Quantity).^2) |> mean
    push!(error, e)
end
error

histogram(error)

p,q = (80., 100.)
D = Curve(get_data(dates[1], hours, "Buy"))
S = Curve(get_data(dates[1], hours, "Sell"))

[𝐏(S ⊕ˢ (p,q),D)] - [𝐏(S,D)]

[𝐏(Supply_Curves[1] ⊕ˢ (p,q), D)]
[𝐏(Supply_Curves[1] ⊕ˢ (p,q), D)] - [𝐏(Supply,D)]

𝐐(Supply, D)
Supply.f.q[end].b
D.f.q[end].b
