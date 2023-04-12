using Pkg
Pkg.activate(".")

# Packages and functions: 
using Distributions, Embeddings, StatsBase

include("SampleScript\\Includes.jl")
include("Mollifiers.jl")
include("functions.jl")


# Get Random Date 
date = rand(Date(2022,1,1):Day(1):Date(2022,11,23))
hour = rand(1:24)

# Get the curves: 
S = Curve(get_data(date, hour, "Sell"))
D = Curve(get_data(date, hour, "Buy"))

# Show point pattern
df = get_data(date, hour, "Sell")

df.Quantity

scatter(df.Quantity, marker = (:circle, 2), alpha = 0.5, label = "Quantities ordered by price")
# scatter(
#     cumsum(df.Quantity), zero(df.Quantity), 
#     alpha = 0.5, label = "Acc. Quantity",
#     framestyle = :origin, 
#     yaxis = false, ylims = (-1,1), 
#     markerstrokewidth = 3, 
#     marker = (:vline, 6),
#     size = (600,200)
# )
scatter(
    df.Quantity, zero(df.Quantity), 
    alpha = 0.5, label = "Quantity, supply bids.",
    framestyle = :origin, 
    yaxis = false, ylims = (-1,1), 
    markerstrokewidth = 3, 
    marker = (:vline, 6),
    size = (600,200)
)

plot(cumsum(sort(df.Quantity)), 1:length(df.Quantity), lt = :steppost) 
plot(1:length(df.Quantity), sort(df.Quantity), lt = :steppost)

clean_df = filter(:Price => !=(-500.0), df)

scatter(clean_df.Quantity, marker = (:circle, 2), alpha = 0.5, label = "Quantities ordered by price")
scatter(
    clean_df.Quantity, zero(clean_df.Quantity), 
    alpha = 0.5, label = "Quantities",
    framestyle = :origin, 
    yaxis = false, ylims = (-1,1), 
    markerstrokewidth = 3, 
    marker = (:vline, 6),
    size = (600,200)
)

scatter(df.Price, marker = (:circle, 2), alpha = 0.5, label = "Prices ordered")
scatter(
    df.Price, zero(df.Price), 
    alpha = 0.5, label = "Price, supply bids.",
    framestyle = :origin, 
    yaxis = false, ylims = (-1,1), 
    markerstrokewidth = 3, 
    marker = (:vline, 6),
    size = (600,200)
)

# spatial plot
p1 = scatter(df.Price, df.Quantity, label = "Data", ylabel = "Quantity")
p2 = scatter(clean_df.Price, clean_df.Quantity, label = "Cleaned data")
plot(p1,p2, layout=(1,2), size=(800,400), xlabel = "Price", marker=(:circle,3))

# Nonparametric estimation of the intensity

from_date = Date(2022,1,1)
to_date = Date(2022,1,5)
dates = from_date:Day(1):to_date
hours = 12
side = "Sell"
point_process = "Price"

# comb_real = CombineRealizations(from_date, to_date, hours, side, point_process)
# price_pp = comb_real.pp
# CumuIntensity = PWLinearEstimator(comb_real)

comb_prices = CombinePriceDF(dates, hours, side)
price_pp = comb_prices.DF.Price
CumuIntensity = PWLinearEstimator((pp = price_pp, k = comb_prices.k, n = comb_prices.n))

plot(price_pp, CumuIntensity, label = "Λ(t)")

scatter!(price_pp, zero(price_pp))
histogram(CumuIntensity.(price_pp))

EstIntensity = DifferentiatePWLinear(price_pp, CumuIntensity)
plot(price_pp, EstIntensity, lt = :steppost, label = "λ(t)")
xlims!(-5,100)

more_unique = unique(round.(price_pp, digits = 0))

EstIntensity2 = DifferentiatePWLinear(more_unique, CumuIntensity)
plot(more_unique, EstIntensity2, lt = :steppost)
xlims!(-5,100)

# Test model

s = CumuIntensity.(price_pp)
s = CumuIntensity.(more_unique)
Δs = [s[1], diff(s)...]
histogram(Δs, normalize = true, label = "Δsₖ")
plot!(t -> (t ≥ 0)*exp(-t), label = "PDF of Standard Exponential")

qqplot(Δs, Exponential(1))
xlabel!("Data")
ylabel!("Exponential Distribution")


# Resample and validate

tₙ = price_pp[price_pp .!== -0.0]
# tₙ = more_unique[more_unique .!== -0.0]

result = OgataThinning(tₙ, EstIntensity)
result.hom_sample
result.sim

plot(tₙ, EstIntensity, lt = :steppost, label = "λ(t)")
scatter!(result.hom_sample, zero, alpha = 0.3, marker = (:circle,2), markerstrokewidth = 0, label = "Hom. PP")
scatter!(result.sim, zero, marker = (:circle,2), label = "Thinned PP")

scatter(tₙ, zero, alpha = 0.5, markerstrokewidth = 0, framestyle = :origin, label = "Data")
ylims!(-1,1)
scatter!(result.sim, zero, alpha = 0.5, markerstrokewidth = 0, label = "Simulation")

bins = tₙ[begin]:((tₙ[end]-tₙ[begin])/100):tₙ[end]
histogram(tₙ, bins = bins, alpha = 0.7, label = "Data")
histogram!(result.sim, bins = bins, alpha = 0.7, label = "Simulation")

length(result.sim)
length(tₙ)/(length(dates)*length(hours))






# Mollified intensity

ϕᵋ = Mollifier(10.)
F̂ = ϕᵋ(CumuIntensity, price_pp)
f̂ = ∂(ϕᵋ)(CumuIntensity, price_pp)

plot(price_pp, CumuIntensity)
plot!(F̂)
plot!(f̂)

plot(price_pp, f̂)
plot!(more_unique, f̂)

plot(more_unique, f̂, label = "Mollified intensity")
# plot!(price_pp, EstIntensity, alpha = 0.5)
plot!(more_unique, EstIntensity2, alpha = 0.5, label = "Intensity")
xlims!(-5,500)

# Test mollified intensity

s = F̂.(price_pp)
Δs = [s[1], diff(s)...]
histogram(Δs, normalize = true, label = "Δsₖ")
plot!(t -> (t ≥ 0)*exp(-t), label = "PDF of Standard Exponential")

qqplot(Δs, Exponential(1))
xlabel!("Data")
ylabel!("Exponential Distribution")






# Resample and validate mollified intensity

result_moll = OgataThinning(tₙ, f̂)

plot(tₙ, f̂, lt = :steppost, label = "λ(t)")
scatter!(result_moll.hom_sample, zero, alpha = 0.3, marker = (:circle,2), markerstrokewidth = 0, label = "Hom. PP")
scatter!(result_moll.sim, zero, marker = (:circle,2), label = "Thinned PP")

scatter(tₙ, zero, alpha = 0.5, markerstrokewidth = 0, framestyle = :origin, label = "Data")
ylims!(-1,1)
scatter!(result_moll.sim, zero, alpha = 0.5, markerstrokewidth = 0, label = "Simulation")

bins = tₙ[begin]:((tₙ[end]-tₙ[begin])/100):tₙ[end]
histogram(tₙ, bins = bins, alpha = 0.7, label = "Data")
histogram!(result_moll.sim, bins = bins, alpha = 0.7, label = "Simulation")

length(result_moll.sim)
length(tₙ)/(length(dates)*length(hours))






# Transform from intensity to density function

a = rand(Uniform(-1,2), 20)
x = Float64[1:20...]
l = LinearEmbedding(a,x)
plot(x, l)
l = MakePos(l)
l = RestrictFun(l)
x2 = [-5:0.1:25...]
plot!(x2, l)


a, x = l.f, l.x


∫(F::LinearEmbedding; l = 0., u = Inf, ε = 0.0) = begin
    f, x = F.f, F.x 

    a, b = extrema(x)

    x̂ = filter(x -> max(l,a) ≤ x ≤ min(u,b + ε), x)  
    f̂ = map(x -> F(x), x̂)

    I = zero(f̂)
    for i in 2:(length(x̂))
        I[i] = I[i-1] + 0.5*(f̂[i] + f̂[i-1])*(x̂[i] - x̂[i-1])
    end

    return Embeddings.LinearEmbedding(I, x̂)
end

@time ∫(f̂; l = 0, u = 10000.)

f̂
MakePos(f̂)

b = ∫(f̂; ε = 50.)(price_pp[end])

plot(price_pp, f̂.(price_pp)/b)










# Combine all dates in one DataFrame

data = DataFrame([[],[],[],[],[],[]], ["Price", "Quantity", "Side", "Date", "Hour", "Curve"])
for date = Date(2022,1,1):Day(1):Date(2022,01,30)
    for hour in 1:24
        append!(data, get_data(date, hour, "Sell"))
        append!(data, get_data(date, hour, "Buy"))
    end
end

data

test = transform(groupby(data, [:Date, :Hour, :Side]), :Quantity .=> cumsum)

# grouped_data = groupby(data, [:Date, :Hour, :Side])

test |>
    filter(:Side => ==("Sell"))

@df test plot(
    :Price,
    :Quantity_cumsum,
    group = (:Date, :Hour),
    legend = false,
    marker = (:circle,2), 
    markerstrokewidth = 0,
    alpha = 0.1
)

plot(data[!,"Price"], data[!,"Quantity"], marker = (:circle,2))

all_dates = Date(2022,1,1):Day(1):Date(2022,11,23)

plot_all_s_curves = function ()
    plot()
    for i in 1:24
        plot!(
            [0], [0], 
            color = i, 
            alpha = 0.7, 
            label = string("Hour ", i), 
            marker = (:circle,2), 
            markerstrokewidth = 0
        )
    end
    xlabel!("Quantity")
    ylabel!("Price")
    for i in 1:length(all_dates)
        for j in 1:24
            S = Curve(get_data(all_dates[i], j, "Sell"))
            plot!(
                S, 
                color = j, 
                alpha = 0.1, 
                marker = (:circle,1),
                markerstrokewidth = 0, 
                label = false
            )
        end
    end
    current()
end

plot_test = plot_all_s_curves()





# Plot the curves 
plot_1 = plot(S, color = 1, label = "Supply", legend = :outerright, marker = (:circle,3))
plot!(D, color = 2, label = "Demand", marker = (:circle,3))
scatter!([𝐐(S,D)],[𝐏(S,D)], label = "Clearing", color = :red)
xlabel!("Quantity")
ylabel!("Price")

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

# Plot the impact
𝐩 = LinRange(-1000, 𝐏(S,D)*1.3, 101)
𝐪 = LinRange(0, 10000, 101)
Z = [impact(S,D, (p,q)) for p ∈ 𝐩, q ∈ 𝐪]
plot_3 = heatmap(𝐪,𝐩, Z)
ylabel!("Price of Bid")
xlabel!("Quantity of Bid")
title!("Impact on Day-Ahead Price (Before EUPHEMIA)")
hline!([𝐏(S,D)], label = "Unimpacted Price", ls = :dash, color = :black)

