using Pkg
Pkg.activate(".")

# Packages and functions: 
using Distributions, Embeddings, StatsBase, LaTeXStrings
using EmpiricalCopulas, Chain, DataFramesMeta, ForwardDiff, Interpolations, BivariateCopulas

include("SampleScript\\Includes.jl")
include("Mollifiers.jl")
include("functions.jl")
include("DiscretizedDistributions.jl")


# Get Random Date 
date = Date(2022,1,1)
hour = 12
# date = rand(Date(2022,1,1):Day(1):Date(2022,11,23))
# hour = rand(1:24)

# Get the curves: 
S = Curve(get_data(date, hour, "Sell"))
D = Curve(get_data(date, hour, "Buy"))

# Show point pattern
raw_df = get_data(date, hour, "Sell")

lab = "Supply bids"
ylab = "Quantity"
xlab = "Price"

plot1 = scatter(raw_df.Price, raw_df.Quantity, label = lab, ylabel = ylab, xlabel = xlab)

unique_df = @chain raw_df begin
    @subset @byrow :Quantity != 0
    groupby([:Price, :Curve])
    combine(:Quantity => sum => :Quantity)
end

plot2 = scatter(unique_df.Price, unique_df.Quantity, label = lab, ylabel = ylab, xlabel = xlab)

df = @chain unique_df begin
    @subset @byrow :Price != -500.0
end

plot3 = scatter(df.Price, df.Quantity, label = lab, ylabel = ylab, xlabel = xlab)

# spatial plot
p1 = scatter(raw_df.Price, raw_df.Quantity, label = string(lab, " raw"))
p2 = scatter(unique_df.Price, unique_df.Quantity, label = string(lab, " unique"))
p3 = scatter(df.Price, df.Quantity, xlabel = xlab, label = string(lab, " cleaned"))
plot4 = plot(p1,p2,p3, layout=(3,1), size=(600,600), ylabel = ylab, marker=(:circle,4))

p4 = scatter(unique_df.Price, unique_df.Quantity, label="")
p5 = scatter(df.Price, df.Quantity, label="", xlabel=xlab)
plot5 = plot(p4,p5, layout=(2,1), size=(600,600), ylabel = ylab, marker=(:circle,4))

plot6 = plot(xlabel = xlab, ylabel = ylab)
scatter!(raw_df.Price, raw_df.Quantity, label = string(lab, " raw"), marker = (:circle, 3), alpha = 1)
scatter!(unique_df.Price, unique_df.Quantity, label = string(lab, " unique"), marker = (:circle, 3), alpha = 1)
scatter!(df.Price, df.Quantity, label = string(lab, " cleaned"), marker = (:circle, 3), alpha = 0.5)
ylims!(-25,800)

savefig(plot5, "Figures/supply_bids.pdf")

# Show point pattern price

df.Price
scatter(raw_df.Price, marker = (:circle, 3), alpha = 0.5, label = "Raw prices ordered")
scatter!(df.Price, marker = (:circle, 3), alpha = 0.5, label = "Prices ordered")
plot6 = scatter(df.Price, marker = (:circle, 3), alpha = 0.5, label = "Prices ordered")
plot7 = scatter(
    df.Price, zero, 
    alpha = 0.5, label = "Prices",
    framestyle = :zerolines,
    yaxis = false, ylims = (-0.8,1.2), 
    markerstrokewidth = 10, 
    marker = (:vline, 10),
    size = (600,180),
    bottom_margin=5mm,
    top_margin=2mm,
    left_margin=-2mm,
    right_margin=3mm
)
savefig(plot7, "Figures/prices_pp.pdf")

# Show point pattern quantity

df.Quantity
scatter(df.Quantity, marker = (:circle, 2), alpha = 0.5, label = "Quantities ordered by price")
scatter(
    df.Quantity, zero, 
    alpha = 0.5, label = "Quantity, supply bids.",
    framestyle = :origin, 
    yaxis = false, ylims = (-1,1), 
    markerstrokewidth = 3, 
    marker = (:vline, 6),
    size = (600,200)
)


# Nonparametric estimation of the intensity

from_date = Date(2022,1,1)
to_date = Date(2022,1,1)
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
EstIntensity = DifferentiatePWLinear(price_pp, CumuIntensity)

p6 = plot(price_pp, CumuIntensity, label=L"\hat\Lambda(p)")
scatter!(price_pp, zero, color=1)

plot(price_pp, EstIntensity, lt = :steppost, label=L"\hat\lambda(p)")
p7 = plot(EstIntensity, xlims=extrema(price_pp).*1.05, label=L"\hat\lambda(p)")
xlabel!("Price")
# scatter!(price_pp, zero), color=1)
# plot!(xlim=(-5,50), ylim=(-1,10))

plot8 = plot(p6,p7, layout=(2,1), size=(600,600), ylabel="Intensity")

# more_unique = unique(round.(price_pp, digits = 0))

# EstIntensity2 = DifferentiatePWLinear(more_unique, CumuIntensity)
# plot(more_unique, EstIntensity2, lt = :steppost)
# xlims!(-5,100)

# Test model

s = CumuIntensity.(price_pp)
# s = CumuIntensity.(more_unique)
Δs = [s[1], diff(s)...]
histogram(Δs, normalize = true, label = "Δsₖ")
plot!(t -> (t ≥ 0)*exp(-t), label = "PDF of Standard Exponential")

qqplot(Δs, Exponential(1))
xlabel!("Data")
ylabel!("Exponential Distribution")


# Resample and validate

tₙ = price_pp[price_pp .!== -0.0]
# tₙ = more_unique[more_unique .!== -0.0]

ogata = OgataThinning(tₙ, EstIntensity)
ogata.hom_sample
ogata.sim

plot(EstIntensity, xlims=extrema(price_pp).*1.05, label=L"\hat\lambda(p)")
scatter!(ogata.hom_sample, zero, alpha = 0.3, marker = (:circle,3), markerstrokewidth = 0, label = "Hom. PP")
scatter!(ogata.sim, zero, marker = (:circle,3), label = "Thinned PP")

scatter(tₙ, zero, alpha = 0.5, markerstrokewidth = 0, framestyle = :origin, label = "Data")
ylims!(-1,1)
scatter!(ogata.sim, zero, alpha = 0.5, markerstrokewidth = 0, label = "Simulation")

bins = tₙ[begin]:((tₙ[end]-tₙ[begin])/100):tₙ[end]
histogram(tₙ, bins = bins, alpha = 0.7, label = "Data")
histogram!(ogata.sim, bins = bins, alpha = 0.7, label = "Simulation")

length(ogata.sim)
length(tₙ)/(length(dates)*length(hours))






# Mollified intensity

ϕᵋ = Mollifier(10.)
MollCumuIntensity = ϕᵋ(CumuIntensity, price_pp)
MollIntensity = ∂(ϕᵋ)(CumuIntensity, price_pp)

plot(price_pp, CumuIntensity)
plot!(MollCumuIntensity)
plot!(MollIntensity)

plot(EstIntensity, xlims=extrema(price_pp).*1.05, label=L"\hat\lambda(p)")
plot!(MollIntensity)
xlims!(-5,500)

# Test mollified intensity

s = MollCumuIntensity.(price_pp)
Δs = [s[1], diff(s)...]
histogram(Δs, normalize = true, label = "Δsₖ")
plot!(t -> (t ≥ 0)*exp(-t), label = "PDF of Standard Exponential")

qqplot(Δs, Exponential(1))
xlabel!("Data")
ylabel!("Exponential Distribution")






# Resample and validate mollified intensity

ogata_moll = OgataThinning(tₙ, MollIntensity)

plot(tₙ, MollIntensity, label=L"\hat\lambda(p)")
scatter!(ogata_moll.hom_sample, zero, alpha = 0.3, marker = (:circle,3), markerstrokewidth = 0, label = "Hom. PP")
scatter!(ogata_moll.sim, zero, marker = (:circle,3), label = "Thinned PP")

scatter(tₙ, zero, alpha = 0.5, markerstrokewidth = 0, framestyle = :origin, label = "Data")
ylims!(-1,1)
scatter!(ogata_moll.sim, zero, alpha = 0.5, markerstrokewidth = 0, label = "Simulation")

bins = tₙ[begin]:((tₙ[end]-tₙ[begin])/100):tₙ[end]
histogram(tₙ, bins = bins, alpha = 0.7, label = "Data")
histogram!(ogata_moll.sim, bins = bins, alpha = 0.7, label = "Simulation")

length(ogata_moll.sim)
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


# ∫(F::LinearEmbedding; l = 0., u = Inf, ε = 0.0) = begin
#     f, x = F.f, F.x 

#     a, b = extrema(x)

#     x̂ = filter(x -> max(l,a) ≤ x ≤ min(u,b + ε), x)  
#     f̂ = map(x -> F(x), x̂)

#     I = zero(f̂)
#     for i in 2:(length(x̂))
#         I[i] = I[i-1] + 0.5*(f̂[i] + f̂[i-1])*(x̂[i] - x̂[i-1])
#     end

#     return Embeddings.LinearEmbedding(I, x̂)
# end

# @time ∫(f̂; l = 0, u = 10000.)

# b = ∫(f̂; ε = 50.)(price_pp[end])
# plot(price_pp, f̂.(price_pp)/b)


MollIntensity
empirical_pdf = MollIntensity |> MakePos |> RestrictFun
plot(price_pp, empirical_pdf)





















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

