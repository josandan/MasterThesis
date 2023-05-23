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

# side = "Buy"
# comb_prices = CombineDemandDF(dates, hours, side)

comb_prices = CombinePriceDF(dates, hours, side)
price_pp = comb_prices.DF.Price
CumuIntensity = PWLinearEstimator((pp = price_pp, k = comb_prices.k, n = comb_prices.n))
EstIntensity = DifferentiatePWLinear(price_pp, CumuIntensity)

p6 = plot(price_pp, CumuIntensity, label=L"\hat\Lambda(p)", color=1)
scatter!(price_pp, zero, color=1, label = L"p_1,\ldots,p_n", marker = (:circle, 4))
# scatter!(
#     df.Price, zero, color=1,
#     alpha = 0.5, label = L"p_1,\ldots,p_n", 
#     markerstrokewidth = 1, 
#     marker = (:vline, 5)
# )

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

p8 = plot(price_pp, CumuIntensity, label=L"\hat\Lambda(p)", legend=:right)
plot!(MollCumuIntensity, label=L"\hat\Lambda(p)*\varphi_\varepsilon")
scatter!(price_pp, zero, color=1, label = L"p_1,\ldots,p_n", marker = (:circle, 3), markerstrokewidth=0.5)
p9 = plot(EstIntensity, xlims=extrema(price_pp).*1.05, label=L"\hat\lambda(p)")
plot!(MollIntensity, label=L"\partial(\hat\lambda(p)*\varphi_\varepsilon)")
xlabel!("Price")
# xlims!(-5,500)

plot9 = plot(p8,p9, layout=(2,1), size=(600,600), ylabel="Intensity")

# savefig(plot9, "Figures/price_intensity.pdf")

# Test mollified intensity

s = MollCumuIntensity.(price_pp)
Δs = [s[1], diff(s)...]
h20 = histogram(Δs, normalize = true, label=L"\Delta \tau")
plot!(t -> (t ≥ 0)*exp(-t), label=L"PDF\,\, of\,\, Exp(1)")
plot!(title=L"Testing\,\, \hat\Lambda(p)*\varphi_\epsilon",title_align=:left)
q20 = qqplot(Δs, Exponential(1), xlabel="Δτ", ylabel="Exp(1)", markerstrokewidth=0.5)
plot20 = plot(h20,q20, layout=(1,2), size=(700,300), bottom_margin=3mm)

# savefig(plot20, "Figures/price_residual_analysis.pdf")



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

