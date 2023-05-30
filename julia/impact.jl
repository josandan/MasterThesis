using Pkg
Pkg.activate(".")

# Packages and functions: 
include("Includes.jl")
include("Mollifiers.jl")
include("DiscretizedDistributions.jl")
include("BetaCopula.jl")
include("functions.jl")

MeanChange = DataFrame(
    Day=Int64[],
    Q100=Float64[], Q200=Float64[], Q300=Float64[], Q400=Float64[], Q500=Float64[], 
    Q600=Float64[], Q700=Float64[], Q800=Float64[], Q900=Float64[], Q1000=Float64[]
)
UBChange = DataFrame(
    Day=Int64[],
    Q100=Float64[], Q200=Float64[], Q300=Float64[], Q400=Float64[], Q500=Float64[], 
    Q600=Float64[], Q700=Float64[], Q800=Float64[], Q900=Float64[], Q1000=Float64[]
)
LBChange = DataFrame(
    Day=Int64[],
    Q100=Float64[], Q200=Float64[], Q300=Float64[], Q400=Float64[], Q500=Float64[], 
    Q600=Float64[], Q700=Float64[], Q800=Float64[], Q900=Float64[], Q1000=Float64[]
)
Executed = DataFrame(
    Day=Int64[],
    Q100=Float64[], Q200=Float64[], Q300=Float64[], Q400=Float64[], Q500=Float64[], 
    Q600=Float64[], Q700=Float64[], Q800=Float64[], Q900=Float64[], Q1000=Float64[]
)

P = Progress(24)
for i in 1:24
# i = 24

    from_date = Date(2022,1,1)
    to_date = Date(2022,1,1)
    dates = from_date:Day(1):to_date
    hours = i
    mollifier_tolerance = 10.
    fixed_bid = comb_prices.fixed_bid
    # S = Curve(push!(copy(bid_df), fixed_bid[1,:]))
    D = GetOtherCurve(dates, hours, "Buy")
    # Supply = Curve(bid_df)
    # 𝐛 = bids(Supply)
    # n = convert(Int, round(comb_prices.n/comb_prices.k))

    comb_prices = CombinePriceDF(dates, hours, "Sell")
    price_pp = comb_prices.DF.Price
    quantity_pp = comb_prices.DF.Quantity

    F_p = GetDensity(price_pp, comb_prices.k, comb_prices.n, mollifier_tolerance)
    G_q = ecdf(quantity_pp)
    # F_p⁻¹ = InverseDensity(cdf(F_p, price_pp), price_pp)
    G_q⁻¹ = x -> quantile(quantity_pp, x)

    U = cdf(F_p, price_pp)
    V = G_q(quantity_pp)
    C = BetaCopula(U, V)
    C = DiscretizedCopula{:PDF}(C, 300)

    Λₚ⁻¹ = GetIntensities(price_pp, comb_prices.k, comb_prices.n, mollifier_tolerance).Λ⁻¹

    simulated_curves = map(1:1000) do _
        # ogata = OgataThinning(price_pp, λₚ).sim
        X̂ = SimulateByInversion(Λₚ⁻¹, price_pp[end])
        Ŷ = SimulateQuantity(X̂, F_p, G_q⁻¹, C)
    
        push!(X̂, fixed_bid.Price[1])
        push!(Ŷ, fixed_bid.Quantity[1])
        DataFrame(:Price => X̂, :Quantity => Ŷ, :Curve => :Supply) |> Curve
    end

    # p_range = -50.:50:150 |> collect
    p = -50.
    q_range = 100.:100:1000 |> collect

    MeanVec = Float64[]
    UBVec = Float64[]
    LBVec = Float64[]
    ExeVec = Float64[]
    for q_index in eachindex(q_range)
        Price_Change = Float64[]
        Bid_executed = Bool[]
        for i = 1:1000
            pc = 𝐏(simulated_curves[i] ⊕ˢ (p,q_range[q_index]),D) - 𝐏(simulated_curves[i],D)
            push!(Price_Change, pc)
            be = p ≤ 𝐏(simulated_curves[i] ⊕ˢ (p,q_range[q_index]),D)
            push!(Bid_executed, be)
        end
        push!(MeanVec, mean(Price_Change))
        push!(UBVec, quantile(Price_Change, 0.975))
        push!(LBVec, quantile(Price_Change, 0.025))
        push!(ExeVec, mean(Bid_executed))
    end

    push!(MeanChange, [i,MeanVec...])
    push!(UBChange, [i,UBVec...])
    push!(LBChange, [i,LBVec...])
    push!(Executed, [i,ExeVec...])

    next!(P)

end

SMeanChange
SUBChange
SLBChange
SExecuted

HMeanChange
HUBChange
HLBChange
HExecuted

MeanChange = SMeanChange
UBChange = SUBChange
LBChange = SLBChange

plot_pc_jan = plot(ylabel="Price change", color_palette = palette(:lightrainbow, 10))
for i in 100:100:1000
    plot!(MeanChange[:,string("Q",i)], label=string(i,"MW"))
end
plot!(bottom_margin=-2mm)
plot_pc_bounds_jan = plot(xlabel="Day in January", ylabel="Price change", color_palette = palette(:lightrainbow, 10))
for i in 1:10
    plot!(MeanChange[:,string("Q",i*100)], label=string(i*100,"MW"), color=i)
    plot!(UBChange[:,string("Q",i*100)], label="", alpha=0.4, color=i, linestyle=:dash)
    plot!(LBChange[:,string("Q",i*100)], label="", alpha=0.4, color=i, linestyle=:dash)
end
plot!(top_margin=-2mm)

# savefig(plot_pc_jan, "Figures/price_change_january.pdf")
# savefig(plot_pc_bounds_jan, "Figures/price_change_bounds_january.pdf")
plot41 = plot(plot_pc_jan, plot_pc_bounds_jan, layout=(2,1), size=(650,650))
# savefig(plot41, "Figures/price_change_jan_2.pdf")


plot_pc_hours = plot(ylabel="Price change", color_palette = palette(:lightrainbow, 10))
for i in 100:100:1000
    plot!(MeanChange[:,string("Q",i)], label=string(i,"MW"))
end
plot!(bottom_margin=-2mm)
plot_pc_bounds_hours = plot(xlabel="Hour", ylabel="Price change", color_palette = palette(:lightrainbow, 10))
for i in 1:10
    plot!(MeanChange[:,string("Q",i*100)], label=string(i*100,"MW"), color=i)
    plot!(UBChange[:,string("Q",i*100)], label="", alpha=0.4, color=i, linestyle=:dash)
    plot!(LBChange[:,string("Q",i*100)], label="", alpha=0.4, color=i, linestyle=:dash)
end
plot!(top_margin=-2mm)

# savefig(plot_pc_hours, "Figures/price_change_hours.pdf")
# savefig(plot_pc_bounds_hours, "Figures/price_change_bounds_hours.pdf")
plot40 = plot(plot_pc_hours, plot_pc_bounds_hours, layout=(2,1), size=(650,650))
# savefig(plot40, "Figures/price_change_hours_2.pdf")
