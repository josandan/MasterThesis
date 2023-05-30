using Pkg
Pkg.activate(".")

# Packages and functions: 
include("Includes.jl")
include("Mollifiers.jl")
include("DiscretizedDistributions.jl")
include("BetaCopula.jl")
include("functions.jl")

NBids = DataFrame(
    Day=Int64[],
    True=Int64[], AvgOgata=Float64[], AvgInv=Float64[], MedOgata=Float64[], MedInv=Float64[], 
    LbOgata=Float64[], UbOgata=Float64[], LbInv=Float64[], UbInv=Float64[]
)
TotalQuantity = DataFrame(
    Day=Int64[],
    True=Float64[], AvgOgata=Float64[], AvgInv=Float64[], MedOgata=Float64[], MedInv=Float64[], 
    LbOgata=Float64[], UbOgata=Float64[], LbInv=Float64[], UbInv=Float64[]
)

P = Progress(31)
for i in 1:31
    from_date = Date(2022,1,i)
    to_date = Date(2022,1,i)
    dates = from_date:Day(1):to_date
    hours = 12
    mollifier_tolerance = 10.

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

    Intensities = GetIntensities(price_pp, comb_prices.k, comb_prices.n, mollifier_tolerance)
    λₚ = Intensities.λ
    Λₚ⁻¹ = Intensities.Λ⁻¹

    ogata_length = Int64[]
    inv_length = Int64[]
    ogata_sum_q = Float64[]
    inv_sum_q = Float64[]
    for _ in 1:1000
        ogata_p = OgataThinning(price_pp, λₚ).sim
        inv_p = SimulateByInversion(Λₚ⁻¹, price_pp[end])
        push!(ogata_length, length(ogata_p))
        push!(inv_length, length(inv_p))
        ogata_q = SimulateQuantity(ogata_p, F_p, G_q⁻¹, C)
        inv_q = SimulateQuantity(inv_p, F_p, G_q⁻¹, C)
        push!(ogata_sum_q, sum(ogata_q))
        push!(inv_sum_q, sum(inv_q))
    end

    mbo = mean(ogata_length)
    mbi = mean(inv_length)
    medbo = median(ogata_length)
    medbi = median(inv_length)
    lbbo = quantile(ogata_length, 0.025)
    ubbo = quantile(ogata_length, 0.975)
    lbbi = quantile(inv_length, 0.025)
    ubbi = quantile(inv_length, 0.975)
    trueb = length(price_pp)/(length(dates)*length(hours))

    push!(NBids, (i,trueb,mbo,mbi,medbo,medbi,lbbo,ubbo,lbbi,ubbi))

    mqo = mean(ogata_sum_q)
    mqi = mean(inv_sum_q)
    medqo = median(ogata_sum_q)
    medqi = median(inv_sum_q)
    lbqo = quantile(ogata_sum_q, 0.025)
    ubqo = quantile(ogata_sum_q, 0.975)
    lbqi = quantile(inv_sum_q, 0.025)
    ubqi = quantile(inv_sum_q, 0.975)
    trueq = sum(quantity_pp)

    push!(TotalQuantity, (i,trueq,mqo,mqi,medqo,medqi,lbqo,ubqo,lbqi,ubqi))

    next!(P)

end


NBids
TotalQuantity

p1 = @df NBids begin
    plot(xlabel="Day in January", ylabel="Number of bids", title="Ogata", titleposition=:left)
    plot!(:True, label="True", bottom_margin=-2mm)
    plot!(:AvgOgata, label="Avg.")
    plot!(:MedOgata, label="Median", color=4)
    plot!(:LbOgata, label="Lower bound", color=3)
    plot!(:UbOgata, label="Upper bound", color=3)
end
p2 = @df NBids begin
    plot(xlabel="Day in January", ylabel="Number of bids", title="Inversion", titleposition=:left)
    plot!(:True, label="True", top_margin=-2mm)
    plot!(:AvgInv, label="Avg.")
    plot!(:MedInv, label="Median", color=4)
    plot!(:LbInv, label="Lower bound", color=3)
    plot!(:UbInv, label="Upper bound", color=3)
end

plot_nbids = plot(p1,p2, layout=(2,1), size=(700,700))

p3 = @df TotalQuantity begin
    plot(xlabel="Day in January", ylabel="Total supplied quantity", title="Ogata", titleposition=:left)
    plot!(:True, label="True", bottom_margin=-2mm)
    plot!(:AvgOgata, label="Avg.")
    plot!(:MedOgata, label="Median", color=4)
    plot!(:LbOgata, label="Lower bound", color=3)
    plot!(:UbOgata, label="Upper bound", color=3)
end
p4 = @df TotalQuantity begin
    plot(xlabel="Day in January", ylabel="Total supplied quantity", title="Inversion", titleposition=:left)
    plot!(:True, label="True", top_margin=-2mm)
    plot!(:AvgInv, label="Avg.")
    plot!(:MedInv, label="Median", color=4)
    plot!(:LbInv, label="Lower bound", color=3)
    plot!(:UbInv, label="Upper bound", color=3)
end

plot_totquant = plot(p3,p4, layout=(2,1), size=(700,700))

# savefig(plot_nbids, "Figures/sim_january_nbids.pdf")
# savefig(plot_totquant, "Figures/sim_january_totquantity.pdf")



results

saveres

latexify(round.(results, digits=1), env=:table)

table = DataFrame(
Date = string.(1:31,"-1-2022"),
TrueN = results.TrueNbids,
MeanNOgata = string.(results.AvgNBidsOgata, "a(", results.StdNbidsOgata, ")"),
MeanNInv = string.(results.AvgNBidsInv, "a(", results.StdNbidsInv, ")"),
TrueQ = results.TrueQuantity,
MeanQOgata = string.(results.AvgQuantityOgata, "a(", results.StdQuantityOgata, ")"),
MeanQInv = string.(results.AvgQuantityInv, "a(", results.StdQuantityInv, ")")
)

latexify(table[1:5,:], env=:table)


Off = @chain results begin
    @transform OffNOgata = abs.(:TrueNbids - :AvgNBidsOgata)
    @transform OffNInv = abs.(:TrueNbids - :AvgNBidsInv)
    @transform OffQOgata = abs.(:TrueQuantity - :AvgQuantityOgata)
    @transform OffQInv = abs.(:TrueQuantity - :AvgQuantityInv)
    @select(:OffNOgata, :OffNInv, :OffQOgata, :OffQInv)
end
mean.(eachcol(Off))
