using Pkg
Pkg.activate(".")

# Packages and functions: 
using Distributions, Embeddings, StatsBase, Polynomials
using EmpiricalCopulas, Chain, DataFramesMeta, ForwardDiff, Interpolations, BivariateCopulas
include("SampleScript\\Includes.jl")
include("functions.jl")

# Get first date hour 12 
date = Date(2022,1,1)
hour = 12

# Show point pattern
ylab = "Quantity"
xlab = "Price"

raw_df_s = get_data(date, hour, "Sell")
raw_df_d = get_data(date, hour, "Buy")

unique_df_s = get_unique_df(raw_df_s)
unique_df_d = get_unique_df(raw_df_d)

df_s = @chain unique_df_s begin
    @subset @byrow :Price != -500.0
end
df_d = @chain unique_df_d begin
    @subset @byrow :Quantity <= 5000.0
end

# p1 = scatter(raw_df_s.Price, raw_df_s.Quantity, label = string(lab, " raw"))
# p2 = scatter(unique_df_s.Price, unique_df_s.Quantity, label = string(lab, " unique"))
# p3 = scatter(df_s.Price, df_s.Quantity, xlabel = xlab, label = string(lab, " cleaned"))
# plot1 = plot(p1,p2,p3, layout=(3,1), size=(600,600), ylabel = ylab, marker=(:circle,4))

p4 = scatter(unique_df_s.Price, unique_df_s.Quantity, label="Supply", color=2)
p5 = scatter(df_s.Price, df_s.Quantity, label="Supply", xlabel=xlab, color=2)
plot2 = plot(p4,p5, layout=(2,1), size=(600,600), ylabel = ylab, marker=(:circle,4))

p4_d = scatter(unique_df_d.Price, unique_df_d.Quantity, label="Demand", color=1)
p5_d = scatter(df_d.Price, df_d.Quantity, label="Demand", xlabel=xlab, color=1)
plot3 = plot(p4,p4_d,p5,p5_d, layout=(2,2), size=(800,600), ylabel=ylab, xlabel=xlab, markerstrokewidth=0.5)

# plot4 = plot(xlabel = xlab, ylabel = ylab)
# scatter!(raw_df.Price, raw_df.Quantity, label = string(lab, " raw"), marker = (:circle, 3), alpha = 1)
# scatter!(unique_df.Price, unique_df.Quantity, label = string(lab, " unique"), marker = (:circle, 3), alpha = 1)
# scatter!(df.Price, df.Quantity, label = string(lab, " cleaned"), marker = (:circle, 3), alpha = 0.5)
# ylims!(-25,800)

# savefig(plot3, "Figures/supply_and_demand_bids.pdf")

# Show point pattern price

df_s.Price
scatter(raw_df_s.Price, marker = (:circle, 3), alpha = 0.5, label = "Raw prices ordered")
scatter!(df_s.Price, marker = (:circle, 3), alpha = 0.5, label = "Prices ordered")
plot5 = scatter(
    df_s.Price, zero, 
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
# savefig(plot5, "Figures/prices_pp.pdf")

# Show point pattern quantity
df_s.Quantity
scatter(df_s.Quantity, marker = (:circle, 2), alpha = 0.5, label = "Quantities ordered by price")
scatter(
    df_s.Quantity, zero, 
    alpha = 0.5, label = "Quantity, supply bids.",
    framestyle = :origin, 
    yaxis = false, ylims = (-1,1), 
    markerstrokewidth = 3, 
    marker = (:vline, 6),
    size = (600,200)
)



# ===========================================
# All dates and hour 12

from_date = Date(2022,1,1)
to_date = Date(2022,10,23)
dates = from_date:Day(1):to_date
hours = 12

# ===========================================
# Supply

@time data_s = CombineDF(dates, hours, "Sell")

DF_s = @chain data_s begin
    @subset @byrow :Quantity != 0
    groupby([:Date, :Hour, :Price, :Curve])
    combine(:Quantity => sum => :Quantity)
end

# Investigating the bids with very high quantities

min_price_df = @chain DF_s begin
    @subset @byrow :Price == -500.0
end

df = @chain DF_s begin
    @subset @byrow :Price != -500.0
end

scatter(min_price_df.Quantity, min_price_df.Price, alpha = 0.6, markerstrokewidth = 0)

# Histogram: quantities, hour 12
h1 = histogram(min_price_df.Quantity, bins=50, color=2, label = "Supply bids (q > 5000)")
xticks!(13000:9000:40000)
ylabel!("Frequency")
xlabel!("Quantity")
h2 = histogram(df.Quantity, bins=200, color=2, label = "Supply bids (q ≤ 5000)")
ylabel!("Frequency")
xlabel!("Quantity")
plot6 = plot(h1,h2, layout=(2,1), size=(600,500))

# savefig(plot6, "Figures/hist_quantity_supply.pdf")

min_price_df.Quantity |> extrema
df.Quantity |> extrema

# Plot: Total quantity pr day with price -500€, hour 12
plot7 = scatter(dates, min_price_df.Quantity, label = "Supply bids (q > 5000)", color=2)
ylabel!("Quantity")
xlabel!("Date")
p = Polynomials.fit([1.:length(min_price_df.Quantity)...], min_price_df.Quantity, 2)
plot!(dates, p.([1.:length(min_price_df.Quantity)...]), lw=2, ls=:dash, color=3, label = "2nd order polynomial")
p = Polynomials.fit([1.:length(min_price_df.Quantity)...], min_price_df.Quantity, 3)
plot!(dates, p.([1.:length(min_price_df.Quantity)...]), lw=2, ls=:dash, color=4, label = "3rd order polynomial")

# savefig(plot7, "Figures/high_quantity_supply.pdf")

# Investigating the number of bids

n_bids = @chain df begin
    groupby([:Date, :Hour, :Curve])
    combine(nrow => :n)
end
extrema(n_bids.n)

p1 = scatter(dates, n_bids.n, label = "Supply", color=2, markerstrokewidth=0.5)
xlabel!("Date")
ylabel!("Number of bids")
p2 = histogram(n_bids.n, bins = 120:10:410, label = "Supply", color=2)
xlabel!("Number of bids")
ylabel!("Frequency")
plot8 = plot(p1,p2, layout=(2,1), size=(600,600))

# savefig(plot8, "Figures/number_of_supply_bids.pdf")

# Investigating the total supplied quantity // market volume

market_quantity = @chain df begin
    groupby([:Date, :Hour, :Curve])
    combine(:Quantity => sum => :MarketQuantity)
end
extrema(market_quantity.MarketQuantity)

p3 = scatter(dates, market_quantity.MarketQuantity, label="Supply",color=2, markerstrokewidth=0.5)
xlabel!("Date")
ylabel!("Market quantity")
p4 = histogram(market_quantity.MarketQuantity, bins = 5000:1000:35000, label="Supply",color=2)
xticks!(5000:10000:35000)
xlabel!("Market quantity")
ylabel!("Frequency")
plot9 = plot(p3,p4, layout=(2,1), size=(600,600))

# savefig(plot9, "Figures/market_supply_quantity.pdf")

# Investigating max prices

max_prices = @chain df begin
    groupby([:Date, :Hour, :Curve])
    transform(:Price => maximum => :MaxPrice)
    @subset @byrow :Price == :MaxPrice
end
extrema(max_prices.Price)

scatter(dates, max_prices.Price)

# ===============================================
# Demand

@time data_d = CombineDF(dates, hours, "Buy")

DF_d = @chain data_d begin
    @subset @byrow :Quantity != 0
    groupby([:Date, :Hour, :Price, :Curve])
    combine(:Quantity => sum => :Quantity)
end

# Investigating the bids with very high quantities

high_q_df = @chain DF_d begin
    @subset @byrow :Quantity > 5000.0
end

# df_d = @chain DF_d begin
#     @subset @byrow :Price != 3000.0
#     @subset @byrow :Price != 4000.0
# end
df_d = @chain DF_d begin
    @subset @byrow :Quantity <= 5000.0
end

scatter(high_q_df.Quantity, high_q_df.Price, markerstrokewidth = 0.5)
scatter(dates, high_q_df.Price, markerstrokewidth = 0.5)

# Histogram: quantities, hour 12
h3 = histogram(high_q_df.Quantity, bins=50, color=1, label = "Demand bids (q > 5000)")
xticks!(14000:7000:40000)
ylabel!("Frequency")
xlabel!("Quantity")
h4 = histogram(df_d.Quantity, bins=200, color=1, label = "Demand bids (q ≤ 5000)")
ylabel!("Frequency")
xlabel!("Quantity")
plot10 = plot(h3,h4, layout=(2,1), size=(600,600))

savefig(plot10, "Figures/hist_quantity_demand.pdf")

high_q_df.Quantity |> extrema
df_d.Quantity |> extrema

# Plot: Total quantity pr day with price -500€, hour 12
plot11 = scatter(dates, high_q_df.Quantity, label = "Demand bids (q > 5000)")
ylabel!("Quantity")
xlabel!("Date")
p = Polynomials.fit([1.:length(high_q_df.Quantity)...], high_q_df.Quantity, 2)
plot!(dates, p.([1.:length(high_q_df.Quantity)...]), lw=2, ls=:dash, color=3, label = "2nd order polynomial")
p = Polynomials.fit([1.:length(high_q_df.Quantity)...], high_q_df.Quantity, 3)
plot!(dates, p.([1.:length(high_q_df.Quantity)...]), lw=2, ls=:dash, color=4, label = "3rd order polynomial")

savefig(plot11, "Figures/high_quantity_demand.pdf")

# Investigating the number of bids

n_bids_d = @chain df_d begin
    groupby([:Date, :Hour, :Curve])
    combine(nrow => :n)
end
extrema(n_bids_d.n)

p5 = scatter(dates, n_bids_d.n, label = "Demand", markerstrokewidth=0.5)
xlabel!("Date")
ylabel!("Number of bids")
p6 = histogram(n_bids_d.n, bins = 190:10:540, label = "Demand")
xlabel!("Number of bids")
ylabel!("Frequency")
plot12 = plot(p5,p6, layout=(2,1), size=(600,600))

savefig(plot12, "Figures/number_of_demand_bids.pdf")

# Investigating the total supplied quantity // market volume

market_quantity_d = @chain df_d begin
    groupby([:Date, :Hour, :Curve])
    combine(:Quantity => sum => :MarketQuantity)
end

p7 = scatter(dates, market_quantity_d.MarketQuantity, label = "Demand", markerstrokewidth=0.5)
xlabel!("Date")
ylabel!("Market quantity")
p8 = histogram(market_quantity_d.MarketQuantity, bins = 12000:1000:31000, label = "Demand")
xlabel!("Market quantity")
ylabel!("Frequency")
plot13 = plot(p7,p8, layout=(2,1), size=(600,600))

savefig(plot13, "Figures/market_demand_quantity.pdf")



# ====================================




plot20 = plot(p1,p5,p2,p6, layout=(2,2), size=(800,600))
extrema(n_bids.n)
extrema(n_bids_d.n)

plot21 = plot(p3,p7,p4,p8, layout=(2,2), size=(800,600))
extrema(market_quantity.MarketQuantity)
extrema(market_quantity_d.MarketQuantity)

plot22 = plot(h1,h3,h2,h4, layout=(2,2), size=(800,600))

savefig(plot20, "Figures/number_of_bids.pdf")
savefig(plot21, "Figures/market_quantity.pdf")
savefig(plot22, "Figures/hist_quantity.pdf")
