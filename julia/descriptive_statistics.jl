using Pkg
Pkg.activate(".")

# Packages and functions: 
using Distributions, Embeddings, StatsBase, Polynomials

include("SampleScript\\Includes.jl")

# Get Random Date 
date = rand(Date(2022,1,1):Day(1):Date(2022,11,23))
hour = rand(1:24)

df = get_data(date, hour, "Sell")

# All dates and hour 12

from_date = Date(2022,1,1)
to_date = Date(2022,10,23)
dates = from_date:Day(1):to_date
hours = 12
side = "Sell"

@time data = CombineDF(dates, hours, side)

DF = @chain data begin
    @subset @byrow :Quantity != 0
    groupby([:Date, :Hour, :Price, :Curve])
    combine(:Quantity => sum => :Quantity)
end

# Investigating the bids with very high quantities

min_price_df = @chain DF begin
    @subset @byrow :Price == -500.0
end

df = @chain DF begin
    @subset @byrow :Price != -500.0
end

scatter(min_price_df.Quantity, min_price_df.Price, alpha = 0.6, markerstrokewidth = 0)

# Histogram: quantities, hour 12
h1 = histogram(min_price_df.Quantity, bins = 50, label = "Bids priced at -500 euro")
ylabel!("Frequency")
xlabel!("Quantity")
h2 = histogram(df.Quantity, bins = 200, label = "All other bids")
ylabel!("Frequency")
xlabel!("Quantity")
plot1 = plot(h1,h2, layout=(2,1), size=(600,600))

savefig(plot1, "Figures/hist_quantity.pdf")

min_price_df.Quantity |> extrema
df.Quantity |> extrema

# Plot: Total quantity pr day with price -500€, hour 12
plot2 = scatter(dates, min_price_df.Quantity, label = "Data")
ylabel!("Quantity")
xlabel!("Date")
p = Polynomials.fit([1.:length(min_price_df.Quantity)...], min_price_df.Quantity, 2)
plot!(dates, p.([1.:length(min_price_df.Quantity)...]), lw = 2, ls = :dash, label = "2nd order polynomial")
p = Polynomials.fit([1.:length(min_price_df.Quantity)...], min_price_df.Quantity, 3)
plot!(dates, p.([1.:length(min_price_df.Quantity)...]), lw = 2, ls = :dash, label = "3rd order polynomial")

savefig(plot2, "Figures/minus500_quantity.pdf")

# Cov matrix

α = 0.95

Σ = [var(df[:,:Price]) cov(df[:,:Price], df[:,:Quantity]);
     cov(df[:,:Price], df[:,:Quantity]) var(df[:,:Quantity])]

μ = [mean(df[:,:Price]), mean(df[:,:Quantity])]

scatter(df[:,:Price], df[:,:Quantity])
covellipse!(μ, 10 * Σ)

# Investigating the number of bids

n_bids = @chain df begin
    groupby([:Date, :Hour, :Curve])
    combine(nrow => :n)
end
extrema(n_bids.n)

p1 = scatter(dates, n_bids.n, label = "")
xlabel!("Date")
ylabel!("Number of bids")
p2 = histogram(n_bids.n, bins = 120:10:410, label = "")
xlabel!("Number of bids")
ylabel!("Frequency")
plot3 = plot(p1,p2, layout=(2,1), size=(600,600))

savefig(plot3, "Figures/number_of_bids.pdf")

# Investigating the total supplied quantity // market volume

market_quantity = @chain df begin
    groupby([:Date, :Hour, :Curve])
    combine(:Quantity => sum => :MarketQuantity)
end
extrema(market_quantity.MarketQuantity)

p3 = scatter(dates, market_quantity.MarketQuantity, label = "")
xlabel!("Date")
ylabel!("Volume")
p4 = histogram(market_quantity.MarketQuantity, bins = 5000:1000:35000, label = "")
xlabel!("Volume")
ylabel!("Frequency")
plot4 = plot(p3,p4, layout=(2,1), size=(600,600))

savefig(plot4, "Figures/market_quantity.pdf")

