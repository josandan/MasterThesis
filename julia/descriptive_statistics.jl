using Pkg
Pkg.activate(".")

# Packages and functions: 
using Distributions, Embeddings, StatsBase

include("c:\\Code\\julia_DayAheadEmbedding\\SampleScript\\Includes.jl")

# Get Random Date 
date = rand(Date(2022,1,1):Day(1):Date(2022,11,23))
hour = rand(1:24)

df = get_data(date, hour, "Sell")

# All dates and hours

from_date = Date(2022,1,1)
to_date = Date(2022,10,23)
dates = from_date:Day(1):to_date
hours = 12
side = "Sell"


data = DataFrame([[],[],[],[],[],[]], ["Price", "Quantity", "Side", "Date", "Hour", "Curve"])
@time for date = dates
    for hour in hours
        df = get_data(date, hour, side)
        max_quantities = filter(:Price => ==(-500.0), df)
        append!(data, max_quantities)
    end
end

filter(:Price => ==(-500.0), df)
filter(:Quantity => >=(5000), df)

scatter(data.Quantity, data.Price, alpha = 0.6, markerstrokewidth = 0)

# Histogram: Bids priced at -500€ for hour 12
histogram(data.Quantity, bins = 50, label = "")
ylabel!("Frequency")
xlabel!("Quantity")

# Histogram: Total quantity pr day with price -500€, hour 12
sum_pr_day = combine(groupby(data, :Date), :Quantity => sum)
histogram(sum_pr_day.Quantity_sum, bins = 50, label = "")
ylabel!("Frequency")
xlabel!("Quantity")

# Plot: Total quantity pr day with price -500€, hour 12
scatter(Date.(sum_pr_day.Date, dateformat"y-m-d"), sum_pr_day.Quantity_sum, label = "")
ylabel!("Quantity")
xlabel!("Date")



data
data.Quantity |> maximum
occurrences = rle(data.Hour)
occurrences[2] |> maximum



α = 0.95

Σ = [var(df[:,:Price]) cov(df[:,:Price], df[:,:Quantity]);
     cov(df[:,:Price], df[:,:Quantity]) var(df[:,:Quantity])]

μ = [mean(df[:,:Price]), mean(df[:,:Quantity])]


scatter(df[:,:Price], df[:,:Quantity])
covellipse!(μ, 10 * Σ)