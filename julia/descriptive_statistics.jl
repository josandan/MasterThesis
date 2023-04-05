using Pkg
Pkg.activate(".")

# Packages and functions: 
using Distributions, Embeddings, StatsBase, Polynomials

include("SampleScript\\Includes.jl")

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

# Investigating the bids with very high quantities

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
scatter(dates, sum_pr_day.Quantity_sum, label = "Data")
ylabel!("Quantity")
xlabel!("Date")
p = Polynomials.fit([1.:length(sum_pr_day.Quantity_sum)...], sum_pr_day.Quantity_sum, 2)
plot!(dates, p.([1.:length(sum_pr_day.Quantity_sum)...]), lw = 2, ls = :dash, label = "2nd order polynomial")
p = Polynomials.fit([1.:length(sum_pr_day.Quantity_sum)...], sum_pr_day.Quantity_sum, 3)
plot!(dates, p.([1.:length(sum_pr_day.Quantity_sum)...]), lw = 2, ls = :dash, label = "3rd order polynomial")

data.Quantity |> maximum

# Cov matrix

α = 0.95

Σ = [var(df[:,:Price]) cov(df[:,:Price], df[:,:Quantity]);
     cov(df[:,:Price], df[:,:Quantity]) var(df[:,:Quantity])]

μ = [mean(df[:,:Price]), mean(df[:,:Quantity])]

scatter(df[:,:Price], df[:,:Quantity])
covellipse!(μ, 10 * Σ)

# Investigating the number of bids

n_bids = Float64[]
@time for date = dates
    for hour in hours
        df = get_data(date, hour, side)
        n = filter(:Price => !=(-500.0), df) |> nrow
        append!(n_bids, n)
    end
end

scatter(dates, n_bids, label = "")
xlabel!("Date")
ylabel!("Number of bids")

histogram(n_bids, bins = 175:25:450, label = "")
xlabel!("Number of bids")
ylabel!("Frequency")

# Investigating the total supplied quantity // market volume

total_quantity = Float64[]
@time for date = dates
    for hour in hours
        df = get_data(date, hour, side)
        tq = filter(:Price => !=(-500.0), df).Quantity |> sum
        append!(total_quantity, tq)
    end
end

scatter(dates, total_quantity, label = "")
xlabel!("Date")
ylabel!("Volume")

histogram(total_quantity, label = "")
xlabel!("Volume")
ylabel!("Frequency")
