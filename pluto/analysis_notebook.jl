### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ f4652d00-c988-11ed-3df1-fb3e0012c517
begin 
	using Pkg; Pkg.activate()

	using DayAheadEmbedding
	using Dates
	using SQLite
	using DataFrames
	using Distributions
	using Plots
	using Measures
	using StatsPlots
	using ProgressMeter
	using PlutoUI
	
	import SQLite:
	    DB,
	    load!

	import SQLite.DBInterface:
	    execute
	
	function get_data(date::Date, hour::Int, side::String)
	    @assert 1 ≤ hour ≤ 24
	    @assert side ∈ ["Buy","Sell"]
	
	    db = DB("../TestData/DayAhead.db")
	
	    execute(db, "SELECT * FROM '$(string(date)*"_"*string(hour)*"_"*side)'") |> DataFrame
	end
	
	PlutoUI.TableOfContents()

end

# ╔═╡ ff2bfabb-bda4-4378-ab37-05c4aee38b6f
md"""
# Introduction

We start by choosing a date and hour randomly, and plot the supply and demand curves.

"""

# ╔═╡ 7b832df5-ea0a-42f3-937f-ec5bd85feb25
begin
	# Get Random Date 
	date = rand(Date(2022,1,1):Day(1):Date(2022,11,23))
	hour = rand(1:24)
	
	# Get the curves: 
	S = Curve(get_data(date, hour, "Sell"))
	D = Curve(get_data(date, hour, "Buy"))
end

# ╔═╡ bdfe9195-16aa-4336-9926-875bf96b9a95
begin
	# Plot curves:
	plot(S, color = 1, label = "Supply", legend = :outerright)
	plot!(D, color = 2, label = "Demand")
	scatter!([𝐐(S,D)],[𝐏(S,D)], label = "Clearing", color = :red)
	xlabel!("Quantity")
	ylabel!("Price")
end

# ╔═╡ 2fbaf81c-a6e0-4ffe-96b1-da5c8d6da425
md"""

We know focus on the supply curve and the underlying supply bids. We start by plotting the bids ordered after price.

"""

# ╔═╡ fc90bc23-e5f7-401d-b602-cb8ce5decc56
begin
	df = get_data(date, hour, "Sell")
	df.Quantity
end

# ╔═╡ 90f2222a-7966-4826-a4e4-3b99523c34a8
begin
	scatter(df.Quantity, marker = (:circle, 2), alpha = 0.5, label = "Supply quantities ordered by price") # the supply quantities ordered after price
	xlabel!("Number")
	ylabel!("Quantity")
end

# ╔═╡ c19d2d0c-a778-4952-9077-d7dbb9596192
md"""

On the plot above, we see that there is one bid of very high quantity and a lot of bids of low and close to zero quantities. 

Then we plot the univariate point pattern of the supply bids quantities. 

"""

# ╔═╡ a7e6a6cd-4285-457e-bfa3-ab0d7170b4be
begin
	scatter(
	    df.Quantity, zero(df.Quantity), 
	    alpha = 0.5, label = "Quantity, supply bids.",
	    framestyle = :origin, 
	    yaxis = false, ylims = (-1,1), 
	    markerstrokewidth = 3, 
	    marker = (:vline, 6),
	    size = (600,200)
	)
end

# ╔═╡ 50647689-860a-4697-950e-f6a0e32bb00b
md"""

The plot above shows the same pattern with almost all of the points being close to zero, meaning that most quantities are quite small. To look further into this, we plot a histogram of the points.

"""

# ╔═╡ 2446d50f-9ac1-4e02-b2b7-36302da34133
begin
	bins = 0:(maximum(df.Quantity)/1000):maximum(df.Quantity)
	histogram(df.Quantity, linewidth = 0, label = "Quantity, supply bids.", bins = bins)
end

# ╔═╡ 2752e69f-a508-4174-b61b-44da93567598
let
	histogram(df.Quantity, linewidth = 0, label = "Quantity, supply bids.", xlims = (-10,300), bins = bins)
end

# ╔═╡ 904a8288-57ed-4d21-8713-dddb1c8eac17
md"""
## Descriptive statistics


"""

# ╔═╡ 29047410-d70a-4940-8380-0297c347cb60
md"""
## Cleaning of data

By removing the bids that are priced at -500€, we obtain the following point pattern for the quantities of the supply bids.

"""

# ╔═╡ c2744cc0-ab6e-4b83-8428-0fdc915411f9
begin
	clean_df = filter(:Price => !=(-500.0), df)

	scatter(clean_df.Quantity, marker = (:circle, 2), alpha = 0.5, label = "Supply quantities ordered by price")
	xlabel!("Number")
	ylabel!("Quantity")
end

# ╔═╡ fac631ca-0a6f-4a4f-bae8-bfa436c7ffb0
let
	scatter(
	    clean_df.Quantity, zero(clean_df.Quantity), 
	    alpha = 0.5, label = "Quantity, supply bids.",
	    framestyle = :origin, 
	    yaxis = false, ylims = (-1,1), 
	    markerstrokewidth = 3, 
	    marker = (:vline, 6),
	    size = (600,200)
	)
end

# ╔═╡ c9703c70-002e-465e-9ca6-d037dc6f0128
md"""

# Nonparametric estimation of the intensity

We will apply the procedure of Leemis [1991] that is a nonparametric estimation procedure of the cumulative intensity function. For this estimation it is required that the points are ordered and multiplicity of points should be handled. For now, we only keep one point for each case of tied values to avoid multiplicity of points.

The following function is for combining different realization, meaning that it takes dates, hours and side as inputs and returns a tuple consisting of the combined point process vector (pp), the number of realizations (k) and the total number of points (n).

"""

# ╔═╡ 710f0385-7454-45b0-93d3-c5038f4fa4b1
function CombineRealizations(from_date::Date, to_date::Date, hours::Union{AbstractVector{N}, Int64}, side::String) where {N}
    dates = from_date:Day(1):to_date

    data = Float64[]
    for date = dates
        for hour in hours
			df = get_data(date, hour, side)
            cleaned_df = filter(:Price => !=(-500.0), df)
            append!(data, cleaned_df.Quantity)
        end
    end

    k = length(dates) * length(hours)
    n = length(data)
    pp = unique(round.(data, digits = 4)) |> sort

    return (pp = pp, k = k, n = n)
end

# ╔═╡ 62bda101-84b9-4207-8635-80acabe5a2a9
md"""

We implement the estimation procedure:

"""

# ╔═╡ 6f8ef284-3ecd-4d1e-b6cc-cfd2e50b0278
function PWLinearEstimator(comb_realization::NamedTuple{(:pp, :k, :n), Tuple{Vector{Float64}, Int64, Int64}})
    pp = comb_realization.pp
    k = comb_realization.k
    n = comb_realization.n

    t -> begin
        pp = sort(pp)
        i = searchsortedlast(pp, t)
        if i == 0
            return 0
        elseif i >= length(pp)
            return ((i-1)*n)/((n+1)*k) + (n*(pp[i]-pp[i-1]))/((n+1)*k*(pp[i]-pp[i-1]))
        else
            return (i*n)/((n+1)*k) + (n*(t-pp[i]))/((n+1)*k*(pp[i+1]-pp[i]))
        end
    end
end

# ╔═╡ 3c6d9b5c-e4fc-489d-9934-40baa044f946
md"""

We use the first six hours of the first three days of the data set, and estimate the cumulative intensity function of that data. The estimated cumulative intensity function is plotted below.

"""

# ╔═╡ 96e676d3-a820-49fc-8002-91f67836f178
begin
	from_date = Date(2022,1,1)
	to_date = Date(2022,1,1)
	hours = 12
	side = "Sell"
	
	comb_real = CombineRealizations(from_date, to_date, hours, side)
	CumuIntensity = PWLinearEstimator(comb_real)
	
	plot(comb_real.pp, CumuIntensity, label = "")
	xlabel!("Quantity, q")
	ylabel!("Λ(q)")
end

# ╔═╡ 0cc6c9c5-7462-4671-b626-84651f6834fd
md"""

The cumulative intensity function estimated in this way is a piece-wise linear function. Hence, the intensity function is piece-wise constant and can be computed by the following function.

"""

# ╔═╡ 549e9f5c-79fa-45e4-8020-4b9a97af0413
function DifferentiatePWLinear(pp::Vector{Float64}, PWLinearFun::Function)
    t -> begin
        sort!(pp)
        i = searchsortedlast(pp,t)
        if i == 0
            return 0
        elseif i >= length(pp)
            return 0
        else
            return (PWLinearFun(pp[i+1]) - PWLinearFun(pp[i])) / (pp[i+1] - pp[i])
        end
    end
end

# ╔═╡ 238f1a77-0285-4fa6-b714-504777e48370
md"""

By use of the above function, we can compute the intensity function, which is plotted below.

"""

# ╔═╡ 9a39f749-eb54-41df-b1c7-8d75e40d5f23
begin
	Est_Intensity = DifferentiatePWLinear(comb_real.pp, CumuIntensity)
	plot(comb_real.pp, Est_Intensity, lt = :steppost, label = "")
	xlabel!("Quantity, q")
	ylabel!("λ(q)")
end

# ╔═╡ d2607758-24ed-4b3d-8d0a-0c4968088a1b
let
	plot(comb_real.pp, Est_Intensity, lt = :steppost, label = "")
	xlims!(-5,100)
	xlabel!("Quantity, q")
	ylabel!("λ(q)")
end

# ╔═╡ 35800923-b37b-42f0-aff3-15cd33ff4df5
md"""
The intensity fluctuates a lot due to having many points on a fine grid.

The smallest amount that can be traded on the EPEX day-ahead electricity market in Germany is 1 MWh. Therefore, we can coarse the grid to such that the smallest change in quantity is 1 unit. This will also get rid of some the fluctuations in the intensity function.

"""

# ╔═╡ 9c7d981b-55d1-4592-ad1b-364ed3172f1a
more_unique = unique(round.(comb_real.pp, digits = 0))

# ╔═╡ d710143d-c57e-4258-bcee-7599e547d754
begin
	Est_Intensity2 = DifferentiatePWLinear(more_unique, CumuIntensity)
	plot(more_unique, Est_Intensity2, lt = :steppost, label = "")
	xlabel!("Quantity, q")
	ylabel!("λ(q)")
end

# ╔═╡ fcb36368-c61c-416d-b308-07ff265a576e
let
	plot(more_unique, Est_Intensity2, lt = :steppost, label = "")
	xlims!(-5,200)
	xlabel!("Quantity, q")
	ylabel!("λ(q)")
end

# ╔═╡ 246dafae-4e72-4f31-bb95-f3aa291c1731
md"""

## Testing the model

Let's test the point process model, which we do by testing the estimated intensity function.

"""

# ╔═╡ 7f035c3d-25d1-45dd-bf37-ea84e221e198
# let
# 	s = CumuIntensity.(more_unique)
# 	Δs = [s[1], diff(s)...]
# 	histogram(Δs, normalize = true, label = "Δsₖ")
# 	plot!(t -> (t ≥ 0)*exp(-t), label = "PDF of Standard Exponential")
# end

# ╔═╡ 196e63d5-7e61-47ba-a6b9-7338f1288546
# let
# 	s = CumuIntensity.(more_unique)
# 	Δs = [s[1], diff(s)...]
# 	qqplot(Δs, Exponential(1))
# 	xlabel!("Data")
# 	ylabel!("Exponential Distribution")
# end

# ╔═╡ 48de2203-a9a4-4b07-85bc-cdd618d5a14d
md"""
# Resample and validate

Now we want to resample from the estimated intensity function. We will do that by using Ogata's thinning algorithm. For the Ogata's thinning algorithm we need an algorithm for simulating homogenuous Poisson processes. Both are implemented below.

"""

# ╔═╡ 08083494-9371-4baf-98ad-1c3404359cb7
function SimulateHomPoisson(T::Real, λ::Real)
	t = 0
	σ = Float64[]

	while t < T
		Δσₖ = first(rand(Exponential(1/λ), 1))
		t = t + Δσₖ
		if t ≤ T
			σₖ = t
			push!(σ, σₖ)
		end
	end

	return σ
end

# ╔═╡ 6f81e81a-4747-4b58-9b21-d4ad6d055d49
function OgataThinning(tₙ::AbstractVector{R}, intensity_func::Function) where {R}
	σ = Float64[]
	T = tₙ[end]
	max_intensity = maximum(intensity_func.(tₙ))

	hom_sample = SimulateHomPoisson(T, max_intensity)

	for i in 1:length(hom_sample)
		u = rand(Uniform(0,1))
		if u <= intensity_func(hom_sample[i])/max_intensity
			push!(σ, hom_sample[i])
		end
	end
	return (hom_sample = hom_sample, sim = σ)
end

# ╔═╡ 46188f37-312d-48fa-86f8-3f1999612b29
begin
	tₙ = more_unique
	result = OgataThinning(tₙ, Est_Intensity2)
end

# ╔═╡ f69e16ca-4961-460d-bbb4-5a03e5ec3e19
let
	plot(tₙ, Est_Intensity, lt = :steppost, label = "λ(t)")
	scatter!(result.hom_sample, zero, alpha = 0.2, marker = (:circle,2), markerstrokewidth = 0, label = "Hom. PP")
	scatter!(result.sim, zero, marker = (:circle,2), label = "Thinned PP")
end

# ╔═╡ ce5cad4c-4904-4404-a2b0-4df790f50c26
let
	scatter(tₙ, zero, alpha = 0.5, markerstrokewidth = 0, framestyle = :origin, label = "Data")
	ylims!(-1,1)
	scatter!(result.sim, zero, alpha = 0.5, markerstrokewidth = 0, label = "Simulation")
end

# ╔═╡ 15cd3191-6f96-4016-8521-318f166ca05a
let
	bins = 0:(tₙ[end]/100):tₙ[end]
	histogram(tₙ, bins = bins, alpha = 0.5, label = "Data")
	histogram!(result.sim, bins = bins, alpha = 0.5, label = "Simulation")
end

# ╔═╡ e32e2abb-84af-4137-9e71-7e73266ed16a
length(result.sim)

# ╔═╡ c706b1e8-58b9-4536-8cfe-68f5d8631b28
length(tₙ)

# ╔═╡ 5dd187d8-e7df-462a-b161-adbbd3c7dd19
round(sum(result.sim), digits = 0)

# ╔═╡ 4d4a7cb2-07c5-4d55-960b-9dac0561c5fe
sum(tₙ)

# ╔═╡ a85190fa-6424-49e9-997e-5628aa1fa4a4


# ╔═╡ Cell order:
# ╠═f4652d00-c988-11ed-3df1-fb3e0012c517
# ╟─ff2bfabb-bda4-4378-ab37-05c4aee38b6f
# ╠═7b832df5-ea0a-42f3-937f-ec5bd85feb25
# ╟─bdfe9195-16aa-4336-9926-875bf96b9a95
# ╟─2fbaf81c-a6e0-4ffe-96b1-da5c8d6da425
# ╠═fc90bc23-e5f7-401d-b602-cb8ce5decc56
# ╠═90f2222a-7966-4826-a4e4-3b99523c34a8
# ╟─c19d2d0c-a778-4952-9077-d7dbb9596192
# ╠═a7e6a6cd-4285-457e-bfa3-ab0d7170b4be
# ╟─50647689-860a-4697-950e-f6a0e32bb00b
# ╠═2446d50f-9ac1-4e02-b2b7-36302da34133
# ╠═2752e69f-a508-4174-b61b-44da93567598
# ╠═904a8288-57ed-4d21-8713-dddb1c8eac17
# ╟─29047410-d70a-4940-8380-0297c347cb60
# ╠═c2744cc0-ab6e-4b83-8428-0fdc915411f9
# ╠═fac631ca-0a6f-4a4f-bae8-bfa436c7ffb0
# ╟─c9703c70-002e-465e-9ca6-d037dc6f0128
# ╠═710f0385-7454-45b0-93d3-c5038f4fa4b1
# ╟─62bda101-84b9-4207-8635-80acabe5a2a9
# ╠═6f8ef284-3ecd-4d1e-b6cc-cfd2e50b0278
# ╟─3c6d9b5c-e4fc-489d-9934-40baa044f946
# ╠═96e676d3-a820-49fc-8002-91f67836f178
# ╟─0cc6c9c5-7462-4671-b626-84651f6834fd
# ╠═549e9f5c-79fa-45e4-8020-4b9a97af0413
# ╟─238f1a77-0285-4fa6-b714-504777e48370
# ╠═9a39f749-eb54-41df-b1c7-8d75e40d5f23
# ╠═d2607758-24ed-4b3d-8d0a-0c4968088a1b
# ╟─35800923-b37b-42f0-aff3-15cd33ff4df5
# ╠═9c7d981b-55d1-4592-ad1b-364ed3172f1a
# ╠═d710143d-c57e-4258-bcee-7599e547d754
# ╠═fcb36368-c61c-416d-b308-07ff265a576e
# ╠═246dafae-4e72-4f31-bb95-f3aa291c1731
# ╠═7f035c3d-25d1-45dd-bf37-ea84e221e198
# ╠═196e63d5-7e61-47ba-a6b9-7338f1288546
# ╠═48de2203-a9a4-4b07-85bc-cdd618d5a14d
# ╠═08083494-9371-4baf-98ad-1c3404359cb7
# ╠═6f81e81a-4747-4b58-9b21-d4ad6d055d49
# ╠═46188f37-312d-48fa-86f8-3f1999612b29
# ╠═f69e16ca-4961-460d-bbb4-5a03e5ec3e19
# ╠═ce5cad4c-4904-4404-a2b0-4df790f50c26
# ╠═15cd3191-6f96-4016-8521-318f166ca05a
# ╠═e32e2abb-84af-4137-9e71-7e73266ed16a
# ╠═c706b1e8-58b9-4536-8cfe-68f5d8631b28
# ╠═5dd187d8-e7df-462a-b161-adbbd3c7dd19
# ╠═4d4a7cb2-07c5-4d55-960b-9dac0561c5fe
# ╠═a85190fa-6424-49e9-997e-5628aa1fa4a4
