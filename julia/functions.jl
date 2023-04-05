# functions

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






function MakePos(l::LinearEmbedding)     
    a = l.f.inner
    x = l.x.inner

    a_rle = rle((a.<0))
    rle_below_0 = a_rle[2][findall(a_rle[1].==1)]

    if (a[1] < 0)
        i = rle_below_0[1]
        α = (a[i+1]-a[i]) / (x[i+1]-x[i])
        x₀ = -a[i]/α + x[i]
        x[i] = x₀
        a[1:i] .= 0
        if (i == 1) insert!(x,1,1); insert!(a,1,0) end
    end
    if (a[end] < 0)
        i = length(x) - rle_below_0[end] + 1
        α = (a[i]-a[i-1]) / (x[i]-x[i-1])
        x₀ = -a[i-1]/α + x[i-1]
        x[i] = x₀
        a[i:end] .= 0
    end
    below_0 = findall(a.<0)
    while (!isempty(below_0))
        a_rle = rle((a.<0))
        rle_below_0 = a_rle[2][findall(a_rle[1].==1)]
        i = below_0[1]
        k = rle_below_0[1]
        αₗ = (a[i] - a[i-1]) / (x[i] - x[i-1])
        xₗ = -a[i-1]/αₗ + x[i-1]
        j = i + k - 1
        αᵣ = (a[j+1] - a[j]) / (x[j+1] - x[j])
        xᵣ = -a[j]/αᵣ + x[j]
        x[i] = xₗ
        if (i == j) 
            insert!(x,i+1,xᵣ); insert!(a,i+1,0) 
        else
            x[j] = xᵣ
        end
        a[i:j] .= 0
        below_0 = findall(a.<0)
    end
    l = LinearEmbedding(a,x)
    return l
end

