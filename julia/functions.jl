# functions

get_unique_df = function (raw_df)
    @chain raw_df begin
        @subset @byrow :Quantity != 0
        groupby([:Price, :Curve])
        combine(:Quantity => sum => :Quantity)
    end
end

function CombineRealizations(
    from_date::Date, to_date::Date, 
    hours::Union{AbstractVector{N}, Int64}, 
    side::String,
    point_process::String = "Quantity"
) where {N}
    dates = from_date:Day(1):to_date

    data = Float64[]
    for date = dates
        for hour in hours
            df = get_data(date, hour, side)
            cleaned_df = filter(:Price => !=(-500.0), df)
            append!(data, cleaned_df[:,point_process])
        end
    end

    k = length(dates) * length(hours)
    n = length(data)
    pp = data |> sort
    # pp = unique(round.(data, digits = 4)) |> sort

    return (pp = pp, k = k, n = n)
end

function CombineDF(
    dates::StepRange{Date, Day}, 
    hours::Union{AbstractVector{N}, Int64}, 
    side::String
) where {N}
    data = DataFrame(Price = Float64[], Quantity = Float64[], Side = [], Date = [], Hour = [], Curve = [])
    for date = dates
        for hour in hours
            df = get_data(date, hour, side)
            append!(data, df)
        end
    end

    return data
end

function CombinePriceDF(
    dates::StepRange{Date, Day}, 
    hours::Union{AbstractVector{N}, Int64}, 
    side::String
) where {N}
    DF = CombineDF(dates, hours, side)

    cleaned_df = @chain DF begin
        @subset @byrow :Price != -500.0
        @subset @byrow :Quantity != 0
    end

    grouped_df = @chain cleaned_df begin
        groupby([:Price, :Curve])
        combine(:Quantity => sum => :Quantity)
    end

    k = length(dates) * length(hours)
    n = k==1 ? nrow(grouped_df) : nrow(cleaned_df)

    return (DF = sort(grouped_df, :Price), k = k, n = n)
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

function SimulateHomPoisson(t₀::Real, T::Real, λ::Real)
	t = t₀
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

function LewisThinning(tₙ::AbstractVector{R}, intensity_func::Function) where {R}
	σ = Float64[]
	t₀ = min(tₙ[begin],0)
    T = tₙ[end]
	max_intensity = maximum(intensity_func.(tₙ))

	hom_sample = SimulateHomPoisson(t₀, T, max_intensity)

	for i in 1:length(hom_sample)
		u = rand(Uniform(0,1))
		if u <= intensity_func(hom_sample[i])/max_intensity
			push!(σ, hom_sample[i])
		end
	end
	return (hom_sample = hom_sample, sim = σ)
end

@views function makechunks(X::AbstractVector, n::Integer)
    c = length(X) ÷ n
    return [X[1+c*k:(k == n-1 ? end : c*k+c)] for k = 0:n-1]
end

# @views function makechunks(X::AbstractVector, n::Integer)
#     X_splits = [X[begin]:((X[end]-X[begin])/n):X[end]...]
#     X_left = X_splits[1:n]
#     X_right = X_splits[2:(n+1)]
#     X_right[n] = X_right[n]+1
#     return [X[(X.>=X_left[i]) .& (X.<X_right[i])] for i = 1:n]
# end

function OgataThinning(tₙ::AbstractVector{R}, intensity_func::Function) where {R}
	tₙ_chunks = makechunks(tₙ, 100)
    σ = Float64[]
    hom_sample = Float64[]
    # max_int = Float64[]

    for k in 1:length(tₙ_chunks)
        if !isempty(tₙ_chunks[k])
            max_intensity = maximum(intensity_func.(tₙ_chunks[k]))
            # push!(max_int, max_intensity)
            hom_sample_chunk = SimulateHomPoisson(tₙ_chunks[k][begin], tₙ_chunks[k][end], max_intensity)
            append!(hom_sample, hom_sample_chunk)
        
            for i in 1:length(hom_sample_chunk)
                u = rand(Uniform(0,1))
                if u <= intensity_func(hom_sample_chunk[i])/max_intensity
                    push!(σ, hom_sample_chunk[i])
                end
            end 
        end
    end
	return (hom_sample = hom_sample, sim = σ)
end




function MakePos(l::LinearEmbedding)     
    a, x = l.f.inner, l.x.inner

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
        if (i == length(x)) push!(x, length(x)); push!(a, 0) end
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

function RestrictFun(l::LinearEmbedding)
    eps = 1e-4
    f, x = l.f.inner, l.x.inner
    insert!(x, 1, x[begin]-eps)
    insert!(f, 1, 0)
    insert!(x, 1, x[begin]-eps)
    insert!(f, 1, 0)

    push!(x, x[end]+eps)
    push!(f, 0)
    push!(x, x[end]+eps)
    push!(f, 0)
    return LinearEmbedding(f,x)
end

function GetDensity(point_pattern::Vector{Float64}, k::Real, n::Real, mollifier_tolerance::Real)
    CumuIntensity = PWLinearEstimator((pp = point_pattern, k = k, n = n))
    # EstIntensity = DifferentiatePWLinear(point_pattern, CumuIntensity)
    ϕᵋ = Mollifier(mollifier_tolerance)
    MollCumuIntensity = ϕᵋ(CumuIntensity, point_pattern)
    # MollIntensity = ∂(ϕᵋ)(CumuIntensity, point_pattern)

    M = MollCumuIntensity.x[end]
    B = MollCumuIntensity.x[begin]
    x = MollCumuIntensity.x
    # Fₓ = map(x -> (MollCumuIntensity(x)-MollCumuIntensity(B))/(MollCumuIntensity(M)-MollCumuIntensity(B)), MollCumuIntensity.x)
    Fₓ = map(x -> (MollCumuIntensity(x))/(MollCumuIntensity(M)), MollCumuIntensity.x)

    return DiscretizedDistribution(LinearEmbedding(Fₓ,x))
end

# ((Λₙ[i+1] - Λₙ[i])/(tₙ[i+1] - tₙ[i])) * (p - tₙ[i]) + Λₙ[i]

function InverseDensity(Λₙ::AbstractVector{R}, tₙ::AbstractVector{Q}) where {R,Q}
	p -> begin
        i = searchsortedfirst(Λₙ, p)

        if i == 1
            return 0
        elseif i == length(tₙ) + 1
            # a = (Λₙ[i-1] - Λₙ[i-2])/(tₙ[i-1] - tₙ[i-2])
            # b = Λₙ[i-2]
            # c = tₙ[i-2]
            # return (Λₙ[i-1] + a*c - b)/a
            return Inf
        elseif Λₙ[i] == Λₙ[i-1]
            return Λₙ[i-1]
        else
            a = (Λₙ[i] - Λₙ[i-1])/(tₙ[i] - tₙ[i-1])
            b = Λₙ[i-1]
            c = tₙ[i-1]
            return (p + a*c - b)/a
        end
	end
end

# function InverseDensity(Λₙ::AbstractVector{R}, tₙ::AbstractVector{Q}) where {R,Q}
# 	p -> begin
#         i = searchsortedlast(Λₙ, p)

#         if i == 0
#             return 0
#         elseif i == length(tₙ)
#             a = (Λₙ[i] - Λₙ[i-1])/(tₙ[i] - tₙ[i-1])
#             b = Λₙ[i-1]
#             c = tₙ[i-1]
#             return (Λₙ[i] + a*c - b)/a
#         elseif Λₙ[i] == Λₙ[i+1]
#             return Λₙ[i]
#         else
#             a = (Λₙ[i+1] - Λₙ[i])/(tₙ[i+1] - tₙ[i])
#             b = Λₙ[i]
#             c = tₙ[i]
#             return (p + a*c - b)/a
#         end
# 	end
# end

function SimulateByInversion(Λ⁻¹::Function, T::Real)
    uᵢ = 0
    # uₙ = Float64[]
    tₙ = Float64[]
    t = 0
    while t ≤ T
        u = rand(Exponential(1), 1)[1]
        uᵢ = uᵢ + u
        # uₙ = push!(uₙ, uᵢ)
        t = Λ⁻¹(uᵢ)
        t ≤ T && push!(tₙ, t)
    end
    return tₙ
end
