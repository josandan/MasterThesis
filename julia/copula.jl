using Pkg
Pkg.activate(".")

# Packages and functions: 
using Distributions, Embeddings, StatsBase
using EmpiricalCopulas, Chain, DataFramesMeta

# using KernelDensity

include("SampleScript\\Includes.jl")

# get_index() = begin 
#     db = DB("./julia/TestData/DayAhead.db")
#     execute(db, "SELECT * FROM 'index'") |> DataFrame 
# end 

# index = get_index()

# Get Supply and Demand 
data = get_data(Date(2022,1,1),1,"Buy") 
data = @chain data begin 
    groupby([:Price,:Curve]) 
    combine(:Quantity => sum => :Quantity) 
    @subset @byrow :Quantity != 0 
end

Demand = Curve(data) 
𝐛 = bids(Demand)

# Plot Buy Bids 
scatter(𝐛, label = "") 
xlabel!("Price") 
ylabel!("Quantity")

bid_dataframe = DataFrame(:Price => [b[1] for b in 𝐛],:Quantity => [b[2] for b in 𝐛])

X,Y = bid_dataframe[:,:Price], bid_dataframe[:,:Quantity]

F̂ = ecdf(bid_dataframe[:,:Price]) 
Ĝ = ecdf(bid_dataframe[:,:Quantity])

plot(LinRange(minimum(Y), maximum(Y), 1001),x -> Ĝ(x))

F̂⁻¹ = x -> quantile(bid_dataframe[:,:Price], x)
Ĝ⁻¹ = x -> quantile(bid_dataframe[:,:Quantity], x)

# Fit Copula 
U = F̂.(bid_dataframe[:,:Price]) 
V = Ĝ.(bid_dataframe[:,:Quantity]) 
C = BernsteinCopula(U,V)

gradient(C, [0.5,0.5])

# Plot
Z = [pdf(C,[u,v]) for u in LinRange(0.01,0.99,101), v in LinRange(0.01,0.99,101)]
heatmap(Z)

W = rand(C, 227)' 
Û = W[:,1] 
V̂ = W[:,2] 

X̂ = F̂⁻¹.(Û) 
Ŷ = Ĝ⁻¹.(V̂) 

scatter(X,Y) 
scatter!(X̂,Ŷ)

Demand₀ = DataFrame(:Price => X̂, :Quantity => Ŷ, :Curve => :Demand) |> Curve

plot(Demand, color = 1)
plot!(Demand₀, color = 2)
