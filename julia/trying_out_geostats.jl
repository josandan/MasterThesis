### Example of Geostats usage

using GeoStats
using Plots
using GeoStatsPlots

# list of properties with coordinates
props = (Z=[1.,0.,1.],)
coord = [(25.,25.), (50.,75.), (75.,50.)]

# georeference data
𝒟 = georef(props, coord)

# estimation domain
𝒢 = CartesianGrid(100, 100)

# estimation problem
problem = EstimationProblem(𝒟, 𝒢, :Z)

# choose a solver from the list of solvers
solver = Kriging(
  :Z => (variogram=GaussianVariogram(range=35.),)
)

# solve the problem
solution = solve(problem, solver)

# plot the solution
contourf(solution, clabels=true)

#______________________________________________

### Another example

𝒟 = georef((Z=[10sin(i/10) + j for i in 1:100, j in 1:200],))

𝒮 = sample(𝒟, 500)

p1 = hscatter(𝒮, :Z, lag=0)
p2 = hscatter(𝒮, :Z, lag=20)
p3 = hscatter(𝒮, :Z, lag=40)
p4 = hscatter(𝒮, :Z, lag=60)

plot(p1, p2, p3, p4)

#______________________________________________


using GeoStats
using Plots
using GeoStatsPlots
using PointPatterns

p = PoissonProcess(0.1)
b = Box((0.,0.), (100.,100.))
s = rand(p, b, 1)
plot(s[1])

p₁ = BinomialProcess(50)
p₂ = BinomialProcess(50)
p  = p₁ ∪ p₂ # 100 points

s = rand(p, b, 2)

plot(plot(s[1]), plot(s[2]))

ishomogeneous(p)

p = PoissonProcess()
s = rand(p, b, 1)
plot(s[1])




