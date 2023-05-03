
# Small example with nonparametric estimator of the intensity

T = 1
n = 10

plot(0.5+0.1*Normal(0,1), xlims=(0,1))
tₙ = rand(0.5+0.1*Normal(0,T), n) |> sort
plot(0.15*Exponential(1), xlims=(0,1))
tₙ = rand(0.15*Exponential(1), n) |> sort

scatter!(tₙ, zero, xlims=(0,1))

input = (pp = tₙ, k = 1, n = n)
Λ_est = PWLinearEstimator(input)
λ_est = DifferentiatePWLinear(tₙ, Λ_est)

plot!(tₙ, 0.02*λ_est.(tₙ), lt=:steppost)
ylims!(-0.2,5)

ϕᵋ = Mollifier(0.2)
MollIntensity = ∂(ϕᵋ)(Λ_est, tₙ)

plot!(tₙ, 0.05*MollIntensity.(tₙ))


# Showing that the estimated intensity is pw constant and the estimated cumulative intensity is pw linear

tₙ = rand(Uniform(0,T), n) |> sort
tₙ = tₙ |> sort

input = (pp = tₙ, k = 1, n = n)
Λ_est = PWLinearEstimator(input)
λ_est = DifferentiatePWLinear(tₙ, Λ_est)

plot1 = scatter(xlabel="Location")
plot!(tₙ, Λ_est, label = L"\hat{\Lambda}(t)", color = 2, marker = (:circle,3), markerstrokewidth=0)
plot!(Λ_est, label = "", color = 2)
plot!([-1], [-1], label = L"\hat{\lambda}(t)", color = 3, marker = (:circle,3), markerstrokewidth=0)
plot!(tₙ, λ_est, xlims=(-0.1,1.1), label = "", color = 3, marker = (:circle,3), markerstrokewidth=0, lt=:steppost)
plot!(λ_est, label = "", color = 3)
scatter!(tₙ, zero, label = L"t", color = 1)

savefig(plot1, "Figures/plot1.pdf")

