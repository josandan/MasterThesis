# Example simulation
using QuadGK

T = 100
λ = x -> 0≤x≤T ? 1+sin(0.1*x) : 0
# Λ = x -> x>T ? T-10*cos(0.1*T) : max(0, x-10*cos(0.1*x))
Λ = x -> x≤T ? quadgk(λ, -100,x)[1] : quadgk(λ, -100,100)[1]
tₙ = [0.1:0.1:T...]
Λ.(tₙ)
Λ⁻¹ = InverseDensity(Λ.(tₙ), tₙ)

p1 = plot(λ, xlim=(-2.5,T+10), label=L"\lambda", color = 1)
p2 = plot(Λ, xlim=(-2.5,T+10), label=L"\Lambda", color = 2)
plot!(Λ⁻¹, label=L"\Lambda^{-1}", color = 3, xlabel="t")
plot1 = plot(p1, p2, layout=(2,1), size=(600,600), ylabel="Intensity")

# savefig(plot1, "Figures/ex_sim_intensity.pdf")

sim_inv = SimulateByInversion(Λ⁻¹, T)
plot2 = plot(λ, xlim=(-2.5,T+2.5), ylim=(-0.1,2.1), label=L"\lambda", xlabel="t", ylabel="Intensity")
scatter!(sim_inv, zero, label=L"t_1,\ldots,t_n", color = 3)

# savefig(plot2, "Figures/ex_sim_inversion.pdf")

ogata = OgataThinning(tₙ, λ)
p3 = plot(λ, xlim=(-2.5,T+2.5), ylim=(-0.1,2.1), label=L"\lambda", xlabel="t", ylabel="Intensity")
# evt piecewise linear λₘ here
scatter!(ogata.hom_sample, zero, markerstrokewidth=0.7, label=L"u_1,\ldots,u_m")
scatter!(ogata.sim, zero, color = 3, marker=:star5, markerstrokewidth=0.3, label=L"t_1,\ldots,t_n")

lewis = LewisThinning(tₙ, λ)
p4 = plot(λ, xlim=(-2.5,T+2.5), ylim=(-0.1,2.1), label=L"\lambda", ylabel="Intensity")
# hline!([2], label=L"\lambda_m", linestyle=:dash)
scatter!(lewis.hom_sample, zero, color = 2, markerstrokewidth=0.7, label=L"u_1,\ldots,u_m")
scatter!(lewis.sim, zero, color = 3, marker=:star5, markerstrokewidth=0.3, label=L"t_1,\ldots,t_n")

plot3 = plot(p4, p3, layout=(2,1), size=(600,600))

# savefig(plot3, "Figures/ex_sim_ogata.pdf")

print(string(
    "Number of kept points are ", length(ogata.sim), "/", length(ogata.hom_sample), 
    " and ", length(lewis.sim), "/", length(lewis.hom_sample), 
    " for Ogata and Lewis, respectively."
))

Λ_hat_inv = PWLinearEstimator((pp = sim_inv, k = 1, n = length(sim_inv)))
Λ_hat_lewis = PWLinearEstimator((pp = lewis.sim, k = 1, n = length(lewis.sim)))
Λ_hat_ogata = PWLinearEstimator((pp = ogata.sim, k = 1, n = length(ogata.sim)))

plot4 = plot(Λ, xlim=(-2.5,T+2.5), label=L"\Lambda", xlabel="t", ylabel="Intensity")
plot!(Λ_hat_inv, label=L"\hat\Lambda \,(Inversion)")
plot!(Λ_hat_lewis, label=L"\hat\Lambda \,(Lewis)")
plot!(Λ_hat_ogata, label=L"\hat\Lambda \,(Ogata)")

# savefig(plot4, "Figures/ex_estim_cumu_int.pdf")

λ_hat_inv = DifferentiatePWLinear(sim_inv, Λ_hat_inv)
λ_hat_lewis = DifferentiatePWLinear(lewis.sim, Λ_hat_lewis)
λ_hat_ogata = DifferentiatePWLinear(ogata.sim, Λ_hat_ogata)

p5 = plot(λ, xlim=(-2.5,T+2.5), label=L"\lambda", ylabel="Intensity")
plot!(λ_hat_inv, label=L"\hat\lambda \,(Inversion)")
plot!(λ_hat_lewis, label=L"\hat\lambda \,(Lewis)")
plot!(λ_hat_ogata, label=L"\hat\lambda \,(Ogata)")

p6 = plot(λ, xlim=(-2.5,T+2.5), label="", ylabel="Intensity", ylim=(-0.1,3.1), xlabel="t")
plot!(λ_hat_inv, label="")
plot!(λ_hat_lewis, label="")
plot!(λ_hat_ogata, label="")

plot5 = plot(p5,p6, layout=(2,1), size=(600,600))

# savefig(plot5, "Figures/ex_estim_int.pdf")

ϕᵋ = Mollifier(10.)
Λ_moll_inv = ϕᵋ(Λ_hat_inv, sim_inv)
Λ_moll_lewis = ϕᵋ(Λ_hat_lewis, lewis.sim)
Λ_moll_ogata = ϕᵋ(Λ_hat_ogata, ogata.sim)
λ_moll_inv = ∂(ϕᵋ)(Λ_hat_inv, sim_inv)
λ_moll_lewis = ∂(ϕᵋ)(Λ_hat_lewis, lewis.sim)
λ_moll_ogata = ∂(ϕᵋ)(Λ_hat_ogata, ogata.sim)

plot6 = plot(Λ, xlim=(-2.5,T+2.5), label=L"\Lambda", ylabel="Intensity", xlabel="t")
plot!(Λ_moll_inv, label=L"\hat\Lambda*\varphi_\epsilon \,(Inversion)")
plot!(Λ_moll_lewis, label=L"\hat\Lambda*\varphi_\epsilon \,(Lewis)")
plot!(Λ_moll_ogata, label=L"\hat\Lambda*\varphi_\epsilon \,(Ogata)")

# savefig(plot6, "Figures/ex_moll_cumu_int.pdf")

plot7 = plot(λ, xlim=(-2.5,T+2.5), label=L"\lambda", ylabel="Intensity", xlabel="t")
plot!(λ_moll_inv, label=L"\partial (\hat\Lambda*\varphi_\epsilon) \,(Inversion)")
plot!(λ_moll_lewis, label=L"\partial (\hat\Lambda*\varphi_\epsilon) \,(Lewis)")
plot!(λ_moll_ogata, label=L"\partial (\hat\Lambda*\varphi_\epsilon) \,(Ogata)")
plot!(legend=:topleft, ylim=(-0.1,3.3))

# savefig(plot7, "Figures/ex_moll_int.pdf")


# =====================================================
# Residual analysis

# true intensity with simulated points
s = Λ.(sim_inv)
Δs = [s[1], diff(s)...]
h7 = histogram(Δs, normalize = true, label=L"\Delta u", title=L"\mathrm{True}\,\, \Lambda",title_align=:left)
plot!(t -> (t ≥ 0)*exp(-t), label=L"\mathrm{PDF\,\, of}\,\, Exp(1)")
q7 = qqplot(Δs, Exponential(1), xlabel="Δu", ylabel="Exp(1)", markerstrokewidth=0.5)

# estimated intensities
s1 = Λ_hat_inv.(sim_inv)
s2 = Λ_hat_lewis.(lewis.sim)
s3 = Λ_hat_ogata.(ogata.sim)
Δs1 = [s1[1], diff(s1)...]
Δs2 = [s2[1], diff(s2)...]
Δs3 = [s3[1], diff(s3)...]
h1 = histogram(Δs1, normalize = true, label=L"\Delta u", title=L"\mathrm{Testing}\,\, \hat\Lambda \,\,\mathrm{(Inversion)}",title_align=:left)
xticks!(LinRange(minimum(Δs1),maximum(Δs1),3))
plot!(t -> (t ≥ 0)*exp(-t), label=L"\mathrm{PDF\,\, of }\,\, Exp(1)")
h2 = histogram(Δs2, normalize = true, label=L"\Delta u", title=L"\mathrm{Testing}\,\, \hat\Lambda \,\,\mathrm{(Lewis)}",title_align=:left)
xticks!(LinRange(minimum(Δs2),maximum(Δs2),3))
plot!(t -> (t ≥ 0)*exp(-t), label=L"\mathrm{PDF\,\, of }\,\, Exp(1)")
h3 = histogram(Δs3, normalize = true, label=L"\Delta u", title=L"\mathrm{Testing}\,\, \hat\Lambda \,\,\mathrm{(Ogata)}",title_align=:left)
xticks!(LinRange(minimum(Δs3),maximum(Δs3),3))
plot!(t -> (t ≥ 0)*exp(-t), label=L"\mathrm{PDF\,\, of }\,\, Exp(1)")

q1 = qqplot(Δs1, Exponential(1), xlabel="Δu", ylabel="Exp(1)", markerstrokewidth=0.5)
xticks!(LinRange(minimum(Δs1),maximum(Δs1),3))
q2 = qqplot(Δs2, Exponential(1), xlabel="Δu", ylabel="Exp(1)", markerstrokewidth=0.5)
xticks!(LinRange(minimum(Δs2),maximum(Δs2),3))
q3 = qqplot(Δs3, Exponential(1), xlabel="Δu", ylabel="Exp(1)", markerstrokewidth=0.5)
xticks!(LinRange(minimum(Δs3),maximum(Δs3),3))

plot8 = plot(h1,q1,h2,q2,h3,q3, layout=(3,2), size=(800,800), right_margin=1cm)

savefig(plot8, "Figures/ex_sim_res_analysis.pdf")

s1 = Λ_moll_inv.(sim_inv)
s2 = Λ_moll_lewis.(lewis.sim)
s3 = Λ_moll_ogata.(ogata.sim)
Δs1 = [s1[1], diff(s1)...]
Δs2 = [s2[1], diff(s2)...]
Δs3 = [s3[1], diff(s3)...]
h4 = histogram(Δs1, normalize = true, label=L"\Delta u", title=L"\mathrm{Testing}\,\, \hat\Lambda*\varphi_\epsilon \,\,\mathrm{(Inversion)}",title_align=:left)
plot!(t -> (t ≥ 0)*exp(-t), label=L"\mathrm{PDF\,\, of }\,\, Exp(1)")
h5 = histogram(Δs2, normalize = true, label=L"\Delta u", title=L"\mathrm{Testing }\,\, \hat\Lambda*\varphi_\epsilon \,\,\mathrm{(Lewis)}",title_align=:left)
plot!(t -> (t ≥ 0)*exp(-t), label=L"\mathrm{PDF\,\, of }\,\, Exp(1)")
h6 = histogram(Δs3, normalize = true, label=L"\Delta u", title=L"\mathrm{Testing }\,\, \hat\Lambda*\varphi_\epsilon \,\,\mathrm{(Ogata)}",title_align=:left)
plot!(t -> (t ≥ 0)*exp(-t), label=L"\mathrm{PDF\,\, of }\,\, Exp(1)")

q4 = qqplot(Δs1, Exponential(1), xlabel="Δu", ylabel="Exp(1)", markerstrokewidth=0.5)
q5 = qqplot(Δs2, Exponential(1), xlabel="Δu", ylabel="Exp(1)", markerstrokewidth=0.5)
q6 = qqplot(Δs3, Exponential(1), xlabel="Δu", ylabel="Exp(1)", markerstrokewidth=0.5)

plot9 = plot(h4,q4,h5,q5,h6,q6, layout=(3,2), size=(800,800))

savefig(plot9, "Figures/ex_moll_res_analysis.pdf")

