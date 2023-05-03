# Example in introduction

date = Date(2022,1,2)
hour = 12

df_s = get_data(date, hour, "Sell")
df_d = get_data(date, hour, "Buy")

S = Curve(get_data(date, hour, "Sell"))
D = Curve(get_data(date, hour, "Buy"))

p1 = scatter(df_s.Price, df_s.Quantity, label = "Supply bids", xlabel = "", ylabel = "Quantity", color=2)
p2 = scatter(df_d.Price, df_d.Quantity, label = "Demand bids", xlabel = "Price", ylabel = "Quantity", color=1)
plot1 = plot(p1,p2, layout=(2,1), size=(600,600), marker=(:circle,4))

scatter(df_s.Price, df_s.Quantity, label = "Supply bids", xlabel = "Price", ylabel = "Quantity", alpha = 0.5, markerstrokewidth = 0, marker=(:circle,3))
scatter!(df_d.Price, df_d.Quantity, label = "Demand bids", alpha = 0.5, markerstrokewidth = 0, marker=(:circle,3))

savefig(plot1, "Figures/ex_bids.pdf")

plot2 = plot(size=(600,400), legend = :outertopright)
plot!(S, color = 2, label = "Supply", marker = (:circle,3), markerstrokewidth = 0)
plot!(D, color = 1, label = "Demand", marker = (:circle,3), markerstrokewidth = 0)
scatter!([𝐐(S,D)],[𝐏(S,D)], label = "Clearing", color = :red)
xlabel!("Quantity")
ylabel!("Price")

savefig(plot2, "Figures/ex_curves.pdf")

p,q = (-300., 1000.)
plot3 = plot(size=(600,400), legend = :outertopright)
plot!(S, color = 2, label = "Supply before", alpha = 0.25, marker = (:circle,3), markerstrokewidth = 0)
plot!(S ⊕ˢ (p,q), color = 3, label = "Supply after", alpha = 0.25, marker = (:circle,3), markerstrokewidth = 0)
plot!(D, color = 1, label = "Demand", alpha = 0.25, marker = (:circle,3), markerstrokewidth = 0)
scatter!([𝐐(S,D)],[𝐏(S,D)], label = "Clearing before", color = :red)
scatter!([𝐐(S ⊕ˢ (p,q),D)],[𝐏(S ⊕ˢ (p,q),D)], label = "Clearing after", color = :yellow)
xlabel!("Quantity")
ylabel!("Price")

savefig(plot3, "Figures/ex_impact.pdf")

[𝐏(S ⊕ˢ (p,q),D)] - [𝐏(S,D)]
# 2000 MW -30, 1000 MW -7.3, 500 MW -4.1, 100 MW 0, 10 MW 0

p,q = (80., 1000.)

result = Float64[]
for q in 1.:2000
    r = 𝐏(S ⊕ˢ (p,q),D) - 𝐏(S,D)
    push!(result,r)
end
result
plot4 = plot(1.:2000, -result, label="")
xlabel!("Bid quantity")
ylabel!("Change in price")

savefig(plot4, "Figures/ex_price_changes.pdf")

plot5 = plot(xlabel="Bid quantity", ylabel="Change in price", color_palette = palette(:lightrainbow, 10))
for i = 1:10
    date = Date(2022,1,i)
    S = Curve(get_data(date, hour, "Sell"))
    D = Curve(get_data(date, hour, "Buy"))

    result = Float64[]
    for q in 1.:2000
        r = 𝐏(S ⊕ˢ (p,q),D) - 𝐏(S,D)
        push!(result,r)
    end
    plot!(1.:2000, -result, label=string(date))
end
plot!()

savefig(plot5, "Figures/ex_price_changes_10d.pdf")

palette([:orange, :red, :purple, :blue], 10)
palette(:lightrainbow, 10)
