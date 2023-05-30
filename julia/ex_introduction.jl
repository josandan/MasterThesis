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

𝐏(S ⊕ˢ (p,q),D) - 𝐏(S,D)
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

p = -300.
hours = 12
PriceChange = DataFrame(
    Day=Int64[],
    Q100=Float64[], Q200=Float64[], Q300=Float64[], Q400=Float64[], Q500=Float64[], 
    Q600=Float64[], Q700=Float64[], Q800=Float64[], Q900=Float64[], Q1000=Float64[]
)
for i = 1:31
    date = Date(2022,1,i)
    S = Curve(get_data(date, hours, "Sell"))
    D = Curve(get_data(date, hours, "Buy"))

    result = Float64[]
    for q in 100.:100:1000
        r = 𝐏(S ⊕ˢ (p,q),D) - 𝐏(S,D)
        push!(result,r)
    end
    push!(PriceChange, [i,result...])
end

plot5 = plot(xlabel="Day in January", ylabel="Change in price", color_palette = palette(:lightrainbow, 10))
for i in 1:10
    plot!(PriceChange[:,string("Q",i*100)], label=string(i*100,"MW"), color=i)
end
plot!()

savefig(plot5, "Figures/ex_price_changes_jan.pdf")



PriceChange = DataFrame(
    Hour=Int64[],
    Q100=Float64[], Q200=Float64[], Q300=Float64[], Q400=Float64[], Q500=Float64[], 
    Q600=Float64[], Q700=Float64[], Q800=Float64[], Q900=Float64[], Q1000=Float64[]
)
for hour = 1:24
    date = Date(2022,1,1)
    S = Curve(get_data(date, hour, "Sell"))
    D = Curve(get_data(date, hour, "Buy"))

    result = Float64[]
    for q in 100.:100:1000
        r = 𝐏(S ⊕ˢ (p,q),D) - 𝐏(S,D)
        push!(result,r)
    end
    push!(PriceChange, [hour ,result...])
end

plot6 = plot(xlabel="Hour", ylabel="Change in price", color_palette = palette(:lightrainbow, 10))
for i in 1:10
    plot!(PriceChange[:,string("Q",i*100)], label=string(i*100,"MW"), color=i)
end
plot!()

savefig(plot6, "Figures/ex_price_changes_hours.pdf")



palette([:orange, :red, :purple, :blue], 10)
palette(:lightrainbow, 10)
