
@recipe function f(C::Curve{S,R,U}) where {S, R <: Real, U <: Real}
    P = [[p,p] for p in C.f.p]
    Q = Vector{U}[]
    n = length(P) 

    for i in 1:n
        push!(Q, [C.f.q[i].a, C.f.q[i].b])
    end

    for i in eachindex(P)
        if i == 1
            @series begin 
                seriestype := :path
                seriescolor --> :black
                label --> ""
                Q[i], P[i]
            end
        else
            @series begin
                seriestype := :path
                seriescolor --> :black
                label := ""
                Q[i], P[i]
            end
        end
    end

    
    @series begin
        seriestype := :scatter
        seriescolor --> :black
        label := ""
        [Q[i][2] for i in eachindex(Q)], [P[i][2] for i in eachindex(P)]
    end
    
end