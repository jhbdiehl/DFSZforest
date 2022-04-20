using Pkg
Pkg.activate("..")

println("Hello there!")

using Symbolics
using Combinatorics
using LinearAlgebra

# For runtime optimization
using BenchmarkTools

using FileIO, JLD2

using StatsPlots
using Plots

include("./helpers.jl")

@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 



model = [u1, u2, u1, d1, d2, d1, l1, l1, l1]
nH = length(unique(model))
quads, multis = get_quads(model)
quads
notu1d1 = isequal.(unique(model),u1) .+ isequal.(unique(model),d1) .== false

isequal.(unique(model)[notu1d1], u2)

# Quads exactly like in mathematica for specific n=5 model

# Comes from Combinatorics package
# Equivalent to Mathematica Subsets
subsets = collect(powerset(quads, nH-2,nH-2));
println("Number of subsets: "*string(length(subsets)))

# Below would be way faster when specifying data type for Array{Any}. However "SymbolicUtils.Sym{Real, Base.ImmutableDict{DataType, Any}}" for some reason doesnt work!
# Equivalent to Mathematica Solve
sol = Array{Num}(undef, nH-2, length(subsets))
errors = []
infties = []
# difference to mathematica: does not accept solutions like {u2->u1, d2->l1} (no third condition)
@time Threads.@threads for i in range(1, size(sol)[2])
    try
        sol[:,i] = Symbolics.solve_for(subsets[i], unique(model)[notu1d1])
    catch
        append!(errors, i)
        sol[:,i] = [NaN for i in 1:nH-2]
    end
end
length(errors)

AnomalyRatio = Array{Num}(undef,size(sol)[2]);
@time Threads.@threads for i in range(1,size(sol)[2])
    try
        mnew = deepcopy(model)
        for j in 1:length(unique(model))-2
            mnew[isequal.(model,unique(model)[notu1d1][j])] .= sol[j,i]
        end
        AnomalyRatio[i] = EoverN(mnew[1:3], mnew[4:6], mnew[7:9])
    catch
        append!(infties, i)
        AnomalyRatio[i] = Inf
    end
end
length(infties)
print("Number of models (no errors): "*string(length(subsets)-length(errors)-length(infties)))

AnomalyRatio = substitute(AnomalyRatio, Dict([var => 1//1 for var in unique(model)]))
AnomalyRatio = Symbolics.value.(AnomalyRatio)
sol = sol[:,isfinite.(AnomalyRatio)]
subsets = subsets[isfinite.(AnomalyRatio)]
AnomalyRatio = AnomalyRatio[isfinite.(AnomalyRatio)]
checks(AnomalyRatio)

sol
subsets

AnomalyRatio[AnomalyRatio .== 8//3]

data = Dict(
    "AR" => AnomalyRatio,
    "sol" => sol,
    "subsets" => subsets,
    "errors" => sort(vcat(errors, infties)),
    "quads" => quads
)

@time FileIO.save("../data/test.jld2", data)
#print("Histogram time!")
@time histogram(AnomalyRatio, bins=-10:0.1:10)

"""
    Do the whole thing numerically
"""

using CUDA
using BenchmarkTools
using StatsBase


quadmat = zeros((length(unique(model)),length(quads)))
for j in 1:length(quads)
    vars = []
    for i in 1:length(unique(model))
        t = substitute(quads[j], Dict(unique(model)[i] => 1))
        t = substitute(t, Dict([var => 0 for var in unique(model)]))
        append!(vars, Symbolics.value(t))
    end
    quadmat[:,j] = vars
end

quadvecvec = [(quadmat[:,j][notu1d1], -sum(quadmat[:,j]) + sum(quadmat[:,j][notu1d1])) for j in 1:size(quadmat)[2]]

sol = [0,0,0]
k = 0
rawLES = collect(powerset(quadvecvec, nH-2, nH-2))
@time for i in 1:size(rawLES)[1]
    A = CuArray(mapreduce(permutedims, vcat, [rawLES[i][j][1] for j in 1:size(rawLES[1])[1]]))
    b = CuArray([rawLES[i][j][2] for j in 1:size(rawLES[1])[1]])
    try
        sol = hcat(sol, A\b)
        k += 1
    catch
        nothing
    end
end
sol = sol[:,2:end]


k = 0
sol = [0,0,0]
rawLES = collect(powerset(quadvecvec, nH-2, nH-2))
@time for i in 1:size(rawLES)[1]
    A = mapreduce(permutedims, vcat, [rawLES[i][j][1] for j in 1:size(rawLES[1])[1]])
    b = [rawLES[i][j][2] for j in 1:size(rawLES[1])[1]]
    try
        sol = hcat(sol, A\b)
        k += 1
    catch
        nothing
    end
end
sol = sol[:,2:end]

EoN = 2/3 .+ 2 .* (2 .+sol[1,:] .+ 3 .* sol[3,:]) ./ (2 .+ sol[1,:] .+ 2 .+sol[2,:])


h = fit(Histogram, EoN, -100.0000001:0.1:100)
h2 = fit(Histogram, AnomalyRatio, -100.0000001:0.1:100)

dif = h2.weights .- h.weights

histogram(sol, bins=-10:1:10)


#FileIO.load("./data/test.jld2", "AR")

#checked for n=3, n=4 (typical),  n=5 (first model)
#checked quads for [u1, u2, u1, d1, d2, d1, l1, l1, l1]