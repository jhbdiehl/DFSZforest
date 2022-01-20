println("Hello there!")

using Symbolics
using Combinatorics
using LinearAlgebra

# For runtime optimization
using BenchmarkTools

using StatsPlots
using Plots

@variables u1::Int u2::Int d1::Int l1::Int d2::Int

"""
    Careful! This ignores different quadrilinears leading to the same condition! This is maybe not the behavior we want!
    (Also this function is incredibly dirty, holy cow!)
"""
function get_quads(model)
    vars = unique(model)
    try
        try
            us = [u1, u2, u3]
        catch
            us = [u1, u2]
        end
    catch
        us = [u1]
    end
    doublets = collect(powerset(vars, 2,2))
    diffs = [length(setdiff(doublets[i], us)) for i in 1:length(doublets)]
    sig = isodd.(diffs) .* 2 .- 1 # lifehack! Maps BitVector (0 for false, 1 for true) to (-1 for false and 1 for true)
    sig = hcat(sig, ones(size(sig)))
    doublets = permutedims(hcat(doublets...))
    doublets .*=sig
    doublets = [doublets[i,:] for i in 1:size(doublets, 1)]
    allquads = collect(powerset(doublets, 2,2))
    allquads = permutedims(hcat(allquads...))
    allquads_conj = deepcopy(allquads)
    allquads_conj[:,1] .*= -1
    allquads_conj
    allquads = vcat(allquads, allquads_conj)
    allquads = [sum(sum(allquads[i,:])) for i in 1:size(allquads)[1]]
    ut = tril(permutedims(hcat([isequal.(allquads, -allquads[i]) for i in 1:length(allquads)]...)),-1) # make sure to remove duplicates like a = -b
    fut = findall(ut)
    for i in 1:length(fut)
        allquads[fut[i][1]] = allquads[fut[i][2]]
    end
    uniquads = unique(allquads) # remove all duplicates
    fillones = substitute(uniquads, Dict([var => 1 for var in vars]))
    fillones = (sign.(fillones) .== 0) + sign.(fillones) # make sign return 1 if element = 0
    fillones = (fillones .== 1) .* 2 .- 1 # weird hack to make sure the array contains only Ints, because unique() does not work with negative floats and symbols
    uniquads .*= fillones
    quads = uniquads[isequal.(uniquads, d1+u1) .!= 1] # effectively is not equal...
    quads .*= (length.(string.(quads)) .<= 7) .+ 1 # weird way of multiplying all quads by 2 if they only contain two different higgs. This is so sum(abs(higgs)) == 4
end

using PyCall

@time begin
py"""
import numpy as np
import sympy as sy

u1, u2, d1, e1 = sy.symbols('u1 u2 d1 e1')

def condition1(u1, u2, d1, e1):
    return u1 + u2 - d1 - e1

def condition2(u1, u2, d1, e1):
    return np.abs(u1) + np.abs(u2) + np.abs(d1) + np.abs(e1) - 4

res = sy.nonlinsolve([condition1(u1, u2, d1, e1), condition2(u1, u2, d1, e1), u1+u2+e1+d1 >= 0], [u1, u2, d1, e1])
print(res)
"""
end


using NLsolve

# 1: u1, 2: u2, 3:d1, 4:e1
function f!(Q, x)
    Q[1] = x[3] + x[4] - x[2]

end


Symbolics.solve_for([u1+d1 ~ 2, -2u1+3d1+3*e1-u2 ~ 0], [e1, u2])

quads = [-2u1 + 2u2, d1-l1-u1+u2, -d1+l1-u1+u2, d1-u1+ 2u2, l1-u1+ 2u2, 2d1-2*l1, 2d1-l1+u2, u2+2*l1-d1, 2d1+ 2u2, d1+l1+ 2u2, 2u2+2*l1, 2d1-l1+u1, u1-d1+ 2*l1, 2d1+u1+u2, d1+l1+u1+u2, u1+2*l1+u2, d1+ 2u1 - u2, l1+ 2u1-u2, d1+l1+2u1, 2u1+2*l1]

#@time begin

function EoverN(ups, downs, leptons)
    2//3 + 2 * simplify_fractions((sum(ups)+sum(leptons)) //(sum(ups)+sum(downs)))
end

function checks(nH::Int64, model::Vector{Num})
    length(model) == 9 ? nothing : throw(AssertionError("model does not specify couplings to all SM fermions"))
    lm = length(unique(model))
    lm == nH ? nothing : throw(AssertionError("model does not specify nH different hidings, but instead $lm"))
end

nH = 5
model = [u1, u2, u1, d1, d2, d1, l1, l1, l1]
checks(nH, model) #Stop user right there if definition is nonsense.
quads = get_quads(model)

notu1d1 = isequal.(unique(model),u1) .+ isequal.(unique(model),d1) .== false


# Comes from Combinatorics package
# Equivalent to Mathematica Subsets
subsets = collect(powerset(quads, nH-2,nH-2));
println(length(subsets))

# Below would be way faster when specifying data type for Array{Any}. However "SymbolicUtils.Sym{Real, Base.ImmutableDict{DataType, Any}}" for some reason doesnt work!
# Equivalent to Mathematica Solve
sol = Array{Num}(undef, nH-2, length(subsets))
errors = 0
infties = 0
for i in range(1, length(subsets))
    try
        sol[:,i] = Symbolics.solve_for(subsets[i], unique(model)[notu1d1])
    catch
        errors += 1
        sol[:,i] = [NaN for i in 1:nH-2]
    end
end
errors

AnomalyRatio = Array{Num}(undef,size(sol)[2]);
for i in range(1,size(sol)[2])
    try 
        mnew = deepcopy(model)
        for j in 1:length(unique(model))-2
            mnew[isequal.(model,unique(model)[notu1d1][j])] .= sol[j,i]
        end
        AnomalyRatio[i] = EoverN(mnew[1:3], mnew[4:6], mnew[7:9])
    catch
        infties += 1
        AnomalyRatio[i] = Inf
    end
end
infties
print(length(subsets)-errors-infties)

ss = deepcopy(AnomalyRatio)

AnomalyRatio = substitute(AnomalyRatio, Dict([var => 1 for var in unique(model)]))

AnomalyRatio = Symbolics.value.(AnomalyRatio)
AnomalyRatio = AnomalyRatio[AnomalyRatio .!= NaN]
AnomalyRatio = AnomalyRatio[isfinite.(AnomalyRatio)]
AnomalyRatio
minimum(AnomalyRatio)
sum(isnan.(AnomalyRatio))
#end
print("lol")
@time histogram(AnomalyRatio, bins=range(-10, stop = 20, length = 25))