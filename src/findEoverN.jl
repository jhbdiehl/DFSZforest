println("Hello there!")

using Symbolics
using Combinatorics

# For runtime optimization
using BenchmarkTools

using StatsPlots
using Plots

@variables u1::Int u2::Int d1::Int l1::Int


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

nH = 4
model = [u1, u2, u1, d1, d1, d1, l1, l1, l1]
checks(nH, model) #Stop user right there if definition is nonsense.

notu1d1 = isequal.(unique(model),u1) .+ isequal.(unique(model),d1) .== false

# Comes from Combinatorics package
# Equivalent to Mathematica Subsets
subsets = collect(powerset(quads, nH-2,nH-2));

# Below would be way faster when specifying data type for Array{Any}. However "SymbolicUtils.Sym{Real, Base.ImmutableDict{DataType, Any}}" for some reason doesnt work!
# Equivalent to Mathematica Solve
sol = Array{Num}(undef, nH-2, length(subsets))
errors = 0
infties = 0
for i in range(1, length(subsets))
    try
        sol[:,i] = Symbolics.solve_for(subsets[i], [l1, u2])
    catch
        errors += 1
        sol[:,i] = [NaN for i in 1:nH-2]
    end
end
errors

sol

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

AnomalyRatio = substitute(AnomalyRatio, Dict([var => 1 for var in unique(model)]))

AnomalyRatio = Symbolics.value.(AnomalyRatio)

#end

histogram(AnomalyRatio)