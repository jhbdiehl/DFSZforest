println("Hello there!")

using Symbolics
using Combinatorics

# For runtime optimization
using BenchmarkTools

using StatsPlots
using Plots

@variables u1::Int u2::Int d1::Int e1::Int

#Todo: implement path to nD4Final here!

nD4Final = [-2u1 + 2u2, d1-e1-u1+u2, -d1+e1-u1+u2, d1-u1+ 2u2, e1-u1+ 2u2, 2d1-2*e1, 2d1-e1+u2, u2+2*e1-d1, 2d1+ 2u2, d1+e1+ 2u2, 2u2+2*e1, 2d1-e1+u1, u1-d1+ 2*e1, 2d1+u1+u2, d1+e1+u1+u2, u1+2*e1+u2, d1+ 2u1 - u2, e1+ 2u1-u2, d1+e1+2u1, 2u1+2*e1]

@time begin

# Comes from Combinatorics package
# Equivalent to Mathematica Subsets
nD4Subsets = collect(powerset(nD4Final, 2,2));

# Below would be way faster when specifying data type for Array{Any}. However "SymbolicUtils.Sym{Real, Base.ImmutableDict{DataType, Any}}" for some reason doesnt work!
# Equivalent to Mathematica Solve
nD4sol = Array{Num}(undef, 2, length(nD4Subsets))
errors = 0
infties = 0
for i in range(1, length(nD4Subsets))
    try
        nD4sol[:,i] = Symbolics.solve_for(nD4Subsets[i], [e1,u2])
    catch
        errors += 1
        nD4sol[:,i] = [NaN, NaN]
    end
end


nD4AnomalyRatio = Array{Num}(undef,size(nD4sol)[2]);
for i in range(1,size(nD4sol)[2])
    try nD4AnomalyRatio[i] = EoverN(u1, nD4sol[2,i], d1, nD4sol[1,i])
    catch
        infties += 1
        nD4AnomalyRatio[i] = Inf
    end
end
nD4AnomalyRatio = substitute(nD4AnomalyRatio, Dict([u1 => 1, u2 => 1, d1 => 1, e1 => 1]))

function EoverN(u1, u2, d1, e1)
    2//3 + 2 * simplify_fractions((2u1 +u2+3*e1) //(2u1+u2+3d1))
end

nD4AnomalyRatio = Symbolics.value.(nD4AnomalyRatio)

end

histogram(nD4AnomalyRatio)