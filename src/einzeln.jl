
using BenchmarkTools
using Symbolics
using Combinatorics
using LinearAlgebra
using StaticArrays
using IterTools, StatsBase
using FileIO, JLD2

include("./drawer.jl")
include("./helpers.jl")

@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 


# 58 cores produces segmentation fault.
# 56 cores doesnt.
model = [u1, u2, u1, d1, d2, d1, l1, l1, l1]
nsamps=nothing


################################################
function do_calculation(model, nsamps)
    un = unique(model)
    nH = length(un)
    quads = get_quads(model)
    tot = binomial(length(quads),nH-2)

    nsamps = checks(tot, nsamps, un)
    as, bs = get_numquads(quads, un, nH)
    myEoN = get_EoNfunc(model)
    drawer = construct_draw(nH)

    chunk = 10^7
    @time begin
    if typeof(nsamps) == Int
        chunk = minimum([chunk, nsamps])
        for j in 1:Int(round(nsamps/chunk))
            @time begin
            @info "You are in iteration $j out of $(Int(round(nsamps/chunk)))"
            sol_array = make_sol_array()
            Threads.@threads for i in rand(1:tot, chunk)
                id = drawer(i)
                sol = mysolve(as, bs, id)
                calculate_EoverN(myEoN, sol, sol_array)
            end
            save_AR(model, sol_array, j)
            end
        end
    elseif typeof(nsamps) == Nothing
        chunk = minimum([chunk, tot])
        for j in 1:Int(round(tot/chunk))
            @time begin
            @info "You are in iteration $j out of $(Int(round(tot/chunk)))"
            sol_array = make_sol_array()
            @time Threads.@threads for i in ((j-1)*chunk + 1):minimum([j*chunk, tot])
                id = drawer(i)
                sol = mysolve(as, bs, id)
                calculate_EoverN(myEoN, sol, sol_array)
            end
            save_AR(model, sol_array, j)
            end
        end
    end
    end
end

do_calculation(model, nsamps)

####################################################

#end

#A = [2.0 0.0 1.0 -1.0; 2.0 1.0 0.0 0.0; -2.0 -1.0 0.0 -3.0]
#inv(A[:,1:3])

#B = [1.0 0.0 0.0 -3.0; 2.0 0.0 0.0 -2.0; 2.0 0.0 1.0 -1.0]
#inv(B[:,1:3])
#tt = [2.0 0.0 1.0; 3.0 0.0 0.0; 2.0 0.0 7.0]
#inv(tt)
#@benchmark rank(tt)
#det(tt) == 0


#=
using StatsPlots
AR = read_AR("n6_u1_u2_u1_d1_d2_d1_l1_l2_l1")

hist = _make_hist(AR; bins=-10:0.001:13)
using Plots
plot()
plot_hist!(hist, yaxis=(:linear))
=#