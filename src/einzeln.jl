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
model = [u1, u2, u3, d1, d2, d3, l1, l2, l1]
nsamps=10^7


################################################
function do_calculation(model, nsamps)
    un = unique(model)
    nH = length(un)
    quads, multis = get_quads(model)
    tot = binomial(length(quads),nH-2)

    nsamps = checks(tot, nsamps, un)
    as, bs = get_numquads(quads, un, nH)
    myEoN = get_EoNfunc(model)
    uds = setdiff(model, [u1,d1,l1,l2,l3])
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
                sol = mysolve(as, bs, multis, id)
                calculate_EoverN(myEoN, uds, sol, sol_array)
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
                sol = mysolve(as, bs, multis, id)
                calculate_EoverN(myEoN, uds, sol, sol_array)
            end
            save_AR(model, sol_array, j)
            end
        end
    end
    end
end

do_calculation(model, nsamps)



un = unique(model)
nH = length(un)
quads, multis = get_quads(model)
tot = binomial(length(quads),nH-2)
multis
maximum(multis)
nsamps = checks(tot, nsamps, un)
as, bs = get_numquads(quads, un, nH)
myEoN = get_EoNfunc(model)
uds = setdiff(model, [u1,d1,l1,l2,l3])
drawer = construct_draw(nH)

sol_array = make_sol_array()

@time id = drawer.(1:10^6)
@benchmark sol = mysolve(as, bs, multis, id)
@benchmark calculate_EoverN(myEoN, uds, sol, sol_array)

sol_array

@time save_AR(model, sol_array, 3)

####################################################

#A = [2.0 0.0 1.0 -1.0; 2.0 1.0 0.0 0.0; -2.0 -1.0 0.0 -3.0]
#inv(A[:,1:3])

#B = [1.0 0.0 0.0 -3.0; 2.0 0.0 0.0 -2.0; 2.0 0.0 1.0 -1.0]
#inv(B[:,1:3])
#tt = [2.0 0.0 1.0; 3.0 0.0 0.0; 2.0 0.0 7.0]
#inv(tt)
#@benchmark rank(tt)
#det(tt) == 0

=
using StatsPlots
using Plots
AR = read_AR("n9_u1_u2_u3_d1_d2_d3_l1_l2_l3")
collect(keys(AR))[collect(values(AR)) .== maximum(values(AR))]
sum(values(AR))
tmpmat = permutedims(hcat(collect(keys(AR))...))
Ns = tmpmat[:,2]
ARs = tmpmat[:,1]
Vs = collect(values(AR))
AR1 = ARs[-0.001 .< Ns .< 0.001]
V1 = Vs[-0.001 .< Ns .< 0.001]
Ns[-0.001 .< Ns .< 0.001]
AR1cm = countmap(AR1, V1)
histAR1 = _make_hist(AR1cm; bins=-10:0.0001:13)
plot()
plot_hist!(histAR1, yaxis=(:linear))

histogram(Ns, bins=-10:0.01:10)
vline!([0,1,2,3,4])

histAR = _make_hist(AR; bins=-10:0.0001:13)
plot()
plot_hist!(histAR, yaxis=(:linear))
vline!([1.924, 2/3, 8/3], ls=:dash)

GAG = get_gag(AR, 40e-6)
histGAG = _make_hist(GAG; bins=-16.5:0.001:-12.5)
plot()
plot_hist!(histGAG, yaxis=(:linear))
vline!([gaγγ.([2/3,8/3], Ref(fa(40e-6)))], ls=:dash)

maximum(keys(AR))

AR2 = read_AR("n9_u1_u2_u3_d1_d2_d3_l1_l2_l3_u11_d11")
histAR2 = _make_hist(AR2; bins=-10:0.0001:13)
plot_hist!(histAR2, yaxis=(:linear))
a = histAR2[1].weights .- histAR[1].weights
plot(a)
GAG2 = get_gag(AR2, 40e-6)
histGAG2 = _make_hist(GAG2; bins=-16.5:0.001:-12.5)
plot()
plot_hist!(histGAG2, yaxis=(:linear))

=#