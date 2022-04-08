
#using CUDA
using BenchmarkTools
using Symbolics
using Combinatorics
using LinearAlgebra
#using CUDA.CUSPARSE
#using SparseArrays
#using StatsPlots
using StaticArrays
using IterTools, StatsBase
using FileIO, JLD2

include("./drawer.jl")
include("./helpers.jl")

@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 


#histogram!(AR, bins=-10:0.1:13)

function get_numquads(quads, un, nH)
    notu1d1 = isequal.(un,u1) .+ isequal.(un,d1) .== false
    numquads = zeros(Int8,length(quads), nH-2)

    u1d1dict = Dict{Num, Int8}(un .=> 0)
    u1d1dict[u1] = 1
    u1d1dict[d1] = 1
    bs = Vector{Int8}(-1 .* Symbolics.value.(substitute.(quads, (u1d1dict,))))

    for i in 1:nH-2
        mydict = Dict(un .=> 0.0)
        mydict[un[notu1d1][i]] = 1.0
        numquads[:,i] = Symbolics.value.(substitute.(quads, (mydict,)))
    end
    return _vecvec(SMatrix{length(quads), nH-2, Int8, (nH-2)*length(quads)}(numquads)), bs
end

function mysolve(as, bs, idxs)
    A = copy(hcat(as[idxs]...)')
    b = bs[idxs]
	if _issolvable(A)
	    return A \ b
    else
		return SVector{length(b),Float64}(NaN for i in 1:length(b))
	end
end

function _issolvable(mat)
    #rank(mat) == size(mat)[1]
    ≈(det(mat), 0.0; atol=1e-10) == false # tolerance chosen arbitrarily. Make sure this fits estimated numerical errors!
end

function calculate_EoverN(myEoN, sol, sol_array)
    AnomalyRatio = myEoN(sol...)
    if AnomalyRatio < 1e3 && AnomalyRatio > -1e3
        push!(sol_array[Threads.threadid()],AnomalyRatio)
    end
end

function make_sol_array()
    solution_data = Vector{Vector{Float64}}()
    for i in 1:Threads.nthreads()
        push!(solution_data,Float64[])
    end
    return solution_data
end





# 58 cores produces segmentation fault.
# 56 cores doesnt.
model = [u1, u2, u3, d1, d2, d3, l1, l2, l3]
nsamps=1 * 10^9


################################################
#function do_calculation(model, nsamps)
un = unique(model)
nH = length(un)
quads = get_quads(model)
tot = binomial(length(quads),nH-2)

nsamps = checks(tot, nsamps)
as, bs = get_numquads(quads, un, nH)
myEoN = get_EoNfunc(model)
drawer = construct_draw(nH)

chunk = 10^7
@time begin
if typeof(nsamps) == Int
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
@info "Maximal E/N: $(maximum(abs.(AR)))"
histogram(AR, bins=-1000:0.1:1000)

sum(collect(values(AR)))


hist = _make_hist(AR; bins=-10:0.001:13)
using Plots
plot()
plot_hist!(hist, yaxis=(:linear))

function _make_hist(cm; bins=-10:0.0001:13)
    hi = fit(Histogram, collect(keys(cm)), Weights(collect(values(cm))),bins)#, nbins=5000)
    m = (hi.edges[1] .+ (hi.edges[1][2] - hi.edges[1][1])/2.)[1:end-1]
    #hi.weights = Vector{Real}(hi.weights ./ sum(hi.weights))
    return hi, m
end

function plot_hist!(hist; yaxis=(:log, [0.00000001, :auto]), kwargs...)
    plot!(hist[2], hist[1].weights ./ sum(hist[1].weights); yaxis=yaxis, kwargs...)
    #plot!(hist[2], _gauss_prediction(gauss,hist[1].edges[1]); kwargs..., label="")
end
=#