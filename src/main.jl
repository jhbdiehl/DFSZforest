using LinearAlgebra, StaticArrays
using Random, Statistics, StatsBase
using BenchmarkTools
using Base.Threads
using Plots
using Random: default_rng
using Symbolics
using Combinatorics
using FileIO, JLD2


include("./drawer.jl")
include("./helpers.jl")

@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 

# 58 cores produces segmentation fault.
# 56 cores doesnt.

function make_idx_bnc(N)
    for (i,t) in enumerate(TupIter{N}())
        i > 10_000 && break
        @assert i == tupidx(t)
    end
    bnc = binom_cache(N, 1000);
    idxarr = similar([1],N);
    return idxarr, bnc
end

function parallel_randeqn_solve_proc!(
    proc_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer},
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer},
    tot::Int
) where N

    idxarr, bnc = make_idx_bnc(N)

    @threads for i in eachindex(proc_rs)
        @inbounds begin
            # When using drawer need to allocate indices. Is there a way around?
            idxs_i =  myidxtup!(idxarr, bnc, rand(1:tot), Val(N))#rand_idxs(default_rng(), eachindex(as), Val(N)) # drawer(13131) #
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = myEoN(r)
            rs_ws[i] = prod(ws[idxs_i])
        end
    end
end

function parallel_alleqn_solve_proc!(
    proc_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer},
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer},
    tot::Int
) where N

    idxarr, bnc = make_idx_bnc(N)

    @threads for i in eachindex(proc_rs)
        @inbounds begin
            # When using drawer need to allocate indices. Is there a way around?
            idxs_i =  myidxtup!(idxarr, bnc, i, Val(N))#rand_idxs(default_rng(), eachindex(as), Val(N)) # drawer(13131) #
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = myEoN(r)
            rs_ws[i] = prod(ws[idxs_i])
        end
    end
end

model = [u1, u2, u3, d1, d2, d3, l1, l1, l3]

using Printf
@time begin
    un = unique(model)
    nH = length(un)
    quads, multis = get_quads(model)
    tot = binomial(length(quads),nH-2)
    @printf "Your model-group has %.3E models \n" tot
    @time as, bs = get_numquads(quads, un, nH)
    myEoN = get_EoNfunc(model)
end


# Save all ARs for a specific model
#proc_rs = similar(bs, tot)
#rs_ws = similar(multis, length(proc_rs))
#@time parallel_alleqn_solve_proc!(proc_rs, rs_ws, as, bs, multis, tot)
#@time save_AR(model, proc_rs, rs_ws, 1)

#Save samples
chunk = 10^7
m=10
@time for i in 1:m
    @info "Computing round $i of $m"
    proc_rs = similar(bs, chunk)
    rs_ws = similar(multis, length(proc_rs))
    @time parallel_randeqn_solve_proc!(proc_rs, rs_ws, as, bs, multis, tot)
    @time save_AR(model, proc_rs, rs_ws, i)
end

# Read and plot data
tt = FileIO.load("./data/n8_u1_u2_u3_d1_d2_d3_l1_l1_l3.jld2", "ARs")
plot(tt, lt = :stepbins, xrange=(-10,13))