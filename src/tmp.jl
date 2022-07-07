using LinearAlgebra, StaticArrays
using Random, Statistics, StatsBase
using BenchmarkTools
using Base.Threads
using Plots
using Random: default_rng
using Symbolics
using Combinatorics
using FileIO, JLD2, HDF5
using LaTeXStrings, Printf


include("./drawer.jl")
include("./helpers.jl")

@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 

# 58 cores produces segmentation fault.
# 56 cores doesnt.

#dataset = "test"

"""
    Functionality:
        run for only one n value. run for specified bilinears. run all models or run cheaper version with assumed multiplicities.
        run
"""


compute_equivalent_theories=false
ns=[4]

if compute_equivalent_theories
    a, ms = generate_all_models()
else
    a, ms = generate_unique_models()
end

if ns == :all
    mya = a
else
    mya = a[length.(unique.(a)) .∈ Ref(ns)]
    ms = ms[length.(unique.(a)) .∈ Ref(ns)]
end

mya

wslists = Array{Any}(undef, length(mya))
rslists = Array{Any}(undef, length(mya))
EoNlists = Array{Any}(undef, length(mya))
ttlists = Array{Any}(undef, length(mya))
dataset="test"
@time for (k, model) in enumerate(mya)
    println(model)
    bilins = get_bilins(model)
    myEoN = get_EoNfunc(model)
    myN = get_Nfunc(model)
    quads, multis = get_quads(model)
    nH = length(unique(model))
    EoNlist = []
    rslist = []
    tt = []
    wslist = []
    for (i,bilin) in enumerate(bilins)
        if nH < 9
            terms = vcat(quads, bilins[i+1:end])
            mult = vcat([1,1], multis, ones(Int64, length(bilins[i+1:end])))
        else
            terms = vcat(quads, bilins[1:end .!= i])
            mult = vcat([1,1], multis, ones(Int64, length(bilins)-1))
        end
        terms = vcat(orthogonality(model), bilin, terms)
        as, bs = get_numquads(terms, unique(model), nH)
        tot = binomial(length(terms)-2,nH-2)
        if i == 1
            tt = terms
        end
        proc_rs = similar(as, tot)
        EoN_rs = similar(bs, tot)
        rs_ws = similar(mult, length(proc_rs))
        parallel_alleqn_solve_proc_fullsol!(proc_rs, EoN_rs, rs_ws, as, bs, mult, tot, myEoN, myN)
        #save_full_NEW(model, proc_rs, EoN_rs, rs_ws, 1; folder=dataset*"/n"*string(nH)*"/", bilin=bilin, ms=ms[1])
        good_idxs = findall(!isnan, EoN_rs)
        EoN_rs = EoN_rs[good_idxs]
        rs_ws = rs_ws[good_idxs]
        append!(EoNlist, EoN_rs[-1e10 .< EoN_rs .< 1e10])
        append!(wslist, rs_ws[-1e10 .< EoN_rs .< 1e10])
        good_proc_rs = proc_rs[good_idxs,:]
        append!(rslist, proc_rs)
    end
    EoNlists[k] = EoNlist
    rslists[k] = rslist
    ttlists[k] = tt
    wslists[k] = wslist
end


check_symmetry(EoNlists; wslists=wslists)

