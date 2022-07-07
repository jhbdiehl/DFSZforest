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
    myEoN = get_EoNfunc_NEW(model)
    myN = get_Nfunc(model)
    quads, multis = get_quads_NEW(model)
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
        as, bs = get_numquads_NEW(terms, unique(model), nH)
        tot = binomial(length(terms)-2,nH-2)
        if i == 1
            tt = terms
        end
        proc_rs = similar(as, tot)
        EoN_rs = similar(bs, tot)
        rs_ws = similar(mult, length(proc_rs))
        parallel_alleqn_solve_NEW_proc_fullsol!(proc_rs, EoN_rs, rs_ws, as, bs, mult, tot, myEoN, myN)
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



eon1 = round.(EoNlists[7][-1e10 .< EoNlists[7] .< 1e10], digits=4)
eon2 = round.(EoNlists[9][-1e10 .< EoNlists[9] .< 1e10], digits=4)
eon3 = round.(EoNlists[30][-1e10 .< EoNlists[30] .< 1e10], digits=6)
histogram(eon3, bins=-10+5/3:0.01:10+5/3)


EoN = [round.(EoNlist, digits=6) for EoNlist in EoNlists]
ws = Int.(vcat(wslists...))
EoN = vcat(EoN...)
cEoN = countmap(EoN, ws)
sEoN = Dict(round.(-1 .* keys(cEoN) .+ 3.33333333, digits=4) .=> values(cEoN))
EoN = [round.(EoNlist, digits=4) for EoNlist in EoNlists]
EoN = vcat(EoN...)
cEoN = countmap(EoN, ws)
mEoN = mergewith( (x,y) -> x-y, countmap(EoN), sEoN)
[mEoN[x] for x in keys(mEoN) if abs(x) .== 0.0] #.-5/3
sum(abs.(values(mEoN))) / 2


fulleon = vcat(EoNlists...)
histogram(EoN, bins=-10+5/3:0.1:10+5/3, weights=ws)

[mEoN[x] for x in keys(mEoN) if -0.001 < x < 0.001]

@time for (k, model) in enumerate(mya) #a[length.(unique.(a)) .== 8]
    bilins = unique(sort.(collect(combinations(model,2)), by=x->Symbol(x)))
    bilins = bilins[length.(unique.(bilins)) .== 2]
    tott = 0
    for (i, bilin) in enumerate(bilins)
        valp1, valp2 = bilinvals(bilin)
        un = unique(model)
        nH = length(un)
        quads, multis = get_quads(model; p1=bilin[1], p2=bilin[2], valp1=valp1, valp2=valp2)
        if exactly_one_bilinear == true
            terms = quads
        elseif exactly_one_bilinear == false
            if length(un) <= sample_n_gt
                bi = bilinsum.(bilins[i+1:end])
            else # If you sample all bilins with equal nr of samples, excluding models you already calculated would lead to bias effect!
                bi = bilinsum.(bilins[1:end])
            end
            terms = vcat(quads,bi)
            multis = vcat(multis, bilin_weights * ones(Int64, size(bi)...))
        end
        as, bs = get_numquads(terms, un, nH; p1=bilin[1], p2=bilin[2], valp1=valp1, valp2=valp2)
        tot = binomial(length(terms),nH-2)
        tott += tot
        myEoN = get_EoNfunc(model; p1=bilin[1], p2=bilin[2], valp1=valp1, valp2=valp2)
        if length(un) <= sample_n_gt && full_solution == false
            # Save all ARs for a specific model
            nothing
        elseif length(un) <= sample_n_gt && full_solution == true
            @time begin
                proc_rs = similar(as, tot)
                EoN_rs = similar(bs, tot)
                rs_ws = similar(multis, length(proc_rs))
                parallel_alleqn_solve_proc_fullsol!(proc_rs, EoN_rs, rs_ws, as, bs, multis, tot, myEoN)
                save_full(model, proc_rs, EoN_rs, rs_ws, 1; folder=dataset*"/n"*string(nH)*"/", bilin=bilin, valp1=valp1, valp2=valp2, ms=ms[k])
            end
        elseif length(un) > sample_n_gt && full_solution == true
            nothing
        else
            nothing
        end
    end
    @info ""
    @info "Finished computing $model !\n"
    @printf "I had to compute %.3E models!\n" tott
    @info ""
end
