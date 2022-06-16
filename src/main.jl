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
function runDFSZ(dataset;sample_n_gt=6, sample_log_nr_mods=9, ns =:all, compute_equivalent_theories=false, exactly_one_bilinear=false, bilin_weights::Integer=1, full_solution=false)

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
                @time begin
                proc_rs = similar(bs, tot)
                rs_ws = similar(multis, length(proc_rs))
                parallel_alleqn_solve_proc!(proc_rs, rs_ws, as, bs, multis, tot, myEoN)
                save_AR(model, proc_rs, rs_ws, 1; folder=dataset*"/n"*string(nH)*"/", bilin=bilin, ms=ms[k])
                end
            elseif length(un) <= sample_n_gt && full_solution == true
                @time begin
                    proc_rs = similar(as, tot)
                    EoN_rs = similar(bs, tot)
                    rs_ws = similar(multis, length(proc_rs))
                    parallel_alleqn_solve_proc_fullsol!(proc_rs, EoN_rs, rs_ws, as, bs, multis, tot, myEoN)
                    save_full(model, proc_rs, EoN_rs, rs_ws, 1; folder=dataset*"/n"*string(nH)*"/", bilin=bilin, valp1=valp1, valp2=valp2, ms=ms[k])
                end
            else
                chunk = 10^(sample_log_nr_mods - 1)
                m=10
                @time for i in 1:m
                    @info "Computing round $i of $m"
                    proc_rs = similar(bs, chunk)
                    rs_ws = similar(multis, length(proc_rs))
                    parallel_randeqn_solve_proc!(proc_rs, rs_ws, as, bs, multis, tot, myEoN)
                    save_AR(model, proc_rs, rs_ws, i; folder=dataset*"/n"*string(nH)*"/", bilin=bilin, ms=ms[k])
                end
            end
        end
        @info ""
        @info "Finished computing $model !\n"
        @printf "I had to compute %.3E models!\n" tott
        @info ""
    end
end


@time runDFSZ("220616-nbilin_fullsol"; sample_n_gt=6, sample_log_nr_mods=9, ns=[3,4,5], exactly_one_bilinear=false, compute_equivalent_theories=true, full_solution=true)

#fid = h5open("./data/DFSZ_models/220607-fullsolutionsv1/n6/full_n6.h5")
#cc = read(fid["3n6_u1_u2_u3_d1_d1_d2_l1_l1_l1"]["bl=u1-u2"]["N"])
#close(fid)

@info ""
@info ""
@info "Start plotting!!!"
@info ""
@info ""




# This would be the procedure to actually calculate really all models!
#=
for model in with_replacement_combinations([u1, d1, l1],9)
    un = unique(model)
    nH = length(un)
    if nH < 3
        nothing
    else
        quads, multis = get_quads(model)
        tot = binomial(length(quads),nH-2)
        @printf "Your model-group has %.3E models \n" tot
        @time as, bs = get_numquads(quads, un, nH)
        myEoN = get_EoNfunc(model)


        # Save all ARs for a specific model
        proc_rs = similar(bs, tot)
        rs_ws = similar(multis, length(proc_rs))
        @time parallel_alleqn_solve_proc!(proc_rs, rs_ws, as, bs, multis, tot)
        @time save_AR(model, proc_rs, rs_ws, 1; folder="n"*string(nH)*"/")
    end
end
=#


# Plot data
#=
for k in 5:5
    @info "$k"
    savefolder = "./plots/"*dataset*"/ARs/n"*string(k)
    mkpath(savefolder)
    folders = readdir("./data/DFSZ_models/"*dataset*"/n"*string(k); join=true)
    #hist_list = similar(files, Any)
    for folder in folders
        files = readdir(folder)
        @time for (i, file) in enumerate(files)
            println(file)
            model = fname2model(file)
            bilin = fname2bilin(file)
            tt = read_AR(model; folder=dataset*"/n"*string(k), bilin=bilin)
            tt = normalize(tt; mode=:probability)
            p1 = plot(tt, lt=:stepbins, label="", title="$(file[1:end-5])", 
                xlabel="E/N", ylabel="Probability", xrange=(5/3-9,5/3+9),
                bottom_margin=2Plots.mm, legend=:topright,
                size=(400,300), lw=2)
            savefig(p1, savefolder*"/$(file[1:end-5])_ARs.pdf")
            #gagh = gag_histogram(tt; mode=:probability)
            #hist_list[i] = gagh
        end
    end
    #gaghs[k-2] = merge(hist_list...)
end
=#

#=
@time for k in 5:5
    @info "$k"
    savefolder = "./plots/"*dataset*"/ARs/n"*string(k)
    mkpath(savefolder)
    folders = readdir("./data/DFSZ_models/"*dataset*"/n"*string(k); join=true)
    
    tttot_list = similar(folders, Any)
    for (j, folder) in enumerate(folders)
        m = fname2m(folder)
        files = readdir(folder)
        hist_list = similar(files, Any)
        @time for (i, file) in enumerate(files)
            model = fname2model(file)
            bilin = fname2bilin(file)
            tt = read_AR(model; folder=dataset*"/n"*string(k), bilin=bilin, m=m)
            #tt = normalize(tt; mode=:probability)
            #gagh = gag_histogram(tt; mode=:probability)
            hist_list[i] = tt
        end
        tttot = hist_list[3]#merge(hist_list...)
        #tttot = normalize(tttot; mode=:probability)
        #tttot.weights .*= parse(Int,string(m))
        tttot_list[j] = tttot
        type = split(folder, "/")[end]
        println(type)
        p2 = plot(tttot, lt=:stepbins, label="", title="$(type)", 
        xlabel="E/N", ylabel="Probability", xrange=(5/3-9,5/3+9),
        bottom_margin=2Plots.mm, legend=:topright,
        size=(400,300), lw=2)
        #savefig(p2, savefolder*"/full_$(type)_ARs.pdf")
    end
    global ttall=tttot_list[6]#ttall = merge(tttot_list...)
    #ttall = normalize(ttall; mode=:probability)
    p1 = plot(ttall, lt=:stepbins, label="", title="n$(k)", 
    xlabel="E/N", ylabel="Probability", xrange=(5/3-9,5/3+9),
    bottom_margin=2Plots.mm, legend=:topright,
    size=(400,300), lw=2)
    #savefig(p1, savefolder*"/full_n$(k)_ARs.pdf")
    #gaghs[k-2] = merge(hist_list...)
end
=#

#ttall.edges[1][1:end-1][ttall.weights .> 0]
#ttall.weights[ttall.weights .> 0]