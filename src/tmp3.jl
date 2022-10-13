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




function runDFSZ_saveall(dataset, model_vec; model_multiplicity=ones(length(model_vec)), log_calculate_all_models_if_below=8, sample_models_factor=0.01)
    @time for (k, model) in enumerate(model_vec)
        println(model)
        bilins = get_bilins(model, same_χH_one_theory=false)
        myEoN = get_EoNfunc(model)
        myN = get_Nfunc(model)
        quads, multis = get_quads(model)
        nD = length(unique(model))
        tot_min = binomial(length(quads), nD-2)
        @info "Total number of models to be calculated above $(length(bilins)*tot_min)"
        
        @time for (i,bilin) in enumerate(bilins)
            #if tot_min <= 10^log_calculate_all_models_if_below
            terms = vcat(quads, bilins[i+1:end])
            mult = vcat([1,1] , multis, ones(Int64, length(bilins[i+1:end]))) #this changes relative bilin to quadrilinear probability
            #else
            #    @info "Sampling route!"
            #    terms = vcat(quads, bilins[1:end .!= i])
            #    mult = vcat([1,1], multis, ones(Int64, length(bilins)-1)) #this changes relative bilin to quadrilinear probability
            #end
            terms = vcat(orthogonality(model), bilin, terms)
            as, bs = get_numquads(terms, unique(model), nD)
            tot = binomial(length(terms)-2,nD-2)
            
            if tot_min <= 10^log_calculate_all_models_if_below
                proc_rs = similar(as, tot)
                EoN_rs = similar(bs, tot)
                rs_ws = similar(mult, length(proc_rs))
                myterms = fill([u1 for i in 1:nD],tot)
                sols = fill([0 for i in 1:nD],tot)
                parallel_alleqn_solve_proc_fullsol!(proc_rs, EoN_rs, rs_ws, as, bs, mult, terms, myterms, sols, tot, myEoN, myN)
            else
                factor=sample_models_factor
                chunk = Int(round(factor * tot))
                #chunk = 10^(log_sample_models-1)
                @info "Too many models to be calculated. Proceeding with calculation in chunks."
                @info "I will compute $(chunk*1) models in total."
                for l in 1:1
                    @info "Calculating chuck $l of 1"
                    proc_rs = similar(as, chunk)
                    EoN_rs = similar(bs, chunk)
                    rs_ws = similar(mult, length(proc_rs))
                    myterms = fill([u1 for i in 1:nD],tot)
                    sols = fill([0 for i in 1:nD],tot)
                    @time parallel_randeqn_solve_proc_fullsol!(proc_rs, EoN_rs, rs_ws, as, bs, mult, terms, myterms, sols, tot, myEoN, myN)
                end
            end

            stringterms = [string.(term) for term in myterms]
            stringsols = ["=" .*string.(sol) .*"s" for sol in sols]
            stringsols = permutedims(hcat(stringsols...))
            stringterms = permutedims(hcat(stringterms...))
            stringterms .*= stringsols
            #println(stringterms)

            # save_full somehow is slow if not precompiled! Let runDFSZ_saveall run first for e.g. n=3, which is inexpensive and later for higher ns. Somehow this makes a difference!
            save_full(model, proc_rs, EoN_rs, rs_ws, stringterms, 1; folder=dataset, bilin=bilin, ms=mult, model_multiplicity=model_multiplicity[k], full=(tot_min <= 10^log_calculate_all_models_if_below))
        end
    end
end

models, mmultis = model_list(nD=[6])
folder = "220916-serioustxtfiles"
#cmlists = runDFSZ_saveall(folder, models; model_multiplicity=mmultis, log_calculate_all_models_if_below=2, sample_models_factor=0.01)
h5totxt("samples_n6"; folder="./data/DFSZ_models/"*folder*"/")