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



"""
    _EoNft(countmap_χH_results, function_calculating_N, function_calculating_Anomaly_Ratio)

Converts charges to anomaly ratio and cleans up a bit (if N ≈ 0, then E/N is either very large, or completely undefined, if also E ≈ 0).
"""
function _EoNft(crs, myN, myEoN)
    crsN = myN.(collect(keys(crs)))
    crsEoN = myEoN.(collect(keys(crs)))
    EoNf = ifelse.(-0.000000001 .< crsN .< 0.000000001, NaN, crsEoN)
end

"""
    runDFSZ(dataset, model_vec; model_multiplicity=ones(length(model_vec)), log_calculate_all_models_if_below=8, log_sample_models=7, same_χH_one_theory=true)

Calculate only aggregate anomaly ratios of specific models and save them into a JLD2 file in the folder dataset. This function cannot be used to reconstruct possible Higgs charges belonging to a specific anomaly ratio. The probability of having a specific quadrilinear term and a specific bilinear term in the potential are assumed to be equally likely. 

# Arguments:
- `dataset::String`: Folder name to save data to.
- `model_vec::Vector{Vector{Num}}`: Vector containing symbolic representations of all models.
- `model_multiplicity::Vector{Any}`: Vector of integers tracking how often specific models can arise.
- `log_calculate_all_models_if_below::Int`: Calculate all possible models, if for only one bilinear the number of possible combinations for the potential is below 10^ this value. To calculate nD=6 fully you shold set this to 8 or above. Calculating nD=5 shold take order of seconds, nD=6 order of minutes, can't recommend to go to nD=7...
- `sample_models_factor::Float`: When sampling, get `sample_models_factor * total_#_of_models` samples for each specific bilinear.
- `same_χH_one_theory::Bool`: When true, multiplicities down to the level of solutions for the charges are ignored. Aka, if the solution is the same, then we can add up the terms to form a new, unique potential that gives this solution.

# Returns
- `EoN_countmaps::Vector{Dicts}`: For each model in model_vec a countmap of (rounded) anomaly ratios is provided.
"""
function runDFSZ(dataset, model_vec; model_multiplicity=ones(length(model_vec)), log_calculate_all_models_if_below=8, sample_models_factor=0.01, same_χH_one_theory=true, NDW1=false, filename="EoNs")
    EoN_countmaps = Array{Any}(undef, length(model_vec))
    @time for (k, model) in enumerate(model_vec)
        println(model)
        bilins = get_bilins(model; same_χH_one_theory=false)
        myEoN = get_EoNfunc(model)
        myN = get_Nfunc(model)
        quads, multis = get_quads(model)
        nD = length(unique(model))
        rslist = Array{Any}(undef, length(bilins))
        tot_min = binomial(length(quads), nD-2)
        @info "Total number of models to be calculated above $(length(bilins)*tot_min)"
        
        @time for (i,bilin) in enumerate(bilins)
            #if tot_min <= 10^log_calculate_all_models_if_below
            if same_χH_one_theory
                terms = quads # since all quads can be constructed out of two bilinears, if you dont care about multiplicities, you can assume only one bilinear wlg
                mult = vcat([1,1], multis)
            else
                terms = vcat(quads, bilins[i+1:end])
                mult = vcat([1,1], multis, ones(Int64, length(bilins[i+1:end]))) #this changes relative bilin to quadrilinear probability
            end
            #else
            #    @info "Sampling route!"
            #    terms = vcat(quads, bilins[1:end .!= i])
            #    mult = vcat([1,1], multis, ones(Int64, length(bilins)-1)) #this changes relative bilin to quadrilinear probability
            #end
            terms = vcat(orthogonality(model), bilin, terms)
            as, bs = get_numquads(terms, unique(model), nD)
            tot = binomial(length(terms)-2,nD-2)
            println(tot)
            
            if tot_min <= 10^log_calculate_all_models_if_below
                proc_rs = similar(as, tot)
                EoN_rs = similar(bs, tot)
                rs_ws = similar(mult, length(proc_rs))
                parallel_alleqn_solve_proc!(proc_rs, rs_ws, as, bs, mult, tot)
                rslist[i] = countmap(proc_rs, rs_ws)
            else
                factor= sample_models_factor
                chunk = Int(round(factor * tot))
                #chunk = 10^(log_sample_models-1)
                rslist[i] = countmap([])
                @info "Too many models to be calculated. Proceeding with calculation in chunks."
                @info "I will compute $(chunk*1) models in total. (model $k, bilin $i)"
                for l in 1:1
                    @info "Calculating chuck $l of 1"
                    proc_rs = similar(as, chunk)
                    EoN_rs = similar(bs, chunk)
                    rs_ws = similar(mult, length(proc_rs))
                    @time parallel_randeqn_solve_proc!(proc_rs, rs_ws, as, bs, mult, tot)
                    dummy = countmap(proc_rs, rs_ws)
                    rslist[i] = mergewith(+, rslist[i], dummy)
                end
            end
            
        end

        crs = mergewith(+,rslist...)
        

        # This adds all of the negative solutions from hermitian conjugated bilinears!
        if same_χH_one_theory
            crsk = collect(keys(crs))
            #crskm = -1 .* crsk
            #crskmgood0 = replace.(crskm, -0.0 => 0.0)
            crskgood0 = replace.(crsk, -0.0 => 0.0)
            #crsm = countmap(crskmgood0, collect(values(crs)))
            crs = countmap(crskgood0, collect(values(crs)))
            #crs = mergewith(+, crs, crsm)
        else
            crsk = collect(keys(crs))
            crskgood0 = replace.(crsk, -0.0 => 0.0)
            crs = countmap(crskgood0, collect(values(crs)))
        end
        #crs = collect(keys(crsi2))


        #crs = convert(Vector{Vector{Float64}}, crs)
        #@time alltmps = Vector{Vector{Float64}}()
        #@time for cr in crs
        #    tmp = similar(cr)
        #    tmp[findall(!iszero,cr)] = cr[findall(!iszero,cr)]
        #    tmp[findall(iszero,cr)] = cr[findall(iszero,cr)].^2
        #    alltmps = vcat(alltmps, [tmp])
        #end
        #crs = permutedims(hcat(crs...))
        #crs[findall(iszero,crs)] = crs[findall(iszero,crs)].^2 # -0.0 and 0.0 is the same thing!
        #crs = collect(eachrow(crs))
        #crs = countmap(crs, collect(values(crsi2)))

        #crs = mergewith(+,rslist...)

        #@time Threads.@threads for i in 1:1 # This loop is necessary, otherwise usage of myN generated function from above leads to world age error when multithreading
        #EoNf = _EoNft(crs, myN, myEoN)
        crsN = Base.invokelatest.(myN, collect(keys(crs)))
        crsEoN = Base.invokelatest.(myEoN, collect(keys(crs)))
        #EoNf = ifelse.(-0.000000001 .< crsN .< 0.000000001, NaN, crsEoN)
        if NDW1
            truN = similar(crsN)
            truN[isnan.(crsN)] .= NaN
            truN[isnan.(crsN) .== 0] .= numerator.(rationalize.(2 .* crsN[isnan.(crsN) .== 0], tol=0.0001)) ./ 2
            EoNf = ifelse.(0.49999 .< abs.(truN) .< 0.50001, crsEoN, NaN)
        else
            EoNf = ifelse.(-0.000000001 .< crsN .< 0.000000001, NaN, crsEoN)
        end

        if same_χH_one_theory
            cEoN = countmap(round.(EoNf, digits=6), ones(length(collect(values(crs)))) .* model_multiplicity[k]) #model_multiplicity[k] ./ collect(values(crs)))#
        else
            cEoN = countmap(round.(EoNf, digits=6), collect(values(crs)) .* model_multiplicity[k])
        end
        
        #clean_countmap!(cEoN)
        #EoN_countmaps[k] = crs#cEoN
        save_EoN(model, cEoN; folder=dataset, filename=filename)
        #end

    end
    #return EoN_countmaps
end

# Suitable sample_models_factor for not too extensive runtimes (Order of hours)
# n=7 -> 0.001
# n=8 -> 0.000001
# n=9 -> 0.000000001
models, multis = model_list(nD=[7])
for x in parse.(Int64, ARGS)
    @info "Computing model $x of $ARGS ."
    dataset = "test"
    filename = "7EoNs_nomulti_NDW1_1"
    println(filename)
    runDFSZ(dataset, [models[x]]; model_multiplicity=[multis[x]], log_calculate_all_models_if_below=20, sample_models_factor=0.01, same_χH_one_theory=true, NDW1=true, filename=filename)
end