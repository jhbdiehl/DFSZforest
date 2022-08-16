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
    model_list([, nH=:all, compute_equivalent_theories=false])    

Construct models that are to be analysed in a manner that the run function can read. If you want to relax specific conditions (like considering models such as [u1,u2,u3,u1,u1,u1,u1,u1,u1]) you will have to change the generating functions `generate_all_models()` or `generate_unique_models()`.

# Arguments
- `nH`: Number of unique higgs particles in your model. Use either `:all` or list of integers. `:all` is equal to `[3,4,5,6,7,8,9]`.
- `compute_equivalent_theories::Bool`: If `false` e.g. [u1,u1,u3,d1,d1,d1,l1,l1,l1] and [u1,u3,u1,d1,d1,d1,l1,l1,l1] will be treated separately. If `true` only one will be calculated and a note will be made that this solution happens three times (in this example).

# Returns
- `model_vec::Vector{Vector{Num}}`: Vector containing symbolic representations of all models that you set up for calculation.
- `model_multiplicity::Vector{Any}`: Vector of integers tracking how often specific models can arise. Equal to `ones(length(model_vec))`, if `compute_equivalent_theories == true`.

"""
function model_list(;nH=:all, compute_equivalent_theories=false)
    if compute_equivalent_theories
        a, ms = generate_all_models()
    else
        a, ms = generate_unique_models()
    end

    if nH == :all
        model_vec = a
    else
        model_vec = a[length.(unique.(a)) .∈ Ref(nH)]
        model_multiplicity = ms[length.(unique.(a)) .∈ Ref(nH)]
    end
    return model_vec, model_multiplicity
end

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
- `log_calculate_all_models_if_below::Int`: Calculate all possible models, if for only one bilinear the number of possible combinations for the potential is below 10^ this value. To calculate nH=6 fully you shold set this to 8 or above. Calculating nH=5 shold take order of seconds, nH=6 order of minutes, can't recommend to go to nH=7...
- `log_sample_models::Int`: When sampling, get `10^log_sample_models` this many samples for one specific bilinear. (I.e. if a specific model has 9 possible bilinears, `9 × 10^log_sample_models` are calculated.)
- `same_χH_one_theory::Bool`: When true, multiplicities down to the level of solutions for the charges are ignored. Aka, if the solution is the same, then we can add up the terms to form a new, unique potential that gives this solution.

# Returns
- `EoN_countmaps::Vector{Dicts}`: For each model in model_vec a countmap of (rounded) anomaly ratios is provided.
"""
function runDFSZ(dataset, model_vec; model_multiplicity=ones(length(model_vec)), log_calculate_all_models_if_below=8, log_sample_models=7, same_χH_one_theory=true)
    EoN_countmaps = Array{Any}(undef, length(model_vec))
    @time for (k, model) in enumerate(model_vec)
        println(model)
        bilins = get_bilins(model)
        myEoN = get_EoNfunc(model)
        myN = get_Nfunc(model)
        quads, multis = get_quads(model)
        nH = length(unique(model))
        rslist = Array{Any}(undef, length(bilins))
        tot_min = binomial(length(quads), nH-2)
        @info "Total number of models to be calculated above $(length(bilins)*tot_min)"
        
        @time for (i,bilin) in enumerate(bilins)
            if tot_min <= 10^log_calculate_all_models_if_below
                terms = vcat(quads, bilins[i+1:end])
                mult = vcat([1,1], multis, ones(Int64, length(bilins[i+1:end]))) #this changes relative bilin to quadrilinear probability
            else
                @info "Sampling route!"
                terms = vcat(quads, bilins[1:end .!= i])
                mult = vcat([1,1], multis, ones(Int64, length(bilins)-1)) #this changes relative bilin to quadrilinear probability
            end
            terms = vcat(orthogonality(model), bilin, terms)
            as, bs = get_numquads(terms, unique(model), nH)
            tot = binomial(length(terms)-2,nH-2)
            println(tot)
            
            if tot_min <= 10^log_calculate_all_models_if_below
                proc_rs = similar(as, tot)
                EoN_rs = similar(bs, tot)
                rs_ws = similar(mult, length(proc_rs))
                parallel_alleqn_solve_proc!(proc_rs, rs_ws, as, bs, mult, tot)
                rslist[i] = countmap(proc_rs, rs_ws)
            else
                chunk = 10^(log_sample_models-1)
                rslist[i] = countmap([])
                @info "Too many models to be calculated. Proceeding with calculation in chunks."
                for l in 1:10
                    @info "Calculating chuck $l of 10"
                    proc_rs = similar(as, chunk)
                    EoN_rs = similar(bs, chunk)
                    rs_ws = similar(mult, length(proc_rs))
                    @time parallel_randeqn_solve_proc!(proc_rs, rs_ws, as, bs, mult, tot)
                    dummy = countmap(proc_rs, rs_ws)
                    rslist[i] = mergewith(+, rslist[i], dummy)
                end
            end
            
        end
        Threads.@threads for i in 1:1 # This loop is necessary, otherwise usage of myN generated function from above leads to world age error when multithreading
            crsi = mergewith(+,rslist...)
            crs = collect(keys(crsi))
            crs = convert(Vector{Vector{Float64}}, crs)
            crs = permutedims(hcat(crs...))
            crs[findall(iszero,crs)] = crs[findall(iszero,crs)].^2 # -0.0 and 0.0 is the same thing!
            crs = collect(eachrow(crs))
            crs = countmap(crs, collect(values(crsi)))
            EoNf = _EoNft(crs, myN, myEoN)

            if same_χH_one_theory
                cEoN = countmap(round.(EoNf, digits=6), ones(length(collect(values(crs)))) .* model_multiplicity[k])
            else
                cEoN = countmap(round.(EoNf, digits=6), collect(values(crs)) .* model_multiplicity[k])
            end
            
            clean_countmap!(cEoN)
            EoN_countmaps[k] = cEoN
            save_EoN(model, cEoN; folder=dataset)
        end

    end
    return EoN_countmaps
end



models, multis = model_list(nH=[5])
cmlists = runDFSZ("test", models; model_multiplicity=multis, log_calculate_all_models_if_below=20, log_sample_models=6, same_χH_one_theory=false)

# Check if countmap is symmetric around 5/3. If not something went wrong
cEoN = mergewith((x,y)->x+y, cmlists...)
clean_countmap!(cEoN)
tt = check_symmetry(Dict{Float64,Float64}(cEoN)) #if this does not say true you are in trouble
sum(abs.(values(tt)))/2
sum(values(totEoN))
# Plot total
h1 = _make_hist(cEoN)
plot()
plot_hist!(h1, yaxis=(:identity, [0.000000001,:auto]))

function plot_EoN(EoN_countmap)
    h1 = _make_hist(EoN_countmap)
    plot()
    plot_hist!(h1, yaxis=(:identity, [0.000000001,:auto]))
end

function read_EoN(dataset, model;specifier="")
    tpath = "./data/DFSZ_models/"*dataset
    name = model2string(model)
    nH = length(unique(model))

    FileIO.load(tpath*"/EoNs.jld2", string(nH)*"/"*name)
end

totEoN = countmap([])
for model in models
    e1=read_EoN("test", model; specifier="")
    totEoN = mergewith(+,totEoN,e1)
end
plot_EoN(totEoN)

#=
=#


"""
    runDFSZ_saveall(dataset, model_vec; log_calculate_all_models_if_below=8, log_sample_models=7)

Calculate full solutions to all LES and store them in H5 files so potential terms, charges and anomaly ratios can be related to one another. This function is slow and takes vast storage space. Do not ever use this for computing EoN distributions above nH=5!

# Arguments:
See runDFSZ().
"""
function runDFSZ_saveall(dataset, model_vec; log_calculate_all_models_if_below=8, log_sample_models=7)
    @time for (k, model) in enumerate(model_vec)
        println(model)
        bilins = get_bilins(model)
        myEoN = get_EoNfunc(model)
        myN = get_Nfunc(model)
        quads, multis = get_quads(model)
        nH = length(unique(model))
        tot_min = binomial(length(quads), nH-2)
        @info "Total number of models to be calculated above $(length(bilins)*tot_min)"
        
        @time for (i,bilin) in enumerate(bilins)
            if tot_min <= 10^log_calculate_all_models_if_below
                terms = vcat(quads, bilins[i+1:end])
                mult = vcat([1,1], multis, ones(Int64, length(bilins[i+1:end]))) #this changes relative bilin to quadrilinear probability
            else
                @info "Sampling route!"
                terms = vcat(quads, bilins[1:end .!= i])
                mult = vcat([1,1], multis, ones(Int64, length(bilins)-1)) #this changes relative bilin to quadrilinear probability
            end
            terms = vcat(orthogonality(model), bilin, terms)
            as, bs = get_numquads(terms, unique(model), nH)
            tot = binomial(length(terms)-2,nH-2)
            
            if tot_min <= 10^log_calculate_all_models_if_below
                proc_rs = similar(as, tot)
                EoN_rs = similar(bs, tot)
                rs_ws = similar(mult, length(proc_rs))
                myterms = fill([u1 for i in 1:nH],tot)
                parallel_alleqn_solve_proc_fullsol!(proc_rs, EoN_rs, rs_ws, as, bs, mult, terms, myterms, tot, myEoN, myN)
            else
                chunk = 10^(log_sample_models-1)
                @info "Too many models to be calculated. Proceeding with calculation of only a part of the models."
                @info "Calculating chuck 1 of 1"
                proc_rs = similar(as, chunk)
                EoN_rs = similar(bs, chunk)
                rs_ws = similar(mult, length(proc_rs))
                myterms = fill([u1 for i in 1:nH],tot)
                @time parallel_randeqn_solve_proc_fullsol!(proc_rs, EoN_rs, rs_ws, as, bs, mult, terms, myterms, tot, myEoN, myN)
            end
            stringterms = [string.(term) for term in myterms]
            stringterms = permutedims(hcat(stringterms...))

            save_full(model, proc_rs, EoN_rs, rs_ws, stringterms, 1; folder=dataset, bilin=bilin, ms=mult)
            
        end
    end
end


models, multis = model_list(nH=[4])
cmlists = runDFSZ_saveall("test", models; model_multiplicity=multis, log_calculate_all_models_if_below=20, log_sample_models=6, same_χH_one_theory=false)

fid = h5open("./data/DFSZ_models/test/full_n4.h5")

read(fid["n4_u1_u1_u1_d1_d1_d1_l1_l1_l3"]["bl=-d1 - u1"]["terms"])