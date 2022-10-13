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
        EoN_countmaps[k] = crs#cEoN
        save_EoN(model, cEoN; folder=dataset, filename=filename)
        #end

    end
    return EoN_countmaps
end

# Suitable sample_models_factor for not too extensive runtimes (Order of hours)
# n=7 -> 0.001
# n=8 -> 0.000001
# n=9 -> 0.000000001
models, multis = model_list(nD=[9])
models
for x in parse.(Int64, ARGS)
    @info "Computing model $x of $ARGS ."
    dataset = "220824-seriousruns"
    filename = "9EoNs_full_nomulti_sample"
    println(filename)
    cmlists = runDFSZ(dataset, [models[x]]; model_multiplicity=[multis[x]], log_calculate_all_models_if_below=2, sample_models_factor=0.000000001, same_χH_one_theory=true, filename=filename)
end

=
cmlists = runDFSZ(dataset, models; model_multiplicity=multis, log_calculate_all_models_if_below=2, log_sample_models=6, same_χH_one_theory=false, filename="6EoNs_samplesample_lowmulti")
cmlists = runDFSZ(dataset, models; model_multiplicity=multis, log_calculate_all_models_if_below=20, log_sample_models=7, same_χH_one_theory=true, filename="6EoNs_full_highmulti")
cmlists = runDFSZ(dataset, models; model_multiplicity=multis, log_calculate_all_models_if_below=2, log_sample_models=7, same_χH_one_theory=true, filename="6EoNs_sample_highmulti")
cmlists[1]
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

function read_EoN(dataset, models; specifier="EoNs")
    totEoN = countmap([])
    for model in models
        tpath = "./data/DFSZ_models/"*dataset
        name = model2string(model)
        nD = length(unique(model))

        e1 = FileIO.load(tpath*"/"*specifier*".jld2", string(nD)*"/"*name)
        totEoN= mergewith(+,totEoN,e1)
    end
    #v = collect(values(totEoN)) ./ sum(values(totEoN))
    #totEoN = countmap(collect(keys(totEoN)), v)
    return totEoN
end

crs1 = read_EoN(dataset, models; specifier="6EoNs_full_nomultitest")
k1 = collect(keys(crs1))
p1 = sort(collect(values(crs1))[sum.([isnan.(k) for k in k1]) .== 0], rev=true)
plot(p1, xscale=:log10, yscale=:log10)
crs2 = read_EoN(dataset, models; specifier="6EoNs_full_nomulti")
k2 = collect(keys(crs2))
tt = k2[isnan.(k2) .== 0]
tt[argmax(abs.(tt .- 1.92))] .* 3
EoNs[argmax(abs.(EoNs .- 1.92))] * 3
p2 = sort(collect(values(crs2))[sum.([isnan.(k) for k in k2]) .== 0], rev=true)
plot!(p2, xscale=:log10, yscale=:log10)
sum(p1)
unique(p1) ./ p1[1]
p1
sum(p2)
p2[35:60]
unique(p2) ./ p2[1]
histogram(p1, bins=0:1:100)
a1 = fit(Histogram, p1, nbins=maximum(p1))
a1.edges
a1.weights
p1
scatter(a1.edges[1][1:end-1], a1.weights ./ sum(p1), xscale=:log10, yaxis=(:log10, [0.00000001, :auto]))
a2 = fit(Histogram, p2, nbins=maximum(p2))
scatter!(a2.edges[1][1:end-1], a2.weights ./ sum(p2), xscale=:log10, yaxis=(:log10, [0.00000001, :auto]))

x1 = 0:1/length(p1):1
plot(x1[2:end],p1, xscale=:log10, yscale=:log10)

x2 = x1[2]:(1-x1[2])/length(p2):1
plot!(x2[2:end], p2, xscale=:log10, yscale=:log10)

e1=read_EoN("220824-seriousruns", models; specifier="4EoNs_full_nomulti")
e2=read_EoN("220816_sample_smalltest", models; specifier="5EoNs_sample_001")
e3=read_EoN("test", models; specifier="6EoNs_sample_01_modified")
e3=read_EoN("test", models; specifier="6EoNs_sample_05_modified")
e5=read_EoN("220816_sample_smalltest", models; specifier="6EoNs_newsample_highmultitest")
sve5 = sum(values(e5))
e1
totEoN = mergewith(-, e1, e2)
maximum(collect(values(totEoN)) ./ sum(abs.(values(totEoN))))
h1[1].edges[1][1:end-1]
h1[1].weights
h1 = _make_hist(e1)
plot()
scatter(h1[1].edges[1][1:end-1], aa, yaxis=(:log10, [0.000001,:auto]))
h2 = _make_hist(e2)
scatter!(h2[1].edges[1][1:end-1], n2weights, yaxis=(:log10, [0.0001,:auto]), xaxis=(:identity, [0.5,3]))



aa = h1[1].weights ./ sum(h1[1].weights)
minimum(aa[h1[1].weights .!= 0])
h2[1].weights
h1[1].weights

myweights = deepcopy(h2[1].weights)
normweights = myweights ./ sum(myweights)
normweights[normweights .!= 0.0] .+= 0.0#0025 #.+= minimum(aa[h1[1].weights .!= 0])
n2weights = normweights ./ sum(normweights)
plot(aa .- n2weights)

plot_hist!(h2, yaxis=(:log10, [0.000000001,:auto]))
h3 = _make_hist(e3)
plot_hist!(h3, yaxis=(:identity, [0.000000001,:auto]))
h4 = _make_hist(e4)
plot_hist!(h4, yaxis=(:identity, [0.000000001,:auto]))
h5 = _make_hist(e5)
plot_hist!(h5, yaxis=(:identity, [0.000000001,:auto]))

h4[1].weights[h4[1].weights .!= 0]
h2[1].weights[h2[1].weights .!= 0]
diff = h3[1].weights - h2[1].weights
plot(diff)
diff[diff .!= 0]
sum(diff)
sum(abs.(diff)) ./ length(diff[diff .!= 0])
n = sve5 / length(e5)
1/sqrt(n)

xlims!((0,2.5))
#=
=#


"""
    runDFSZ_saveall(dataset, model_vec; log_calculate_all_models_if_below=8, log_sample_models=7)

Calculate full solutions to all LES and store them in H5 files so potential terms, charges and anomaly ratios can be related to one another. This function is slow and takes vast storage space. Do not ever use this for computing EoN distributions above nD=5!

# Arguments:
See runDFSZ().
"""
function runDFSZ_saveall(dataset, model_vec; model_multiplicity=ones(length(model_vec)), log_calculate_all_models_if_below=8, sample_models_factor=0.01)
    @time for (k, model) in enumerate(model_vec)
        println(model)
        bilins = get_bilins(model)
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

models, mmultis = model_list(nD=[5])
cmlists = runDFSZ_saveall("221007-seriousnewtxtfiles", models; model_multiplicity=mmultis, log_calculate_all_models_if_below=20, sample_models_factor=0.001)

@time h5totxt("full_n4", folder="./data/DFSZ_models/221007-seriousnewtxtfiles/")

models, mmultis = model_list(nD=[4])
runDFSZ("220824-seriousruns", models; model_multiplicity=mmultis, log_calculate_all_models_if_below=20, sample_models_factor=0.001, NDW1=true, filename="6EoNs_nomulti_NDW1")

e1=read_EoN("220824-seriousruns", models; specifier="EoNs_nomulti_NDW1")

ke1 = collect(keys(e1))
ke1[isnan.(ke1) .== 0]
ve1 = collect(values(e1))
sum(ve1[isnan.(ke1) .== 0])/sum(ve1)

plot()
plot_EoN(e1)


fid = h5open("./data/DFSZ_models/221007-seriousnewtxtfiles/full_n4.h5")

#read(fid["n6_u1_u2_u3_d1_d1_d3_l1_l1_l1"]["bl=-d1 - u1"]["EoN"])

#EoNs = []
mychi = Matrix{Float64}(undef,0,10) 
#Ns = []

for (i, a) in enumerate(fid)
    k = read(a)
    if i != 1
        println(keys(fid)[i])
        for tuple in k#collect(values(k))
            dat = tuple[2]
            for j in 1:dat["model_multiplicity"]
                mychi = vcat(mychi, round.(dat["Chis"], digits=5))
            end

            #append!(EoNs, round.(dat["EoN"], digits=5))
            #append!(Ns, round.(dat["N"], digits=5))
                #write(io, model*"   "*mmult*"      "*rpad(dat["multis"][j],3)*"                 "*lpad(round(dat["EoN"][j], digits=3),8)*" "*lpad(Int(round(dat["E"][j])),5)*lpad(Int(round(dat["N"][j])),5)*(*(lpad.(round.(dat["Chis"][j,:],digits=3),9)...))*" "*(*(lpad.(filter.(x -> !isspace(x), dat["terms"][j,:]),18)...))*"\n")
        end
    end
end

tt4 = collect(eachrow(mychi))#unique(eachrow(mychi))
aEoN = [EoverN(tt[1:9]...) for tt in tt4]
aN = [N(tt[1:9]...) for tt in tt4]
aEoN[0.499 .<= abs.(aN) .<= 0.501]

close(fid)
=#
