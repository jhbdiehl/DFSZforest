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




models, multis = model_list(nD=[7])
dataset = "220824-seriousruns"
cmlists = runDFSZ(dataset, models; model_multiplicity=multis, log_calculate_all_models_if_below=20, sample_models_factor=1, same_χH_one_theory=true, filename="6EoNs_full_nomultitest")
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
    #plot()
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

crs1 = Dict{Float64, Int64}(crs1)
check_symmetry(crs1)
plot_EoN(crs1)
models, multis = model_list(nD=[7])
crs1 = read_EoN("220824-seriousruns", models; specifier="7EoNs_full_nomulti")
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

models
e1 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti.jld2", "7/n7_u1_u2_u3_d1_d1_d3_l1_l1_l3")
e1[0.0]
e1[-0.0]
e2 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti2.jld2", "7/n7_u1_u1_u3_d1_d2_d3_l1_l1_l3")
sum(values(e2))
sum(values(e1))
ee = mergewith(+, e1,e2)
plot_EoN(ee)
check_symmetry(ee)
e1[3.333333] - e2[0]
e1[0.0]
e1[-0.0]
e2[0]
e2[-0.0]

scatter(a1.edges[1][1:end-1], a1.weights ./ sum(p1), xscale=:log10, yaxis=(:log10, [0.00000001, :auto]))
a2 = fit(Histogram, p2, nbins=maximum(p2))
scatter!(a2.edges[1][1:end-1], a2.weights ./ sum(p2), xscale=:log10, yaxis=(:log10, [0.00000001, :auto]))

x1 = 0:1/length(p1):1
plot(x1[2:end],p1, xscale=:log10, yscale=:log10)

x2 = x1[2]:(1-x1[2])/length(p2):1
plot!(x2[2:end], p2, xscale=:log10, yscale=:log10)

models, multis = model_list(nD=[9])
e1=read_EoN("220824-seriousruns", models; specifier="9EoNs_full_nomulti_sample")
e2=read_EoN("220830-n5-tests", models; specifier="5EoNs_onebilinwithhc")
e3=read_EoN("220830-n5-tests", models; specifier="5EoNs_onebilinnohc")
plot()
plot_EoN(e3)
plot_EoN(e2)
plot_EoN(e1)
check_symmetry(Dict{Float64,Int64}(e1))
check_symmetry(Dict{Float64,Int64}(e2))
check_symmetry(Dict{Float64,Int64}(e3))

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

minimum(collect(keys(e1))[isnan.(collect(keys(e1))) .== 0])

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
function runDFSZ_saveall(dataset, model_vec; model_multiplicity=ones(length(model_vec)), log_calculate_all_models_if_below=8, sample_models_factor=0.01, same_χH_one_theory=true)
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
            println(tot)
            
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

models, mmultis = model_list(nD=[3])
cmlists = runDFSZ("220824-seriousruns", models; model_multiplicity=mmultis, log_calculate_all_models_if_below=20, sample_models_factor=0.001)

@time h5totxt("full_n3"; folder="./data/DFSZ_models/test3/")


fid = h5open("./data/DFSZ_models/221007-seriousnewtxtfiles/full_n4.h5")

#read(fid["n6_u1_u2_u3_d1_d1_d3_l1_l1_l1"]["bl=-d1 - u1"]["EoN"])

EoNs = []
mychi = Matrix{Float64}(undef,0,10)
tot = 0 

for (i, a) in enumerate(fid)
    k = read(a)
    if i != 1
        println(keys(fid)[i])
        for tuple in k#collect(values(k))
            dat = tuple[2]
            tot += sum(dat["model_multiplicity"] .* dat["multis"])
            mychi = vcat(mychi, round.(dat["Chis"], digits=5))
            append!(EoNs, round.(dat["EoN"], digits=5))
                #write(io, model*"   "*mmult*"      "*rpad(dat["multis"][j],3)*"                 "*lpad(round(dat["EoN"][j], digits=3),8)*" "*lpad(Int(round(dat["E"][j])),5)*lpad(Int(round(dat["N"][j])),5)*(*(lpad.(round.(dat["Chis"][j,:],digits=3),9)...))*" "*(*(lpad.(filter.(x -> !isspace(x), dat["terms"][j,:]),18)...))*"\n")
        end
    end
end
tot
mychi
EoNs
tt2 = unique(eachrow(mychi))

length(EoNs[abs.(EoNs .- 1.92) .< 0.04]) ./ length(EoNs)

tt3 = unique(EoNs)
minimum(tt3)
close(fid)
countmap(EoNs)
dd = Dict{Float64, Int64}(countmap(EoNs))


ddmbwh = get_dd(h5open("./data/DFSZ_models/220830-n5-tests/full_n5_multiplebilinswithhc.h5"))
dd1bwh = get_dd(h5open("./data/DFSZ_models/220830-n5-tests/full_n5_onebilinwithhc.h5"))
dd1bnh = get_dd(h5open("./data/DFSZ_models/220830-n5-tests/full_n5_onebilinnohc.h5"))

plot()
plot_EoN(ddmbwh)
plot_EoN(dd1bwh)
plot_EoN(dd1bnh)

check_symmetry(dd)

EoNs
reshape(EoNs, (8068,10))
EoNs[argmax(abs.(EoNs .- 1.92))] * 3

close(fid)


models, mmultis = model_list(nD=[6])
cmlists = howmany("test", models; model_multiplicity=mmultis, log_calculate_all_models_if_below=20, sample_models_factor=0.001)

0.0e1 + 1180524705248944737

function howmany(dataset, model_vec; model_multiplicity=ones(length(model_vec)), log_calculate_all_models_if_below=8, sample_models_factor=7, same_χH_one_theory=true, filename="EoNs")
    EoN_countmaps = Array{Any}(undef, length(model_vec))
    sumtot = 0
    @time for (k, model) in enumerate(model_vec)
        println(model)
        bilins = get_bilins(model)
        myEoN = get_EoNfunc(model)
        myN = get_Nfunc(model)
        quads, multis = get_quads(model)
        nD = length(unique(model))
        rslist = Array{Any}(undef, length(bilins))
        tot_min = binomial(length(quads), nD-2)
        #@info "Total number of models to be calculated above $(length(bilins)*tot_min)"
        
        @time for (i,bilin) in enumerate(bilins)
            #if tot_min <= 10^log_calculate_all_models_if_below
            terms = vcat(quads, bilins[i+1:end])
            #else
            #    @info "Sampling route!"
            #    terms = vcat(quads, bilins[1:end .!= i])
            #    mult = vcat([1,1], multis, ones(Int64, length(bilins)-1)) #this changes relative bilin to quadrilinear probability
            #end
            terms = vcat(orthogonality(model), bilin, terms)
            tot = binomial(length(terms)-2,nD-2)
            #println(tot)
            sumtot += tot
        end
    end
    println(sumtot)
end

models, mmultis = model_list(nD=[4], compute_equivalent_theories=true)
a = read_EoN("220824-seriousruns", models, specifier="4EoNs_full_nomulti_equivalents")

tt = collect(keys(a))
maximum(tt[isnan.(tt) .==0])






a,m = get_quads([u1,u1,u1,d1,d1,d1,l1,l1,l1])
m
aa = get_bilins([u1,u1,u1,d1,d1,d1,l1,l1,l1])



models, mmultis = model_list(nD=[5])
dataset="test"
model_vec=models
model_multiplicity=mmultis
log_calculate_all_models_if_below=8
sample_models_factor=7
same_χH_one_theory=true
filename="EoNs"
EoN_countmaps = Array{Any}(undef, length(model_vec))
@time for (k, model) in enumerate(model_vec)
    println(model)
    bilins = get_bilins(model)
    myEoN = get_EoNfunc(model)
    myN = get_Nfunc(model)
    quads, multis = get_quads(model)
    nD = length(unique(model))
    rslist = Array{Any}(undef, length(bilins))
    tot_min = binomial(length(quads), nD-2)
    @info "Total number of models to be calculated above $(length(bilins)*tot_min)"
    @time for (i,bilin) in enumerate(bilins)
        #if tot_min <= 10^log_calculate_all_models_if_below
        terms = vcat(quads, bilins[i+1:end])
        mult = vcat([1,1], multis, ones(Int64, length(bilins[i+1:end]))) #this changes relative bilin to quadrilinear probability
        terms = vcat(orthogonality(model), bilin, terms)
        as, bs = get_numquads(terms, unique(model), nD)
        t1 = deepcopy(as)
        for n in 1:length(unique(model))-2
            for sp in multiset_combinations(t1[3:end], n)
                if rank(permutedims(hcat(t1[1], t1[2], sp...))) != n+2
                    for j in n:length(sp)
                        t1 = t1[t1 .!= Ref(sp[j])]
                    end
                end
            end
        end
        println(length(t1))
    end
        #=
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
            @info "I will compute $(chunk*10) models in total."
            for l in 1:1
                @info "Calculating chuck $l of 10"
                proc_rs = similar(as, chunk)
                EoN_rs = similar(bs, chunk)
                rs_ws = similar(mult, length(proc_rs))
                @time parallel_randeqn_solve_proc!(proc_rs, rs_ws, as, bs, mult, tot)
                dummy = countmap(proc_rs, rs_ws)
                rslist[i] = mergewith(+, rslist[i], dummy)
            end
        end
        =#
        
end
    #=
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
            cEoN = countmap(round.(EoNf, digits=6), ones(length(collect(values(crs)))) .* model_multiplicity[k]) #model_multiplicity[k] ./ collect(values(crs)))#
        else
            cEoN = countmap(round.(EoNf, digits=6), collect(values(crs)) .* model_multiplicity[k])
        end

        clean_countmap!(cEoN)
        EoN_countmaps[k] = cEoN
        save_EoN(model, cEoN; folder=dataset, filename=filename)
    end
    =#
t1 = [[1, 1, 2, 3],[1, 2, 2, 2],[3 ,2, 3, 3],[2,2,4,6],[2,3,4,5],[4,7,2,5],[2,9,11,13]]
sy = multiset_combinations(t1,2)
for s in sy
    println(rank(permutedims(hcat(s...))))
end
rank(permutedims(hcat(sy...)))
for n in 2:2
    for sp in multiset_combinations(t1, n)
        if rank(permutedims(hcat(sp...))) != n
            for i in n:length(sp)
                t1 = t1[t1 .!= Ref(sp[i])]
            end
        end
    end
    println(t1)
end
length(collect(multiset_combinations(t1, 4)))

t = [[-1.0, -1.0, 1.0, 1.0], [0.0, 2.0, 0.0, 2.0], [2.0, -1.0, 0.0, 1.0], [1.0, -1.0, -1.0, 1.0], [0.0, 1.0, 2.0, -1.0], [0.0, 2.0, 1.0, 1.0], [-2.0, 2.0, 0.0, 0.0], [-1.0, 1.0, -1.0, 1.0], [2.0, 0.0, 2.0, 0.0], [2.0, -1.0, 1.0, 0.0], [0.0, 2.0, 2.0, 0.0], [1.0, 1.0, 2.0, 0.0], [0.0, 0.0, -2.0, 2.0], [-1.0, 2.0, 1.0, 0.0], [-1.0, 2.0, 0.0, 1.0], [2.0, 0.0, 0.0, 2.0], [2.0, 0.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0], [0.0, 1.0, -1.0, 2.0], [1.0, 0.0, -1.0, 2.0], [1.0, 1.0, 0.0, 2.0], [1.0, 0.0, 2.0, -1.0], [1.0, -1.0, 0.0, 0.0], [1.0, 0.0, 1.0, 0.0], [1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 1.0, 0.0], [0.0, 1.0, 0.0, 1.0], [0.0, 0.0, 1.0, -1.0], [-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, -1.0, 0.0], [-1.0, 0.0, 0.0, -1.0], [0.0, -1.0, -1.0, 0.0], [0.0, -1.0, 0.0, -1.0], [0.0, 0.0, -1.0, 1.0]]




mymod = [u1,u1,u1,d1,d1,d1,l1,l1,l1]
a = get_bilins(mymod)
b, mm = get_quads(mymod)
b
a

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

models, mmultis = model_list(nD=[5])
cmlists = runDFSZ_saveall("220916-serioustxtfiles", models; model_multiplicity=mmultis, log_calculate_all_models_if_below=20, sample_models_factor=0.001)
h5totxt("full_n4"; folder="./data/DFSZ_models/220916-serioustxtfiles/")


n=6
models, mmultis = model_list(nD=[n])
crs = read_EoN("220824-seriousruns", models; specifier=string(n)*"EoNs_full_nomulti")
crsv = abs.(collect(values(crs)))
crsk = collect(keys(crs))
sum(crsv[(!isnan).(crsk)])
@info sort(crsk)
sum(collect(values(crs))[abs.(collect(keys(crs)) .- 1.92) .< 0.04]) / sum(collect(values(crs)))

plot()
plot_EoN(crs)
xlims!((2/3,8/3))
ylims!((0.000, 0.015))

NaN < 0.04

e1 = FileIO.load("./data/DFSZ_models/220824-seriousruns/4EoNs_full_nomulti.jld2", "4/n4_u1_u1_u1_d1_d1_d1_l1_l1_l3")
e2 = FileIO.load("./data/DFSZ_models/220824-seriousruns/4EoNs_full_nomulti.jld2", "4/n4_u1_u1_u3_d1_d1_d1_l1_l1_l1")
e3 = FileIO.load("./data/DFSZ_models/220824-seriousruns/4EoNs_full_nomulti.jld2", "4/n4_u1_u1_u1_d1_d1_d3_l1_l1_l1")

ee = [e1,e2,e3]

eek = vcat(collect.(keys.(ee))...)
eev = vcat(collect.(values.(ee))...)
eev[(!isnan).(eek)]
sum(sum(eev[(!isnan).(eek)]))/3

e1 = FileIO.load("./data/DFSZ_models/220824-seriousruns/6EoNs_full_nomulti.jld2", "6/n6_u1_u1_u3_d1_d2_d3_l1_l1_l1")
e1 = FileIO.load("./data/DFSZ_models/220824-seriousruns/6EoNs_full_nomulti.jld2", "6/n6_u1_u1_u3_d1_d1_d3_l1_l1_l3")
e1 = FileIO.load("./data/DFSZ_models/220824-seriousruns/6EoNs_full_nomulti.jld2", "6/n6_u1_u1_u1_d1_d1_d3_l1_l2_l3")
e1 = FileIO.load("./data/DFSZ_models/220824-seriousruns/6EoNs_full_nomulti.jld2", "6/n6_u1_u2_u3_d1_d1_d1_l1_l1_l3")
e1 = FileIO.load("./data/DFSZ_models/220824-seriousruns/6EoNs_full_nomulti.jld2", "6/n6_u1_u1_u1_d1_d2_d3_l1_l1_l3")
e1 = FileIO.load("./data/DFSZ_models/220824-seriousruns/6EoNs_full_nomulti.jld2", "6/n6_u1_u1_u3_d1_d1_d1_l1_l2_l3")
e1 = FileIO.load("./data/DFSZ_models/220824-seriousruns/6EoNs_full_nomulti.jld2", "6/n6_u1_u2_u3_d1_d1_d3_l1_l1_l1")



model = [u1,u1,u1,d1,d1,d1,l1,l1,l3]

q,m = get_quads(model)

Nf = get_Nfunc(model)
EoNf = get_EoNfunc(model)

vals = [-1/3, -1/3, -1/3, -1, -1, -1, -1/3, -1/3, 1]

EoNf(vals)

Nf(vals)

aa = runDFSZ("test", [model])
aa