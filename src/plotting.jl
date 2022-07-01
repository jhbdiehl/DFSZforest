using StaticArrays
using Symbolics
using LaTeXStrings#, Plots
using FileIO
using LinearAlgebra
using Combinatorics
using HDF5

import PyPlot
const plt = PyPlot

include("helpers.jl")
include("ksvz.jl")
@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 

#model = fname2model(readdir("./data/DFSZ_models/preliminary2/n3")[1])
#fname=model2string(model)
#lab = fname[1:11]*"\n    "*fname[12:20]*"\n    "*fname[21:end]

function plot_AR(myAR, dataset, folder; ec="maroon", lw=2, alpha=1.0)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.stairs(myAR.weights, myAR.edges[1], ec=ec, lw=lw, alpha=alpha)
    plt.xlim([5/3-20,5/3+20])
    #plt.yscale("log")
    #plt.ylim([1e-5,1])
    plt.xlabel("Anomaly Ratio E/N")
    plt.ylabel("Nr of models")
    #plt.title(L"All DFSZ-like $n_H = 9$ models")
    #plt.axvline(5/3, ls=":", color="grey")
    #plt.axvline(2/3, ls=":", color="k")
    #plt.axvline(8/3, ls=":", color="k")
    #plt.text(4,1000, "symmetry axis 5/3", color="grey")
    #plt.text(-5,2000, "DFSZ-II", color="k")
    #plt.text(4,2000, "DFSZ-I", color="k")
    plt.savefig("plots/"*dataset*"/"*folder*".pdf")
end

function read_data(dataset; ns=collect(3:9), do_plot=nothing)

    gaghs = similar(ns, Any)
    ARhs = similar(ns, Any)
    @time for (k,n) in enumerate(ns)
        @info "$n"
        savefolder = "./plots/"*dataset*"/ARs/n"*string(n)
        mkpath(savefolder)
        folders = readdir("./data/DFSZ_models/"*dataset*"/n"*string(n); join=true)
        folders = folders[occursin.(".h5",folders) .== 0] # Throw away output of full_solution=true
        
        ARtot_list = similar(folders, Any)
        gagtot_list = similar(folders, Any)
        for (j, folder) in enumerate(folders)
            println(folder)
            m = fname2m(folder)
            files = readdir(folder)
            hist1_list = similar(files, Any)
            hist2_list = similar(files, Any)
            @time for (i, file) in enumerate(files)
                model = fname2model(file)
                bilin = fname2bilin(file)
                ARh = read_AR(model; folder=dataset*"/n"*string(n), bilin=bilin, m=m)
                #ARh = normalize(ARh; mode=:probability)
                gagh = gag_histogram(ARh; mode=:probability, edges=-16.5:0.001:-12)
                gagh = normalize(gagh; mode=:probability)
                #gagh = gag_histogram(tt; mode=:probability)
                hist1_list[i] = ARh
                hist2_list[i] = gagh
                if do_plot == :bilinear || do_plot == :all
                    plot_AR(ARh, dataset, "ARs/n$n/"*split(file, ".")[1])
                end
            end
            ARtot = merge(hist1_list...)
            gagtot = merge(hist2_list...)
            gagtot = normalize(gagtot; mode=:probability)
            gagtot.weights .*= parse(Int,string(m))
            ARtot.weights .*= parse(Int,string(m))
            ARtot_list[j] = ARtot
            gagtot_list[j] = gagtot
            if do_plot == :yukawamodel || do_plot == :all
                plot_AR(ARtot, dataset, "ARs/n$n/"*split(folder,"/")[end])
            end
        end
        ARall = merge(ARtot_list...)
        gagall = merge(gagtot_list...)
        gagall = normalize(gagall; mode=:probability)
        gaghs[k] = gagall
        ARhs[k] = ARall
        if do_plot == :full || do_plot == :all
            plot_AR(ARall, dataset, "ARs/n$n/n$n"*"_fullAR")
        end
    end
    return ARhs, gaghs
end

Arhs, gaghs = read_data("test"; ns=[5], do_plot=nothing)

function limit(frac, H)
    H.edges[1][1:end-1][cdf(H) .> frac][end]
end


# Read KSVZ anomaly ratios from Plakkots data.
KSVZ_ARs, KSVZgag, n_dw = ksvz("all"; edges=-50:0.01:50)
KSVZcdf = cdf(KSVZ_ARs)
KSVZcdf[5153]
KSVZ_ARs.edges[1][5153]

ARhs, gaghs = read_data("220615-nbilinears")
@time gagcdfs = cdf.(gaghs)
@time ARcdfs = cdf.(ARhs)


gagh = normalize(KSVZgag, mode=:probability)
allgagcdf = cdf(gagh)
allgagcdf[912]
gagh.edges[1][912]

gaghs = normalize.(gaghs, mode=:probability)
allgagcdf = cdf(merge(gaghs...))
allgagcdf[366]
gaghs[7].edges[1][365]

ARcdfs[7][100167]
ARhs[7].edges[1][100167]
i = 7
sum(ARhs[i].weights .* (ARhs[i].edges[1][1:end-1] .+ 0.005)) / sum(ARhs[i].weights)
sum(KSVZ_ARs.weights .* (KSVZ_ARs.edges[1][1:end-1] .+ 0.005)) / sum(KSVZ_ARs.weights)

# make a histogram of just the abs(EoN - 1.92) part to feed to limit plot
Arr = rescale_histogram(merge(ARhs...))
limit(0.68, Arr)


c1 = "mediumseagreen"
c2 = "maroon"

c1a = 1.0
c2a = 0.4


#=
# Plot all n=4 models in one histogram
k=4
myAR = ARhs[2]#normalize(ARhs[k-2]; mode=:probability)
fig, ax = plt.subplots(figsize=(6, 4))
ax.stairs(myAR.weights, myAR.edges[1], ec=c2, lw=2, alpha=c2a)
plt.xlim([5/3-20,5/3+20])
plt.yscale("log")
plt.ylim([1e0,5e3])
plt.xlabel("Anomaly Ratio E/N")
plt.ylabel("Number of models")
plt.title(L"All DFSZ-like $n_H = $"*string(k)*" models")
plt.axvline(5/3, ls=":", color="grey")
plt.axvline(2/3, ls=":", color="k")
plt.axvline(8/3, ls=":", color="k")
plt.text(4,1000, "symmetry axis 5/3", color="grey")
plt.text(-5,2000, "DFSZ-II", color="k")
plt.text(4,2000, "DFSZ-I", color="k")
plt.savefig("plots/220615-nbilinears/PDFn"*string(k)*"tot.pdf")
#################################################
=#

#=
# Plot all n=4 models separately
k = 4
@info "$k"
fold = "220615-nbilinears_addendum/n"*string(k)
folders = readdir("./data/DFSZ_models/"*fold, join=true)

ARtot_list = similar(folders, Any)
for (j, folder) in enumerate(folders)
    m = fname2m(folder)
    files = readdir(folder)
    hist1_list = similar(files, Any)
    @time for (i, file) in enumerate(files)
        model = fname2model(file)
        bilin = fname2bilin(file)
        ARh = read_AR(model; folder="220615-nbilinears_addendum/n"*string(k), bilin=bilin, m=m)
        #ARh = normalize(ARh; mode=:probability)
        hist1_list[i] = ARh
    end
    ARtot = merge(hist1_list...)
    ARtot.weights .*= parse(Int,string(m))
    ARtot_list[j] = ARtot
end


slist= ["up", "charm", "top", "down", "strange", "bottom", "electron", "muon", "tau"]
a = collect(with_replacement_combinations(1:3, 2))
b = collect(permutations(1:3,2))
ab = sort(unique(vcat(a,b)))

fig, ax = plt.subplots(3,3,figsize=(8, 6), sharex=true, sharey=true)
i = 1
for h in ARtot_list[end:-1:1]
    #h = normalize(h; mode=:probability)
    ax[ab[i][1],ab[i][2]].set_xlim([5/3-20,5/3+20])
    ax[ab[i][1],ab[i][2]].set_yscale("log")
    #ax[ab[i][1],ab[i][2]].set_ylim([1e-3,3e-1])
    ax[3,ab[i][2]].set_xlabel("Anomaly Ratio E/N")
    ax[ab[i][1],1].set_ylabel("Number of models")
    ax[ab[i][1],ab[i][2]].stairs(h.weights, h.edges[1], lw=2, ec=c2, alpha=0.4)
    ax[ab[i][1],ab[i][2]].set_title(L"..."*slist[i])
    i += 1
end
plt.suptitle(L"DFSZ-like $n_H = 4$ models, special coupling to...")
plt.savefig("plots/220615-nbilinears/PDFn4all.pdf")
#################################################
=#

#=
# comparison plot E/N cdf of DFSZ vs KSVZ
#myAR = normalize(merge(ARhs...); mode=:probability)
mygag = normalize(merge(gaghs...); mode=:probability)
gagcdf = cdf(merge(gaghs...))
#KSVZAR = normalize(KSVZ_ARs; mode=:probability)
KSVZcdf = cdf(KSVZgag)

fig, ax = plt.subplots(figsize=(6, 4))
ax.step(KSVZgag.edges[1][2:end], KSVZcdf; where="pre", lw=3, color=c1, label="KSVZ-like (all)")
ax.step(mygag.edges[1][2:end], gagcdf; where="pre", lw=3, color=c2, alpha=0.4, label="DFSZ-like (all)")
plt.xlim([-16.5,-12])
plt.legend(loc="best")
plt.xlabel(L"\log{g_{a \gamma}[\mathrm{GeV}^{-1}]}\: \: \mathrm{at}\: \: m_a = 40\; \mu\mathrm{eV}")
plt.ylabel("Probability")
plt.savefig("plots/220615-nbilinears/CDFcompare.pdf")
#################################################
=#

#=
# Make comparison plot E/N pdf of DFSZ vs KSVZ
myAR = normalize(merge(ARhs...); mode=:probability)
KSVZAR = normalize(KSVZ_ARs; mode=:probability)

fig, ax = plt.subplots(figsize=(7, 4))
ax.stairs(KSVZAR.weights, KSVZAR.edges[1], linewidth=4000, lw=1.5, ec="mediumseagreen", alpha=1, label="KSVZ-like (all)")
ax.stairs(myAR.weights, myAR.edges[1], linewidth=4000, lw=1.5, ec="maroon", alpha=0.4, label="DFSZ-like (all)")
#ax.step(myAR.edges[1][2:end], myAR.weights)
plt.xlim([-50,50])
plt.yscale("log")
plt.ylim([1e-5,3e-1])
plt.xlabel("Anomaly Ratio E/N")
plt.ylabel("Probability")
plt.legend(loc="best")
plt.savefig("plots/220615-nbilinears/PDFcompare.pdf")
#################################################
=#

#=
# Make zoomed in comparison plot E/N pdf of DFSZ vs KSVZ
myAR = normalize(merge(ARhs...); mode=:probability)
KSVZAR = normalize(KSVZ_ARs; mode=:probability)

fig, ax = plt.subplots(figsize=(7, 4))
ax.stairs(KSVZAR.weights, KSVZAR.edges[1], linewidth=4000, fill=true, fc="mediumseagreen", alpha=1, label="KSVZ-like (all)")
ax.stairs(myAR.weights, myAR.edges[1], linewidth=4000, fill=true, fc="maroon", alpha=0.4, label="DFSZ-like (all)")
#ax.step(myAR.edges[1][2:end], myAR.weights)
plt.xlim([-0,3])
plt.yscale("log")
plt.ylim([1e-5,3e-1])
plt.xlabel("Anomaly Ratio E/N")
plt.ylabel("Probability")
plt.legend(loc="lower center")
plt.savefig("plots/220615-nbilinears/PDFcompareZoom.pdf")
#################################################
=#







#=
using HDF5



=
fid

fold = fid["3n4_u1_u1_u1_d1_d2_d2_l1_l1_l1"]
println(fold)
cclist = Array{Float64}(undef, 0, 10)
EoNlist = []
for bil in fold
    cc = read(bil["Chis"])
    cclist = vcat(cclist, cc)
    EoN = read(bil["EoN"])
    append!(EoNlist, EoN)
end

myEoNu = EoNlist[unique(i -> cclist[i,:], 1:size(cclist)[1])]

unique(i -> cclist[i,:], 1:size(cclist)[1])
unique(i -> round.(cclist[i,:],digits=4), 1:size(cclist)[1])

using StatsPlots
histogram(EoNlistd[-1e10 .< EoNlistd .< 1e10].-5/3,bins=-20:0.1:20)
histogram(-EoNlistu[-1e10 .< EoNlistu .< 1e10].+5/3,bins=-20:0.1:20)

myEoNlistu = []
for ccs in eachrow(cclistu)
    append!(myEoNlistu, EoverN.(ccs[1:9]...))
end
myEoNlistu .≈ EoNlistu

sortslices(cclistd, by=i->EoverN(i[1:9]...), dims=1)
sortslices(cclistu, by=i->EoverN(i[1:9]...), dims=1)

sort(fEoNlistd, by=i->abs(i))

cclistd[unique(i -> round.(cclistd[i,:], digits=4), 1:size(cclistd)[1]),:]
cclistu[unique(i -> round.(cclistu[i,:], digits=4), 1:size(cclistu)[1]),:]
sortslices(cclistu, dims=1)
bcclistd = deepcopy(cclistd)
bcclistd[:,1:3] = -cclistd[:,4:6]
bcclistd[:,4:6] = -cclistd[:,1:3]
sortslices(bcclistd, dims=1)
fEoNlistd = EoNlistd[-1e10 .< EoNlistd .< 1e10]
mybool = sort(EoNlistd[-1e10 .< EoNlistd .< 1e10].-5/3) .≈ sort(-EoNlistu[-1e10 .< EoNlistu .< 1e10].+5/3)
sort(EoNlistd[-1e10 .< EoNlistd .< 1e10].-5/3)[mybool .== 0]
cclist[unique(i -> cclist[i,:], 1:size(cclist)[1]), :]
myEoNd
myEoNu

EoNlist

myEoN = myEoN[-1e10 .< myEoN .< 1e10]
append!(allEoN, myEoN)
=#
fid = h5open("./data/DFSZ_models/test/n5/full_n5.h5")
for fold in fid
    println(split(split(string(fold), "/")[2], " ")[1])
end
close(fid)
nlist=[5]

function read_full_data(dataset; ns=collect(3:9), do_plot=nothing)
    gaghs = similar(ns, Any)
    Arhs = similar(ns, Any)
    for (i, n) in enumerate(ns)
        fid = h5open("./data/DFSZ_models/"*dataset*"/n$n/full_n$n.h5")
        savefolder = "./plots/"*dataset*"/AR_nomulti/n$n"
        mkpath(savefolder)

        allEoN = []
        allweights = Real[]
        for fold in fid
            if occursin("Chis order: u1u2u3d1d2d3l1l2l3s", string(fold))
                nothing
            else
                @time begin
                    println(fold)
                    cclist = Array{Float64}(undef, 0, 10)
                    EoNlist = []
                    for bil in fold
                        cc = read(bil["Chis"])
                        cclist = vcat(cclist, cc)
                        EoN = read(bil["EoN"])
                        append!(EoNlist, EoN)
                    end
                    mult = parse(Int64,split(split(string(fold),"/")[2],"n")[1])
                    myEoN = EoNlist[unique(i -> round.(cclist[i,:],digits=4), 1:size(cclist)[1])]
                    myEoN = myEoN[-1e10 .< myEoN .< 1e10] 
                    if do_plot == :all || do_plot == :yukawamodel
                        ARhyuk = fit(Histogram, myEoN, FrequencyWeights(mult .* ones(size(myEoN))), -50:0.01:50)
                        plot_AR(ARhyuk, dataset, "AR_nomulti/n$n/"*split(split(string(fold), "/")[2], " ")[1])
                    end
                    append!(allEoN, myEoN)
                    append!(allweights, mult .* ones(size(myEoN)))
                end
            end
        end

        ARh = fit(Histogram, allEoN, FrequencyWeights(allweights), -50:0.01:50)
        Arhs[i] = ARh
        gagh = gag_histogram(ARh; mode=:probability, edges=-16.5:0.001:-12)
        gagh = normalize(gagh; mode=:probability)
        gaghs[i] = gagh
        if do_plot == :all || do_plot == :full
            plot_AR(ARh, dataset, "AR_nomulti/n$n")
        end
        close(fid)
    end
    return Arhs, gaghs
end

noMarhs, noMgaghs = read_full_data("test", ns=[5], do_plot=nothing)

nomgaghall = merge(nomgaghlist...)
nomgaghall = merge(nomgaghlist[1:3]...)

fig, ax = plt.subplots(figsize=(7, 4))
ax.stairs(ARhs[2].weights, ARhs[2].edges[1], linewidth=4000, lw=1.5, ec="mediumseagreen", alpha=1, label="all n=4 models")
ax.stairs(nomARh4.weights, nomARh4.edges[1], linewidth=4000, lw=1.5, ec="maroon", alpha=1, label="different PQ charges")
#ax.step(myAR.edges[1][2:end], myAR.weights)
plt.xlim([5/3-20,5/3+20])
plt.yscale("log")
plt.ylim([1e0,4e3])
plt.xlabel("Anomaly Ratio E/N")
plt.ylabel("Number of models")
plt.legend(loc="best")
plt.savefig("plots/220616-nbilin_fullsol/PDFvsmulti4.pdf")



nomgaghlist[1]
gagcdf = cdf(nomgaghlist[1])


fig, ax = plt.subplots(figsize=(6, 4))
ax.step(KSVZgag.edges[1][2:end], cdf(KSVZgag); where="pre", lw=3, color=c1, label="KSVZ-like (all)")
ax.step(nomgaghlist[1].edges[1][2:end], gagcdf; where="pre", lw=3, color=c2, alpha=0.4, label="DFSZ-like (n=5)")
plt.xlim([-16.5,-12])
plt.legend(loc="best")
plt.xlabel(L"\log{g_{a \gamma}[\mathrm{GeV}^{-1}]}\: \: \mathrm{at}\: \: m_a = 40 \; \mu\mathrm{eV}")
plt.ylabel("Probability")
plt.savefig("plots/220616-nbilin_fullsol/CDFcompare5.pdf")

=#


function my_read_full_data(dataset; ns=collect(3:9), do_plot=nothing)
    gaghs = similar(ns, Any)
    Arhs = similar(ns, Any)
    for (i, n) in enumerate(ns)
        fid = h5open("./data/DFSZ_models/"*dataset*"/n$n/full_n$n.h5")
        savefolder = "./plots/"*dataset*"/AR_nomulti/n$n"
        mkpath(savefolder)

        allEoN = []
        allweights = Real[]
        for fold in fid
            if occursin("Chis order: u1u2u3d1d2d3l1l2l3s", string(fold))
                nothing
            else
                @time begin
                    println(fold)
                    cclist = Array{Float64}(undef, 0, 10)
                    EoNlist = []
                    for bil in fold
                        cc = read(bil["Chis"])
                        cclist = vcat(cclist, cc)
                        EoN = read(bil["EoN"])
                        append!(EoNlist, EoN)
                    end
                    mult = parse(Int64,split(split(string(fold),"/")[2],"n")[1])
                    myEoN = EoNlist[unique(i -> round.(cclist[i,:],digits=4), 1:size(cclist)[1])]
                    myEoN = myEoN[-1e10 .< myEoN .< 1e10] 
                    println(length(myEoN))
                    if do_plot == :all || do_plot == :yukawamodel
                        ARhyuk = fit(Histogram, myEoN, FrequencyWeights(mult .* ones(size(myEoN))), -50:0.01:50)
                        plot_AR(ARhyuk, dataset, "AR_nomulti/n$n/"*split(split(string(fold), "/")[2], " ")[1])
                    end
                    append!(allEoN, myEoN)
                    append!(allweights, mult .* ones(size(myEoN)))
                end
            end
        end

        ARh = fit(Histogram, allEoN, FrequencyWeights(allweights), -50:0.01:50)
        Arhs[i] = ARh
        gagh = gag_histogram(ARh; mode=:probability, edges=-16.5:0.001:-12)
        gagh = normalize(gagh; mode=:probability)
        gaghs[i] = gagh
        if do_plot == :all || do_plot == :full
            plot_AR(ARh, dataset, "AR_nomulti/n$n")
        end
        close(fid)
    end
    return Arhs, gaghs
end

ARh5, gagh5 = my_read_full_data("test", ns=[5])


fid = h5open("./data/DFSZ_models/220616-nbilin_fullsol/n5/full_n5.h5");
for fold in fid
    if occursin("Chis order: u1u2u3d1d2d3l1l2l3s", string(fold))
        nothing
    else
        for bil in fold
            println(bil)
        end
    end
end
close(fid)
EoNlists = Array{Any}(undef, length(fid)-1)
cclists = Array{Any}(undef, length(fid)-1)

#=
io = open("./data/DFSZ_models/test/foremi.txt", "w")

write(io, "# yukawa \t bilinear \t maxEoN \t corresponding charges (u1u2u3d1d2d3l1l2l3s) \t minEoN \t corresponding charges (u1u2u3d1d2d3l1l2l3s) \n")


for (i, fold) in enumerate(fid)
    if occursin("Chis order: u1u2u3d1d2d3l1l2l3s", string(fold))
        nothing
    else
        @time begin
            for bil in fold
                println(bil)
                cc = read(bil["Chis"])
                EoN = read(bil["EoN"])
                k = -1e10 .< EoN .< 1e10
                EoN = EoN[k]
                cc = cc[k,:]
                println(cc[1:5,:])
                maxEoN = round(maximum(EoN); digits=2)
                minEoN = round(minimum(EoN); digits=2)
                mincc = cc[findall(EoN .== minimum(EoN)),:]
                maxcc = cc[findall(EoN .== maximum(EoN)),:]
                f = split(split(string(fold),"/")[2], " ")[1]
                b = split(split(string(bil),"/")[3], " ")[1]
                write(io, "$f \t $b \t $maxEoN \t $maxcc \t $minEoN \t $mincc \n")
                println()
                println(size(cc))
            end
        end
    end
end

close(io)
=#

fid = h5open("./data/DFSZ_models/test/n4/full_n4.h5");
EoNlists = Array{Any}(undef, length(fid)-1)
cclists = Array{Any}(undef, length(fid)-1)
for (i, fold) in enumerate(fid)
    if occursin("Chis order: u1u2u3d1d2d3l1l2l3s", string(fold))
        nothing
    else
        @time begin
            println(fold)
            cclist = Array{Float64}(undef, 0, 10)
            EoNlist = []
            for bil in fold
                println(bil)
                cc = read(bil["Chis"])
                cclist = vcat(cclist, cc)
                EoN = read(bil["EoN"])
                append!(EoNlist, EoN)
            end
            mult = parse(Int64,split(split(string(fold),"/")[2],"n")[1])
            myEoN = EoNlist[unique(i -> round.(cclist[i,:],digits=4), 1:size(cclist)[1])]
            myEoN = myEoN[-1e10 .< myEoN .< 1e10] 
            EoNlists[i] = EoNlist
            println(size(cclist))
            cclists[i] = cclist
            println(length(EoNlist))
        end
    end
end
close(fid)
EoNlists[1][ -1e10 .< EoNlists[1] .< 1e10]
EoNlists2[1][ -1e10 .< EoNlists2[1] .< 1e10]

histogram(EoNlists[1])
histogram(EoNlists2[1])

cclists[3]

rslists = Array{Any}(undef, 10-1)
EoNlists2 = Array{Any}(undef, 10-1)
aslists = Array{Any}(undef, 10-1)
bslists = Array{Any}(undef, 10-1)
quadlists = Array{Any}(undef, 10-1)
chislists = Array{Any}(undef, 10-1)
a, ms = generate_all_models()
mya = a[length.(unique.(a)) .∈ Ref([4])]
#mya = [[u2, u1, u1, d1, d1, d1, l1, l1, l1],  [u1, u2, u2, d1, d1, d1, l1, l1, l1]]
odd = (1,1)
even = (1,-1)

@time for (k, model) in enumerate(mya) #a[length.(unique.(a)) .== 8]
    bilins = unique(sort.(collect(combinations(model,2)), by=x->Symbol(x)))
    bilins = bilins[length.(unique.(bilins)) .== 2]
    tott = 0
    EoNlist = []
    rslist= []
    quadlist=[]
    aslist=[]
    bslist=[]
    cclist = Array{Float64}(undef, 0, 10)
    println(model)
    for (i, bilin) in enumerate(bilins)
        for (odd, even) in [((1,1),(1,-1)), ((-1,-1),(-1,1))]
            if i <10
                println(bilin)
                valp1, valp2 = bilinvals(bilin; odd=odd, even=even)
                un = unique(model)
                nH = length(un)
                quads, multis = get_quads(model; p1=bilin[1], p2=bilin[2], valp1=valp1, valp2=valp2)
                append!(quadlist, quads)
                if length(un) >= 7
                    bi = bilinsum.(bilins[i+1:end]; odd=odd, even=even)
                    bi1 = bilinsum.(bilins[i+1:end]; odd=odd, even=even)
                    bi2 = bilinsum.(bilins[i+1:end]; odd=-1 .* odd, even=-1 .*even)
                else # If you sample all bilins with equal nr of samples, excluding models you already calculated would lead to bias effect!
                    bi = bilinsum.(bilins[1:end .!= i]; odd=odd, even=even)
                    bi1 = bilinsum.(bilins[1:end .!= i]; odd=odd, even=even)
                    bi2 = bilinsum.(bilins[1:end .!= i]; odd=-1 .* odd, even=-1 .*even)
                end
                #println(bi)
                #terms = vcat(quads,bi)
                terms = vcat(quads,bi1,bi2)
                multis = vcat(multis, 1 * ones(Int64, size(vcat(bi1,bi2))...))
                
                as, bs = get_numquads(terms, un, nH; p1=bilin[1], p2=bilin[2], valp1=valp1, valp2=valp2)
                append!(aslist, as)
                append!(bslist, bs)
                tot = binomial(length(terms),nH-2)
                tott += tot
                myEoN = get_EoNfunc(model; p1=bilin[1], p2=bilin[2], valp1=valp1, valp2=valp2)
                proc_rs = similar(as, tot)
                EoN_rs = similar(bs, tot)
                rs_ws = similar(multis, length(proc_rs))
                parallel_alleqn_solve_proc_fullsol!(proc_rs, EoN_rs, rs_ws, as, bs, multis, tot, myEoN)
                #save_full(model, proc_rs, EoN_rs, rs_ws, 1; folder="test/n"*string(nH)*"/", bilin=bilin, valp1=valp1, valp2=valp2, ms=ms[k])
                good_idxs = findall(!isnan, EoN_rs)
                EoN_rs = EoN_rs[good_idxs]
                good_proc_rs = proc_rs[good_idxs,:]
                chi_s = ( abs(valp1) + abs(valp2) ) / 2
                mydict = Dict{Num, Float64}(un .=> 0)
                mydict[bilin[1]]= valp1
                mydict[bilin[2]]= valp2
                notp1p2 = isequal.(un,bilin[1]) .+ isequal.(un,bilin[2]) .== false

                Chis = Matrix{Float64}(undef, length(EoN_rs), 10)
                for i in 1:length(rs_ws[good_idxs])
                    for (j, higgs) in enumerate(un[notp1p2])
                        mydict[higgs] = good_proc_rs[i][j]
                    end
                    chivec = Symbolics.value.(substitute.(model, (mydict,)))
                    append!(chivec, chi_s)
                    Chis[i,:] .= chivec #Construct matrix with chi values, last one is charge of the singlet (always 1)
                end
                append!(EoNlist, EoN_rs)
                append!(rslist, good_proc_rs)
                cclist = vcat(cclist, Chis)
            else
                nothing
            end
        end
    end
    EoNlists2[k] = EoNlist
    rslists[k] = rslist
    quadlists[k] = quadlist
    aslists[k] = aslist
    bslists[k] = bslist
    chislists[k] = cclist
end

for (odd, even) in [((1,1),(1,-1)), ((-1,-1),(-1,1))]
    println(odd,even)
end
model = [u1,u1,u2,d1,d1,d1,l1,l1,l1]
bilins = unique(sort.(collect(combinations(model,2)), by=x->Symbol(x)))
bilins = bilins[length.(unique.(bilins)) .== 2]
b1 = bilinsum.(bilins; odd=(1,1), even=(1,-1))
b2 = bilinsum.(bilins; odd=(-1,-1), even=(-1,1))
vcat(b1, b2)

fid = h5open("./data/DFSZ_models/test/n4/full_n4.h5");
EoNlists = Array{Any}(undef, length(fid)-1)
cclists = Array{Any}(undef, length(fid)-1)
for (i, fold) in enumerate(fid)
    if occursin("Chis order: u1u2u3d1d2d3l1l2l3s", string(fold))
        nothing
    else
        @time begin
            println(fold)
            cclist = Array{Float64}(undef, 0, 10)
            EoNlist = []
            for bil in fold
                println(bil)
                cc = read(bil["Chis"])
                cclist = vcat(cclist, cc)
                EoN = read(bil["EoN"])
                append!(EoNlist, EoN)
            end
            mult = parse(Int64,split(split(string(fold),"/")[2],"n")[1])
            #myEoN = EoNlist[unique(i -> round.(cclist[i,:],digits=4), 1:size(cclist)[1])]
            #myEoN = myEoN[-1e10 .< myEoN .< 1e10] 
            EoNlists[i] = EoNlist
            println(size(cclist))
            cclists[i] = cclist
            println(length(EoNlist))
        end
    end
end
close(fid)

tt = Num[u2 + (2//1)*d1 - l1, (2//1)*l1 + (2//1)*u2, l1 + 2u1 - u2, u2 + 2l1 - d1, 2d1 + 2u1, 2l1 - (2//1)*d1, d1 + l1 + u1 + u2, d1 + 2u2 - u1, l1 + 2u2 - u1, l1 + u2 - d1 - u1, (2//1)*l1 + (2//1)*u1, d1 + l1 + 2u1, u1 + u2 + 2l1, d1 + l1 + 2u2, d1 + 2u1 - u2, l1 + u1 - d1 - u2, (2//1)*d1 + (2//1)*u2, u1 + 2l1 - d1, u1 + u2 + 2d1, u1 + (2//1)*d1 - l1, 2u1, 2u1, 2u2, 2u2, l1 - d1]



cpEoNlists2 = deepcopy(EoNlists2)
cprslists = deepcopy(rslists)
cpaslists = deepcopy(aslists)
cpbslists = deepcopy(bslists)
cpquadlists = deepcopy(quadlists)
cpchislists = deepcopy(chislists)

isequal(cpquadlists[1],quadlists[1])
aslists[1] == cpaslists[1]
bslists[1] == cpbslists[1]
bslist
cpbslists

roundrs1 = [round.(el, digits=4) for el in rslists[1]]
cproundrs1 = [round.(el, digits=4) for el in cprslists[1]]
roundrs3 = [round.(el, digits=4) for el in rslists[3]]
countmap(roundrs1) == countmap(cproundrs1)

c1 = chislists[1]
c2 = cpchislists[1]
c3 = chislists[3]
sum(c1[:,3] .== 1 .&& c1[:,1] .== -1)
countmap(c1[:,1])# == countmap(c1[:,3])
mergewith( (x,y) -> abs(x-y), countmap(c1[:,1]), countmap(c1[:,3]))

c1[-1e10 .< EoNlists2[1] .< 1e10,:]
eon1 = round.(EoNlists2[1][-1e10 .< EoNlists2[1] .< 1e10], digits=4)
eon2 = round.(cpEoNlists2[1][-1e10 .< cpEoNlists2[1] .< 1e10], digits=4)
eon3 = round.(EoNlists2[3][-1e10 .< EoNlists2[3] .< 1e10], digits=4)
countmap(eon3)# == 
countmap(eon1)
eon1 .== eon2
sum(values(mergewith( (x,y) -> abs(x-y), countmap(eon1), countmap(eon3))))
countmap(eon1)

EoNlisttmp = EoNlists[end:-1:1]
EoNlists = EoNlisttmp
a = EoNlists2[1][-1e10 .< EoNlists2[1] .< 1e10]
b = EoNlists2[3][-1e10 .< EoNlists2[3] .< 1e10]

a1 = countmap(round.(a, digits=4))
b1 = countmap(round.(b, digits=4))

mergewith( (x,y) -> abs(x-y), a1, b1)


a = [[1,1] [2,3]]
#diff.([1,2], [3,2])
isequal(a,b)
issetequal(a,b)

using StatsBase

cclists[1][round.(EoNlists2[1], digits=4) .== 0.4167,:]
countmap(round.(EoNlists[1][-1e10 .< EoNlists[1] .< 1e10], digits=4))
sum(round.(EoNlists[1], digits=4) .== 0.4167)
sum(isequal.(rslists[1],Ref(rslists[1][ind][3])))

ind = findall(x -> round.(x,digits=4) .== 0.4167, EoNlists[1])
findall(x -> isequal.(x, Ref(rslists[1][ind][3])), rslists[1])
cclists[1]
rslists[1][round.(EoNlists[1],digits=4) .== 0.4167]
rslists[3][round.(EoNlists[1],digits=4) .== 0.4167]
# Should be equal not only for 1,2 but also for 1,3!
a = countmap(sort.(collect(eachrow(round.(quadlists[1][-1e10 .< quadlists[1] .< 1e10], digits=4)))))
b = countmap(sort.(collect(eachrow(round.(quadlists[3][-1e10 .< quadlists[3] .< 1e10],digits=4)))))
a
a == b

mergewith( (x,y) -> abs(x-y), a, b)

cc1us = sortslices(cclists[1][:, 1:3], dims=1)
cc1ds = sortslices(cclists[1][:, 4:6], dims=1)
cc1ds[cc1ds .!= cc1us]
sum(cclists[7])
sortslices(cclists[2], dims=1)
sortslices(cclists[3], dims=1)
isequal(cclists[1], cclists[4])
EoNlists
issetequal(EoNlists[1], EoNlists[3])

idxarr, bnc = make_NEW_idx_bnc(5)
idxarr
idxs_i =  myNEWidxtup!(idxarr, bnc, 12312, Val(5))