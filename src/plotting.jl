using StaticArrays
using Symbolics
using LaTeXStrings#, Plots
using FileIO
using LinearAlgebra
using Combinatorics

import PyPlot
const plt = PyPlot

include("helpers.jl")
include("ksvz.jl")
@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 

#model = fname2model(readdir("./data/DFSZ_models/preliminary2/n3")[1])
#fname=model2string(model)
#lab = fname[1:11]*"\n    "*fname[12:20]*"\n    "*fname[21:end]

dataset="220615-nbilinears"

function all_data(dataset; ns=:all, do_plot=nothing)
    if ns == :all
        ns = collect(3:9)
    else
        ns = ns
    end

    gaghs = similar(ns, Any)
    ARhs = similar(ns, Any)
    @time for (k,n) in enumerate(ns)
        @info "$n"
        savefolder = "./plots/"*dataset*"/ARs/n"*string(n)
        mkpath(savefolder)
        folders = readdir("./data/DFSZ_models/"*dataset*"/n"*string(n); join=true)
        
        ARtot_list = similar(folders, Any)
        gagtot_list = similar(folders, Any)
        for (j, folder) in enumerate(folders)
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
                if do_plot == :bilinear
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
            if do_plot == :yukawamodel
                plot_AR(ARtot, dataset, "ARs/n$n/"*split(folder,"/")[end])
            end
        end
        ARall = merge(ARtot_list...)
        gagall = merge(gagtot_list...)
        gagall = normalize(gagall; mode=:probability)
        gaghs[k] = gagall
        ARhs[k] = ARall
        if do_plot == :full
            plot_AR(ARall, dataset, "ARs/n$n/n$n"*"_fullAR")
        end
    end
    return ARhs, gaghs
end

Arhs, gaghs = all_data(dataset; ns=:all, do_plot=nothing)

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
.

#=
function all_data(dataset;merge=true)
    gaghs = similar(3:9, Any)
    ARhs = similar(3:9, Any)
    for k in 3:9
        @info "$k"
        fold = dataset*"/n"*string(k)
        folders = readdir("./data/DFSZ_models/"*fold)
        ARtot_list = similar(folders, Any)
        gagtot_list = similar(folders, Any)
        for (j, folder) in enumerate(folders)
            files = readdir("./data/DFSZ_models/"*fold*"/"*folder)
            m = fname2m(folder)
            hist1_list = similar(files, Any)
            hist2_list = similar(files, Any)
            @time for (i, file) in enumerate(files)
                ARh, gagh = get_data(file; folder=fold, m=m)
                hist1_list[i] = ARh
                hist2_list[i] = gagh
            end
            ARtot = hist1_list[1] #merge(hist1_list...)
            gagtot = hist2_list[1] #merge(hist2_list...)
            #tttot = normalize(tttot; mode=:probability)
            #tttot.weights .*= parse(Int,string(m))
            ARtot_list[j] = ARtot
            gagtot_list[j] = gagtot
        end

        ARhs[k-2] = merge(ARtot_list...)
        gaghs[k-2] = merge(gagtot_list...)
    end
    return ARhs, gaghs
end

function get_data(file; folder="/preliminary2/n"*string(4)*"/", m=m)
    model = fname2model(file)
    bilin = fname2bilin(file)
    ARh = read_AR(model; folder=folder, m=m, bilin=bilin)
    #ARh = normalize(ARh; mode=:probability)
    gagh = gag_histogram(ARh; mode=:probability, edges=-16.5:0.001:-12)
    gagh = normalize(gagh; mode=:probability)
    #ARh.weights .*= parse(Int,string(m))
    #gagh.weights .*= parse(Int,string(m))
    return ARh, gagh
end
=#
#=
function init_cdf()
    p1 = plot(label="", title="DFSZ axion model CDF", 
        xlabel=L"g_{a\gamma\gamma} \;\; [\log\;\mathrm{GeV}^{-1}]", ylabel="Probability for bigger gaγγ",
        bottom_margin=2Plots.mm, legend=:topright,
        size=(400,300), lw=2)
    return p1
end

function plot_cdf!(H, xH; kwargs...)
    p1 = plot!(xH.edges[1][1:end-1], H; label="", lw=2, kwargs...)
    return p1
end

function init_pdf()
    p1 = plot(lt=:stepbins, label="", title="DFSZ axion model PDF", 
        xlabel=L"ga\gamma\gamma \;\; [\log\;\mathrm{GeV}^{-1}]", ylabel="Probability",
        bottom_margin=2Plots.mm, legend=:topright,
        size=(400,300), lw=2)
    return p1
end

function plot_pdf!(H; kwargs...)
    p1 = plot!(H; lt=:stepbins, lw=1, label="", kwargs...)
    return p1
end
=#

function limit(frac, H)
    H.edges[1][1:end-1][cdf(H) .> frac][end]
end


# Read KSVZ anomaly ratios from Plakkots data.
KSVZ_ARs, KSVZgag, n_dw = ksvz("all"; edges=-50:0.01:50)
KSVZcdf = cdf(KSVZ_ARs)
KSVZcdf[5153]
KSVZ_ARs.edges[1][5153]

ARhs, gaghs = all_data("220615-nbilinears")
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
plt.xlabel("Anomaly Ratio E/N")
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

nlist=[4,5,6,7]
nomgaghlist = similar(nlist, Any)
for (i, n) in enumerate(nlist)
    fid = h5open("./data/DFSZ_models/220616-nbilin_fullsol/n$n/full_n$n.h5")

    allEoN = []
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
                myEoN = EoNlist[unique(i -> round.(cclist[i,:],digits=4), 1:size(cclist)[1])]
                myEoN = myEoN[-1e10 .< myEoN .< 1e10]
                append!(allEoN, myEoN)
            end
        end
    end

    nomARh = fit(Histogram, allEoN, minimum(allEoN)-0.01:0.01:maximum(allEoN)+0.01)
    gagh = gag_histogram(nomARh; mode=:probability, edges=-16.5:0.001:-12)
    gagh = normalize(gagh; mode=:probability)
    nomgaghlist[i] = gagh
    plot_AR(nomARh, "220616-nbilin_fullsol", "ARs_nomulti_n$n")
    close(fid)
end
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



nomgaghlist[3]
gagcdf = cdf(nomgaghlist[3])


fig, ax = plt.subplots(figsize=(6, 4))
ax.step(KSVZgag.edges[1][2:end], cdf(KSVZgag); where="pre", lw=3, color=c1, label="KSVZ-like (all)")
ax.step(nomgaghlist[3].edges[1][2:end], gagcdf; where="pre", lw=3, color=c2, alpha=0.4, label="DFSZ-like (all)")
plt.xlim([-16.5,-12])
plt.legend(loc="best")
plt.xlabel("Anomaly Ratio E/N")
plt.ylabel("Probability")
plt.savefig("plots/220616-nbilin_fullsol/CDFcompare6.pdf")

=#
