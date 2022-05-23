using StaticArrays
using Symbolics
#using LaTeXStrings, Plots
using FileIO
using LinearAlgebra

import PyPlot
const plt = PyPlot

include("helpers.jl")
include("ksvz.jl")
@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 

#model = fname2model(readdir("./data/DFSZ_models/preliminary2/n3")[1])
#fname=model2string(model)
#lab = fname[1:11]*"\n    "*fname[12:20]*"\n    "*fname[21:end]

function all_data()
    gaghs = similar(3:9, Any)
    ARhs = similar(3:9, Any)
    for k in 3:9
        @info "$k"
        fold = "preliminary2/n"*string(k)
        files = readdir("./data/DFSZ_models/"*fold)
        hist1_list = similar(files, Any)
        hist2_list = similar(files, Any)
        @time for (i, file) in enumerate(files)
            ARh, gagh = get_data(file; folder=fold)
            hist1_list[i] = ARh
            hist2_list[i] = gagh
        end
        ARhs[k-2] = merge(hist1_list...)
        gaghs[k-2] = merge(hist2_list...)
    end
    ARh =  merge(ARhs...)
    gagh = merge(gaghs...)
    return ARhs, gaghs
end

function get_data(file; folder="/preliminary2/n"*string(4)*"/")
    model = fname2model(file)
    ARh = read_AR(model; folder=folder*"/")
    #ARh = normalize(ARh; mode=:probability)
    gagh = gag_histogram(ARh; mode=:probability)
    gagh = normalize(gagh; mode=:probability)
    return ARh, gagh
end

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

function limit(frac, H)
    H.edges[1][1:end-1][cdf(H) .> frac][end]
end


# Read KSVZ anomaly ratios from Plakkots data.
KSVZ_ARs, KSVZgag, n_dw = ksvz("all"; edges=-50:0.01:50)
KSVZcdf = cdf(KSVZ_ARs)

ARhs, gaghs = all_data()
@time gagcdfs = cdf.(gaghs)
@time ARcdfs = cdf.(ARhs)

# make a histogram of just the abs(EoN - 1.92) part to feed to limit plot
Arr = rescale_histogram(merge(ARhs...))
limit(0.68, Arr)

#=
# Make comparison plot E/N pdf of DFSZ vs KSVZ
myAR = normalize(merge(ARhs...); mode=:probability)
KSVZAR = normalize(KSVZ_ARs; mode=:probability)

fig, ax = plt.subplots(figsize=(7, 4))
ax.stairs(KSVZAR.weights, KSVZAR.edges[1], linewidth=4000, lw=1.5, ec="mediumseagreen", alpha=1)
ax.stairs(myAR.weights, myAR.edges[1], linewidth=4000, lw=1.5, ec="maroon", alpha=0.4)
#ax.step(myAR.edges[1][2:end], myAR.weights)
plt.xlim([-50,50])
plt.yscale("log")
plt.ylim([1e-5,3e-1])
plt.xlabel("Anomaly Ratio E/N")
plt.ylabel("Probability")
plt.savefig("plots/preliminary2/PDFcompare.pdf")
#################################################
=#

#=
# Make zoomed in comparison plot E/N pdf of DFSZ vs KSVZ
myAR = normalize(merge(ARhs...); mode=:probability)
KSVZAR = normalize(KSVZ_ARs; mode=:probability)

fig, ax = plt.subplots(figsize=(7, 4))
ax.stairs(KSVZAR.weights, KSVZAR.edges[1], linewidth=4000, fill=true, fc="mediumseagreen", alpha=1)
ax.stairs(myAR.weights, myAR.edges[1], linewidth=4000, fill=true, fc="maroon", alpha=0.4)
#ax.step(myAR.edges[1][2:end], myAR.weights)
plt.xlim([-0,3])
plt.yscale("log")
plt.ylim([1e-5,3e-1])
plt.xlabel("Anomaly Ratio E/N")
plt.ylabel("Probability")
plt.savefig("plots/preliminary2/PDFcompareZoom.pdf")
#################################################
=#
