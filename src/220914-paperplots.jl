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


#########################################################################################################
#
# Plot all n=4 models separately
#
#########################################################################################################

folder = "220824-seriousruns"

hs = similar(1:9, Any)
models, multis = model_list(nD=[4], compute_equivalent_theories=true)
for i in 1:9
    tmp = read_EoN(folder, [models[i]]; specifier="4EoNs_full_nomulti_equivalents")
    hs[i], m = _make_hist(tmp, bins=-20:0.0001:21)
    #hs[i] = normalize(hs[i]; mode=:probability)
end

fig, ax = plt.subplots(3,3,figsize=(8, 6), sharex=true, sharey=true)


slist= ["up", "charm", "top", "down", "strange", "bottom", "electron", "muon", "tau"]
a = collect(with_replacement_combinations(1:3, 2))
b = collect(permutations(1:3,2))
ab = sort(unique(vcat(a,b)))

fig, ax = plt.subplots(3,3,figsize=(8, 6), sharex=true, sharey=true)
i = 1
for h in hs
    #h = normalize(h; mode=:probability)
    ax[ab[i][1],ab[i][2]].set_xlim([5/3-20,5/3+20])
    #ax[ab[i][1],ab[i][2]].set_yscale("log")
    ax[ab[i][1],ab[i][2]].set_ylim([1e0,3e1])
    ax[ab[i][1],ab[i][2]].axvline(5/3, ls="--", lw=0.8, c="gray")
    ax[3,ab[i][2]].set_xlabel("Anomaly Ratio E/N")
    ax[ab[i][1],1].set_ylabel("Number of models")
    ax[ab[i][1],ab[i][2]].stairs(h.weights, h.edges[1], lw=2, ec="maroon", alpha=0.4)
    ax[ab[i][1],ab[i][2]].set_title(L"..."*slist[i])
    i += 1
end
plt.suptitle(L"DFSZ-type $n_D = 4$ models, special coupling to...")
plt.savefig("plots/220914-paperplots/PDFn4all.pdf")







#########################################################################################################
#
# Plot n=3 to n=9
#
#########################################################################################################

folder = "220824-seriousruns"
# Load n=3 to n=7
hs = similar(3:7, Any)
for nD in 3:7
    models, multis = model_list(nD=[nD])
    tmp = read_EoN(folder, models; specifier=string(nD)*"EoNs_full_nomulti")
    hs[nD-2], m = _make_hist(tmp)
    hs[nD-2] = normalize(hs[nD-2]; mode=:probability)
end

# Construct upper and lower limit for n=8 and n=9
models, multis = model_list(nD=[7])
crs7 = read_EoN(folder, models; specifier="7EoNs_full_nomulti")
models, multis = model_list(nD=[6])
crs6 = read_EoN(folder, models; specifier="6EoNs_full_nomulti")

crs8 = forecast(crs6, crs7)
crs9 = forecast(crs7, crs8)

crs8_sample = deepcopy(crs7)
crs9_sample = deepcopy(crs7)

begin
fig, ax = plt.subplots(3,3,figsize=(8, 6), sharex=true, sharey=true)
for a in ax
    #a.axhline(1e-2, ls=":", color="lightgrey", lw=1)
    a.grid(ls=":", axis="y" , lw=1, alpha=0.5, color="darkgrey", zorder=1)
end

ax[2,1].axvline(5/3, ls="--", lw=0.8, color="darkgrey", zorder=2)
ax[2,1].axvline(2/3, ls="--", lw=0.7, color="k", zorder=2, alpha=0.7)
ax[2,1].axvline(8/3, ls="--", lw=0.7, color="k", zorder=2, alpha=0.7)
ax[2,1].text(9/3, 1e-1, L"DFSZ$_2$-I", bbox=Dict("alpha" => 0.8, "pad" => 0.1, "color" => "white", "ec" => "white"))
ax[2,1].text(-21/3, 1e-1, L"DFSZ$_2$-II", bbox=Dict("alpha" => 0.8, "pad" => 0.1, "color" => "white", "ec" => "white"))
ax[2,1].text(11/3, 1.5e-2, "symmetry\n axis 5/3", color="darkgrey", bbox=Dict("alpha" => 0.6, "pad" => 0.1, "color" => "white", "ec" => "white"))


h8, m = _make_hist(crs8)
h8 = normalize(h8; mode=:probability)
ax[3,2].stairs(h8.weights, h8.edges[1], lw=2, ec="mediumseagreen", alpha=0.4)
h8s, m = _make_hist(crs8_sample)
h8s = normalize(h8s; mode=:probability)
ax[3,2].stairs(h8s.weights, h8s.edges[1], lw=2, ec="maroon", alpha=0.4)

h9, m = _make_hist(crs9)
h9 = normalize(h9; mode=:probability)
ax[3,3].stairs(h9.weights, h9.edges[1], lw=2, ec="mediumseagreen", alpha=0.4)
h9s, m = _make_hist(crs8_sample)
h9s = normalize(h9s; mode=:probability)
ax[3,3].stairs(h9s.weights, h9s.edges[1], lw=2, ec="maroon", alpha=0.4)


slist= [L"n_D = 3", L"n_D = 4", L"n_D = 5", L"n_D = 6", L"n_D = 7", L"n_D = 8", L"n_D = 9"]

ax[1,1].stairs(hs[1].weights, hs[1].edges[1], lw=2, ec="maroon", alpha=0.4)
ax[1,1].set_title(L"..."*slist[1])
ax[2,1].stairs(hs[2].weights, hs[2].edges[1], lw=2, ec="maroon", alpha=0.4, zorder=3)
ax[2,1].set_title(L"..."*slist[2])
ax[2,2].stairs(hs[3].weights, hs[3].edges[1], lw=2, ec="maroon", alpha=0.4)
ax[2,2].set_title(L"..."*slist[3])
ax[2,3].stairs(hs[4].weights, hs[4].edges[1], lw=2, ec="maroon", alpha=0.4)
ax[2,3].set_title(L"..."*slist[4])
ax[3,1].stairs(hs[5].weights, hs[5].edges[1], lw=2, ec="maroon", alpha=0.4)
ax[3,1].set_title(L"..."*slist[5])
ax[3,2].set_title(L"..."*slist[6])
ax[3,3].set_title(L"..."*slist[7])



ax[3,1].set_xlabel("Anomaly Ratio E/N")
ax[3,2].set_xlabel("Anomaly Ratio E/N")
ax[3,3].set_xlabel("Anomaly Ratio E/N")
ax[1,1].set_ylabel("Probability")
ax[2,1].set_ylabel("Probability")
ax[3,1].set_ylabel("Probability")
ax[1,2].set_xlim([-10,13])
#ax[ab[i][1],ab[i][2]].axhline([2e-2])
ax[1,2].set_yscale("log")
ax[1,2].set_ylim([1e-4,3e-1])
fig.delaxes(ax[1,2])
fig.delaxes(ax[1,3])

plt.suptitle("DFSZ-type models with ...")
plt.savefig("plots/220914-paperplots/PDFn3to9.pdf")
end






#########################################################################################################
#
# Plot n=3 to n=9 for minimal potentials
#
#########################################################################################################

folder = "220824-seriousruns"
# Load n=3 to n=7
hs = similar(3:7, Any)
for nD in 3:7
    models, multis = model_list(nD=[nD])
    if nD < 7 
        tmp = read_EoN(folder, models; specifier=string(nD)*"EoNs_full_wmulti")
    else
        tmp = read_EoN(folder, models; specifier=string(nD)*"EoNs_full_wmulti_sample")
    end
    hs[nD-2], m = _make_hist(tmp)
    hs[nD-2] = normalize(hs[nD-2]; mode=:probability)
end

# Construct upper and lower limit for n=8 and n=9
models, multis = model_list(nD=[8])
crs8 = read_EoN(folder, models; specifier="8EoNs_full_wmulti_sample")
models, multis = model_list(nD=[9])
crs9 = read_EoN(folder, models; specifier="9EoNs_full_wmulti_sample")


begin
fig, ax = plt.subplots(3,3,figsize=(8, 6), sharex=true, sharey=true)
for a in ax
    #a.axhline(1e-2, ls=":", color="lightgrey", lw=1)
    a.grid(ls=":", axis="y" , lw=1, alpha=0.5, color="darkgrey")
end

h8, m = _make_hist(crs8)
h8 = normalize(h8; mode=:probability)
ax[3,2].stairs(h8.weights, h8.edges[1], lw=2, ec="maroon", alpha=0.4)


h9, m = _make_hist(crs9)
h9 = normalize(h9; mode=:probability)
ax[3,3].stairs(h9.weights, h9.edges[1], lw=2, ec="maroon", alpha=0.4)


slist= [L"n_D = 3", L"n_D = 4", L"n_D = 5", L"n_D = 6", L"n_D = 7", L"n_D = 8", L"n_D = 9"]

ax[1,1].stairs(hs[1].weights, hs[1].edges[1], lw=2, ec="maroon", alpha=0.4)
ax[1,1].set_title(L"..."*slist[1])
ax[2,1].stairs(hs[2].weights, hs[2].edges[1], lw=2, ec="maroon", alpha=0.4)
ax[2,1].set_title(L"..."*slist[2])
ax[2,2].stairs(hs[3].weights, hs[3].edges[1], lw=2, ec="maroon", alpha=0.4)
ax[2,2].set_title(L"..."*slist[3])
ax[2,3].stairs(hs[4].weights, hs[4].edges[1], lw=2, ec="maroon", alpha=0.4)
ax[2,3].set_title(L"..."*slist[4])
ax[3,1].stairs(hs[5].weights, hs[5].edges[1], lw=2, ec="maroon", alpha=0.4)
ax[3,1].set_title(L"..."*slist[5])
ax[3,2].set_title(L"..."*slist[6])
ax[3,3].set_title(L"..."*slist[7])



ax[3,1].set_xlabel("Anomaly Ratio E/N")
ax[3,2].set_xlabel("Anomaly Ratio E/N")
ax[3,3].set_xlabel("Anomaly Ratio E/N")
ax[1,1].set_ylabel("Probability")
ax[2,1].set_ylabel("Probability")
ax[3,1].set_ylabel("Probability")
ax[1,2].set_xlim([-10,13])
#ax[ab[i][1],ab[i][2]].axhline([2e-2])
ax[1,2].set_yscale("log")
ax[1,2].set_ylim([1e-4,3e-1])
fig.delaxes(ax[1,2])
fig.delaxes(ax[1,3])

plt.suptitle("DFSZ-type models with ...")
plt.savefig("plots/220914-paperplots/PDFn3to9-minimalpotentials.pdf")
end









#########################################################################################################
#
# Make comparison plot E/N pdf of DFSZ vs KSVZ
#
#########################################################################################################

KSVZ_ARs, KSVZgag, n_dw = ksvz("all"; edges=-50:0.001:50)
KSVZcdf = cdf(KSVZ_ARs)
KSVZAR = normalize(KSVZ_ARs; mode=:probability)

KSVZ_ARs

folder = "220824-seriousruns"

# Construct upper and lower limit for n=8 and n=9
models, multis = model_list(nD=[3])
crs3 = read_EoN(folder, models; specifier="3EoNs_full_nomulti")
models, multis = model_list(nD=[4])
crs4 = read_EoN(folder, models; specifier="4EoNs_full_nomulti")
models, multis = model_list(nD=[5])
crs5 = read_EoN(folder, models; specifier="5EoNs_full_nomulti")
models, multis = model_list(nD=[6])
crs6 = read_EoN(folder, models; specifier="6EoNs_full_nomulti")
models, multis = model_list(nD=[7])
crs7 = read_EoN(folder, models; specifier="7EoNs_full_nomulti")


crs8 = forecast(crs6, crs7)
crs9 = forecast(crs7, crs8)

crs8_sample = deepcopy(crs7)
crs9_sample = deepcopy(crs7)


crs3 = normalize_cm(crs3)
crs4 = normalize_cm(crs4)
crs5 = normalize_cm(crs5)
crs6 = normalize_cm(crs6)
crs7 = normalize_cm(crs7)
crs8 = normalize_cm(crs8)
crs9 = normalize_cm(crs9)
crs8_sample = normalize_cm(crs8_sample)
crs9_sample = normalize_cm(crs9_sample)

crs8

crs_lower = merge(+, crs3, crs4, crs5, crs6, crs7, crs8, crs9)
crs_upper = merge(+, crs3, crs4, crs5, crs6, crs7, crs8_sample, crs9_sample)

hl, m = _make_hist(crs_lower, bins=-50:0.0001:50)
hl = normalize(hl; mode=:probability)
hu, m = _make_hist(crs_upper, bins=-50:0.0001:50)
hu = normalize(hu; mode=:probability)


fig, ax = plt.subplots(figsize=(7, 4))
ax.stairs(KSVZAR.weights, KSVZAR.edges[1], linewidth=4000, lw=1.5, ec="mediumseagreen", alpha=1, label="KSVZ-type (all)")
ax.stairs(hu.weights, hu.edges[1], linewidth=4000, lw=1.5, ec="maroon", alpha=0.3, label="DFSZ-type (all)")
ax.stairs(hl.weights, hl.edges[1], linewidth=4000, lw=1.5, ec="maroon", alpha=0.3, label="")

#ax.step(myAR.edges[1][2:end], myAR.weights)
plt.xlim([-50,50])
plt.yscale("log")
plt.ylim([1e-5,3e-1])
plt.xlabel("Anomaly Ratio E/N")
plt.ylabel("Probability")
plt.legend(loc="best")
plt.savefig("plots/test/PDFcompare.pdf")



#########################################################################################################
#
# Make comparison plot E/N pdf of DFSZ vs KSVZ (zoomed)
#
#########################################################################################################
fig, ax = plt.subplots(figsize=(7, 4))
ax.stairs(KSVZAR.weights, KSVZAR.edges[1],  lw=1.5, ec="mediumseagreen", color="mediumseagreen", alpha=1, fill=true, label="KSVZ-type (all)")
ax.stairs(hu.weights, hu.edges[1], lw=1.5, ec="maroon", alpha=0.3, label="DFSZ-type (all)")
ax.stairs(hl.weights, hl.edges[1], lw=1.5, ec="maroon", alpha=0.3, label="")

plt.yscale("log")
plt.ylim([1e-5,3e-1])
plt.xlabel("Anomaly Ratio E/N")
plt.ylabel("Probability")
plt.legend(loc="best")
plt.xlim([-0,3])
plt.savefig("plots/220914-paperplots/PDFcompareZoom.pdf")




#########################################################################################################
#
# Make comparison plot E/N cdf of DFSZ vs KSVZ
#
#########################################################################################################

KSVZ_ARs, KSVZgag, n_dw = ksvz("all"; edges=-500:0.0001:500, Cagdist=true)
KSVZcdf = cdf(KSVZgag)

10^(KSVZgag.edges[1][findfirst(KSVZcdf .< 0.68)-1]) * 2 * π * fa(40e-6) / αem()
10^(KSVZgag.edges[1][findfirst(KSVZcdf .< 0.95)-1]) * 2 * π * fa(40e-6) / αem()

hgag = gag_histogram(hl, edges=-17:0.000001:-12, Cagdist=true)
cdfgag = cdf(hgag)

limit_diff(0.68, hgag, KSVZgag)

gaγγ(8/3, fa(40e-6))

10^( (1.92 - (8/3))) * 2 * π * fa(40e-6) / αem()
10^(hgag.edges[1][findfirst(cdfgag .< 0.95)-1]) * 2 * π * fa(40e-6) / αem()
hgag.edges[1]

fig, ax = plt.subplots(figsize=(4.8, 3.5))

plt.fill_between(x=hgag.edges[1][2:end], y1=cdfgag, where=(hgag.edges[1][2:end] .> hgag.edges[1][findfirst(cdfgag .< 0.95)-1]), alpha=0.4, color="darkgrey")
plt.fill_between(x=KSVZgag.edges[1][2:end], y1=KSVZcdf, where=(KSVZgag.edges[1][2:end] .> KSVZgag.edges[1][findfirst(KSVZcdf .< 0.95)-1]), alpha=0.4, color="darkgrey")
plt.fill_between(x=hgag.edges[1][2:end], y1=cdfgag, where=(hgag.edges[1][2:end] .> hgag.edges[1][findfirst(cdfgag .< 0.68)-1]), alpha=0.8, color="darkgrey")
plt.fill_between(x=KSVZgag.edges[1][2:end], y1=KSVZcdf, where=(KSVZgag.edges[1][2:end] .> KSVZgag.edges[1][findfirst(KSVZcdf .< 0.68)-1]), alpha=0.8, color="darkgrey")

ax.set_ylim([-0.05,1.05])
ymin, ymax = ax.get_ybound()
yrange = ymax - ymin
ax.axvline(hgag.edges[1][findfirst(cdfgag .< 0.68)], (0- ymin)/yrange, (cdfgag[findfirst(cdfgag .< 0.68)-1] - ymin)/yrange, ls="--", c="maroon", alpha=0.4)
ax.axvline(KSVZgag.edges[1][findfirst(KSVZcdf .< 0.68)], (0- ymin)/yrange, (KSVZcdf[findfirst(KSVZcdf .< 0.68)-1] - ymin)/yrange, ls="--", c="mediumseagreen")
ax.axvline(hgag.edges[1][findfirst(cdfgag .< 0.95)], (0- ymin)/yrange, (cdfgag[findfirst(cdfgag .< 0.95)-1] - ymin)/yrange, ls="--", c="maroon", alpha=0.4)
ax.axvline(KSVZgag.edges[1][findfirst(KSVZcdf .< 0.95)], (0- ymin)/yrange, (KSVZcdf[findfirst(KSVZcdf .< 0.95)-1] - ymin)/yrange, ls="--", c="mediumseagreen")

ax.axvline(gaγγ(8/3, fa(40e-6)), ls=":", c="k")
ax.axvline(gaγγ(2/3, fa(40e-6)),  ls=":", c="k")
ax.text(-13.9, 0.97, L"DFSZ$_2$-II", bbox=Dict("alpha" => 0.0, "pad" => 0.1, "color" => "white", "ec" => "white"))
ax.text(-14.95, 0.97, L"DFSZ$_2$-I", bbox=Dict("alpha" => 0.0, "pad" => 0.1, "color" => "white", "ec" => "white"))


ax.step(KSVZgag.edges[1][2:end], KSVZcdf; where="pre", lw=3, color="mediumseagreen", label="KSVZ-type (all)")
ax.step(hgag.edges[1][2:end], cdfgag; where="pre", lw=3, color="maroon", alpha=0.4, label="DFSZ-type (all)")
#plt.fill_between(x=KSVZgag.edges[1][2:end], y1=KSVZcdf, where=(KSVZgag.edges[1][2:end] .> KSVZgag.edges[1][findfirst(KSVZcdf .< 0.68)-1]), alpha=0.9, color="darkgrey")
#plt.fill_between(x=KSVZgag.edges[1][2:end], y1=KSVZcdf, where=(KSVZgag.edges[1][2:end] .> KSVZgag.edges[1][findfirst(KSVZcdf .< 0.95)-1]), alpha=0.6, color="lightgrey")
#plt.fill_between(x=hgag.edges[1][2:end], y1=cdfgag, where=(hgag.edges[1][2:end] .> hgag.edges[1][findfirst(cdfgag .< 0.68)-1]), alpha=0.9, color="darkgrey")
#plt.fill_between(x=hgag.edges[1][2:end], y1=cdfgag, where=(hgag.edges[1][2:end] .> hgag.edges[1][findfirst(cdfgag .< 0.95)-1]), alpha=0.6, color="lightgrey")
plt.xlim([-16.5,-12])
plt.legend(loc="lower left")
plt.xlabel(L"\log{g_{a \gamma}[\mathrm{GeV}^{-1}]}\: \: \mathrm{at}\: \: m_a = 40\; \mu\mathrm{eV}")
plt.ylabel("Probability")
plt.savefig("plots/220914-paperplots/CDFcompare-sample.pdf", bbox_inches="tight")




#########################################################################################################
#
# n=3 to n=9 cdf plot
#
#########################################################################################################

cdfs = similar(3:7, Any)
for nD in 3:7
    @info nD
    models, multis = model_list(nD=[nD])
    tmp = read_EoN(folder, models; specifier=string(nD)*"EoNs_full_nomulti")
    hs, m = _make_hist(tmp; bins=-100:0.0001:100)
    hs = normalize(hs; mode=:probability)
    hgag = gag_histogram(hs, edges=-17:0.001:-12, Cagdist=true)
    cdfs[nD-2] = cdf(hgag)  
end
cdfs

slist= [L"n_D = 3", L"n_D = 4", L"n_D = 5", L"n_D = 6", L"n_D = 7"]
clist = ["mediumseagreen", "yellowgreen", "darkkhaki", "darkgoldenrod", "maroon"]
alphalist = [1.0, 0.8, 0.8, 0.8, 0.4]

fig, ax = plt.subplots(figsize=(4.8, 3.5))
for (i, cdff) in enumerate(cdfs)
    ax.step(hgag.edges[1][2:end], cdff; where="pre", alpha=alphalist[i], c=clist[i], lw=3, label=slist[i])
end
plt.xlim([-16.5,-12])
plt.legend(loc="best")
plt.xlabel(L"\log{g_{a \gamma}[\mathrm{GeV}^{-1}]}\: \: \mathrm{at}\: \: m_a = 40\; \mu\mathrm{eV}")
plt.ylabel("Probability")
plt.savefig("plots/220914-paperplots/CDFn3to9.pdf", bbox_inches="tight")


#########################################################################################################
#
# n=3 to n=9 cdf plot
#
#########################################################################################################

using Distributions


models, mmultis = model_list(nD=[4])
e1=read_EoN("220824-seriousruns", models; specifier="4EoNs_full_nomulti")
hs, m = _make_hist(e1; bins=-500:0.0001:500)
hs = normalize(hs; mode=:probability)
hgag4 = gag_histogram(hs, edges=-17:0.001:-12, Cagdist=true)
cdf4 = cdf(hgag4) 

e1=read_EoN("220824-seriousruns", models; specifier="4EoNs_nomulti_NDW1")
hs, m = _make_hist(e1; bins=-500:0.0001:500)
hs = normalize(hs; mode=:probability)
hgag4dw = gag_histogram(hs, edges=-17:0.001:-12, Cagdist=true)
cdf4dw = cdf(hgag4dw)


models, mmultis = model_list(nD=[5])
e1=read_EoN("220824-seriousruns", models; specifier="5EoNs_full_nomulti")
hs, m = _make_hist(e1; bins=-500:0.0001:500)
hs = normalize(hs; mode=:probability)
hgag5 = gag_histogram(hs, edges=-17:0.001:-12, Cagdist=true)
cdf5 = cdf(hgag5) 

e1=read_EoN("220824-seriousruns", models; specifier="5EoNs_nomulti_NDW1")
hs, m = _make_hist(e1; bins=-500:0.0001:500)
hs = normalize(hs; mode=:probability)
hgag5dw = gag_histogram(hs, edges=-17:0.001:-12, Cagdist=true)
cdf5dw = cdf(hgag5dw)


models, mmultis = model_list(nD=[6])
e1=read_EoN("220824-seriousruns", models; specifier="6EoNs_full_nomulti")
hs, m = _make_hist(e1; bins=-500:0.0001:500)
hs = normalize(hs; mode=:probability)
hgag6 = gag_histogram(hs, edges=-17:0.001:-12, Cagdist=true)
cdf6 = cdf(hgag6) 

e1=read_EoN("220824-seriousruns", models; specifier="6EoNs_nomulti_NDW1")
hs, m = _make_hist(e1; bins=-500:0.0001:500)
hs = normalize(hs; mode=:probability)
hgag6dw = gag_histogram(hs, edges=-17:0.001:-12, Cagdist=true)
cdf6dw = cdf(hgag6dw)


fig, ax = plt.subplots(figsize=(4.8, 3.5))
ax.step(hgag4.edges[1][2:end], cdf4; where="pre", alpha=1, c="mediumseagreen", lw=3, label=L"n_D=4\; \mathrm{all}\; N_{DW}")
ax.step(hgag4dw.edges[1][2:end], cdf4dw; where="pre", alpha=0.8, c="darkgoldenrod", lw=3, label=L"n_D=4\; \mathrm{w/}\; N_{DW} = 1")
ax.step(hgag5.edges[1][2:end], cdf5; where="pre", alpha=0.8, c="yellowgreen", lw=3, label=L"n_D=5\; \mathrm{all}\; N_{DW}")
ax.step(hgag5dw.edges[1][2:end], cdf5dw; where="pre", alpha=0.8, c="gold", lw=3, label=L"n_D=5\; \mathrm{w/}\; N_{DW} = 1")
ax.step(hgag6.edges[1][2:end], cdf6; where="pre", alpha=0.4, c="olive", lw=3, label=L"n_D=6\; \mathrm{all}\; N_{DW}")
ax.step(hgag6dw.edges[1][2:end], cdf6dw; where="pre", alpha=0.4, c="maroon", lw=3, label=L"n_D=6\; \mathrm{w/}\; N_{DW} = 1")
plt.xlim([-16.5,-12])
plt.legend(loc="best")
plt.xlabel(L"\log{g_{a \gamma}[\mathrm{GeV}^{-1}]}\: \: \mathrm{at}\: \: m_a = 40\; \mu\mathrm{eV}")
plt.ylabel("Probability")
plt.savefig("plots/220914-paperplots/CDFdw.pdf", bbox_inches="tight")







#########################################################################################################
#
# Generate the numbers in the overview table
#
#########################################################################################################

# #Veb 
# comes from calculation with formula given in paper.

# unique solutions
fid = h5open("./data/DFSZ_models/221007-seriousnewtxtfiles/full_n5.h5")

mychi = Matrix{Float64}(undef,0,10) 
for (i, a) in enumerate(fid)
    k = read(a)
    if i != 1
        println(keys(fid)[i])
        for tuple in k#collect(values(k))
            dat = tuple[2]
            mychi = vcat(mychi, round.(dat["Chis"], digits=5))
        end
    end
end

tt4 = unique(eachrow(mychi))

close(fid)
# for n=6 comes from 220916-serioustxtfiles, but number of solutions didn't change between those two iterations.

# unique EoNs
n = 7
models, mmultis = model_list(nD=[n])
e1=read_EoN("220824-seriousruns", models; specifier=string(n)*"EoNs_full_nomulti")

ke1 = collect(keys(e1))
ve1 = collect(values(e1))
ve1 = ve1[isnan.(ke1) .== 0]
ke1 = ke1[isnan.(ke1) .== 0]


#E/N hcat
minimum(ke1) * 3

# % photophobic
sum(ve1[1.88 .< ke1 .< 1.96]) / sum(ve1)

# % N_DW = 1
n = 5
models, mmultis = model_list(nD=[n])
e1=read_EoN("220824-seriousruns", models; specifier=string(n)*"EoNs_full_nomulti")
e2=read_EoN("220824-seriousruns", models; specifier=string(n)*"EoNs_nomulti_NDW1")

ke1 = collect(keys(e1))
ve1 = collect(values(e1))
ve1 = ve1[isnan.(ke1) .== 0]
ke1 = ke1[isnan.(ke1) .== 0]

ke2 = collect(keys(e2))
ve2 = collect(values(e2))
ve2 = ve2[isnan.(ke2) .== 0]
ke2 = ke2[isnan.(ke2) .== 0]

sum(ve2) / sum(ve1)



#########################################################################################################
#
# Generate axion bands
#
#########################################################################################################

function limit(lim, gaghist)
    mycdf = cdf(gaghist)
    10^(gaghist.edges[1][findfirst(mycdf .< lim)-1]) * 2 * π * fa(40e-6) / αem() #+ 1.92
end

limits = [0.68, 0.95, 0.16, 0.84, 0.025, 0.975]
folder = "220824-seriousruns"

# DiLuzios add one field band:
EoNmin = 5/3
EoNmax = 44/3
diLmin = 1.92 - 1.6666666667
diLmax = 13.66666666667 - 1.92


### Our DFSZ results

# Construct upper and lower limit for n=8 and n=9
models, multis = model_list(nD=[3])
crs3 = read_EoN(folder, models; specifier="3EoNs_full_nomulti")
models, multis = model_list(nD=[4])
crs4 = read_EoN(folder, models; specifier="4EoNs_full_nomulti")
models, multis = model_list(nD=[5])
crs5 = read_EoN(folder, models; specifier="5EoNs_full_nomulti")
models, multis = model_list(nD=[6])
crs6 = read_EoN(folder, models; specifier="6EoNs_full_nomulti")
models, multis = model_list(nD=[7])
crs7 = read_EoN(folder, models; specifier="7EoNs_full_nomulti")


crs8 = forecast(crs6, crs7)
crs9 = forecast(crs7, crs8)

crs8_sample = deepcopy(crs7)
crs9_sample = deepcopy(crs7)

crs3 = normalize_cm(crs3)
crs4 = normalize_cm(crs4)
crs5 = normalize_cm(crs5)
crs6 = normalize_cm(crs6)
crs7 = normalize_cm(crs7)
crs8 = normalize_cm(crs8)
crs9 = normalize_cm(crs9)
crs8_sample = normalize_cm(crs8_sample)
crs9_sample = normalize_cm(crs9_sample)

crs_lower = merge(+, crs3, crs4, crs5, crs6, crs7, crs8, crs9)
crs_upper = merge(+, crs3, crs4, crs5, crs6, crs7, crs8_sample, crs9_sample)

hl, m = _make_hist(crs_lower, bins=-500:0.0001:500)
hl = normalize(hl; mode=:probability)
hu, m = _make_hist(crs_upper, bins=-500:0.0001:500)
hu = normalize(hu; mode=:probability)

hgag = gag_histogram(hu, edges=-17:0.000001:-12, Cagdist=true)

# E/N = 
limit(0.68, hgag)
limit(0.95, hgag)
limit(0.16, hgag)
limit(0.84, hgag)
limit(0.975, hgag)
limit(0.025, hgag)

hgagvec = limit.(limits, Ref(hgag))

### Plakkot KSVZ results

KSVZ_ARs, KSVZgag, n_dw = ksvz("all"; edges=-500:0.0001:500, Cagdist=true)

# E/N = 
limit(0.68, KSVZgag)
limit(0.95, KSVZgag)
limit(0.16, KSVZgag)
limit(0.84, KSVZgag)
limit(0.975, KSVZgag)
limit(0.025, KSVZgag)

KSVZgagvec = limit.(limits, Ref(KSVZgag))

### Combined results

e_n, e_n_counts, n_dw = py"plakkot"("histogram_all_LP_allowed_models")
e_n = convert(Vector{Float64}, e_n)
n_dw = convert(Vector{Float64}, n_dw)
ARcm = countmap(e_n, e_n_counts)
KSVZar = normalize_cm(ARcm)
DFSZar = normalize_cm(crs_upper)

totAR = merge(+, DFSZar, KSVZar)

hu, m = _make_hist(totAR, bins=-500:0.0001:500)
hu = normalize(hu; mode=:probability)

totgag = gag_histogram(hu, edges=-17:0.000001:-12, Cagdist=true)

# E/N = 
limit(0.68, totgag)
limit(0.95, totgag)
limit(0.16, totgag)
limit(0.84, totgag)
limit(0.975, totgag)
limit(0.025, totgag)

totgagvec = limit.(limits, Ref(totgag))


### Combined NDW=1 results

models, multis = model_list(nD=[4])
crs4 = read_EoN(folder, models; specifier="4EoNs_nomulti_NDW1")
models, multis = model_list(nD=[5])
crs5 = read_EoN(folder, models; specifier="5EoNs_nomulti_NDW1")
models, multis = model_list(nD=[6])
crs6 = read_EoN(folder, models; specifier="6EoNs_nomulti_NDW1")

crs4 = normalize_cm(crs4)
crs5 = normalize_cm(crs5)
crs6 = normalize_cm(crs6)
crsndw1 = merge(+, crs4, crs5, crs6)


e_n, e_n_counts, n_dw = py"plakkot"("histogram_all_LP_allowed_models")
e_n = convert(Vector{Float64}, e_n)
n_dw = convert(Vector{Float64}, n_dw)
ARcm = countmap(e_n[n_dw .== 1], e_n_counts[n_dw .== 1])
KSVZar = normalize_cm(ARcm)
DFSZar = normalize_cm(crsndw1)

totAR = merge(+, DFSZar, KSVZar)


hu, m = _make_hist(totAR, bins=-500:0.0001:500)
hu = normalize(hu; mode=:probability)

hdw1gag = gag_histogram(hu, edges=-17:0.000001:-12, Cagdist=true)

# E/N = 
limit(0.68, hdw1gag)
limit(0.95, hdw1gag)
limit(0.16, hdw1gag)
limit(0.84, hdw1gag)
limit(0.975, hdw1gag)
limit(0.025, hdw1gag)

hdw1gagvec = limit.(limits, Ref(hdw1gag))



fold="./data/DFSZ_models/221007-seriousnewtxtfiles/"
dataset="AxionBands"

io = open(fold*dataset*".txt", "w")
    write(io, "# Detailed explanation of different scenarios:\n")
    write(io, "# diLuzio : arXiv [1705.05370], classic axion band between E/N = 44/3 and E/N = 5/3.\n")
    write(io, "#           Contains all KSVZ models with one additional quark and all DFSZ models with one additional Higgs doublet.\n")
    write(io, "#           Not a band in a statistical sense, rather a region of the parameterspace where several DISCREET models lie. Use when interested in minimal extensions.\n")
    write(io, "# KSVZ    : arXiv [2107.12378], bands and one-sided limits derived from Plakkot&Hoof.\n")
    write(io, "#           Contains complete set of preferred KSVZ models with up to nine additional quark, assuming equal probability for all models.\n")
    write(io, "#           Use when specifically interested in KSVZ models beyond minimal extensions.\n")
    write(io, "# DFSZ    : this work, bands and one-sided limits.\n")
    write(io, "#           Contains DFSZ models with up to nine higgs doublets, extrapolating beyond DFSZ₇, assuming equal probability for all numbers of Higgs doublets.\n")
    write(io, "#           Contrary to KSVZ case no phenomenological selection criteria were applied.\n")
    write(io, "#           Use when specifically interested in DFSZ models beyond minimal extensions.\n")
    write(io, "# Combined: bands and one-sided limits for KSVZ and DFSZ case combined assuming equal probability for a DFSZ axion and for a KSVZ axion.\n")
    write(io, "#           Use if you have no preference for DFSZ or KSVZ and want to consider beyond minimal extensions.\n")
    write(io, "# N_DW=1  : bands and one-sided limits for KSVZ and DFSZ case combined assuming equal probability for a DFSZ axion and for a KSVZ axion, only considering models with domain wall number of one.\n")
    write(io, "#           DFSZ case only includes DFSZ₂ to DFSZ₆. No extrapolation beyond.\n")
    write(io, "#           Use if you have no preference for DFSZ or KSVZ, want to consider beyond minimal extensions, are interested in the post-inflationary scenario and believe that the domain wall problem is indeed a problem.\n")
    write(io, "#\n")
    write(io, "# All values given are |E/N - 1.92| corresponding to the presented limits.\n")
    write(io, "# All limits given are one-sided. For a central 68% axion band use 0.16 and 0.84 limits, or a central 95% axion band use 0.975 and 0.025.\n")
    write(io, "#\n")
    write(io, "#            minimal value      maximal value\n")
    write(io, "diLuzio         "*rpad(round(diLmin,digits=4),15)*rpad(round(diLmax,digits=4),5)*"\n")
    write(io, "#\n")
    write(io, "# Limit:      0.68     0.95     0.16     0.84    0.025    0.975\n")
    write(io, "KSVZ     "*([lpad(round(val,digits=4),9) for val in KSVZgagvec]...)*"\n")
    write(io, "DFSZ     "*([lpad(round(val,digits=4),9) for val in hgagvec]...)*"\n")
    write(io, "Combined "*([lpad(round(val,digits=4),9) for val in totgagvec]...)*"\n")
    write(io, "N_DW=1   "*([lpad(round(val,digits=4),9) for val in hdw1gagvec]...)*"\n")
    close(io)