using StaticArrays
using Symbolics
using StatsPlots
using Distributions
using HDF5
using StatsBase

import PyPlot
const plt = PyPlot


include("helpers.jl")

@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 

# More different possible values for charges -> less peaked structure
# Restrict to positive charges -> very similar to additive models of Plakkot Hoof
# KSVZ introduces more variation in E and N values due to quark representations not fixed
# If mean(charges) = 0 this ensures median(EoNs) = 5/3. However for n=5 it seems to magically work without!

n = 1000000

charges = rand(-5:10:5, n, 9)

EoNd = EoNhist(charges,n)
EoNd = EoNd .+rand(Normal(0.0, 0.04), length(EoNd))
plot()
histogram(EoNd, bins=-0:0.0001:5)

charges = rand(Uniform(-5,5), n, 9)

EoNc = EoNhist(charges,n)

histogram(rand(Normal(0.0, 0.04), n))
EoNc = EoNc .+rand(Normal(0.0, 0.04), length(EoNc))

histogram(EoNc, bins=-0:0.00001:5)

mean(EoN)

hd, md = _make_hist(countmap(EoNd), bins=-100:0.0001:100)#, bins=-18:0.001:-12)
hc, mc = _make_hist(countmap(EoNc), bins=-100:0.0001:100)

hug = gag_histogram(hc; edges=-17:0.0001:-12, Cagdist=true)
hlg = gag_histogram(hd; edges=-17:0.0001:-12, Cagdist=true)

cdfc = cdf(hug)
cdfd = cdf(hlg)
hc = hlg
hd = hug

limit_diff(0.68, hlg, hug)

limit_diff(0.68, hc, hd)

ac = 10^hc.edges[1][findfirst(cdfc .< 0.95)]
ad = 10^hd.edges[1][findfirst(cdfd .< 0.95)]
2 * abs(ac-ad) / (ac+ad)
cdfd[findfirst(cdfd .< 0.68) - 1]


fig, ax = plt.subplots(figsize=(4.8, 3.5))
plt.fill_between(x=hd.edges[1][2:end], y1=cdfd, where=(hd.edges[1][2:end] .> hd.edges[1][findfirst(cdfd .< 0.95)-1]), alpha=0.4, color="darkgrey")
plt.fill_between(x=hc.edges[1][2:end], y1=cdfc, where=(hc.edges[1][2:end] .> hc.edges[1][findfirst(cdfc .< 0.95)-1]), alpha=0.4, color="darkgrey")
plt.fill_between(x=hd.edges[1][2:end], y1=cdfd, where=(hd.edges[1][2:end] .> hd.edges[1][findfirst(cdfd .< 0.68)-1]), alpha=0.8, color="darkgrey")
plt.fill_between(x=hc.edges[1][2:end], y1=cdfc, where=(hc.edges[1][2:end] .> hc.edges[1][findfirst(cdfc .< 0.68)-1]), alpha=0.8, color="darkgrey")

ax.set_ylim([-0.05,1.05])
ymin, ymax = ax.get_ybound()
yrange = ymax - ymin
ax.axvline(hd.edges[1][findfirst(cdfd .< 0.68)], (0- ymin)/yrange, (cdfd[findfirst(cdfd .< 0.68)-1] - ymin)/yrange, ls="--", c="maroon")
ax.axvline(hc.edges[1][findfirst(cdfc .< 0.68)], (0- ymin)/yrange, (cdfc[findfirst(cdfc .< 0.68)-1] - ymin)/yrange, ls="--", c="mediumseagreen")
ax.axvline(hd.edges[1][findfirst(cdfd .< 0.95)], (0- ymin)/yrange, (cdfd[findfirst(cdfd .< 0.95)-1] - ymin)/yrange, ls="--", c="maroon")
ax.axvline(hc.edges[1][findfirst(cdfc .< 0.95)], (0- ymin)/yrange, (cdfc[findfirst(cdfc .< 0.95)-1] - ymin)/yrange, ls="--", c="mediumseagreen")

ax.step(hc.edges[1][2:end], cdfc; where="pre", lw=3, color="mediumseagreen", label="continuous")
ax.step(hd.edges[1][2:end], cdfd; where="pre", lw=3, color="maroon", alpha=1, label="discreet")

plt.xlim([-16.5,-12])
plt.legend(loc="best")
plt.xlabel(L"\log{g_{a \gamma}[\mathrm{GeV}^{-1}]}\: \: \mathrm{at}\: \: m_a = 40\; \mu\mathrm{eV}")
plt.ylabel("Probability")
plt.savefig("plots/220914-paperplots/CDFcomparedvc-sample.pdf", bbox_inches="tight")



function EoNhist(charges, n)
    EoNs = Array{Float64}(undef, n)
    Ns = Array{Float64}(undef, n)

    for i in eachindex(EoNs)
        EoNs[i] = EoverN(charges[i,:]...)
        Ns[i] = N(charges[i,:]...)
    end

    EoNs = EoNs[abs.(Ns) .> 0.0000001]
    @info "Median value: $(median(EoNs))"

    h, m = _make_hist(countmap(EoNs))
    return EoNs
end


c1 = "mediumseagreen"
c2 = "maroon"

fig, ax = plt.subplots(4,2,figsize=(5, 7), sharex=true)

models_list = [rand(-1:0.00001:1, n, 9), rand(-10:0.5:10, n, 9), rand(-5:1:5, n, 9), rand(0:1:10, n, 9)]

for (ind, models) in enumerate(models_list)
    h = EoNhist(models, n)
    #mh, mm = _make_hist(countmap(vec(models)))
    #mh = normalize(mh, mode=:probability)
    #h = normalize(h, mode=:pdf)

    ax[ind,2].axvline(x=5/3, c="gray", ls="--", lw=1)
    ax[ind,1].hist(vec(models), lw=2, ec=c2, alpha=1, bins=500, range=(-10,10))
    #ax[ind,2].stairs(h.weights, h.edges[1], lw=2, ec=c2, alpha=0.4)
    ax[ind,2].hist(h, lw=1, ec=c2, alpha=1, bins=5000, range=(-10,10))
    ax[ind,2].set_yscale("log")
    ax[ind,2].set_ylim([1e1,6e4])
    ax[ind,2].set_yticklabels([]) 
    ax[ind,1].set_ylabel("Probability (a.u.)")
    ax[ind,1].set_yticklabels([]) 
    
    #ax[ind,1].set_xlim([-10,10])
    #ax[ind,2].set_xlim([-10,10])
end

#ax[1,1].set_yscale("log")
#ax[1,1].set_ylim([8e-8,2e-5])
ax[1,1].set_xlim([-10.5,10.5])
#ax[1,2].hist(rand(Normal(0,2), 9*n), bins=500, range=(-10,10))
ax[length(models_list),2].set_xlabel("Distribution of Anomaly Ratios")
ax[length(models_list),1].set_xlabel("Distribution of Charges")
#plt.suptitle("Influence of charge distribution on E/N")
plt.savefig("plots/220901-EoNtests/EoNtest.pdf", bbox_inches="tight")