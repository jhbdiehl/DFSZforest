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


models, multis = model_list(nD=[5])
crs1 = read_EoN("220824-seriousruns", models; specifier="5EoNs_full_nomulti")

plot()
plot_EoN(crs1)

models, multis = model_list(nD=[6])
crs2 = read_EoN("220824-seriousruns", models; specifier="6EoNs_full_nomulti")
plot_EoN(crs2)

models, multis = model_list(nD=[7])

e1 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti1_2.jld2", "7/n7_u1_u2_u3_d1_d2_d3_l1_l1_l1")
e2 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti1_2.jld2", "7/n7_u1_u2_u3_d1_d1_d3_l1_l1_l3")
e3 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti3.jld2", "7/n7_u1_u1_u3_d1_d2_d3_l1_l1_l3")
e4 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti4.jld2", "7/n7_u1_u2_u3_d1_d1_d1_l1_l2_l3")
e5 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti5.jld2", "7/n7_u1_u1_u3_d1_d1_d3_l1_l2_l3")
#for e in [e1,e2,e3,e4]
#    e = Dict(round.(collect(keys(e)), digits=6) .=> collect(values(e)))
#end
crs3 = merge(+, e1,e2,e3,e4, e5)
models, multis = model_list(nD=[8])
crs4 = read_EoN("220824-seriousruns", models; specifier="8EoNs_full_nomulti_sample_short")
plot_EoN(crs3)
plot_EoN(crs4)

@info "Note very important: merge(-, a,b) does not take negative values if entry in b is not in a! Therefore only subtract if all entries of b are also in a. Otherwise catastrophe!"
crs1s = Dict(collect(keys(crs1)) .=> collect(values(crs1)) ./ sum(collect(values(crs1))))
crs2s = Dict(collect(keys(crs2)) .=> collect(values(crs2)) ./ sum(collect(values(crs2))))
diff1ss = merge(-, crs2s, crs1s)

crs3s = Dict(collect(keys(crs3)) .=> collect(values(crs3)) ./ sum(collect(values(crs3))))
diff2ss = merge(-, crs3s, crs2s)

crs4s = Dict(collect(keys(crs4)) .=> collect(values(crs4)) ./ sum(collect(values(crs4))))

plot()
plot_EoN(diff1ss)
plot_EoN(diff2ss)

hi1 = _make_hist(diff1ss)
hi2 = _make_hist(diff2ss)
plot()
plot!(hi1[2], hi1[1].weights .* -1; yaxis=(:identity, [0.00000001, :auto]))
plot!(hi2[2], hi2[1].weights .* -1; yaxis=(:identity, [-0.0001, :auto]))

plot()
plot_EoN(crs3s)
plot_EoN(f3)
diff1ss
diff1ssm = Dict(collect(keys(diff1ss)) .=> -1. .* collect(values(diff1ss)))
diff2ssm = Dict(collect(keys(diff2ss)) .=> -1. .* collect(values(diff2ss)))
tmp = merge(+, diff1ssm, crs1s)
f3 = merge(-, tmp, diff1ssm, diff1ssm, diff1ssm)




###########################################################################################################################
#
# Show how analogous extrapolation would work in case where sampling is unbiased:
#
###########################################################################################################################




###########################################################################################################################
#
# Use sampling and above considerations to make a lower and upper limit:
#
###########################################################################################################################

import PyPlot
const plt = PyPlot


models, multis = model_list(nD=[5])
crs1 = read_EoN("220824-seriousruns", models; specifier="5EoNs_full_nomulti")

models, multis = model_list(nD=[6])
crs2 = read_EoN("220824-seriousruns", models; specifier="6EoNs_full_nomulti")

e1 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti1_2.jld2", "7/n7_u1_u2_u3_d1_d2_d3_l1_l1_l1")
e2 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti1_2.jld2", "7/n7_u1_u2_u3_d1_d1_d3_l1_l1_l3")
e3 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti3.jld2", "7/n7_u1_u1_u3_d1_d2_d3_l1_l1_l3")
e4 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti4.jld2", "7/n7_u1_u2_u3_d1_d1_d1_l1_l2_l3")
e5 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti5.jld2", "7/n7_u1_u1_u3_d1_d1_d3_l1_l2_l3")
e6 = FileIO.load("./data/DFSZ_models/220824-seriousruns/7EoNs_full_nomulti6.jld2", "7/n7_u1_u1_u1_d1_d2_d3_l1_l2_l3")

#models, multis = model_list(nD=[7])
#cmp = [e1,e2,e3,e4,e5,e6]
#for i in 1:6
#    save_EoN(models[i], cmp[i]; folder="220824-seriousruns", filename="7EoNs_full_nomulti")
#end

crs3 = merge(+, e1,e2,e3,e4,e5, e6)

crs4 = read_EoN("220824-seriousruns", models; specifier="8EoNs_full_nomulti_sample_short")

@info "Note very important: merge(-, a,b) does not take negative values if entry in b is not in a! Therefore only subtract if all entries of b are also in a. Otherwise catastrophe!"
crs1s = Dict(collect(keys(crs1)) .=> collect(values(crs1)) ./ sum(collect(values(crs1))))
crs2s = Dict(collect(keys(crs2)) .=> collect(values(crs2)) ./ sum(collect(values(crs2))))
diff1ss = merge(-, crs2s, crs1s)
crs3s = Dict(collect(keys(crs3)) .=> collect(values(crs3)) ./ sum(collect(values(crs3))))
diff2ss = merge(-, crs3s, crs2s)
crs4s = Dict(collect(keys(crs4)) .=> collect(values(crs4)) ./ sum(collect(values(crs4))))

diff1ssm = Dict(collect(keys(diff1ss)) .=> -1. .* collect(values(diff1ss)))
diff2ssm = Dict(collect(keys(diff2ss)) .=> -1. .* collect(values(diff2ss)))
f4 = merge(-, crs3s, diff2ssm)




hi3 = _make_hist(crs4s)
hi4 = _make_hist(f4)

crs4sn = normalize(hi3[1]; mode=:probability)
f4n = normalize(hi4[1]; mode=:probability)
tt = [crs4sn.weights, f4n.weights]

fig, ax = plt.subplots(figsize=(7, 4))
ax.stairs(elementwise_f(max, tt), f4n.edges[1], linewidth=4000, lw=1.5, ec="lightgrey", alpha=1, label="upper bound")
ax.stairs(elementwise_f(min, tt), crs4sn.edges[1], linewidth=4000, lw=1.5, ec="darkgrey", alpha=1, label="lower bound")
#ax.step(myAR.edges[1][2:end], myAR.weights)
plt.xlim([-10,13])
plt.yscale("log")
plt.ylim([1e-4,1e-1])
plt.xlabel("Anomaly Ratio E/N")
plt.ylabel("Probability")
plt.legend(loc="best")
plt.savefig("plots/220914-paperplots/220913-upperlowerbound.pdf")


hi1 = _make_hist(crs1s)
crs1sn = normalize(hi1[1]; mode=:probability)
crs1gag = gag_histogram(crs1sn)
crs1cdf = cdf(crs1gag)
crs3gag = gag_histogram(crs3sn)
crs3cdf = cdf(crs3gag)
crs4gag = gag_histogram(crs4sn)
crs4cdf = cdf(crs4gag)
f4gag = gag_histogram(f4n)
f4cdf = cdf(f4gag)

fig, ax = plt.subplots(figsize=(6, 4))
ax.step(f4gag.edges[1][2:end], f4cdf; where="pre", lw=3, color="lightgrey", label="extrapolation")
ax.step(crs3gag.edges[1][2:end], crs3cdf; where="pre", lw=3, color="darkgrey", alpha=0.4, label="sampling")
plt.xlim([-16.5,-12])
plt.legend(loc="best")
plt.xlabel(L"\log{g_{a \gamma}[\mathrm{GeV}^{-1}]}\: \: \mathrm{at}\: \: m_a = 40\; \mu\mathrm{eV}")
plt.ylabel("Probability")
plt.savefig("plots/test/220913-upperlowercdf.pdf")

# In the end it doesnt matter!
limit_diff(0.68, f4gag, crs1gag)


