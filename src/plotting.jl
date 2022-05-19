using StaticArrays
using Symbolics
using LaTeXStrings, Plots
using FileIO
using LinearAlgebra

include("helpers.jl")
include("ksvz.jl")
@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 

model = fname2model(readdir("./data/DFSZ_models/preliminary2/n3")[1])
fname=model2string(model)
lab = fname[1:11]*"\n    "*fname[12:20]*"\n    "*fname[21:end]
tt = read_AR(model; folder="/preliminary2/n"*string(3)*"/")
tt = normalize(tt; mode=:probability)

# Read KSVZ anomaly ratios from Plakkots data.
KSVZ_ARs, n_dw = ksvz("all")

function get_data(file)
    model = fname2model(file)
end

gaghs = similar(3:9, Any)
for k in 3:9
    @info "$k"
    files = readdir("./data/DFSZ_models/preliminary2/n"*string(k))
    hist_list = similar(files, Any)
    @time for (i, file) in enumerate(files)
        model = fname2model(file)

        tt = read_AR(model; folder="/preliminary2/n"*string(k)*"/")
        tt = normalize(tt; mode=:probability)
#        p1 = plot(tt, lt=:stepbins, label="", title="$(file[1:end-5])", 
#            xlabel="E/N", ylabel="Probability", xrange=(-10,13),
#            bottom_margin=2Plots.mm, legend=:topright,
#            size=(400,300), lw=2)
#        savefig(p1, "plots/preliminary2/ARs/$(file[1:end-5])_ARs.pdf")
        gagh = gag_histogram(tt; mode=:probability)
        hist_list[i] = gagh
    end
    gaghs[k-2] = merge(hist_list...)
end

gaghs = normalize.(gaghs; mode=:probability)
gagscdf = gag_cdf.(gaghs)

p1 = plot(gaghs[2], lt=:stepbins, label=lab, title="DFSZ axion model PDF", 
    xlabel=L"ga\gamma\gamma \;\; [\log\;\mathrm{GeV}^{-1}]", ylabel="Probability",
    bottom_margin=2Plots.mm, legend=:topright,
    size=(400,300), lw=2)
#savefig(p1, "plots/test.pdf")

plot()
for gagh in gaghs[2:end]
    p2 = plot!(gagh.edges[1][1:end-1], gag_cdf(gagh), label="", title="DFSZ axion model CDF", 
        xlabel=L"ga\gamma\gamma \;\; [\log\;\mathrm{GeV}^{-1}]", ylabel="Probability for bigger gaγγ",
        bottom_margin=2Plots.mm, legend=:topright,
        size=(400,300), lw=2)
end
p3 = plot!()







KSVZ_ARs = fit(Histogram, abs.(e_n .- 1.92), FrequencyWeights(e_n_counts), -1:0.01:20)

KSVZgag = gag_histogram(KSVZ_ARs, mode=:probability, edges=-17:0.001:-12)
KSVZgag = normalize(KSVZgag; mode=:probability)
KSVZcdf = gag_cdf(KSVZ_ARs)

KSVZ_ARs.edges[1][1:end-1][KSVZcdf .> 0.68][end]
KSVZ_ARs.edges[1][1:end-1][KSVZcdf .> 0.95][end]


plot(KSVZ_ARs.edges[1][1:end-1], KSVZcdf, title="Axion model CDF KSVZ vs DFSZ comparison", 
    xlabel=L"ga\gamma\gamma \;\; [\log\;\mathrm{GeV}^{-1}]", ylabel="Probability for bigger gaγγ",
    bottom_margin=2Plots.mm, legend=:bottomleft,
    size=(550,400), lw=3, label="KSVZ-like models (all)")

#plot(normalize(KSVZ_ARs; mode=:probability), label="", title="DFSZ axion model CDF", 
#    xlabel=L"ga\gamma\gamma \;\; [\log\;\mathrm{GeV}^{-1}]", ylabel="Probability for bigger gaγγ",
#    bottom_margin=2Plots.mm, legend=:topright,
#    size=(400,300), lw=2)

gaghs = similar(3:9, Any)
for k in 3:9
    @info "$k"
    files = readdir("./data/DFSZ_models/preliminary2/n"*string(k))
    hist_list = similar(files, Any)
    @time for (i, file) in enumerate(files)
        model = fname2model(file)

        tt = read_AR(model; folder="/preliminary2/n"*string(k)*"/")
        tt = normalize(tt; mode=:probability)
#        p1 = plot(tt, lt=:stepbins, label="", title="$(file[1:end-5])", 
#            xlabel="E/N", ylabel="Probability", xrange=(-10,13),
#            bottom_margin=2Plots.mm, legend=:topright,
#            size=(400,300), lw=2)
#        savefig(p1, "plots/preliminary2/ARs/$(file[1:end-5])_ARs.pdf")
        ttrescale = abs.(collect(tt.edges...) .+ 1/2 * (tt.edges[1][2] - tt.edges[1][1]) .- 1.92)[1:end-1]
        gagh = fit(Histogram, ttrescale, FrequencyWeights(tt.weights), -1:0.01:20)
        #gagh = gag_histogram(tt; mode=:probability, edges=-17:0.001:-12)
        hist_list[i] = gagh
    end
    gaghs[k-2] = merge(hist_list...)
end

gagtot = merge(gaghs...)

gagtot = normalize(gagtot; mode=:probability)
gagtotcdf = gag_cdf(gagtot)

gagtot.edges[1][1:end-1][gagtotcdf .> 0.68][end]
gagtot.edges[1][1:end-1][gagtotcdf .> 0.95][end]

plot!(gagtot.edges[1][1:end-1], gagtotcdf, lw=3, label="DFSZ-like models")
savefig("./plots/preliminary2/CDFcompare.pdf")




KSVZgag = gag_histogram(KSVZ_ARs, mode=:probability, edges=-18:0.001:-12)
KSVZgag = normalize(KSVZgag; mode=:probability)

plot(KSVZgag, lt=:stepbins, label=lab, title="DFSZ axion model PDF", 
    xlabel=L"ga\gamma\gamma \;\; [\log\;\mathrm{GeV}^{-1}]", ylabel="Probability",
    bottom_margin=2Plots.mm, legend=:topright,
    size=(800,600), lw=2)

vline!([2/3, 5/3,8/3])

KSVZcdf = gag_cdf(KSVZgag)
plot(KSVZgag.edges[1][1:end-1], KSVZcdf, label="", title="DFSZ axion model CDF", 
    xlabel=L"ga\gamma\gamma \;\; [\log\;\mathrm{GeV}^{-1}]", ylabel="Probability for bigger gaγγ",
    bottom_margin=2Plots.mm, legend=:topright,
    size=(400,300), lw=2)
