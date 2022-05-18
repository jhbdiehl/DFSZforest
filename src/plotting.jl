using StaticArrays
using Symbolics
using LaTeXStrings, Plots
using FileIO
using LinearAlgebra

include("helpers.jl")
@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 

model = fname2model(readdir("./data/DFSZ_models/preliminary2/n3")[1])
fname=model2string(model)
lab = fname[1:11]*"\n    "*fname[12:20]*"\n    "*fname[21:end]
tt = read_AR(model; folder="/preliminary2/n"*string(3)*"/")
tt = normalize(tt; mode=:probability)






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