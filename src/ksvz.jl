# Manipulate Plakkot Hoof data in the necessary way

using PyCall

py"""
import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction

# Select histogram file to read
def plakkot(str):
    e_num, e_denom, n_num, n_denom, counts = np.genfromtxt('./data/KSVZ-hists/'+str+'.txt',
                                                        delimiter='\t', unpack=True, dtype='i')

    cond = (n_num != 0)
    e_n = np.array([Fraction(en,ed)/Fraction(nn,nd) for en,ed,nn,nd in zip(e_num[cond],e_denom[cond],n_num[cond],n_denom[cond])])
    e_n_counts = counts[cond]
    n_dw = np.array([2*Fraction(nn,nd) for nn,nd in zip(n_num[cond],n_denom[cond])])
    return e_n, e_n_counts, n_dw

"""
e_n, e_n_counts, n_dw = py"plakkot"("histogram_additive_LP_allowed_models")

e_n = convert(Vector{Float64}, e_n)
n_dw = convert(Vector{Float64}, n_dw)

KSVZ_ARs = fit(Histogram, e_n, FrequencyWeights(e_n_counts), -10:0.01:50)

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


    ###############################################################
# Make comparison CDF
e_n, e_n_counts, n_dw = py"plakkot"("histogram_all_LP_allowed_models")

e_n = convert(Vector{Float64}, e_n)
n_dw = convert(Vector{Float64}, n_dw)

KSVZ_ARs = fit(Histogram, e_n, FrequencyWeights(e_n_counts), -50:0.01:50)

KSVZgag = gag_histogram(KSVZ_ARs, mode=:probability, edges=-17:0.001:-12)
KSVZgag = normalize(KSVZgag; mode=:probability)
KSVZcdf = gag_cdf(KSVZgag)

plot(KSVZgag.edges[1][1:end-1], KSVZcdf, title="Axion model CDF KSVZ vs DFSZ comparison", 
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
        gagh = gag_histogram(tt; mode=:probability, edges=-17:0.001:-12)
        hist_list[i] = gagh
    end
    gaghs[k-2] = merge(hist_list...)
end

gagtot = merge(gaghs...)

gagtot = normalize(gagtot; mode=:probability)
gagtotcdf = gag_cdf(gagtot)

plot!(gagtot.edges[1][1:end-1], gagtotcdf, lw=3, label="DFSZ-like models")
savefig("./plots/preliminary2/CDFcompare.pdf")