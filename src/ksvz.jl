# Manipulate Plakkot Hoof data in the necessary way

using PyCall
using Random, Statistics, StatsBase

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

function ksvz(str; edges=-10:0.01:50)
    if str ∈ ["additive", "all", "same_reps"]
        e_n, e_n_counts, n_dw = py"plakkot"("histogram_"*str*"_LP_allowed_models")

        e_n = convert(Vector{Float64}, e_n)
        n_dw = convert(Vector{Float64}, n_dw)
        KSVZ_ARs = fit(Histogram, e_n, FrequencyWeights(e_n_counts), edges)
        return KSVZ_ARs, n_dw
    else
        error("Can only read three different histograms from Plakkot: additive, all and same_reps")
    end
end