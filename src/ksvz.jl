# Manipulate Plakkot Hoof data in the necessary way

using PyCall
using Random, Statistics, StatsBase

py"""
import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction

# Select histogram file to read FILE PATH HAS TO BE CHANGED!!!
def plakkot(str):
    e_num, e_denom, n_num, n_denom, counts = np.genfromtxt('./data/KSVZ-hists/'+str+'.txt',
                                                        delimiter='\t', unpack=True, dtype='i')

    cond = (n_num != 0)
    e_n = np.array([Fraction(en,ed)/Fraction(nn,nd) for en,ed,nn,nd in zip(e_num[cond],e_denom[cond],n_num[cond],n_denom[cond])])
    e_n_counts = counts[cond]
    n_dw = np.array([2*Fraction(nn,nd) for nn,nd in zip(n_num[cond],n_denom[cond])])
    return e_n, e_n_counts, n_dw

"""

function ksvz(str; edges=-10:0.01:50, Cagdist=false)
    if str ∈ ["additive", "all", "same_reps"]
        e_n, e_n_counts, n_dw = py"plakkot"("histogram_"*str*"_LP_allowed_models")
        println(sum(e_n_counts))

        e_n = convert(Vector{Float64}, e_n)
        n_dw = convert(Vector{Float64}, n_dw)
        
        KSVZ_ARs = fit(Histogram, e_n, FrequencyWeights(e_n_counts), edges)
        

        KSVZgag = gag_histogram(KSVZ_ARs, mode=:probability, edges=-17:0.000001:-12, Cagdist=Cagdist)
        KSVZgag = normalize(KSVZgag; mode=:probability)

        return KSVZ_ARs, KSVZgag, n_dw
    else
        error("Can only read from Plakkot: additive_LP_allowed_models, all_LP_allowed_models,  same_reps_LP_allowed_models")
    end
end