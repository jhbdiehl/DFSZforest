# Manipulate Plakkot Hoof data in the necessary way

using PyCall
using Random, Statistics, StatsBase

"""
    Note: Datasets from Plakkot&Hoof are not included in this repository due to copyright and have to be downloaded separately from zenodo: https://doi.org/10.5281/zenodo.5091707
"""

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
        

        KSVZgag = Cag_histogram(KSVZ_ARs, edges=-3:0.05:2.5, Cagdist=Cagdist, mode=:probability)
        KSVZgag = normalize(KSVZgag; mode=:probability)

        return KSVZ_ARs, KSVZgag, n_dw
    else
        error("Can only read from Plakkot: additive_LP_allowed_models, all_LP_allowed_models,  same_reps_LP_allowed_models")
    end
end


py"""
import sys
import numpy as np
import h5py as h5
from collections import Counter
from fractions import Fraction

def get_nq_countmaps(Ndw1=False):
    cclist = []
    for nq in range(1, 29):
        print(nq)
        try:
            with h5.File('./data/KSVZ-hists/catalogue_all_LP_allowed_models.h5', 'r') as file:
                e_num = file['NQ_{:d}/E_numerator'.format(nq)][:]
                n_num = file['NQ_{:d}/N_numerator'.format(nq)][:]
                e_den = file['NQ_{:d}/E_denominator'.format(nq)][:]
                n_den = file['NQ_{:d}/N_denominator'.format(nq)][:]
                models = file['NQ_{:d}/representations'.format(nq)][:]
                try:
                    lps = file['NQ_{}/LP'.format(nq)][:]
                except:
                    pass
                if Ndw1==False:
                    cond = (n_num != 0)
                else:
                    cond = (n_num != 0) & (2 * n_num == n_den)
                e_n_ratios = [Fraction(int(en), int(ed))/Fraction(int(nn), int(nd)) for en, ed, nn, nd in zip(e_num[cond], e_den[cond], n_num[cond], n_den[cond])]
                pref_mod = models[cond]
                cc = Counter([float(e_n_ratio) for e_n_ratio in e_n_ratios])
                cclist.append(cc)

        except FileNotFoundError:
            print("File catalogue_all_LP_allowed_models.h5 doesn't exist! You need to uncompress the catalogue tar.gz files first.")
    
    return cclist
"""
