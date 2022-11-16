# DFSZforest

Compute anomaly ratios for different DFSZ-type models. This code has been used for the paper [*DFSZ axions and where to find them*](https://arxiv.com) by [Johannes Diehl](https://scholar.google.com/citations?hl=en&user=BrHSTFwAAAAJ) by [Emmanouil Koutsangelas](). The data files included can also be found separately on [Zenodo](https://zenodo.org).

## Installation

Open a terminal and go to your favourite directory. Then run

```shell
git clone https://github.com/jhbdiehl/DFSZforest.git
cd DFSZforest
```

The project is written mainly in Julia, so you will need to have this programming language installed. Go [here](https://julialang.org/) for more information. To be able to run it from the terminal maybe have a look [here](https://julialang.org/downloads/platform/). Once you have it installed and are able to call it from the terminal, the only thing you have to do to start the example notebook is run the following command on the terminal. 

```shell
julia runnotebook.jl
```

This may take a while when running it for the first time, since it installs all dependencies as well as a contained jupyter installation for running notebooks.

## Usage

First let's clear up some jargon: There are three types of fermions with three generations each: up-, down and lepton type. Up-type fermions would for example be the up, charm and top quarks. Instead of abbreviating different generations differently, we give a generational index. Here and throughout the code terms like `ui` are used to denote the charge of the Higgs doublet coupling to the up-type quarks. E.g. u3 is the charge of $H_{u_3}$, the Higgs which couples to the top quark and l2 the charge of $H_{l_2}$, the Higgs coupling to the muon.

To write specific models, we use a list of symbols, e.g. `[u1,u1,u3,d1,d1,d1,l1,l1,l1]` is a model with four Higgs doublets ($n_D=4$), where all lepton-type fermions couple to one Higgs, all down-type fermions couple to a second Higgs and the top quark couples to a different Higgs than up and charm quark.

For specifying the Yukawa sector, run
`
  model_list([nH, compute_equivalent_theories])
`
which will output a list of all possible Yukawa sectors with $n_D$ Higgs particles. `compute_equivalent_theories` is a Boolean to determine, if you want to compute theories like `[u1,u1,u3,d1,d1,d1,l1,l1,l1]` and `[u1,u3,u1,d1,d1,d1,l1,l1,l1]` separately or not. These theories are equivalent in the sense, that the formulae for calculating the anomalies are generation independent.

The package is constructed to give you maximal freedom in choosing your theory assumptions:

`model_list()` constructs all models, which obey type invariance, but not generational invariance, i.e. models like `[u1,u1,u3,d1,d1,d1,l1,l1,l1]` are possible, but not models like `[u1,u1,u1,d1,d1,u1,l1,l1,l1]`. If you want to compute other models, you can easily write and include a different generating function that also takes your favourite models into account, as long as every Higgs doublet only couples to one fermion.

`runDFSZ()` is the core of the package. It constructs all possible minimal $V_{explicit breaking}$ potentials with $n_D-1$ terms, computes the solutions for the PQ charges, computes the anomaly ratio from the charges and saves it as a countmap. With bad settings this function can take longer than the age of the universe to execute, we suggest to get familiar with the amount of models generated by `model_list()` before doing anything here. Depending on your choice of keywords there are multiple different execution modes:
1. Treat minimal potentials which lead to the same PQ charges as one theory. To do this, set the keyword `same_$\chi$H_one_theory = true`. In this case fewer computations have to be made, since many potentials can a priori be seen to have the same solution. However it is not possible to do sampling in this case.
2. Treat potentials which lead to the same PQ charges as different theories. To do this, set the keyword `same_$\chi$H_one_theory = false`.
3. When choosing option 2 you can decide to not compute all theories, but to sample representatively to obtain a approximately correct anomaly ratio distribution. To do this, you can specify the logarithm of the maximal number of models in a theory you want to calculate (`log_calculate_all_models_if_below`) and the fraction of models you want to sample (`sample_models_factor`). You probably want to start to sample at $n_D=6$ or $n_D=7$.
4. Additionally there is a function `runDFSZ_saveall()`, which not only saves a countmap of anomaly ratios, but every single minimal potential with corresponding solution for PQ charges, anomalies, anomaly ratio and multiplicities. Obviously this function needs more memory and runtime, run above $n_D=5$ at your own risk.