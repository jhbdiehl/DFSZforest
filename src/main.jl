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

a,m = get_quads([u1,u1,u1,d1,d1,d1,l1,l1,l1])
m
aa = get_bilins([u1,u1,u1,d1,d1,d1,l1,l1,l1])


dataset = "test"
filename = "testEoNs"

# Suitable sample_models_factor for not too extensive runtimes (Order of hours)
# n=7 -> 0.001
# n=8 -> 0.000001
# n=9 -> 0.000000001
models, multis = model_list(nD=[4])
runDFSZ(dataset, models; model_multiplicity=multis, same_χH_one_theory=true, filename=filename)
runDFSZ_saveall(dataset, models; model_multiplicity=mmultis, filename=filename)
h5totxt(filename; folder="./data/DFSZ_models/"*dataset*"/")



models, mmultis = model_list(nD=[4])
e1=read_EoN(dataset, models; specifier=testEoNs)

ke1 = collect(keys(e1))
ve1 = collect(values(e1))
ve1 = ve1[isnan.(ke1) .== 0]
ke1 = ke1[isnan.(ke1) .== 0]

plot()
plot_EoN(e1)

check_symmetry(e1)


EoNs = []
mychi = Matrix{Float64}(undef,0,10)
tot = 0 

for (i, a) in enumerate(fid)
    k = read(a)
    if i != 1
        println(keys(fid)[i])
        for tuple in k#collect(values(k))
            dat = tuple[2]
            tot += sum(dat["model_multiplicity"] .* dat["multis"])
            mychi = vcat(mychi, round.(dat["Chis"], digits=5))
            append!(EoNs, round.(dat["EoN"], digits=5))
                #write(io, model*"   "*mmult*"      "*rpad(dat["multis"][j],3)*"                 "*lpad(round(dat["EoN"][j], digits=3),8)*" "*lpad(Int(round(dat["E"][j])),5)*lpad(Int(round(dat["N"][j])),5)*(*(lpad.(round.(dat["Chis"][j,:],digits=3),9)...))*" "*(*(lpad.(filter.(x -> !isspace(x), dat["terms"][j,:]),18)...))*"\n")
        end
    end
end
tot
mychi
EoNs