using StaticArrays
using Symbolics
using StatsPlots
using Distributions
using HDF5

include("helpers.jl")

@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 

# More different possible values for charges -> less peaked structure
# Restrict to positive charges -> very similar to additive models of Plakkot Hoof
# KSVZ introduces more variation in E and N values due to quark representations not fixed
# If mean(charges) = 0 this ensures median(EoNs) = 5/3. However for n=5 it seems to magically work without!


n = 1000000
models = rand(Uniform(-1.0,1.0), n, 9)
models = rand(2*Poisson(2) - 3.5, n, 9)
models = rand(-0:0.01:100, n, 9)
mean(vec(models))
EoNs = Array{Float64}(undef, n)

for i in eachindex(EoNs)
    EoNs[i] = EoverN(models[i,:]...)
end

EoNs = EoNs[-1e10 .< EoNs .< 1e10]
median(EoNs)

histogram(EoNs, bins=5/3-20:0.01:5/3+20)
vline!([5/3])



nlist=[5]
EoNlist = similar(nlist, Any)
allEoN = []
for (i, n) in enumerate(nlist)
    fid = h5open("./data/DFSZ_models/220616-nbilin_fullsol/n$n/full_n$n.h5")

    for fold in fid
        if occursin("Chis order: u1u2u3d1d2d3l1l2l3s", string(fold))
            nothing
        else
            @time begin
                println(fold)
                global cclist = Array{Float64}(undef, 0, 10)
                EoNlist = []
                for bil in fold
                    cc = read(bil["Chis"])
                    cclist = vcat(cclist, cc)
                    EoN = read(bil["EoN"])
                    append!(EoNlist, EoN)
                end
                myEoN = EoNlist#[unique(i -> round.(cclist[i,:],digits=4), 1:size(cclist)[1])]
                myEoN = myEoN[-1e10 .< myEoN .< 1e10]
                append!(allEoN, myEoN)
            end
        end
    end
end
median(allEoN) - 5/3

EoverN(mean(cclist[:,1:9], dims=1)...)

histogram(vec(cclist[:,1:9]), bins=-6:0.1:6)
histogram(allEoN,bins=-10:0.01:15)

histogram(rand(2*Poisson(2)-3, 100000))