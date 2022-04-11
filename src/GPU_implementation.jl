
#using CUDA
#using BenchmarkTools
using Symbolics
using Combinatorics
using LinearAlgebra
#using CUDA.CUSPARSE
#using SparseArrays
#using StatsPlots
using StaticArrays
using IterTools, StatsBase
using FileIO, JLD2



include("./helpers.jl")

@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 


#histogram!(AR, bins=-10:0.1:13)

function get_EoNfunc(model)
    un = unique(model)
    nH = length(un)
    notu1d1 = isequal.(un,u1) .+ isequal.(un,d1) .== false
    EoN = EoverN(model[1:3], model[4:6], model[7:9])
    EoN = substitute(EoN, Dict(u1=>1,d1=>1))
    myEoNtmp = build_function(EoN, un[notu1d1]...)
    myEoN = eval(myEoNtmp)
end

function solve(un, nH, subs)

    @info "Solving LES parallelized on CPU."

    notu1d1 = isequal.(un,u1) .+ isequal.(un,d1) .== false

    sol = zeros(Float64, size(subs)[1], nH-2)
    Threads.@threads for i in 1:size(subs)[1]
        A = SMatrix{nH-2, nH-2, Float64}(subs[i,1:7,:]')
        b = Float64.(subs[i,8,:])
        if _issolvable(A)
            sol[i,:] = A \ b
        else
            sol[i,:] = [NaN for j in 1:nH-2]
        end
    end
    sol = sol[isnan.(sol[:,1]) .== 0,:]
    return sol
end

function calculate_EoverN(myEoN, sol)

    @info "Calculating Anomaly Ratios from solutions to equation systems."

    AnomalyRatio = Array{Float64}(undef,size(sol)[1]);
    Threads.@threads for i in 1:size(sol)[1]
        AnomalyRatio[i] = myEoN(sol[i,:]...)
    end
    AnomalyRatio = AnomalyRatio[isfinite.(AnomalyRatio)]
    AR = AnomalyRatio[AnomalyRatio .< 1e3 .&& AnomalyRatio .> -1e3]
    return countmap(AR)
end

function calculate_subsets(model; nsamps::Union{Nothing, Int}=nothing)
    un = unique(model)
    nH = length(un)
    notu1d1 = isequal.(un,u1) .+ isequal.(un,d1) .== false
    quads = get_quads(model)
    numquads = zeros(Int8,length(quads), nH-1)
    tot = binomial(length(quads),nH-2)
    
    println("")
    @info " Calculating subsets for n=$nH model!"

    u1d1dict = Dict{Num, Int8}(un .=> 0)
    u1d1dict[u1] = 1
    u1d1dict[d1] = 1
    numquads[:,nH-1] = Vector{Int8}(-1 .* Symbolics.value.(substitute.(quads, (u1d1dict,))))

    for i in 1:nH-2
        mydict = Dict(un .=> 0.0)
        mydict[un[notu1d1][i]] = 1.0
        numquads[:,i] = Symbolics.value.(substitute.(quads, (mydict,)))
    end

    if typeof(nsamps) == Nothing
        tot > 1e6 && @warn "You are trying to calculate $tot subsets. This will take more than $(tot/1e6) Mbyte of RAM and may take a while. If you think this is a bad idea, please abort!"
        subsets = collect( powerset(_vecvec(numquads), nH-2, nH-2))

    elseif typeof(nsamps) == Int
        if nsamps > tot
             @warn "You tried to choose nsamps = $nsamps > total number of subsets = $tot ! I will calculate all models instead and not draw samples."
             subsets = collect(powerset(_vecvec(numquads), nH-2, nH-2))
        else
            nsamps > 1e-3 * tot && @warn "This code does not have an implementation to make sure you don't choose samples twice! nsamps = $nsamps is dangerously close to the total number of subsets ($tot)!"
            #subsets = Vector{Vector{Vector{Int8}}}()
            subsets = zeros(Int8, nsamps, nH-1, nH-2)
            nq = _vecvec(numquads)
            for i in 1:nsamps
                subsets[i,:,:] = hcat(sample(nq, nH-2; replace=false)...)
            end
            subsets = _vecvecvec(subsets)
        end
    end
    return subsets
end


function calculate_subsets!(iter::Union{Iterators.Take, Iterators.Drop, Iterators.Flatten}; steps=100000)
    subsets = collect( Iterators.take(iter, steps))
    global iter = Iterators.drop(iter, steps)
    return subsets
end

function calculate_subsets(nsamps::Int, nH::Int, nq)
    subsets = zeros(Int8, nsamps, nH-1, nH-2)
    Threads.@threads for i in 1:nsamps
        subsets[i,:,:] = hcat(sample(nq, nH-2; replace=false)...)
    end
    #@time subsets = _vecvecvec(subsets)
    return subsets
end


function get_numquads(quads, un, nH)
    notu1d1 = isequal.(un,u1) .+ isequal.(un,d1) .== false
    numquads = zeros(Int8,length(quads), nH-1)

    u1d1dict = Dict{Num, Int8}(un .=> 0)
    u1d1dict[u1] = 1
    u1d1dict[d1] = 1
    numquads[:,nH-1] = Vector{Int8}(-1 .* Symbolics.value.(substitute.(quads, (u1d1dict,))))

    for i in 1:nH-2
        mydict = Dict(un .=> 0.0)
        mydict[un[notu1d1][i]] = 1.0
        numquads[:,i] = Symbolics.value.(substitute.(quads, (mydict,)))
    end
    return _vecvec(SMatrix{length(quads), nH-1, Int8, (nH-1)*length(quads)}(numquads))
end

function checks(tot, nsamps)
    if typeof(nsamps) == Int
        if nsamps > tot
            @warn "You tried to choose nsamps = $nsamps > total number of subsets = $tot ! I will calculate all models instead and not draw samples."
            nsamps = nothing
        elseif 1e-3 * tot < nsamps < tot
            error("This code does not have an implementation to make sure you don't choose samples twice! nsamps = $nsamps is dangerously close to the total number of subsets ($tot)!")
        end
    elseif typeof(nsamps) == Nothing
        tot > 1e7 && @warn "You are trying to calculate $tot subsets. This will take more than $(tot/1e6) Mbyte of RAM and may take a while. If you think this is a bad idea, please abort!"
    end
    return nsamps
end

function _issolvable(mat)
    #rank(mat) == size(mat)[1]
    ≈(det(mat), 0.0; atol=1e-10) == false # tolerance chosen arbitrarily. Make sure this fits estimated numerical errors!
end




# 58 cores produces segmentation fault.
# 56 cores doesnt.
model = [u1, u2, u1, d1, d2, d1, l1, l1, l1]
nsamps=20000000


################################################
#function do_calculation(model, nsamps)
un = unique(model)
nH = length(un)
quads = get_quads(model)
tot = binomial(length(quads),nH-2)

nsamps = checks(tot, nsamps)
nq = get_numquads(quads, un, nH)
myEoN = get_EoNfunc(model)

println("")
@info " Calculating subsets for n=$nH model!"
if typeof(nsamps) == Nothing
    iter = powerset(nq, nH-2, nH-2)
    k = 1
    while k <= 1e100
        subs = calculate_subsets!(iter; steps=10000)
        if subs == []
            break
        end
        @time sol = solve(model, subs)
        @time AR = calculate_EoverN(myEoN, sol)
        @time save_AR(model, AR, i)
        k += 1
    end
elseif typeof(nsamps) == Int
    for j in 1:Int(round(nsamps/10000000))
        @info "You are at iteration $j of $(Int(round(nsamps/10000000)))"
        @time subs = calculate_subsets(10000000, nH, nq)
        @time sol = solve(un, nH, subs)
        @time AR = calculate_EoverN(myEoN, sol)
        @time save_AR(model, AR, 2)
    end
end
####################################################

#end

#A = [2.0 0.0 1.0 -1.0; 2.0 1.0 0.0 0.0; -2.0 -1.0 0.0 -3.0]
#inv(A[:,1:3])

#B = [1.0 0.0 0.0 -3.0; 2.0 0.0 0.0 -2.0; 2.0 0.0 1.0 -1.0]
#inv(B[:,1:3])
#tt = [2.0 0.0 1.0; 3.0 0.0 0.0; 2.0 0.0 7.0]
#inv(tt)
#@benchmark rank(tt)
#det(tt) == 0


#=
using StatsPlots
AR = read_AR("n9_u1_u2_u3_d1_d2_d3_l1_l2_l3")
@info "Maximal E/N: $(maximum(abs.(AR)))"
histogram(AR, bins=-1000:0.1:1000)
@time cm2 = countmap(AR)
sum(collect(values(AR)))

mcm = merge(+,cm, cm2)

if isfile("./data/test.jld2")
end
cm_old = FileIO.load("./data/test2.jld2", "ARs")
cm_new = merge(+,cm_old,cm2)
FileIO.save("./data/test.jld2", Dict("ARs" => cm_new))

D1 = Dict("a" => 1, "b" => 3)
D2 = Dict("a" => 2, "c" => 1)
merge(+,D1,D2)

hist = _make_hist(AR; bins=-10:0.001:13)
using Plots
plot()
plot_hist!(hist, yaxis=(:linear))

function _make_hist(cm; bins=-10:0.0001:13)
    hi = fit(Histogram, collect(keys(cm)), Weights(collect(values(cm))),bins)#, nbins=5000)
    m = (hi.edges[1] .+ (hi.edges[1][2] - hi.edges[1][1])/2.)[1:end-1]
    #hi.weights = Vector{Real}(hi.weights ./ sum(hi.weights))
    return hi, m
end

function plot_hist!(hist; yaxis=(:log, [0.00000001, :auto]), kwargs...)
    plot!(hist[2], hist[1].weights ./ sum(hist[1].weights); yaxis=yaxis, kwargs...)
    #plot!(hist[2], _gauss_prediction(gauss,hist[1].edges[1]); kwargs..., label="")
end
=#