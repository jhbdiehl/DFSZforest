using LinearAlgebra, StaticArrays
using Random, Statistics, StatsBase
using BenchmarkTools
using Base.Threads
using Plots
using Random: default_rng
using Symbolics
using Combinatorics
using FileIO, JLD2

function I_static(::Val{N}, ::Type{T}) where {N,T<:Real}
    convert(SMatrix{N,N,T}, Diagonal(SVector(ntuple(i -> one(T), Val(N)))))
end

function mysolve(as::AbstractVector{<:SVector{N,T}}, bs, idxs) where {N,T<:Real}
    A_true = copy(hcat(as[idxs]...)')
    solveable = !(det(A_true) ≈ 0)
    A_dummy = A_true .* zero(T) .+ I_static(Val(N), T)
    A = ifelse(solveable, A_true, A_dummy)

    b = bs[idxs]
    b_cand = A \ b
    b_nonsolve = SVector(ntuple(i -> convert(T, NaN), Val(N)))
    ifelse(solveable, b_cand, b_nonsolve)
end

function rand_idxs(rng::AbstractRNG, from::AbstractVector{Int}, ::Val{N}) where N
    SVector(ntuple(i -> rand(rng, from), Val(N)))
end

process_solution(r::AbstractVector{<:Real}) = sum(r) ## Dummy

###################### Functions by Johannes ####################################

function next_tup(i, tup::NTuple{N,Int}) where N
    N == 1 && return (tup[1]+1,)
    if tup[1]+1 != tup[2]
        return (tup[1]+1, Base.tail(tup)...)
    else
        return (i, next_tup(i+1, Base.tail(tup))...)
    end
end;

struct TupIter{N} end
function Base.iterate(tp::TupIter{N}) where N
    t = ntuple(i->i, N)
    return (t,t)
end;

function Base.iterate(::TupIter{N}, state) where N
    nx = next_tup(1, state)
    return nx, nx
end;

function tupidx(tup::NTuple{N,Int}) where N
    s = 1
    for i=1:N
        s += binomial(tup[i]-1, i)
    end
    return s
end;

function make_idx_bnc(N)
    for (i,t) in enumerate(TupIter{N}())
        i > 10_000 && break
        @assert i == tupidx(t)
    end
    bnc = binom_cache(N, 1000);
    idxarr = similar([1],N);
    return idxarr, bnc
end

function binom_cache(k, N)
    res = Vector{Vector{Int}}(undef, k)
    for k_ = 1:k
        bn = res[k_] = Vector{Int}(undef, N)
        for n=1:N
            bn[n] = binomial(n, k_)
        end
    end
    res
end;

function myidxtup!(idxarr::Vector{<:Real}, bnc, idx::Int, ::Val{k}) where k
    if k > 1
        for j in k:-1:2
            last = searchsortedlast(bnc[j], idx-1)
            idx -= bnc[j][last]
            idxarr[j] = last+1
        end
    else
        nothing
    end
    idxarr[1] = idx
    SVector(ntuple(i -> idxarr[i], Val(k)))
end

function get_EoNfunc(model)
    un = unique(model)
    nH = length(un)
    notu1d1 = isequal.(un,u1) .+ isequal.(un,d1) .== false
    EoN = EoverN(model...)
    EoN = substitute(EoN, Dict(u1=>1,d1=>1))
    myEoNtmp = build_function(EoN, un[notu1d1])
    myEoN = eval(myEoNtmp)
end

function save_AR(model, proc_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer}, i::Int)
    fname = _model2string(model)
    good_idxs = findall(!isnan, proc_rs)
    good_proc_rs = proc_rs[good_idxs]
    good_rs_ws = rs_ws[good_idxs]

    ARh = fit(Histogram, good_proc_rs, FrequencyWeights(good_rs_ws), -1000:0.01:1000)

    @info "Saving..."

    if isfile("./"*fname*".jld2") && i != 1 # i.e. overwrite if file exists before calculation
        cm_old = FileIO.load("./"*fname*".jld2", "ARs")
        if cm_old.edges == ARh.edges
            cm_new = Histogram(ARh.edges, cm_old.weights + ARh.weights)
        end
    else
        cm_new = ARh
    end
    FileIO.save("./"*fname*".jld2", Dict("ARs" => cm_new))

    @info "Done!"
end

function _model2string(model)
    un = unique(model)
    nH = length(un)

    str = "n"*string(nH)
    for arg in model
        str *= "_"*string(arg) 
    end
    return str
end

function EoverN(u1, u2, u3, d1, d2, d3, l1, l2, l3)
    2/3 + 2 * (u1 + u2 + u3 + l1 + l2 + l3) / (u1 + u2 + u3 + d1 + d2 + d3)
end

function get_quads(model)
    vars = unique(model)
    us = Num[]
    for u in [u1, u2, u3]
        if sum(isequal.(u, vars)) !== 0
            append!(us, u)
        else
            nothing
        end
    end
    doublets = collect(powerset(vars, 2,2))
    diffs = [length(setdiff(doublets[i], us)) for i in 1:length(doublets)]
    sig = isodd.(diffs) .* 2 .- 1 # lifehack! Maps BitVector (0 for false, 1 for true) to (-1 for false and 1 for true)
    sig = hcat(sig, ones(size(sig)))
    doublets = permutedims(hcat(doublets...))
    doublets .*= Int64.(sig).* 2//1 .* 1//2
    doublets = [doublets[i,:] for i in 1:size(doublets, 1)]
    allquads = collect(powerset(doublets, 2,2))
    allquads = permutedims(hcat(allquads...))
    allquads_conj = deepcopy(allquads)
    allquads_conj[:,1] .*= -1//1
    allquads = vcat(allquads, allquads_conj)
    allquads = [sum(sum(allquads[i,:])) for i in 1:size(allquads)[1]]
    ut = tril(permutedims(hcat([isequal.(allquads, -allquads[i]) for i in 1:length(allquads)]...)),-1) # make sure to remove duplicates like a = -b
    fut = findall(ut)
    for i in 1:length(fut)
        allquads[fut[i][1]] = allquads[fut[i][2]]
    end
    uniquads = allquads
    fillones = substitute(uniquads, Dict([var => 1//1 for var in vars]))
    fillones = (sign.(fillones) .== 0) + sign.(fillones) # make sign return 1 if element = 0
    fillones = (fillones .== 1) .* 2//1 .- 1//1 # weird hack to make sure the array contains only Ints, because unique() does not work with negative floats and symbols
    uniquads .*= fillones
    quads = uniquads[isequal.(uniquads, d1+u1) .!= 1//1] # effectively is not equal...
    quads .*= (length.(string.(quads)) .<= 7) .+ 1//1 # weird way of multiplying all quads by 2 if they only contain two different higgs. This is so sum(abs(higgs)) == 4
    quads = countmap(quads)
    return collect(keys(quads)), collect(values(quads))
end

function get_numquads(quads, un, nH)
    notu1d1 = isequal.(un,u1) .+ isequal.(un,d1) .== false
    numquads = zeros(Int8,length(quads), nH-2)

    u1d1dict = Dict{Num, Int8}(un .=> 0)
    u1d1dict[u1] = 1
    u1d1dict[d1] = 1
    bs = SVector{length(quads), Float64}(-1 .* Symbolics.value.(substitute.(quads, (u1d1dict,))))

    for i in 1:nH-2
        mydict = Dict(un .=> 0.0)
        mydict[un[notu1d1][i]] = 1.0
        numquads[:,i] = Symbolics.value.(substitute.(quads, (mydict,)))
    end
    as = Vector{SVector{nH-2, Float64}}(_vecvec(numquads))
    return as, bs
end

_vecvec(mat) = [mat[i,:] for i in 1:size(mat,1)]

@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 

############### Johannes Initializer ##############################################

model = [u1, u2, u3, d1, d2, d3, l1, l2, l3]

@time begin
    un = unique(model)
    nH = length(un)
    quads, multis = get_quads(model)
    tot = binomial(length(quads),nH-2)
    as, bs = get_numquads(quads, un, nH)
    myEoN = get_EoNfunc(model)
end

#####################################################################################




proc_rs = similar(bs, 10^8);
rs_ws = similar(multis, length(proc_rs))

function parallel_randeqn_solve_proc!(
    proc_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer},
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer}
) where N

    idxarr, bnc = make_idx_bnc(N)

    @threads for i in eachindex(proc_rs)
        @inbounds begin
            idxs_i = myidxtup!(idxarr, bnc, i, Val(N))#rand_idxs(default_rng(), eachindex(as), Val(N))
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = myEoN(r) # process_solution(r) #
            rs_ws[i] = prod(ws[idxs_i])
        end
    end
end

@time begin parallel_randeqn_solve_proc!(proc_rs, rs_ws, as, bs, multis)
 save_AR(model, proc_rs, rs_ws, 1) end

fname = _model2string(model)
tt = FileIO.load("./"*fname*".jld2", "ARs")
plot(tt, lt = :stepbins, xrange=(-10,13))


########### Benchmarks ##################################################################

idxs_i = rand_idxs(default_rng(), eachindex(as), Val(7))
# Referenz: Gut implementierter equation solver
@btime mysolve($as, $bs, $idxs_i)

r = mysolve(as, bs, idxs_i)
# Similar runtime, same allocations :-)
@btime myEoN(r)
@btime process_solution(r)

# Different runtime, different allocations :-(
@btime myEoN($r)
@btime process_solution($r) 

@btime rand_idxs(default_rng(), eachindex($as), Val(7))
@time idxarr, bnc = make_idx_bnc(7)
# This function takes twice as long as yours, but is still 1/10 of mysolve runtime.
# Should not influence total time by more than 10%, but does so?!
# Probably the threads don't like that they are feeding on the same object? But somehow it seems to work!
@btime myidxtup!($idxarr, $bnc, 512, Val(7))
