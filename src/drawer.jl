# This code was adapted from https://discourse.julialang.org/t/is-there-a-more-efficient-way-to-define-an-iterator-through-ordered-tuples/20819

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

function col_tups(tp::TupIter{N}, num) where N
    res = Vector{NTuple{N, Int}}(undef, num)
    for (i,t) in enumerate(tp)
        res[i]=t
        i == num && break
    end
    res
end;

function tupidx(tup::NTuple{N,Int}) where N
    s = 1
    for i=1:N
        s += binomial(tup[i]-1, i)
    end
    return s
end;

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

function tupidx(bnc, tup::NTuple{N,Int}) where N
    s = 1
    for i=1:N
        (tup[i]>1) && (s += bnc[i][tup[i]-1])
    end
    return s
end;

function idxtup(bnc, idx, ::Val{k}) where k
    k == 1 && return (idx,)
    last = searchsortedlast(bnc[k], idx-1)
    idx -= bnc[k][last]
    return [idxtup(bnc, idx, Val(k-1))...,last+1]
end;

function check_tups(bnc, tp::TupIter{N}, n) where N
    vn = Val(N)
    for (i,t) in enumerate(tp)
        if t !== idxtup(bnc, i, vn)
            @show N, i, t, idxtup(bnc, i, vn)
            error("fail")
        end
        i >= n && break
    end
end;

function construct_draw(nD)
    for (i,t) in enumerate(TupIter{nD-2}())
        i > 10_000 && break
        @assert i == tupidx(t)
    end
    bnc = binom_cache(nD-2, 1000);
    idxarr = similar([1],nD-2);
    ac = i->myidxtup!(idxarr, bnc, i, Val(nD-2));
end

function myidxtup!(idxarr::Vector{<:Real}, bnc, idx::Int, ::Val{k}) where k
    if k > 1
        for j in k:-1:2
            last = searchsortedlast(bnc[j], idx-1)
            idx -= bnc[j][last]
            idxarr[j+2] = last+1
        end
    else
        nothing
    end
    idxarr[3] = idx
    idxarr[2] = 0
    idxarr[1] = -1
    SVector(ntuple(i -> idxarr[i] +2, Val(k+2)))
end

function make_idx_bnc(N)
    for (i,t) in enumerate(TupIter{N}())
        i > 10_000 && break
        @assert i == tupidx(t)
    end
    bnc = binom_cache(N, 1000);
    idxarr = similar([1],N+2);
    return idxarr, bnc
end