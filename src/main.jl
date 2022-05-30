using LinearAlgebra, StaticArrays
using Random, Statistics, StatsBase
using BenchmarkTools
using Base.Threads
using Plots
using Random: default_rng
using Symbolics
using Combinatorics
using FileIO, JLD2
using LaTeXStrings, Printf


include("./drawer.jl")
include("./helpers.jl")

@variables u1::Int u2::Int u3::Int d1::Int d2::Int d3::Int l1::Int l2::Int l3::Int 

# 58 cores produces segmentation fault.
# 56 cores doesnt.

function make_idx_bnc(N)
    for (i,t) in enumerate(TupIter{N}())
        i > 10_000 && break
        @assert i == tupidx(t)
    end
    bnc = binom_cache(N, 1000);
    idxarr = similar([1],N);
    return idxarr, bnc
end

function parallel_randeqn_solve_proc!(
    proc_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer},
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer},
    tot::Int, myEoN
) where N

    idxarr, bnc = make_idx_bnc(N)

    @threads for i in eachindex(proc_rs)
        @inbounds begin
            # When using drawer need to allocate indices. Is there a way around?
            idxs_i =  myidxtup!(idxarr, bnc, rand(1:tot), Val(N))#rand_idxs(default_rng(), eachindex(as), Val(N)) # drawer(13131) #
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = myEoN(r)
            rs_ws[i] = prod(ws[idxs_i])
        end
    end
end

function parallel_alleqn_solve_proc!(
    proc_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer},
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer},
    tot::Int, myEoN
) where N

    idxarr, bnc = make_idx_bnc(N)

    @threads for i in eachindex(proc_rs)
        @inbounds begin
            # When using drawer need to allocate indices. Is there a way around?
            idxs_i =  myidxtup!(idxarr, bnc, i, Val(N))#rand_idxs(default_rng(), eachindex(as), Val(N)) # drawer(13131) #
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = myEoN(r)
            rs_ws[i] = prod(ws[idxs_i])
        end
    end
end


dataset = "220530-mediumrun"

#a = generate_all_models()
b, ms = generate_unique_models()

@time for (k, model) in enumerate(b) #a[length.(unique.(a)) .== 8]
    bilins = unique(sort.(collect(combinations(model,2)), by=x->Symbol(x)))
    goodbilins = bilins[length.(unique.(bilins)) .== 2]
    for bilin in goodbilins
        p1, p2 = bilin[1], bilin[2]
        valp1, valp2 = bilinvals(bilin)
        un = unique(model)
        nH = length(un)
        @time begin
        quads, multis = get_quads(model; p1=p1, p2=p2, valp1=valp1, valp2=valp2)
        tot = binomial(length(quads),nH-2)
        @info ""
        @info "Next model group: $model with bilinear $bilin" 
        @info ""
        @printf "Your model-group has %.3E models \n" tot
        as, bs = get_numquads(quads, un, nH; p1=p1, p2=p2, valp1=valp1, valp2=valp2)
        myEoN = get_EoNfunc(model; p1=p1, p2=p2, valp1=valp1, valp2=valp2)
        end
        if tot <= 10^9
            # Save all ARs for a specific model
            @time begin
            proc_rs = similar(bs, tot)
            rs_ws = similar(multis, length(proc_rs))
            parallel_alleqn_solve_proc!(proc_rs, rs_ws, as, bs, multis, tot, myEoN)
            save_AR(model, proc_rs, rs_ws, 1; folder=dataset*"/n"*string(nH)*"/", bilin=bilin, ms=ms[k])
            end
        else
            chunk = 10^8
            m=10
            @time for i in 1:m
                @info "Computing round $i of $m"
                proc_rs = similar(bs, chunk)
                rs_ws = similar(multis, length(proc_rs))
                parallel_randeqn_solve_proc!(proc_rs, rs_ws, as, bs, multis, tot, myEoN)
                save_AR(model, proc_rs, rs_ws, i; folder=dataset*"/n"*string(nH)*"/", bilin=bilin, ms=ms[k])
            end
        end
    end
end

@info ""
@info ""
@info "Start plotting!!!"
@info ""
@info ""


# This would be the procedure to actually calculate really all models!
#=
for model in with_replacement_combinations([u1, d1, l1],9)
    un = unique(model)
    nH = length(un)
    if nH < 3
        nothing
    else
        quads, multis = get_quads(model)
        tot = binomial(length(quads),nH-2)
        @printf "Your model-group has %.3E models \n" tot
        @time as, bs = get_numquads(quads, un, nH)
        myEoN = get_EoNfunc(model)


        # Save all ARs for a specific model
        proc_rs = similar(bs, tot)
        rs_ws = similar(multis, length(proc_rs))
        @time parallel_alleqn_solve_proc!(proc_rs, rs_ws, as, bs, multis, tot)
        @time save_AR(model, proc_rs, rs_ws, 1; folder="n"*string(nH)*"/")
    end
end
=#


# Plot data
#=
for k in 5:5
    @info "$k"
    savefolder = "./plots/"*dataset*"/ARs/n"*string(k)
    mkpath(savefolder)
    folders = readdir("./data/DFSZ_models/"*dataset*"/n"*string(k); join=true)
    #hist_list = similar(files, Any)
    for folder in folders
        files = readdir(folder)
        @time for (i, file) in enumerate(files)
            println(file)
            model = fname2model(file)
            bilin = fname2bilin(file)
            tt = read_AR(model; folder=dataset*"/n"*string(k), bilin=bilin)
            tt = normalize(tt; mode=:probability)
            p1 = plot(tt, lt=:stepbins, label="", title="$(file[1:end-5])", 
                xlabel="E/N", ylabel="Probability", xrange=(5/3-9,5/3+9),
                bottom_margin=2Plots.mm, legend=:topright,
                size=(400,300), lw=2)
            savefig(p1, savefolder*"/$(file[1:end-5])_ARs.pdf")
            #gagh = gag_histogram(tt; mode=:probability)
            #hist_list[i] = gagh
        end
    end
    #gaghs[k-2] = merge(hist_list...)
end
=#

@time for k in 4:9
    @info "$k"
    savefolder = "./plots/"*dataset*"/ARs/n"*string(k)
    mkpath(savefolder)
    folders = readdir("./data/DFSZ_models/"*dataset*"/n"*string(k); join=true)
    
    tttot_list = similar(folders, Any)
    for (j, folder) in enumerate(folders)
        m = fname2m(folder)
        files = readdir(folder)
        hist_list = similar(files, Any)
        @time for (i, file) in enumerate(files)
            model = fname2model(file)
            bilin = fname2bilin(file)
            tt = read_AR(model; folder=dataset*"/n"*string(k), bilin=bilin, m=m)
            tt = normalize(tt; mode=:probability)
            #gagh = gag_histogram(tt; mode=:probability)
            hist_list[i] = tt
        end
        tttot = merge(hist_list...)
        tttot = normalize(tttot; mode=:probability)
        tttot.weights .*= parse(Int,string(m))
        tttot_list[j] = tttot
        type = split(folder, "/")[end]
        println(type)
        p2 = plot(tttot, lt=:stepbins, label="", title="$(type)", 
        xlabel="E/N", ylabel="Probability", xrange=(5/3-9,5/3+9),
        bottom_margin=2Plots.mm, legend=:topright,
        size=(400,300), lw=2)
        savefig(p2, savefolder*"/full_$(type)_ARs.pdf")
    end
    ttall = merge(tttot_list...)
    ttall = normalize(ttall; mode=:probability)
    p1 = plot(ttall, lt=:stepbins, label="", title="n$(k)", 
    xlabel="E/N", ylabel="Probability", xrange=(5/3-9,5/3+9),
    bottom_margin=2Plots.mm, legend=:topright,
    size=(400,300), lw=2)
    savefig(p1, savefolder*"/full_n$(k)_ARs.pdf")
    #gaghs[k-2] = merge(hist_list...)
end