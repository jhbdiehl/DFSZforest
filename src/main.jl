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

#=
nHs = []
for model in generate_all_models()
    append!(nHs, length(unique(model)))
end
sum(nHs .== 8)
string(nHs[1])
=#

for model in generate_all_models()[111:end]
    un = unique(model)
    nH = length(un)
    @time begin
    quads, multis = get_quads(model)
    tot = binomial(length(quads),nH-2)
    @info ""
    @info "Next model group: $model"
    @info ""
    @printf "Your model-group has %.3E models \n" tot
    as, bs = get_numquads(quads, un, nH)
    myEoN = get_EoNfunc(model)
    end
    if tot <= 10^10
        # Save all ARs for a specific model
        @time begin
        proc_rs = similar(bs, tot)
        rs_ws = similar(multis, length(proc_rs))
        parallel_alleqn_solve_proc!(proc_rs, rs_ws, as, bs, multis, tot, myEoN)
        save_AR(model, proc_rs, rs_ws, 1; folder="n"*string(nH)*"/")
        end
    else
        chunk = 10^8
        m=10
        @time for i in 1:m
            @info "Computing round $i of $m"
            proc_rs = similar(bs, chunk)
            rs_ws = similar(multis, length(proc_rs))
            parallel_randeqn_solve_proc!(proc_rs, rs_ws, as, bs, multis, tot, myEoN)
            save_AR(model, proc_rs, rs_ws, i; folder="n"*string(nH)*"/")
        end
    end
end


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

#Save samples
#=
chunk = 10^7
m=10
@time for i in 1:m
    @info "Computing round $i of $m"
    proc_rs = similar(bs, chunk)
    rs_ws = similar(multis, length(proc_rs))
    @time parallel_randeqn_solve_proc!(proc_rs, rs_ws, as, bs, multis, tot)
    @time save_AR(model, proc_rs, rs_ws, i)
end
=#

# Read and plot data
#=
gaghs = similar(3:9, Any)
for k in 3:9
    @info "$k"
    files = readdir("./data/DFSZ_models/preliminary2/n"*string(k))
    #hist_list = similar(files, Any)
    @time for (i, file) in enumerate(files)
        model = fname2model(file)

        tt = read_AR(model; folder="/preliminary2/n"*string(k)*"/")
        tt = normalize(tt; mode=:probability)
        p1 = plot(tt, lt=:stepbins, label="", title="$(file[1:end-5])", 
            xlabel="E/N", ylabel="Probability", xrange=(-10,13),
            bottom_margin=2Plots.mm, legend=:topright,
            size=(400,300), lw=2)
        savefig(p1, "plots/preliminary2/ARs/$(file[1:end-5])_ARs.pdf")
        #gagh = gag_histogram(tt; mode=:probability)
        #hist_list[i] = gagh
    end
    #gaghs[k-2] = merge(hist_list...)
end

files = readdir("./data/DFSZ_models/preliminary2/n"*string(8))[1]
files[1:end-5]
model = fname2model(files)
tt = read_AR(model; folder="preliminary2/n"*string(8)*"/")
plot(tt, lt=:stepbins, xrange=(-10,13))
gaghs = normalize.(gaghs; mode=:probability)
gagscdf = gag_cdf.(gaghs)

plot()
#for gagh in gaghs[1]
p1 = plot!(gaghs[end], lt=:stepbins, label="", title="DFSZ n=9 axion models PDF", 
    xlabel=L"ga\gamma\gamma \;\; [\log\;\mathrm{GeV}^{-1}]", ylabel="Probability",
    bottom_margin=2Plots.mm, legend=:topright,
    size=(400,300), lw=2)
#end
plot!()
savefig(p1, "plots/preliminary/n9_pdf.pdf")

plot()
for gagh in gaghs[2:end]
    p2 = plot!(gagh.edges[1][1:end-1], gag_cdf(gagh), label="", title="DFSZ axion model CDF", 
        xlabel=L"ga\gamma\gamma \;\; [\log\;\mathrm{GeV}^{-1}]", ylabel="Probability for bigger gaγγ",
        bottom_margin=2Plots.mm, legend=:topright,
        size=(400,300), lw=2)
end
p3 = plot!()
savefig(p3, "plots/preliminary/alln_cdf.pdf")



# example plot
model = fname2model(readdir("./data/DFSZ_models/n3")[1])
fname=model2string(model)
lab = fname[1:11]*"\n    "*fname[12:20]*"\n    "*fname[21:end]
p1 = plot(gagh, lt=:stepbins, label=lab, title="DFSZ axion model PDF", 
    xlabel=L"ga\gamma\gamma \;\; [\log\;\mathrm{GeV}^{-1}]", ylabel="Probability",
    bottom_margin=2Plots.mm, legend=:topright,
    size=(400,300), lw=2)
savefig(p1, "plots/test.pdf")

p2 = plot(gagh.edges[1][1:end-1], gagcdf, label=lab, title="DFSZ axion model CDF", 
    xlabel=L"ga\gamma\gamma \;\; [\log\;\mathrm{GeV}^{-1}]", ylabel="Probability for bigger gaγγ",
    bottom_margin=2Plots.mm, legend=:topright,
    size=(400,300), lw=2)
savefig(p2, "plots/test2.pdf")
=#