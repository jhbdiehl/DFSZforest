"""
"""
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
    
    doublets = get_doublets(vars, us)
    
    allquads = collect(with_replacement_combinations(doublets, 2))
    allquads = permutedims(hcat(allquads...))
    allquads_conj = deepcopy(allquads)
    allquads_conj[:,1] .*= -1
    allquads = vcat(allquads, allquads_conj)
    allquads = [sum(sum(allquads[i,:])) for i in 1:size(allquads)[1]]
    allquads .*= (length.(string.(allquads)) .<= 8) .+ 1 # weird way of multiplying all quads by 2 if they only contain two different higgs. This is so sum(abs(higgs)) == 4
    ut = tril(permutedims(hcat([isequal.(allquads, -allquads[i]) for i in 1:length(allquads)]...)),-1) # make sure to remove duplicates like a = -b
    fut = findall(ut)
    for i in 1:length(fut)
        allquads[fut[i][1]] = allquads[fut[i][2]]
    end
    
    quads = allquads
    fillones = substitute(quads, Dict([var => 1 for var in vars]))
    fillones = (sign.(fillones) .== 0) + sign.(fillones) # make sign return 1 if element = 0
    fillones = (fillones .== 1) .* 2 .- 1 # weird hack to make sure the array contains only Ints, because unique() does not work with negative floats and symbols
    quads .*= fillones
    quads = quads[isequal.(quads, 0) .== 0] # remove quads that dont give sensible conditions
    quads = countmap(quads)
    return collect(keys(quads)), collect(values(quads))
end

function get_doublets(vars, us)
    doublets = collect(powerset(vars, 2,2))
    diffs = [length(setdiff(doublets[i], us)) for i in 1:length(doublets)]
    sig = isodd.(diffs) .* 2 .- 1 # lifehack! Maps BitVector (0 for false, 1 for true) to (-1 for false and 1 for true)
    sig = hcat(sig, ones(size(sig)))
    doublets = permutedims(hcat(doublets...))
    doublets .*= Int64.(sig)#.* 2//1 .* 1//2
    doublets = [doublets[i,:] for i in 1:size(doublets, 1)]
    return doublets
end

function get_bilins(model)
    vars=unique(model)
    bilins = unique(sort.(collect(combinations(model,2)), by=x->Symbol(x)))
    bilins = bilins[length.(unique.(bilins)) .== 2]
    o = [iseven(sum(in.(Symbol.(bilin),Ref(Symbol.([u1,u2,u3]))))) for bilin in bilins]
    bi = [bilins[i][1] .+ (-1).^o[i] .* bilins[i][2] for i in 1:length(bilins)]
    vcat(bi, bi .* -1)
end

function orthogonality(model)
    un = unique(model)
    o = in.(Symbol.(un),Ref(Symbol.([u1,u2,u3])))
    sum(un .* (-1).^o)
end


"""
    Calculate Anomaly Ratio, given PQ charges of up- and down quarks as well as leptons.
"""
function EoverN(u1, u2, u3, d1, d2, d3, l1, l2, l3)
    2/3 + 2 * (u1 + u2 + u3 + l1 + l2 + l3) / (u1 + u2 + u3 + d1 + d2 + d3)
end

function N(u1,u2,u3,d1,d2,d3,l1,l2,l3)
    1/2 * (u1 + u2 + u3 + d1 + d2 + d3)
end

function get_Nfunc(model)
    un = unique(model)
    tN = N(model...)
    myNtmp = build_function(tN, un)
    myN = eval(myNtmp)
end

"""
    Get multiplicities of unique E/N values.
"""
function multiplicities(AnomalyRatio)
    sort!(AnomalyRatio)
    multi = ones(Int64, length(unique(AnomalyRatio)))
    j = 1
    for i in 1:length(AnomalyRatio)-1
        if isequal(AnomalyRatio[i], AnomalyRatio[i+1])
            multi[j] += 1
        else
            j+=1
        end
    end
    return multi
end

function model2string(model)
    un = unique(model)
    nH = length(un)

    str = "n"*string(nH)
    for arg in model
        str *= "_"*string(arg) 
    end
    return str
end

function fname2model(model_str)
    args = split(model_str, "_")
    args[end] = split(args[end], ".")[1]
    symargs = string2num(args)
    return symargs[2:end-1]
end

function fname2m(model_str)
    args = split(model_str, "/")[end]
    args = split(args, "n")[1]
    return args
end

function bilin2string(bilin)
    if bilin === nothing
        return ""
    else
        return "_bl="*string(bilin)
    end
end

function m2string(m)
    if m == NaN
        return ""
    else
        return string(m)
    end
end

function fname2bilin(model_str)
    args = split(model_str, "_")
    args[end] = split(args[end], ".")[1]
    bilin = split(args[end], "=")[2]
    bilin = split(bilin, "-")
    symargs = string2num(bilin)
    return symargs
end

function string2num(varlist)
    symargs = similar(varlist, Num)
    for var in [u1, u2, u3, d1, d2, d3, l1, l2, l3]
        symargs[varlist .== String(Symbol(var))] .= var
    end
    return symargs
end

function clean_countmap!(cmap)
    if any(collect(keys(cmap)) .=== -0.0) && any(collect(keys(cmap)) .=== 0.0)
        cmap[0.0] += cmap[-0.0]
        cmap[-0.0] = 0
    end
    return cmap
end

_vecvec(mat) = [mat[i,:] for i in 1:size(mat,1)]
_vecvecvec(mat) = [[mat[j,:,i] for i in 1:size(mat,3)] for j in 1:size(mat,1)]

function bilinvals(bilin; odd=(1, 1), even=(1, -1))
    bstring = bilin2string(bilin)
    NrOfUs = length(findall( x -> x == 'u', bstring))
    if isodd(NrOfUs)
        return odd
    elseif iseven(NrOfUs)
        return even
    end
end

function bilinsum(bilin; odd=(2,0), even=(-1,1))
    v1, v2 = bilinvals(bilin; odd=odd, even=even)
    return bilin[1] * v1 + bilin[2] * v2
end

"""
    Save Anomaly Ratios E/N. Two modes: :hist saves E/N as histogram data (2e5 datapoints, saves time and space for huge datasets), :all (store every single datapoint separately, may be prohibitive above 1e8?)
"""
function save_AR(model, proc_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer}, i::Int; folder="", bilin=nothing, ms=NaN)
    fname = model2string(model)
    bilinname = bilin2string(bilin)
    good_idxs = findall(!isnan, proc_rs)
    good_proc_rs = proc_rs[good_idxs]
    good_rs_ws = rs_ws[good_idxs]

    ARh = fit(Histogram, good_proc_rs, FrequencyWeights(good_rs_ws), -1000:0.01:1000)

    @info "Saving..."

    if ms==NaN
        savefile = "./data/DFSZ_models/"*folder*"/"*fname*"/"*fname*bilinname*".jld2"
    else
        savefile = "./data/DFSZ_models/"*folder*"/"*string(ms)*fname*"/"*fname*bilinname*".jld2"
    end

    if isfile(savefile) && i != 1 # i.e. overwrite if file exists before calculation
        cm_old = FileIO.load(savefile, "ARs")
        cm_new = merge(cm_old, ARh)
    else
        cm_new = ARh
    end
    FileIO.save(savefile, Dict("ARs" => cm_new))

    @info "Done!"
end

function get_EoNfunc(model)
    un = unique(model)
    EoN = EoverN(model...)
    myEoNtmp = build_function(EoN, un)
    myEoN = eval(myEoNtmp)
end


function checks(tot, nsamps, un)
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

    m = ifelse(isnothing(nsamps), tot, nsamps)

    @info "You are calculating $m models for the $un Higgs model"
    return nsamps
end

function check_symmetry(EoNlists::Array{Any}; wslists=[ones(size(EoNlist)) for EoNlist in EoNlists])
    EoN = [round.(EoNlist, digits=6) for EoNlist in EoNlists]
    ws = Int.(vcat(wslists...))
    EoN = vcat(EoN...)
    cEoN = countmap(EoN, ws)
    clean_countmap!(cEoN)
    sEoN = Dict(round.(-1 .* keys(cEoN) .+ 3.33333333, digits=4) .=> values(cEoN))
    EoN = [round.(EoNlist, digits=4) for EoNlist in EoNlists]
    EoN = vcat(EoN...)
    cEoN = countmap(EoN, ws)
    clean_countmap!(cEoN)
    mEoN = mergewith( (x,y) -> x-y, cEoN, sEoN)
    if sum(abs.(values(mEoN))) != 0
        if sum(abs.(values(mEoN))) == abs(mEoN[0.0]) + abs(mEoN[-0.0])
            return true
        else
            return mEoN
        end
    else
        return true
    end
end

function check_symmetry(cEoN::Dict{<:Real, <:Real})
    clean_countmap!(cEoN)
    sEoN = Dict(round.(-1.0 .* keys(cEoN) .+ 3.33333333, digits=4) .=> values(cEoN))
    cEoN = Dict(round.(collect(keys(cEoN)), digits=4) .=> collect(values(cEoN)))
    mEoN = mergewith( (x,y) -> x-y, cEoN, sEoN)
    if sum(abs.(values(mEoN))) != 0
        if sum(abs.(values(mEoN))) == abs(mEoN[0.0]) + abs(mEoN[-0.0])
            return true
        else
            @warn "Your resulting EoN Histogram is not symmetric around 5/3. This means some (probably numerical) issue arose. Investigate!"
            return mEoN
        end
    else
        return true
    end
end


function read_AR(model; folder="", bilin=nothing, m=NaN)
    fname = model2string(model)
    bilinname=bilin2string(bilin)
    m = m2string(m)
    return FileIO.load("./data/DFSZ_models/"*folder*"/"*m*fname*"/"*fname*bilinname*".jld2", "ARs")
end


function get_numquads(quads, un, nH)
    numquads = zeros(Int8,length(quads), nH)

    for i in 1:nH
        mydict = Dict(un .=> 0.0)
        mydict[un[i]] = 1.0
        numquads[:,i] = Symbolics.value.(substitute.(quads, (mydict,)))
    end
    chi_s = 1
    as = Vector{SVector{nH, Float64}}(_vecvec(numquads))
    bs = SVector{length(quads), Float64}(2 .* chi_s .* (length.(string.(quads)) .<= 8))
    return as, bs
end


function mysolve(as::AbstractVector{<:SVector{N,T}}, bs, idxs) where {N,T<:Real}
    A_true = copy(hcat(as[idxs]...)')
    solveable = !(≈(det(A_true), 0.0; atol=1e-10))
    A_dummy = A_true .* zero(T) .+ I_static(Val(N), T)
    A = ifelse(solveable, A_true, A_dummy)

    b = bs[idxs]
    b_cand = A \ b
    b_nonsolve = SVector(ntuple(i -> convert(T, NaN), Val(N)))
    ifelse(solveable, b_cand, b_nonsolve)
end


function _make_hist(cm; bins=-10:0.0001:13)
    hi = fit(Histogram, collect(keys(cm)), Weights(collect(values(cm))),bins)#, nbins=5000)
    m = (hi.edges[1] .+ (hi.edges[1][2] - hi.edges[1][1])/2.)[1:end-1]
    #hi.weights = Vector{Real}(hi.weights ./ sum(hi.weights))
    return hi, m
end

function plot_hist!(hist; yaxis=(:log, [0.00000001, :auto]), kwargs...)
    plot!(hist[2], hist[1].weights ./ sum(hist[1].weights); yaxis=yaxis, kwargs...)
    #plot!(hist[2], _gauss_prediction(gauss,hist[1].edges[1]); kwargs..., label="")
    #vline!([1.924], ls=:dash, c=:grey)
end

function fa(ma) #::[GeV]
    return 1e12 * (5.7064e-6/ma)
end

αem() = (1.602176634e-19)^2 / (4*pi* 8.8541878128e-12 * 1.054571817e-34 * 299792458.0)

function gaγγ(EoverN, fa) #::[log(GeV^-1)]
    return log10(αem() / (2.0 * pi * fa) * abs(EoverN - 1.924))
end

function rescale_histogram(tt; edges=-50:0.01:50, mode=:pdf)
    ARrs = abs.(collect(tt.edges...) .+ 1/2 * (tt.edges[1][2] - tt.edges[1][1]) .- 1.924)[1:end-1]
    ARrh = fit(Histogram, ARrs, FrequencyWeights(tt.weights), edges)
    if mode ∈ [:pdf, :probability]
        ARrh = normalize(ARrh; mode=mode)
    end
end

function gag_histogram(tt; ma=40e-6, edges=-16:0.001:-12, mode=:pdf)
    gags = gaγγ.(collect(tt.edges...) .+ 1/2 * (tt.edges[1][2] - tt.edges[1][1]), fa(ma))[1:end-1]
    gagh = fit(Histogram, gags, FrequencyWeights(tt.weights), edges)
    if mode ∈ [:pdf, :probability]
        gagh = normalize(gagh; mode=mode)
    end

    return gagh
end


function cdf(Hist) # I want mw[i] =  sum(gagh.weights[i:end])/ s aka cumsum from the back
    return cumsum(Hist.weights[end:-1:1])[end:-1:1] / sum(Hist.weights)
end

_makeoptions(a,b,c) = [
    [a, a, a],
    [a, a, b],
    [a, b, a],
    [a, b, b],
    [a, b, c]
]

function generate_all_models()
    us = _makeoptions(u1,u2,u3)
    ds = _makeoptions(d1,d2,d3)
    ls = _makeoptions(l1,l2,l3)

    tmp_models = collect(Iterators.product(us, ds, ls));
    length(tmp_models)
    model_list = [u1 for i in 1:9]#similar([1 for i in 1:length(tmp_models)], Any)
    for tmp in collect(Iterators.product(us, ds, ls))
        model_list = hcat(model_list, [tmp[1]...,tmp[2]...,tmp[3]...])
    end
    a = _vecvec(model_list[:,2:end]')
    return a, ones(Int64, length(a))
end

function _makeuniqueoptions(a,b,c; us=nothing, ds=nothing)
    if us === nothing && ds === nothing
        return [[a, a, a],
                [a, a, c],
                [a, b, c]]
    elseif us !== nothing && ds === nothing
        if Symbol.(us) == Symbol.([u1, u1, u1])
            opt = [[u1,u1,u1]]
        elseif Symbol.(us) == Symbol.([u1,u1,u2])
            opt = [[u1,u1,u1], [u1,u1,u2]]
        elseif Symbol.(us) == Symbol.([u1,u2,u3])
            opt = [[u1,u1,u1], [u1,u1,u2], [u1,u2,u3]]
        end
        return[[a, a, a],
               [a, a, b],
               [a, b, c], opt...]
    elseif us!== nothing && ds !== nothing
        return[[a, a, a],
               [a, a, b],
               [a, b, c], us..., ds...]
    end
end

function generate_unique_models()
    us = _makeuniqueoptions(u1,u2,u3)
    ds = _makeuniqueoptions(d1,d2,d3)
    ls = _makeuniqueoptions(l1,l2,l3)

    #tmp_models = collect(Iterators.product(us, ds, ls));
    m = [1,3,1]
    multis = collect(Iterators.product(m,m,m))
    #println(tmp_models)
    #length(tmp_models)
    model_list = [u1 for i in 1:9]#similar([1 for i in 1:length(tmp_models)], Any)
    multi_list = []
    for (i, tmp) in enumerate(collect(Iterators.product(us, ds, ls)))
        mult = *(multis[i]...)
        append!(multi_list, mult)
        model_list = hcat(model_list, [tmp[1]...,tmp[2]...,tmp[3]...])
    end
    return _vecvec(model_list[:,2:end]'), multi_list
end

#=
opt = _makeuniqueoptions(u1,u2,u3)

for nn in [d1,d2,d3]
    myopt = similar(opt)
    for (i, o) in enumerate(opt)
        myopt[i] = [o..., d1]
        for mynum in unique(o)
            append!(myopt, Ref([o..., mynum]))
        end
    end
    opt = deepcopy(myopt)
end
opt
unique(opt)

myk = similar(us, 3)
for (i, u) in enumerate(us)
    myk[i] = [u..., d1]
    for mynum in unique(u)
        append!(myk, Ref([u..., mynum]))
    end
end
us
myk

myl = similar(myk)
for (i, k) in enumerate(myk)
    myl[i] = [k..., d2]
    for mynum in unique(k)
        append!(myl, Ref([k..., mynum]))
    end
end
myl
function makeoptions()
    return [
    [u1,u1,u1,u1,u1,u1,u1,u1,u1],
    [u1,u1,u3,u1,u1,u1,u1,u1,u1],
    [u1,u1,u1,u1,u1,d3,u1,u1,u1],
    [u1,u1,u1,u1,u1,u1,u1,u1,l3],
    [u1,u1,u3,u1,u1,u3,u1,u1,u1],
    [u1,u1,u3,u1,u1,u1,u1,u1,u3],
    [u1,u1,u1,u1,u1,d3,u1,u1,d3],
    [u1,u1,u3,u1,u1,u3,u1,u1,u3],

    ]
end

myus = [
    [u1,u1,u1],
    [u1,u1,u3],
    [u1,u2,u3]
]
mymultis = [1,3,1]

function myds(us, multi)
    ds = [
        [d1,d1,d1],
        [d1,d1,d3],
        [d1,d2,d3]
    ]
    if length(unique(us)) == 1
        m = [1,3,3,3]
        append!(ds, [us[1],us[1],us[1]])
        append!(ds, [us[1], us[1], d3])
        append!(ds, [us[1], d2, d3])
        append!(ds, [us[1], d2, d2])
    elseif length(unique(us)) == 2
        m = vcat([1,3,3,3],[1,3,3,3],[1,3,3])
        for u in unique(us)
            append!(ds, [u,u,u])
            append!(ds, [u, u, d3])
            append!(ds, [u, d2, d3])
            append!(ds, [u, d2, d2])
        end
        append!(ds, [unique(us)[1], unique(us)[2], d3])
        append!(ds, [unique(us)[1], unique(us)[2], unique(us)[2]])
        append!(ds, [unique(us)[1], unique(us)[1], unique(us)[2]])

    end
    multis = multi .* vcat([1,3,1], m)

    return ds
end

function generate_typeviolating_models()
    collect(with_replacement_combinations([u1,u2,u3, d1, d2, d3, l1, l2, l3],9))
end


a = generate_typeviolating_models()
my2 = a[length.(unique.(a)) .== 2]

collect(with_replacement_combinations([u1,u2, u3], 3))
my2[132][1] == u1
numarr = [u1,u2,u3,d1,d2,d3,l1,l2,l3]

powarr = [numarr[1:i] for i in 1:length(numarr)]

powarr[2]
collect(powerset(powarr, 9,9))

if Symbol(u1) in Symbol.(numarr)
    println("tt")
end

for (i, mod) in enumerate(my2[172:end])
    println(mod)
    for (j,par) in enumerate(numarr)
        println(mod[j])
        if Symbol(mod[j]) ∉ Symbol.(numarr[1:j])
            my2[i] = substitute(mod, Dict(mod[j] => par))
        else
            nothing
        end
    end
    break
end

my2 = clean_models(my2)

unique(my2)
function clean_models(my2)
    for i in 1:length(my2)
        replaced = []
        println("\n")
        for j in 1:9
            println(Symbol(my2[i][j]), "\t", Symbol.(replaced))
            if Symbol(my2[i][j]) ∉ Symbol.(replaced)
                my2[i] = substitute(my2[i], Dict(my2[i][j] => numarr[j]))
                append!(replaced, numarr[j])
                println(my2[i])
            else
                nothing
            end
        end
    end
    return my2
end



unique(my2)

my2[13][1]
=#
function I_static(::Val{N}, ::Type{T}) where {N,T<:Real}
    convert(SMatrix{N,N,T}, Diagonal(SVector(ntuple(i -> one(T), Val(N)))))
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

function parallel_randeqn_solve_proc!(
    proc_rs::AbstractVector{<:SVector{N,<:Real}}, rs_ws::AbstractVector{<:Integer},
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer},
    tot::Int
) where N

    idxarr, bnc = make_idx_bnc(N-2)

    Threads.@threads for i in eachindex(proc_rs)
        @inbounds begin
            # When using drawer need to allocate indices. Is there a way around?
            idxs_i =  myidxtup!(idxarr, bnc, rand(1:tot), Val(N-2))#rand_idxs(default_rng(), eachindex(as), Val(N)) # drawer(13131) #
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = round.(r, digits=10)
            rs_ws[i] = prod(ws[idxs_i])
        end
    end
end

function parallel_alleqn_solve_proc!(
    proc_rs::AbstractVector{<:SVector{N,<:Real}}, rs_ws::AbstractVector{<:Integer},
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer},
    tot::Int
) where N

    idxarr, bnc = make_idx_bnc(N-2)

    Threads.@threads for i in eachindex(proc_rs)
        @inbounds begin
            # When using drawer need to allocate indices. Is there a way around?
            idxs_i =  myidxtup!(idxarr, bnc, i, Val(N-2))#rand_idxs(default_rng(), eachindex(as), Val(N)) # drawer(13131) #
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = round.(r, digits=10)
            rs_ws[i] = prod(ws[idxs_i])
        end
    end
end

function parallel_randeqn_solve_proc_fullsol!(
    proc_rs::AbstractVector{<:SVector{N,<:Real}}, EoN_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer},
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer}, terms_all::AbstractVector{<:Num}, terms,
    tot::Int, myEoN, myN
) where N

    idxarr, bnc = make_idx_bnc(N-2)

    Threads.@threads for i in eachindex(proc_rs)
        @inbounds begin
            # When using drawer need to allocate indices. Is there a way around?
            idxs_i =  myidxtup!(idxarr, bnc, rand(1:tot), Val(N-2))#rand_idxs(default_rng(), eachindex(as), Val(N)) # drawer(13131) #
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = round.(r, digits=10)

            EoN_rs[i] = ifelse(-0.00000001 .< myN(r) .< 0.00000001, NaN, myEoN(r))
            terms[i] = Vector(terms_all[idxs_i])
            rs_ws[i] = prod(ws[idxs_i])
        end
    end
end

function parallel_alleqn_solve_proc_fullsol!(
    proc_rs::AbstractVector{<:SVector{N,<:Real}}, EoN_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer},
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer}, terms_all::AbstractVector{<:Num}, terms, sols,
    tot::Int, myEoN, myN
) where N

    idxarr, bnc = make_idx_bnc(N-2)

    Threads.@threads for i in eachindex(proc_rs)
        @inbounds begin
            # When using drawer need to allocate indices. Is there a way around?
            idxs_i =  myidxtup!(idxarr, bnc, i, Val(N-2))#rand_idxs(default_rng(), eachindex(as), Val(N)) # drawer(13131) #
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = round.(r, digits=10)

            EoN_rs[i] = ifelse(-0.00000001 .< myN(r) .< 0.00000001, NaN, myEoN(r))
            #if 82 < EoN_rs[i] < 10000
            #    println(EoN_rs[i])
            #    println(idxs_i)
            #end
            terms[i] = Vector(terms_all[idxs_i])
            sols[i] = bs[idxs_i]
            rs_ws[i] = prod(ws[idxs_i])
        end
    end
end


function save_EoN(model, EoN_countmap; folder="", new_file=false, filename="EoNs")
    
    @info "Saving..."

    
    modstring = model2string(model)
    un = unique(model)
    nH = length(un)

    tpath = "./data/DFSZ_models/"*folder
    group = string(nH)*"/"*modstring
    fname = tpath*"/"*filename*".jld2"
    
    if isfile(fname) == false || new_file
        mkpath(tpath)
        FileIO.save(fname, Dict(group => EoN_countmap))
    else
        jldopen(fname, "a+") do file
            file[group] = EoN_countmap
        end
    end
    

    @info "Done!"
end


function save_full(model, proc_rs::AbstractVector{<:SVector{L,<:Real}}, EoN_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer}, myterms, i::Int; folder="", bilin=nothing, ms=NaN, model_multiplicity=NaN) where L
    
    @info "Saving..."

    if i != 1
        error("i has to be equal to 1. Incremental save not yet supported")
    end
    
    fname = model2string(model)
    bilinname = bilin2string(bilin)
    good_idxs = findall(!isnan, EoN_rs)
    good_EoN_rs = EoN_rs[good_idxs]
    good_proc_rs = proc_rs[good_idxs,:]
    good_rs_ws = rs_ws[good_idxs]
    good_myterms = myterms[good_idxs,:]

    chi_s = 1
    un = unique(model)
    nH = length(un)

    mydict = Dict{Num, Float64}(un .=> 0)

    E = similar(good_EoN_rs)
    N = similar(good_EoN_rs)
    Chis = Matrix{Float64}(undef, length(good_EoN_rs), 10)
    for i in 1:length(good_rs_ws)
        for (j, higgs) in enumerate(un)
            mydict[higgs] = good_proc_rs[i][j]
        end
        chivec = Symbolics.value.(substitute.(model, (mydict,)))
        E[i] = 4 * sum(chivec[1:3]) + 1 * sum(chivec[4:6]) + 3 * sum(chivec[7:9])
        N[i] = 3/2 * sum(chivec[1:6])
        append!(chivec, chi_s)
        Chis[i,:] .= chivec #Construct matrix with chi values, last one is charge of the singlet (always 1)
    end

    
    tpath = "./data/DFSZ_models/"*folder*"/"
    group = fname*"/"*bilinname[2:end]
    
    if isfile(tpath*"full_n"*string(nH)*".h5") == false
        mkpath(tpath)
        h5write(tpath*"full_n"*string(nH)*".h5", "Chis order: u1u2u3d1d2d3l1l2l3s", "")
    end

    h5write(tpath*"full_n"*string(nH)*".h5", group*"/Chis",Chis)
    h5write(tpath*"full_n"*string(nH)*".h5", group*"/E",E)
    h5write(tpath*"full_n"*string(nH)*".h5", group*"/N",N)
    h5write(tpath*"full_n"*string(nH)*".h5", group*"/EoN",good_EoN_rs)
    h5write(tpath*"full_n"*string(nH)*".h5", group*"/multis",good_rs_ws)
    h5write(tpath*"full_n"*string(nH)*".h5", group*"/terms",good_myterms)
    if ms !== NaN
        h5write(tpath*"full_n"*string(nH)*".h5", group*"/ms",ms)
    end

    if model_multiplicity !== NaN
        h5write(tpath*"full_n"*string(nH)*".h5", group*"/model_multiplicity",model_multiplicity)
    end

    @info "Done!"
end

function dict_value_shift!(d::Dict{Ti,Tv},shift::Tv) where Ti where Tv <: Signed
    # In version v1.2+ this can use hasfield(Dict,:vals) instead
    @assert isdefined(d,:vals) "If this fails implimentation of Dict has changed"
    @inbounds for n in 1:length(d.vals)
        d.vals[n] *= shift
    end
end

function h5totxt(dataset::String; folder="./data/DFSZ_models/test/")
    fid = h5open(folder*dataset*".h5")
    io = open(folder*dataset*".txt", "w")
    write(io, "# Explanation of the columns:\n")
    write(io, "# - model string consists of number of Higgs particles and indicators which higgs are meant to be equal (i.e. u1_u1_u1 means that the Higgs particle H_u1 couples to all of the up-type quarks and therefore χHu1, χHu2 and χHu3 obviously also have to be equal, because they are actually the charge of the same particle)\n")
    write(io, "# - model multiplicity indicates 'how often' a model is realized. Completely analogous models like u1_u2_u1 and u1_u1_u3 (i.e. two up-type quarks couple to the same Higgs, one to another) are only calculated once. Since there are three options for this specific case, but e.g. for all up-type quarks coupling to the same Higgs there's only one options, this difference in probability between the two models is accounted for by our parameter 'model multiplicity'.\n")
    write(io, "# - terms multiplicity indicates differences in the probability of the specific terms1-9 arising due to a specific potential. This can be ignored if potentials leading to the same solutions for Higgs charges are deemed equivalent. The difference in probability comes from multiple (quadrilinear) potential terms possibly leading to the same equation (e.g. (Hu1 Hu2') (Hu2 Hd1) leads to the same equation as (Hu1 Hl1) (Hl1' Hd1).\n")
    write(io, "# - Anomaly Ratio is E/N\n")
    write(io, "# - electromagnetic anomaly E\n")
    write(io, "# - color anomaly N\n")
    write(io, "# - χHi: Charge of the Higgs particle coupling to the quark i. If i=s, Charge of the Higgs singlet, always set =1 (wlg for calculating the Anomaly Ratio).\n")
    write(io, "# - eqi: Equation number i. At least nH equations are needed to fix all nH Higgs charges. The nomenclature e.g. u1 is here simplified for χHu1. Together all equations give the matrix that is calculated to obtain solutions for the charges.\n")
    write(io, "# \n")
    write(io, "# model           model multiplicity  terms multiplicity  Anomaly Ratio   E    N     χHu1     χHu2     χHu3     χHd1     χHd2     χHd3     χHl1     χHl2     χHl3      χHs               eq1               eq2               eq3               eq4               eq5               eq6               eq7               eq8               eq9\n")
    for (i, a) in enumerate(fid)
        k = read(a)
        if i != 1
            for tuple in k#collect(values(k))
                dat = tuple[2]
                model = keys(fid)[i]
                mmult = rpad(dat["model_multiplicity"],3)
                for j in 1:length(dat["EoN"])
                    write(io, model*"   "*mmult*"      "*rpad(dat["multis"][j],3)*"                 "*lpad(round(dat["EoN"][j], digits=3),8)*" "*lpad(Int(round(dat["E"][j])),5)*lpad(Int(round(dat["N"][j])),5)*(*(lpad.(round.(dat["Chis"][j,:],digits=3),9)...))*(*(lpad.(filter.(x -> !isspace(x), dat["terms"][j,:]),18)...))*"\n")
                end
            end
            println(keys(fid)[i])
            #println(collect(values(k))[1]["model_multiplicity"])
        end
    end
    close(io)
    close(fid)
end