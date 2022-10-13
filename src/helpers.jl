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

function get_bilins(model; same_χH_one_theory=false)
    vars=unique(model)
    bilins = unique(sort.(collect(combinations(model,2)), by=x->Symbol(x)))
    bilins = bilins[length.(unique.(bilins)) .== 2]
    o = [iseven(sum(in.(Symbol.(bilin),Ref(Symbol.([u1,u2,u3]))))) for bilin in bilins]
    bi = [bilins[i][1] .+ (-1).^o[i] .* bilins[i][2] for i in 1:length(bilins)]
    if same_χH_one_theory
        return bi
    else
        return vcat(bi, bi .* -1)
    end
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
    nD = length(un)

    str = "n"*string(nD)
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

function check_symmetry(cEoN)
    clean_countmap!(cEoN)
    sEoN = Dict(round.(-1.0 .* keys(cEoN) .+ 3.33333333, digits=4) .=> values(cEoN))
    clean_countmap!(sEoN)
    cEoN = Dict(round.(collect(keys(cEoN)), digits=4) .=> collect(values(cEoN)))
    mEoN = mergewith( (x,y) -> x-y, cEoN, sEoN)
    clean_countmap!(mEoN)
    if sum(abs.(values(mEoN))) != 0
        @warn "Your resulting EoN Histogram is not symmetric around 5/3. This means some (probably numerical) issue arose. Investigate!"
        prettyreturn = sortslices([collect(keys(mEoN))[collect(values(mEoN)) .!= 0] collect(values(mEoN))[collect(values(mEoN)) .!= 0]], dims=1, by=x->x[1], rev=false)
        return prettyreturn
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


function get_numquads(quads, un, nD)
    numquads = zeros(Int8,length(quads), nD)

    for i in 1:nD
        mydict = Dict(un .=> 0.0)
        mydict[un[i]] = 1.0
        numquads[:,i] = Symbolics.value.(substitute.(quads, (mydict,)))
    end
    chi_s = 1
    as = Vector{SVector{nD, Float64}}(_vecvec(numquads))
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
    return log10(αem() / (2.0 * pi * fa) * abs(EoverN - 1.92))
end


function rescale_histogram(tt; edges=-50:0.01:50, mode=:pdf)
    ARrs = abs.(collect(tt.edges...) .+ 1/2 * (tt.edges[1][2] - tt.edges[1][1]) .- 1.92)[1:end-1]
    ARrh = fit(Histogram, ARrs, FrequencyWeights(tt.weights), edges)
    if mode ∈ [:pdf, :probability]
        ARrh = normalize(ARrh; mode=mode)
    end
end

function gag_histogram(tt; ma=40e-6, edges=-16:0.001:-12, mode=:pdf, Cagdist=false)

    if Cagdist
        ttvec = sample(tt.edges[1][1:end-1] .+ 1/2 * (tt.edges[1][2] - tt.edges[1][1]), FrequencyWeights(tt.weights), 10000000)
        ttvec .+= rand(Normal(0.0,0.04), length(ttvec))
        ttcm = countmap(ttvec)
        tt, tmp = _make_hist(ttcm; bins=-50:0.0001:50)
    end

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
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer}, terms_all::AbstractVector{<:Num}, terms, sols,
    tot::Int, myEoN, myN
) where N

    idxarr, bnc = make_idx_bnc(N-2)

    Threads.@threads for i in eachindex(proc_rs)
        @inbounds begin
            # When using drawer need to allocate indices. Is there a way around?
            idxs_i =  myidxtup!(idxarr, bnc, rand(1:tot), Val(N-2))#rand_idxs(default_rng(), eachindex(as), Val(N)) # drawer(13131) #
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = round.(r, digits=10)

            EoN_rs[i] = ifelse(-0.00000001 .< Base.invokelatest(myN, r) .< 0.00000001, NaN, Base.invokelatest(myEoN, r))
            terms[i] = Vector(terms_all[idxs_i])
            sols[i] = bs[idxs_i]
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

            EoN_rs[i] = ifelse(-0.00000001 .< Base.invokelatest(myN, r) .< 0.00000001, NaN, Base.invokelatest(myEoN, r))
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
    nD = length(un)

    tpath = "./data/DFSZ_models/"*folder
    group = string(nD)*"/"*modstring
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


function save_full(model, proc_rs::AbstractVector{<:SVector{L,<:Real}}, EoN_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer}, myterms, i::Int; folder="", bilin=nothing, ms=NaN, model_multiplicity=NaN, full=true) where L
    
    @info "Saving..."

    if i != 1
        error("i has to be equal to 1. Incremental save not yet supported")
    end
    
    fname = model2string(model)
    bilinname = bilin2string(bilin)
    good_idxs = findall(!isnan, EoN_rs)
    good_EoN_rs = EoN_rs[good_idxs]
    good_proc_rs = proc_rs[good_idxs,:]

    alltmps = Vector{Vector{Float64}}()
    for good_proc_r in good_proc_rs
        tmp = similar(good_proc_r)
        tmp[findall(!iszero,good_proc_r)] = good_proc_r[findall(!iszero,good_proc_r)]
        tmp[findall(iszero,good_proc_r)] = good_proc_r[findall(iszero,good_proc_r)].^2
        alltmps = vcat(alltmps, [tmp])
    end

    good_rs_ws = rs_ws[good_idxs]
    good_myterms = myterms[good_idxs,:]

    chi_s = 1
    un = unique(model)
    nD = length(un)

    mydict = Dict{Num, Float64}(un .=> 0)

    E = similar(good_EoN_rs)
    N = similar(good_EoN_rs)
    Chis = Matrix{Float64}(undef, length(good_EoN_rs), 10)
    for i in 1:length(good_rs_ws)
        for (j, higgs) in enumerate(un)
            mydict[higgs] = alltmps[i][j]
        end
        chivec = Symbolics.value.(substitute.(model, (mydict,)))
        Ntmp = 1/2 * sum(chivec[1:6])
        N_DW = rationalize(2 * Ntmp, tol=0.0001)
        N[i] = numerator(N_DW) / 2
        chi_s = denominator(N_DW)
        chivec .*= chi_s
        E[i] = 4/3 * sum(chivec[1:3]) + 1/3 * sum(chivec[4:6]) + sum(chivec[7:9])
        append!(chivec, chi_s)
        Chis[i,:] .= chivec #Construct matrix with chi values, last one is charge of the singlet (always 1)
    end

    
    tpath = "./data/DFSZ_models/"*folder*"/"
    group = fname*"/"*bilinname[2:end]
    if full
        prefix = "full"
    else
        prefix= "samples"
    end

    if isfile(tpath*prefix*"_n"*string(nD)*".h5") == false
        mkpath(tpath)
        h5write(tpath*prefix*"_n"*string(nD)*".h5", "Chis order: u1u2u3d1d2d3l1l2l3s", "")
    end

    h5write(tpath*prefix*"_n"*string(nD)*".h5", group*"/Chis",Chis)
    h5write(tpath*prefix*"_n"*string(nD)*".h5", group*"/E",E)
    h5write(tpath*prefix*"_n"*string(nD)*".h5", group*"/N",N)
    h5write(tpath*prefix*"_n"*string(nD)*".h5", group*"/EoN",good_EoN_rs)
    h5write(tpath*prefix*"_n"*string(nD)*".h5", group*"/multis",good_rs_ws)
    h5write(tpath*prefix*"_n"*string(nD)*".h5", group*"/terms",good_myterms)
    if ms !== NaN
        h5write(tpath*prefix*"_n"*string(nD)*".h5", group*"/ms",ms)
    end

    if model_multiplicity !== NaN
        h5write(tpath*prefix*"_n"*string(nD)*".h5", group*"/model_multiplicity",model_multiplicity)
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
    write(io, "# - eqi: Equation number i. At least nD equations are needed to fix all nD Higgs charges. The nomenclature e.g. u1 is here simplified for χHu1. Together all equations give the matrix that is calculated to obtain solutions for the charges.\n")
    write(io, "# \n")
    write(io, "# model           model multiplicity  terms multiplicity  Anomaly Ratio   E    N     χHu1     χHu2     χHu3     χHd1     χHd2     χHd3     χHl1     χHl2     χHl3      χHs                eq1               eq2               eq3               eq4               eq5               eq6               eq7               eq8               eq9\n")
    for (i, a) in enumerate(fid)
        k = read(a)
        if i != 1
            for tuple in k#collect(values(k))
                dat = tuple[2]
                model = keys(fid)[i]
                mmult = rpad(dat["model_multiplicity"],3)
                for j in 1:length(dat["EoN"])
                    try
                        write(io, model*"   "*mmult*"      "*rpad(dat["multis"][j],3)*"                 "*lpad(round(dat["EoN"][j], digits=3),8)*" "*lpad(round(dat["E"][j], digits=3),8)*lpad(round(dat["N"][j], digits=1),6)*(*(lpad.(round.(dat["Chis"][j,:],digits=3),9)...))*" "*(*(lpad.(filter.(x -> !isspace(x), dat["terms"][j,:]),18)...))*"\n")
                    catch 
                        write(io, model*"   "*mmult*"      "*rpad(dat["multis"][j],3)*"                 "*lpad(round(dat["EoN"][j], digits=3),8)*" "*lpad(round(dat["E"][j], digits=3),8)*lpad(round(dat["N"][j], digits=1),6)*(*(lpad.(round.(dat["Chis"][j,:],digits=3),9)...))*" "*(*(lpad.(filter.(x -> !isspace(x), dat["terms"][j,:]),18)...))*"\n")
                    end
                end
            end
            println(keys(fid)[i])
            #println(collect(values(k))[1]["model_multiplicity"])
        end
    end
    close(io)
    close(fid)
end



function limit_diff(limit, h1, h2)
    cdf1 = cdf(h1)
    cdf2 = cdf(h2)
    a1 = 10^h1.edges[1][findfirst(cdf1 .< limit)]
    a2 = 10^h2.edges[1][findfirst(cdf2 .< limit)]
    2 * abs(a1-a2) / (a1+a2)
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

function plot_EoN(EoN_countmap; bins=-10:0.0001:13)
    h1 = _make_hist(EoN_countmap; bins=bins)
    #plot()
    plot_hist!(h1, yaxis=(:identity, [0.000000001,:auto]))
end

function read_EoN(dataset, models; specifier="EoNs")
    totEoN = countmap([])
    for model in models
        tpath = "./data/DFSZ_models/"*dataset
        name = model2string(model)
        nD = length(unique(model))

        e1 = FileIO.load(tpath*"/"*specifier*".jld2", string(nD)*"/"*name)
        totEoN= mergewith(+,totEoN,e1)
    end
    #v = collect(values(totEoN)) ./ sum(values(totEoN))
    #totEoN = countmap(collect(keys(totEoN)), v)
    return totEoN
end

elementwise_f(f, v) = reduce((x, y) -> f.(x, y), v)

function forecast(crs1, crs2)
    check = collect(keys(crs1))[(!in).(collect(keys(crs1)), Ref(collect(keys(crs2))))]
    if length(check) > 1 # 1 because NaN is not in the array for comparison with !in
        error("crs2 seems to contain $(length(check)-1) elements which are not in crs1. This leads to catastrophe when merging with -, therefore I cannot let you do that!")
        return check
    end
    crs1s = normalize_cm(crs1)
    crs2s = normalize_cm(crs2)
    diff1ss = merge(-, crs2s, crs1s)
    diff1ssm = Dict(collect(keys(diff1ss)) .=> -1. .* collect(values(diff1ss)))
    crs3 = merge(-, crs2s, diff1ssm)
    return crs3
end

function normalize_cm(cm)
    Dict(collect(keys(cm)) .=> collect(values(cm)) ./ sum(collect(values(cm))))
end


"""
    model_list([, nD=:all, compute_equivalent_theories=false])    

Construct models that are to be analysed in a manner that the run function can read. If you want to relax specific conditions (like considering models such as [u1,u2,u3,u1,u1,u1,u1,u1,u1]) you will have to change the generating functions `generate_all_models()` or `generate_unique_models()`.

# Arguments
- `nD`: Number of unique higgs particles in your model. Use either `:all` or list of integers. `:all` is equal to `[3,4,5,6,7,8,9]`.
- `compute_equivalent_theories::Bool`: If `false` e.g. [u1,u1,u3,d1,d1,d1,l1,l1,l1] and [u1,u3,u1,d1,d1,d1,l1,l1,l1] will be treated separately. If `true` only one will be calculated and a note will be made that this solution happens three times (in this example).

# Returns
- `model_vec::Vector{Vector{Num}}`: Vector containing symbolic representations of all models that you set up for calculation.
- `model_multiplicity::Vector{Any}`: Vector of integers tracking how often specific models can arise. Equal to `ones(length(model_vec))`, if `compute_equivalent_theories == true`.

"""
function model_list(;nD=:all, compute_equivalent_theories=false)
    if compute_equivalent_theories
        a, ms = generate_all_models()
    else
        a, ms = generate_unique_models()
    end

    if nD == :all
        model_vec = a
    else
        model_vec = a[length.(unique.(a)) .∈ Ref(nD)]
        model_multiplicity = ms[length.(unique.(a)) .∈ Ref(nD)]
    end
    return model_vec, model_multiplicity
end

function find_max_EoN(cm)
    keysnonan = collect(keys(cm))[(!isnan).(collect(keys(cm)))]
    keysnonan[argmax(abs.(keysnonan .-1.92))]
end

"""
    _EoNft(countmap_χH_results, function_calculating_N, function_calculating_Anomaly_Ratio)

Converts charges to anomaly ratio and cleans up a bit (if N ≈ 0, then E/N is either very large, or completely undefined, if also E ≈ 0).
"""
function _EoNft(crs, myN, myEoN)
    crsN = myN.(collect(keys(crs)))
    crsEoN = myEoN.(collect(keys(crs)))
    EoNf = ifelse.(-0.000000001 .< crsN .< 0.000000001, NaN, crsEoN)
end

"""
    runDFSZ(dataset, model_vec; model_multiplicity=ones(length(model_vec)), log_calculate_all_models_if_below=8, log_sample_models=7, same_χH_one_theory=true)

Calculate only aggregate anomaly ratios of specific models and save them into a JLD2 file in the folder dataset. This function cannot be used to reconstruct possible Higgs charges belonging to a specific anomaly ratio. The probability of having a specific quadrilinear term and a specific bilinear term in the potential are assumed to be equally likely. 

# Arguments:
- `dataset::String`: Folder name to save data to.
- `model_vec::Vector{Vector{Num}}`: Vector containing symbolic representations of all models.
- `model_multiplicity::Vector{Any}`: Vector of integers tracking how often specific models can arise.
- `log_calculate_all_models_if_below::Int`: Calculate all possible models, if for only one bilinear the number of possible combinations for the potential is below 10^ this value. To calculate nD=6 fully you shold set this to 8 or above. Calculating nD=5 shold take order of seconds, nD=6 order of minutes, can't recommend to go to nD=7...
- `sample_models_factor::Float`: When sampling, get `sample_models_factor * total_#_of_models` samples for each specific bilinear.
- `same_χH_one_theory::Bool`: When true, multiplicities down to the level of solutions for the charges are ignored. Aka, if the solution is the same, then we can add up the terms to form a new, unique potential that gives this solution.
- `NDW1::Bool`: When true only theories with domain wall number equal to 1 are considered. Everything else is set to NaN. When false all domain wall numbers are considered. In either case domain wall number is not given explicitly, use runDFSZ_saveall for that.

# Returns
- `EoN_countmaps::Vector{Dicts}`: For each model in model_vec a countmap of (rounded) anomaly ratios is provided.
"""
function runDFSZ(dataset, model_vec; model_multiplicity=ones(length(model_vec)), log_calculate_all_models_if_below=8, sample_models_factor=0.01, same_χH_one_theory=true, NDW1=false, filename="EoNs")
    EoN_countmaps = Array{Any}(undef, length(model_vec))
    @time for (k, model) in enumerate(model_vec)
        println(model)
        bilins = get_bilins(model; same_χH_one_theory=false)
        myEoN = get_EoNfunc(model)
        myN = get_Nfunc(model)
        quads, multis = get_quads(model)
        nD = length(unique(model))
        rslist = Array{Any}(undef, length(bilins))
        tot_min = binomial(length(quads), nD-2)
        @info "Total number of models to be calculated above $(length(bilins)*tot_min)"
        
        @time for (i,bilin) in enumerate(bilins)
            #if tot_min <= 10^log_calculate_all_models_if_below
            if same_χH_one_theory
                terms = quads # since all quads can be constructed out of two bilinears, if you dont care about multiplicities, you can assume only one bilinear wlg
                mult = vcat([1,1], multis)
            else
                terms = vcat(quads, bilins[i+1:end])
                mult = vcat([1,1], multis, ones(Int64, length(bilins[i+1:end]))) #this changes relative bilin to quadrilinear probability
            end
            #else
            #    @info "Sampling route!"
            #    terms = vcat(quads, bilins[1:end .!= i])
            #    mult = vcat([1,1], multis, ones(Int64, length(bilins)-1)) #this changes relative bilin to quadrilinear probability
            #end
            terms = vcat(orthogonality(model), bilin, terms)
            as, bs = get_numquads(terms, unique(model), nD)
            tot = binomial(length(terms)-2,nD-2)
            println(tot)
            
            if tot_min <= 10^log_calculate_all_models_if_below
                proc_rs = similar(as, tot)
                EoN_rs = similar(bs, tot)
                rs_ws = similar(mult, length(proc_rs))
                parallel_alleqn_solve_proc!(proc_rs, rs_ws, as, bs, mult, tot)
                rslist[i] = countmap(proc_rs, rs_ws)
            else
                factor= sample_models_factor
                chunk = Int(round(factor * tot))
                #chunk = 10^(log_sample_models-1)
                rslist[i] = countmap([])
                @info "Too many models to be calculated. Proceeding with calculation in chunks."
                @info "I will compute $(chunk*1) models in total. (model $k, bilin $i)"
                for l in 1:1
                    @info "Calculating chuck $l of 1"
                    proc_rs = similar(as, chunk)
                    EoN_rs = similar(bs, chunk)
                    rs_ws = similar(mult, length(proc_rs))
                    @time parallel_randeqn_solve_proc!(proc_rs, rs_ws, as, bs, mult, tot)
                    dummy = countmap(proc_rs, rs_ws)
                    rslist[i] = mergewith(+, rslist[i], dummy)
                end
            end
            
        end

        crs = mergewith(+,rslist...)
        

        # This adds all of the negative solutions from hermitian conjugated bilinears!
        if same_χH_one_theory
            crsk = collect(keys(crs))
            #crskm = -1 .* crsk
            #crskmgood0 = replace.(crskm, -0.0 => 0.0)
            crskgood0 = replace.(crsk, -0.0 => 0.0)
            #crsm = countmap(crskmgood0, collect(values(crs)))
            crs = countmap(crskgood0, collect(values(crs)))
            #crs = mergewith(+, crs, crsm)
        else
            crsk = collect(keys(crs))
            crskgood0 = replace.(crsk, -0.0 => 0.0)
            crs = countmap(crskgood0, collect(values(crs)))
        end
        #crs = collect(keys(crsi2))


        #crs = convert(Vector{Vector{Float64}}, crs)
        #@time alltmps = Vector{Vector{Float64}}()
        #@time for cr in crs
        #    tmp = similar(cr)
        #    tmp[findall(!iszero,cr)] = cr[findall(!iszero,cr)]
        #    tmp[findall(iszero,cr)] = cr[findall(iszero,cr)].^2
        #    alltmps = vcat(alltmps, [tmp])
        #end
        #crs = permutedims(hcat(crs...))
        #crs[findall(iszero,crs)] = crs[findall(iszero,crs)].^2 # -0.0 and 0.0 is the same thing!
        #crs = collect(eachrow(crs))
        #crs = countmap(crs, collect(values(crsi2)))

        #crs = mergewith(+,rslist...)

        #@time Threads.@threads for i in 1:1 # This loop is necessary, otherwise usage of myN generated function from above leads to world age error when multithreading
        #EoNf = _EoNft(crs, myN, myEoN)
        crsN = Base.invokelatest.(myN, collect(keys(crs)))
        crsEoN = Base.invokelatest.(myEoN, collect(keys(crs)))
        #EoNf = ifelse.(-0.000000001 .< crsN .< 0.000000001, NaN, crsEoN)
        if NDW1
            truN = similar(crsN)
            truN[isnan.(crsN)] .= NaN
            truN[isnan.(crsN) .== 0] .= numerator.(rationalize.(2 .* crsN[isnan.(crsN) .== 0], tol=0.0001)) ./ 2
            EoNf = ifelse.(0.49999 .< abs.(truN) .< 0.50001, crsEoN, NaN)
        else
            EoNf = ifelse.(-0.000000001 .< crsN .< 0.000000001, NaN, crsEoN)
        end

        if same_χH_one_theory
            cEoN = countmap(round.(EoNf, digits=6), ones(length(collect(values(crs)))) .* model_multiplicity[k]) #model_multiplicity[k] ./ collect(values(crs)))#
        else
            cEoN = countmap(round.(EoNf, digits=6), collect(values(crs)) .* model_multiplicity[k])
        end
        
        #clean_countmap!(cEoN)
        #EoN_countmaps[k] = crs#cEoN
        save_EoN(model, cEoN; folder=dataset, filename=filename)
        #end

    end
    #return EoN_countmaps
end


"""
    runDFSZ_saveall(dataset, model_vec; log_calculate_all_models_if_below=8, log_sample_models=7)

Calculate full solutions to all LES and store them in H5 files so potential terms, charges and anomaly ratios can be related to one another. This function is slow and takes vast storage space. Do not ever use this for computing EoN distributions above nD=5!

# Arguments:
See runDFSZ().
"""
function runDFSZ_saveall(dataset, model_vec; model_multiplicity=ones(length(model_vec)), log_calculate_all_models_if_below=8, sample_models_factor=0.01)
    @time for (k, model) in enumerate(model_vec)
        println(model)
        bilins = get_bilins(model)
        myEoN = get_EoNfunc(model)
        myN = get_Nfunc(model)
        quads, multis = get_quads(model)
        nD = length(unique(model))
        tot_min = binomial(length(quads), nD-2)
        @info "Total number of models to be calculated above $(length(bilins)*tot_min)"
        
        @time for (i,bilin) in enumerate(bilins)
            #if tot_min <= 10^log_calculate_all_models_if_below
            terms = vcat(quads, bilins[i+1:end])
            mult = vcat([1,1] , multis, ones(Int64, length(bilins[i+1:end]))) #this changes relative bilin to quadrilinear probability
            #else
            #    @info "Sampling route!"
            #    terms = vcat(quads, bilins[1:end .!= i])
            #    mult = vcat([1,1], multis, ones(Int64, length(bilins)-1)) #this changes relative bilin to quadrilinear probability
            #end
            terms = vcat(orthogonality(model), bilin, terms)
            as, bs = get_numquads(terms, unique(model), nD)
            tot = binomial(length(terms)-2,nD-2)
            
            if tot_min <= 10^log_calculate_all_models_if_below
                proc_rs = similar(as, tot)
                EoN_rs = similar(bs, tot)
                rs_ws = similar(mult, length(proc_rs))
                myterms = fill([u1 for i in 1:nD],tot)
                sols = fill([0 for i in 1:nD],tot)
                parallel_alleqn_solve_proc_fullsol!(proc_rs, EoN_rs, rs_ws, as, bs, mult, terms, myterms, sols, tot, myEoN, myN)
            else
                factor=sample_models_factor
                chunk = Int(round(factor * tot))
                #chunk = 10^(log_sample_models-1)
                @info "Too many models to be calculated. Proceeding with calculation in chunks."
                @info "I will compute $(chunk*1) models in total."
                for l in 1:1
                    @info "Calculating chuck $l of 1"
                    proc_rs = similar(as, chunk)
                    EoN_rs = similar(bs, chunk)
                    rs_ws = similar(mult, length(proc_rs))
                    myterms = fill([u1 for i in 1:nD],tot)
                    sols = fill([0 for i in 1:nD],tot)
                    @time parallel_randeqn_solve_proc_fullsol!(proc_rs, EoN_rs, rs_ws, as, bs, mult, terms, myterms, sols, tot, myEoN, myN)
                end
            end

            stringterms = [string.(term) for term in myterms]
            stringsols = ["=" .*string.(sol) .*"s" for sol in sols]
            stringsols = permutedims(hcat(stringsols...))
            stringterms = permutedims(hcat(stringterms...))
            stringterms .*= stringsols
            #println(stringterms)

            # save_full somehow is slow if not precompiled! Let runDFSZ_saveall run first for e.g. n=3, which is inexpensive and later for higher ns. Somehow this makes a difference!
            save_full(model, proc_rs, EoN_rs, rs_ws, stringterms, 1; folder=dataset, bilin=bilin, ms=mult, model_multiplicity=model_multiplicity[k], full=(tot_min <= 10^log_calculate_all_models_if_below))
        end
    end
end