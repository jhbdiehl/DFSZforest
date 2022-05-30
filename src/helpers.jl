"""
    Careful! This ignores different quadrilinears leading to the same condition! This is maybe not the behavior we want!
    (Also this function is incredibly dirty, holy cow!)
"""
function get_quads(model; old=false, p1=u1, p2=d1, valp1=1, valp2=1)
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
    
    if old == true
        allquads = collect(powerset(doublets, 2,2)) # This does not take conditions like (ue) (ue) into account, that dont produce new E/N values but should be important for multiplicities!
    else
        allquads = collect(with_replacement_combinations(doublets, 2))
    end
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
    #println(length(allquads))
    #uniquads = unique(allquads, dims=1) # remove all duplicates
    #println(length(uniquads))
    quads = allquads
    fillones = substitute(quads, Dict([var => 1//1 for var in vars]))
    fillones = (sign.(fillones) .== 0) + sign.(fillones) # make sign return 1 if element = 0
    fillones = (fillones .== 1) .* 2//1 .- 1//1 # weird hack to make sure the array contains only Ints, because unique() does not work with negative floats and symbols
    quads .*= fillones
    quads .*= (length.(string.(quads)) .<= 7) .+ 1//1 # weird way of multiplying all quads by 2 if they only contain two different higgs. This is so sum(abs(higgs)) == 4
    if valp1 + valp2 == 2
        quads = quads[isequal.(quads, 2//1 * p1 + 2//1 * p2) .!= 1//1] # effectively is not equal...
    else
        quads = quads[isequal.(quads, 2//1 * p1 - 2//1 * p2) .!= 1//1] # effectively is not equal...
    end
    #quads .*= 2//1 # now all quads are multiplied by two. This ensures type stability and does not affect the endresult!
    #return quads, [1 for i in 1:length(quads)]
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
    doublets .*= Int64.(sig).* 2//1 .* 1//2
    doublets = [doublets[i,:] for i in 1:size(doublets, 1)]
    return doublets
end
"""
    Calculate Anomaly Ratio, given PQ charges of up- and down quarks as well as leptons.
"""
function EoverN(u1, u2, u3, d1, d2, d3, l1, l2, l3)
    2/3 + 2 * (u1 + u2 + u3 + l1 + l2 + l3) / (u1 + u2 + u3 + d1 + d2 + d3)
end

function checks(AR::Vector)
    sum(typeof.(AnomalyRatio) .== Rational{Int64}) == length(AnomalyRatio) ? nothing : error(
    "Your final AnomalyRatio array has not been calculated fully analytically, at least one element does not have type Rational{Int64}.")
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
        return "_bl="*string(bilin[1])*"-"*string(bilin[2])
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

function bilinvals(bilin)
    bstring = bilin2string(bilin)
    NrOfUs = length(findall( x -> x == 'u', bstring))
    if isodd(NrOfUs)
        return 1, 1
    elseif iseven(NrOfUs)
        return 1, -1
    end
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

function get_EoNfunc(model; p1=u1, p2=d1, valp1=1, valp2=1)
    un = unique(model)
    nH = length(un)
    notu1d1 = isequal.(un,p1) .+ isequal.(un,p2) .== false
    EoN = EoverN(model...)
    EoN = substitute(EoN, Dict(p1=>valp1,p2=>valp2))
    myEoNtmp = build_function(EoN, un[notu1d1])
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


function read_AR(model; folder="", bilin=nothing, m=NaN)
    fname = model2string(model)
    bilinname=bilin2string(bilin)
    m = m2string(m)
    return FileIO.load("./data/DFSZ_models/"*folder*"/"*m*fname*"/"*fname*bilinname*".jld2", "ARs")
end

_vecvec(mat) = [mat[i,:] for i in 1:size(mat,1)]
_vecvecvec(mat) = [[mat[j,:,i] for i in 1:size(mat,3)] for j in 1:size(mat,1)]


function get_numquads(quads, un, nH; p1=u1, p2 = d1, valp1=1, valp2=1)
    notu1d1 = isequal.(un,p1) .+ isequal.(un,p2) .== false
    numquads = zeros(Int8,length(quads), nH-2)

    u1d1dict = Dict{Num, Int8}(un .=> 0)
    u1d1dict[p1] = valp1
    u1d1dict[p2] = valp2
    bs = SVector{length(quads), Float64}(-1 .* Symbolics.value.(substitute.(quads, (u1d1dict,))))

    for i in 1:nH-2
        mydict = Dict(un .=> 0.0)
        mydict[un[notu1d1][i]] = 1.0
        numquads[:,i] = Symbolics.value.(substitute.(quads, (mydict,)))
    end
    as = Vector{SVector{nH-2, Float64}}(_vecvec(numquads))
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
    _vecvec(model_list[:,2:end]')
end

_makeuniqueoptions(a,b,c) = [
    [a, a, a],
    [a, a, b],
    [a, b, c]
]

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