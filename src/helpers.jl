"""
    Careful! This ignores different quadrilinears leading to the same condition! This is maybe not the behavior we want! (Not true anymore!)
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
    elseif valp1 + valp2 == 0
        quads = quads[isequal.(quads, 2//1 * p1 - 2//1 * p2) .+ isequal.(quads, 2//1 * p2 - 2//1 * p1) .!= 1//1] # effectively is not equal...
    else
        error("I think your bilinear values for the bilinear you used explicitly are not good!")
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

function bilinvals(bilin; odd=(1, 1), even=(1, -1))
    bstring = bilin2string(bilin)
    NrOfUs = length(findall( x -> x == 'u', bstring))
    if isodd(NrOfUs)
        return odd
    elseif iseven(NrOfUs)
        return even
    end
end

function bilinsum(bilin)
    v1, v2 = bilinvals(bilin)
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
    bs = SVector{length(quads), Float64}(-1 .* Symbolics.value.(substitute.(quads, (u1d1dict,))) .+ (2 .* (length.(string.(quads)) .<= 7))) # i.e. set up les normally and add +2 or +0 depending on if its a bilinear or a quadrilinear

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
    a = _vecvec(model_list[:,2:end]')
    return a, ones(Int8, length(a))
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

=#

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
    idxarr = similar([1],N);
    return idxarr, bnc
end

function parallel_randeqn_solve_proc!(
    proc_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer},
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer},
    tot::Int, myEoN
) where N

    idxarr, bnc = make_idx_bnc(N)

    Threads.@threads for i in eachindex(proc_rs)
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

    Threads.@threads for i in eachindex(proc_rs)
        @inbounds begin
            # When using drawer need to allocate indices. Is there a way around?
            idxs_i =  myidxtup!(idxarr, bnc, i, Val(N))#rand_idxs(default_rng(), eachindex(as), Val(N)) # drawer(13131) #
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = myEoN(r)
            rs_ws[i] = prod(ws[idxs_i])
        end
    end
end

function parallel_randeqn_solve_proc_fullsol!(
    proc_rs::AbstractVector{<:SVector{N,<:Real}}, EoN_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer},
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer},
    tot::Int, myEoN
) where N

    idxarr, bnc = make_idx_bnc(N)

    Threads.@threads for i in eachindex(proc_rs)
        @inbounds begin
            # When using drawer need to allocate indices. Is there a way around?
            idxs_i =  myidxtup!(idxarr, bnc, rand(1:tot), Val(N))#rand_idxs(default_rng(), eachindex(as), Val(N)) # drawer(13131) #
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = r

            EoN_rs[i] = myEoN(r)
            rs_ws[i] = prod(ws[idxs_i])
        end
    end
end

function parallel_alleqn_solve_proc_fullsol!(
    proc_rs::AbstractVector{<:SVector{N,<:Real}}, EoN_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer},
    as::AbstractVector{<:SVector{N,<:Real}}, bs::AbstractVector{<:Real}, ws::AbstractVector{<:Integer},
    tot::Int, myEoN
) where N

    idxarr, bnc = make_idx_bnc(N)

    Threads.@threads for i in eachindex(proc_rs)
        @inbounds begin
            # When using drawer need to allocate indices. Is there a way around?
            idxs_i =  myidxtup!(idxarr, bnc, i, Val(N))#rand_idxs(default_rng(), eachindex(as), Val(N)) # drawer(13131) #
            r = mysolve(as, bs, idxs_i)
            proc_rs[i] = r

            EoN_rs[i] = myEoN(r)
            rs_ws[i] = prod(ws[idxs_i])
        end
    end
end




function save_full(model, proc_rs::AbstractVector{<:SVector{L,<:Real}}, EoN_rs::AbstractVector{<:Real}, rs_ws::AbstractVector{<:Integer}, i::Int; folder="", bilin=nothing, valp1=1, valp2=1, ms=NaN) where L
    
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

    chi_s = ( abs(valp1) + abs(valp2) ) / 2
    un = unique(model)
    nH = length(un)

    mydict = Dict{Num, Float64}(un .=> 0)
    mydict[bilin[1]]= valp1
    mydict[bilin[2]]= valp2

    notp1p2 = isequal.(un,bilin[1]) .+ isequal.(un,bilin[2]) .== false

    E = similar(good_EoN_rs)
    N = similar(good_EoN_rs)
    Chis = Matrix{Float64}(undef, length(good_EoN_rs), 10)
    for i in 1:length(good_rs_ws)
        for (j, higgs) in enumerate(un[notp1p2])
            mydict[higgs] = good_proc_rs[i][j]
        end
        chivec = Symbolics.value.(substitute.(model, (mydict,)))
        E[i] = 4 * sum(chivec[1:3]) + 1 * sum(chivec[4:6]) + 3 * sum(chivec[7:9])
        N[i] = 3/2 * sum(chivec[1:6])
        append!(chivec, chi_s)
        Chis[i,:] .= chivec #Construct matrix with chi values, last one is charge of the singlet (always 1)
    end

    tpath = "./data/DFSZ_models/"*folder
    if ms == NaN
        group = fname*"/"*bilinname[2:end]
    else
        group = string(ms)*fname*"/"*bilinname[2:end]
    end
    
    if isfile(tpath*"full_n"*string(nH)*".h5") == false
        mkpath(tpath)
        h5write(tpath*"full_n"*string(nH)*".h5", "Chis order: u1u2u3d1d2d3l1l2l3s", "")
    end

    h5write(tpath*"full_n"*string(nH)*".h5", group*"/Chis",Chis)
    h5write(tpath*"full_n"*string(nH)*".h5", group*"/E",E)
    h5write(tpath*"full_n"*string(nH)*".h5", group*"/N",N)
    h5write(tpath*"full_n"*string(nH)*".h5", group*"/EoN",good_EoN_rs)
    h5write(tpath*"full_n"*string(nH)*".h5", group*"/multis",good_rs_ws)

    @info "Done!"
end