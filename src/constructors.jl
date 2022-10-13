############################################################################################################
#                                                                                                          #
#                           Construct E, N and corresponding functions.                                    #
#                                                                                                          #
############################################################################################################

"""
    Calculate Anomaly Ratio, given PQ charges of up- and down quarks as well as leptons.
"""
function EoverN(u1, u2, u3, d1, d2, d3, l1, l2, l3)
    2/3 + 2 * (u1 + u2 + u3 + l1 + l2 + l3) / (u1 + u2 + u3 + d1 + d2 + d3)
end

function N(u1,u2,u3,d1,d2,d3,l1,l2,l3)
    1/2 * (u1 + u2 + u3 + d1 + d2 + d3)
end

function get_EoNfunc(model)
    un = unique(model)
    EoN = EoverN(model...)
    myEoNtmp = build_function(EoN, un)
    myEoN = eval(myEoNtmp)
end

function get_Nfunc(model)
    un = unique(model)
    tN = N(model...)
    myNtmp = build_function(tN, un)
    myN = eval(myNtmp)
end


############################################################################################################
#                                                                                                          #
#  Construct conditions from explicit symmetry breaking potential, analytically as well as numerically.    #
#                                                                                                          #
############################################################################################################

"""
    Takes a specific model and spits out all possible quadrilinears together with their multiplicity (given this model).
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

############################################################################################################
#                                                                                                          #
#                            Construct all possible Yukawa sectors.                                        #
#                                                                                                          #
############################################################################################################

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