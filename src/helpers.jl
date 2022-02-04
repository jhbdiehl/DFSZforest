"""
    Careful! This ignores different quadrilinears leading to the same condition! This is maybe not the behavior we want!
    (Also this function is incredibly dirty, holy cow!)
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
    allquads_conj
    allquads = vcat(allquads, allquads_conj)
    allquads = [sum(sum(allquads[i,:])) for i in 1:size(allquads)[1]]
    ut = tril(permutedims(hcat([isequal.(allquads, -allquads[i]) for i in 1:length(allquads)]...)),-1) # make sure to remove duplicates like a = -b
    fut = findall(ut)
    for i in 1:length(fut)
        allquads[fut[i][1]] = allquads[fut[i][2]]
    end
    uniquads = unique(allquads) # remove all duplicates
    fillones = substitute(uniquads, Dict([var => 1//1 for var in vars]))
    fillones = (sign.(fillones) .== 0) + sign.(fillones) # make sign return 1 if element = 0
    fillones = (fillones .== 1) .* 2//1 .- 1//1 # weird hack to make sure the array contains only Ints, because unique() does not work with negative floats and symbols
    uniquads .*= fillones
    quads = uniquads[isequal.(uniquads, d1+u1) .!= 1//1] # effectively is not equal...
    quads .*= (length.(string.(quads)) .<= 7) .+ 1//1 # weird way of multiplying all quads by 2 if they only contain two different higgs. This is so sum(abs(higgs)) == 4
    quads .*= 2//1 # now all quads are multiplied by two. This ensures type stability and does not affect the endresult!
end

"""
    Calculate Anomaly Ratio, given PQ charges of up- and down quarks as well as leptons.
"""
function EoverN(ups, downs, leptons)
    2//3 + 2//1 * (sum(ups)+sum(leptons)) //(sum(ups)+sum(downs))
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