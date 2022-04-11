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
    uniquads = allquads
    fillones = substitute(uniquads, Dict([var => 1//1 for var in vars]))
    fillones = (sign.(fillones) .== 0) + sign.(fillones) # make sign return 1 if element = 0
    fillones = (fillones .== 1) .* 2//1 .- 1//1 # weird hack to make sure the array contains only Ints, because unique() does not work with negative floats and symbols
    uniquads .*= fillones
    quads = uniquads[isequal.(uniquads, d1+u1) .!= 1//1] # effectively is not equal...
    quads .*= (length.(string.(quads)) .<= 7) .+ 1//1 # weird way of multiplying all quads by 2 if they only contain two different higgs. This is so sum(abs(higgs)) == 4
    #quads .*= 2//1 # now all quads are multiplied by two. This ensures type stability and does not affect the endresult!
    #return quads, [1 for i in 1:length(quads)]
    quads = countmap(quads)
    return collect(keys(quads)), [1 for i in 1:length(keys(quads))]#collect(values(quads))
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

function _model2string(model)
    un = unique(model)
    nH = length(un)

    str = "n"*string(nH)
    for arg in model
        str *= "_"*string(arg) 
    end
    return str
end

function _save_subsets(model, subsets, xsubset)
    svdict = Dict()
    fname = _model2string(model)
    svdict["subsets"] = subsets
    svdict["xsubset"] = xsubset
    @info "Saving..."
    FileIO.save("./data/"*fname*".jld2", svdict)
    @info "Done!"
end

"""
    Save Anomaly Ratios E/N. Two modes: :hist saves E/N as histogram data (2e5 datapoints, saves time and space for huge datasets), :all (store every single datapoint separately, may be prohibitive above 1e8?)
"""
function save_AR(model, ARcm, i::Int)
    AR = countmap(vcat(ARcm...))
    fname = _model2string(model)

    @info "Saving..."

    if isfile("./data/"*fname*".jld2") && i != 1 # i.e. overwrite if file exists before calculation
        cm_old = FileIO.load("./data/"*fname*".jld2", "ARs")
        cm_new = merge(+,cm_old,AR)
    else
        cm_new = AR
    end
    FileIO.save("./data/"*fname*".jld2", Dict("ARs" => cm_new))

    @info "Done!"
end

function get_EoNfunc(model)
    un = unique(model)
    nH = length(un)
    notu1d1 = isequal.(un,u1) .+ isequal.(un,d1) .== false
    EoN = EoverN(model[1:3], model[4:6], model[7:9])
    EoN = substitute(EoN, Dict(u1=>1,d1=>1))
    myEoNtmp = build_function(EoN, un[notu1d1]...)
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


function read_AR(file)
    return FileIO.load("./data/"*file*".jld2", "ARs")
end

_vecvec(mat) = [mat[i,:] for i in 1:size(mat,1)]
_vecvecvec(mat) = [[mat[j,:,i] for i in 1:size(mat,3)] for j in 1:size(mat,1)]




function get_numquads(quads, un, nH)
    notu1d1 = isequal.(un,u1) .+ isequal.(un,d1) .== false
    numquads = zeros(Int8,length(quads), nH-2)

    u1d1dict = Dict{Num, Int8}(un .=> 0)
    u1d1dict[u1] = 1
    u1d1dict[d1] = 1
    bs = Vector{Int8}(-1 .* Symbolics.value.(substitute.(quads, (u1d1dict,))))

    for i in 1:nH-2
        mydict = Dict(un .=> 0.0)
        mydict[un[notu1d1][i]] = 1.0
        numquads[:,i] = Symbolics.value.(substitute.(quads, (mydict,)))
    end
    return _vecvec(SMatrix{length(quads), nH-2, Int8, (nH-2)*length(quads)}(numquads)), bs
end

function mysolve(as, bs, multis, idxs)
    A = copy(hcat(as[idxs]...)')
    b = bs[idxs]
    multi = prod(multis[idxs])
	if _issolvable(A)
	    return (A \ b, multi)
    else
		return (SVector{length(b),Float64}(NaN for i in 1:length(b)), multi)
	end
end

function _issolvable(mat)
    #rank(mat) == size(mat)[1]
    ≈(det(mat), 0.0; atol=1e-10) == false # tolerance chosen arbitrarily. Make sure this fits estimated numerical errors!
end

function calculate_EoverN(myEoN, sol, sol_array)
    AnomalyRatio = myEoN(sol[1]...)
    if AnomalyRatio < 1e3 && AnomalyRatio > -1e3
        for i in 1:sol[2]
            push!(sol_array[Threads.threadid()],AnomalyRatio)
        end
    end
end

function make_sol_array()
    solution_data = Vector{Vector{Float64}}()
    for i in 1:Threads.nthreads()
        push!(solution_data,Float64[])
    end
    return solution_data
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
end