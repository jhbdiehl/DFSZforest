############################################################################################################
#                                                                                                          #
#                    Helper functions that didn't fit into another file.                                   #
#                                                                                                          #
############################################################################################################

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

function bilin2string(bilin)
    if bilin === nothing
        return ""
    else
        return "_bl="*string(bilin)
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

_vecvec(mat) = [mat[i,:] for i in 1:size(mat,1)]
_vecvecvec(mat) = [[mat[j,:,i] for i in 1:size(mat,3)] for j in 1:size(mat,1)]

elementwise_f(f, v) = reduce((x, y) -> f.(x, y), v)

function I_static(::Val{N}, ::Type{T}) where {N,T<:Real}
    convert(SMatrix{N,N,T}, Diagonal(SVector(ntuple(i -> one(T), Val(N)))))
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