############################################################################################################
#                                                                                                          #
#  Solve the LSEs from the explicit breaking potential and obtain the PQ charges. The heart of the code.   #
#                                                                                                          #
############################################################################################################
# Special thanks to Oliver Schulz (https://github.com/oschulz), most of the code here belongs to him.

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