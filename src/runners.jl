############################################################################################################
#                                                                                                          #
#                   Essential high level functions to beexecuted by the user.                              #
#                                                                                                          #
############################################################################################################

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
        
        crsN = Base.invokelatest.(myN, collect(keys(crs)))
        crsEoN = Base.invokelatest.(myEoN, collect(keys(crs)))
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

        EoN_countmaps[k] = cEoN
        save_EoN(model, cEoN; folder=dataset, filename=filename)
    end
    return EoN_countmaps
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