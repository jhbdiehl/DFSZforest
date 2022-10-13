############################################################################################################
#                                                                                                          #
#                             Functions used for saving.                                                   #
#                                                                                                          #
############################################################################################################

"""
    Saving function called by runDFSZ.
"""
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

"""
    Saving function called by runDFSZ_saveall.
"""
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

"""
    Take h5 output from save_full and convert it to a human readable textfile. Don't do this if file is large!
"""
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
        end
    end
    close(io)
    close(fid)
end