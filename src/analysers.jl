function fa(ma) #::[GeV]
    return 1e12 * (5.7064e-6/ma)
end

αem() = (1.602176634e-19)^2 / (4*pi* 8.8541878128e-12 * 1.054571817e-34 * 299792458.0)

function gaγγ(EoverN, fa) #::[log(GeV^-1)]
    return log10(αem() / (2.0 * pi * fa) * abs(EoverN - 1.92))
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


function find_max_EoN(cm)
    keysnonan = collect(keys(cm))[(!isnan).(collect(keys(cm)))]
    keysnonan[argmax(abs.(keysnonan .-1.92))]
end

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





function rescale_histogram(tt; edges=-50:0.01:50, mode=:pdf)
    ARrs = abs.(collect(tt.edges...) .+ 1/2 * (tt.edges[1][2] - tt.edges[1][1]) .- 1.92)[1:end-1]
    ARrh = fit(Histogram, ARrs, FrequencyWeights(tt.weights), edges)
    if mode ∈ [:pdf, :probability]
        ARrh = normalize(ARrh; mode=mode)
    end
end

function gag_histogram(tt; ma=40e-6, edges=-16:0.001:-12, mode=:pdf, Cagdist=false)

    if Cagdist
        ttvec = sample(tt.edges[1][1:end-1] .+ 1/2 * (tt.edges[1][2] - tt.edges[1][1]), FrequencyWeights(tt.weights), 100000000)
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

function clean_countmap!(cmap)
    if any(collect(keys(cmap)) .=== -0.0) && any(collect(keys(cmap)) .=== 0.0)
        cmap[0.0] += cmap[-0.0]
        cmap[-0.0] = 0
    end
    return cmap
end