
function normalisation(ms::msJ.MScontainer)
    factor = 100. / ms.basePeakIntensity
    return ms.int .* factor
end

@recipe function f(ms::msJ.MSscan; method = :relative) 
    seriestype --> :path
    seriescolor --> :red
    label --> ""
    xlabel --> "m/z"
    if method == :relative
        y = normalisation(ms)
        ylabel --> "Intensity (%)"
    else
        y = ms.int
        ylabel --> "Intensity (a.u.)"
    end
    ms.mz, y
end
@recipe function g(ms::msJ.MSscans; method = :relative) 
    seriestype --> :path
    seriescolor --> :red
    label --> ""
    xlabel --> "m/z"
    if method == :relative
        y = normalisation(ms)
        ylabel --> "Intensity (%)"
    elseif method == :absolute
        y = ms.int
        ylabel --> "Intensity (a.u.)"
    end
    ms.mz, y
end
