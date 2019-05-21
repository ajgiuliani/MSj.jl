
function normalisation(ms::msJ.MScontainer)
    factor = 100. / ms.basePeakIntensity
    return ms.int .* factor
end

function normalisation(cr::msJ.Chromatogram)
    factor = 100. / cr.maxic
    return cr.ic .* factor
end
function scaling(cr::msJ.Chromatogram)
    return cr.rt ./ 60.
end

@recipe function f(ms::msJ.MSscan; method = :relative) 
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

@recipe function h(cr::msJ.Chromatogram; method = :relative) 
    seriestype  --> :path
    seriescolor --> :blue
    fillrange   --> 0 
    fillalpha   --> 0.3
    label       --> ""
    xlabel      --> "time (mins)"
    if method == :relative
        y = normalisation(cr)
        ylabel  --> "Intensity (%)"
    elseif method == :absolute
        y = cr.ic
        ylabel  --> "Intensity (a.u.)"
    end
    scaling(cr), y
end
