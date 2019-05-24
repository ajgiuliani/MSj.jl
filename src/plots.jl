"""
Plotting module for MScontainer data type (MSscan, MSscans and Chromatogram).
```julia-repl
julia> plot(scans[1])
julia> plot(chr)
```
"""
module plots


using Plots
using RecipesBase


"""
    normalisation(ms::msJ.MScontainer)
Normalization function for plotting mass spectra in relative intensity.
"""
function normalisation(ms::msJ.MScontainer)
    factor = 100. / ms.basePeakIntensity
    return ms.int .* factor
end

"""
    normalisation(cr::msJ.Chromatogram)
Normalization function for plotting chromatograms in raltive intensity.
"""
function normalisation(cr::msJ.Chromatogram)
    factor = 100. / cr.maxic
    return cr.ic .* factor
end

"""
    scaling(cr::msJ.Chromatogram)
Scaling function to display retention times of chromatograms in minutes instead of seconds.
"""
function scaling(cr::msJ.Chromatogram)
    return cr.rt ./ 60.
end


"""
    f(ms::MSscan; method = :relative) 
Allows plotting directly mass spectra MSscan. The defaults relative intensity plotting may be changed by setting method = :absolute.
"""
@recipe function f(ms::MSscan; method = :relative) 
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

"""
    g(ms::msJ.MSscans; method = :relative)
Allows plotting directly mass spectra MSscans. The defaults relative intensity plotting may be changed by setting method = :absolute.
"""
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


"""
    h(cr::msJ.Chromatogram; method = :relative) 
Allows plotting directly chromatograms. The defaults relative intensity plotting may be changed by setting method = :absolute.
"""
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



end
