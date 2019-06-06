"""
Processing functions submodule. 
"""

using Statistics     # used for Perasons correlation calculation
using LsqFit         # used for curve fitting
using DSP            # used for convolution



# User Interface.
# ---------------

export smooth, centroid


# Mass spectra
"""
    smooth(scan::MScontainer; method::MethodType=SG(5, 9))
Smooth the intensity of the input data and returns a similar structure.
# Examples
```julia-repl
julia> smoothed_data = msJ.smooth(scans)
msJ.MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
```
"""
function smooth(scan::MScontainer; method::MethodType=SG(5, 9, 0))
    if method isa msJ.SG
        return savitzky_golay_filtering(scan, method.order, method.window, method.derivative)
    end  
end



"""
    savitzky_golay_filtering(scan::msJ.MScontainer, order::Int, window::Int, deriv::Int)
Savinsky and Golay filtering of mz and int data within the MSscan(s) container.
"""
function savitzky_golay_filtering(scan::MScontainer, order::Int, window::Int, deriv::Int)
    if window % 2 != 1
        return ErrorException("Window has to be an odd number.")
    elseif window < 1
        return ErrorException("window has to be a positive number.")
    elseif window < order + 2
        return ErrorException("window is too small for the order.")
    end
    order_range = range(1, length=(order+1))
    half_window = Int( (window-1) / 2 )

    b = zeros(window, order+1)

    for i = 0:order
        b[:,i+1] = [x for x = -half_window:half_window].^(i)
    end
    
    m = b * LinearAlgebra.pinv(b' * b)
    coefs = m[:,deriv + 1] * factorial(deriv)
    yfirst = scan.int[1]*ones(half_window)
    ylast = scan.int[end]*ones(half_window-1)
    pad = vcat(yfirst, scan.int, ylast)
    y = conv(coefs[end:-1:1], pad)[2 * half_window + 1 : end - 2 * half_window]
    
    basePeakIntensity = ceil(maximum(y))
    basePeakIndex = num2pnt(y, basePeakIntensity)
    basePeakMz = scan.mz[basePeakIndex]
    
    if scan isa MSscan
        return MSscan(scan.num, scan.rt, scan.tic, scan.mz, y, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy)
    elseif scan isa MSscans

        return MSscans(scan.num, scan.rt, scan.tic, scan.mz, y, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, scan.s)
    end
end
 

"""
    centroid(scan::MScontainer; resolution::Symbol = :medium, R::Real = 4500., shape::Symbol = :gauss, threshold::Real = 0.2)
Checks the file extension and calls the right function to load the mass spectra if it exists. Returns an array of individual mass spectra 
# Examples
```julia-repl
julia> reduced_data = find_peaks(scans)
MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
```
"""
function centroid(scan::MScontainer; method::MethodType=TBPD(:gauss, 4500., 0.2) )
    if method isa TBPD
        #return tbpd(scan, method.shape, convert(Float64, method.resolution), convert(Float64,method.threshold))
        ∆mz = 500.0 / method.resolution       # according to ∆mz / mz  = R, we take the value @ m/z 500
        if method.shape == :gauss
            return tbpd(scan, gauss, ∆mz, convert(Float64,method.threshold))
        elseif method.shape == :lorentz
            return tbpd(scan, lorentz, ∆mz, convert(Float64,method.threshold))
        elseif method.shape == :voigt
            return tbpd(scan, voigt, ∆mz, convert(Float64,method.threshold))
        else
            ErrorException("Unsupported peak profile. Use :gauss, :lorentz or :voigt.")
        end

#    elseif method isa SNRA()
#        return snra(scan, method.threshold)
#    elseif method isa CWT()
#        return cwt(scan, method.threshold)
#    else
#        ErrorException("Unsupported method.")
    end
    
end

"""
    gauss(x::Float64, p::Vector{Float64})
Gaussian shape function used by the TBPD method
"""
function gauss(x::Float64, p::Vector{Float64})
    # Gaussian shape function
    # width            = p[1]
    # x0               = p[2]
    # height           = p[3]
    # background level = p[4]
    # model(x, p) = p[4] + p[3] * exp(- ( (x-p[2])/p[1] )^2)
    return  p[4] + p[3] * exp(- ( (x-p[2])/p[1] )^2)
end

"""
    lorentz(x::Float64, p::Vector{Float64})
Cauchy-Lorentz shape function used by the TBPD method
"""
function lorentz(x::Float64, p::Vector{Float64})
    # Lorentzian shape function
    # width            = p[1]
    # x0               = p[2]
    # height           = p[3]
    # background level = p[4]
    # model(x, p) = p[4] + (p[3] / ( p[1] * (x-p[2])^2) )
    return p[4] + p[3]*π*p[1]/(π*p[1] + ( (x - p[2]) / p[1])^2)
end

"""
    voigt(x::Float64, p::Vector{Float64})
Pseudo-Voigt profile function used by the TBPD method
"""
function voigt(x::Float64, p::Vector{Float64})
    # pseudo-voigt profile
    gammaG = p[1] / (2.0 * sqrt(log(2.0)))
    gammaL = p[1] / 2.0
    Gamma = (gammaG^5 + 2.69269 * gammaG^4 * gammaL + 2.42843 * gammaG^3 * gammaL^2 + 4.47163 * gammaG^2 * gammaL^3 + 0.07842 * gammaG * gammaL^4 + gammaL^5)^(1/5)
    eta = 1.36603 *(gammaL / Gamma) - 0.47719 * (gammaL / Gamma)^2 + 0.11116 * (gammaL / Gamma)^3

    L(x,Gam,x0) = (Gam / π) / ((x-x0)^2 + Gam^2)
    G(x,Gam,x0) = exp( -( (x-x0)^2) / (2.0 * Gam^2) ) / Gam * sqrt(2π)
   return  p[4] + p[3]  * ( eta * L(x,Gamma,p[2]) + (1 - eta) * G(x,Gamma,p[2]) ) 
end


"""
    tbpd(scan::msJ.MScontainer, shape::Symbol,  R::Real, thres::Real)
Template based beak detection algorithm
"""
#function tbpd(scan::MScontainer, shape::Symbol,  R::Real, thres::Real)   #template based peak detection
function tbpd(scan::MScontainer, model::Function,  ∆mz::Real, thres::Real)   #template based peak detection
    box = num2pnt(scan.mz, scan.mz[1]+0.4) - 1        # taking a box of 0.5 width m/z
    correlation = zeros(length(scan.mz))
    maxi = maximum(scan.int)
    val = 0.0
    for i = 1:1:length(scan.mz)-box
        level = scan.int[i]
        if level >=  maxi * thres / 100. 
            bkg = 0.0
            p0 = [∆mz, scan.mz[i], level, bkg]
            ydata = [model(el, p0) for el in scan.mz[i:i+box]]
            val = Statistics.cor(scan.int[i:i+box], ydata)
        else
            val = 0.0
        end
        if val >= 0.62
            correlation[i] = val
        else
            correlation[i] = 0.0
        end
    end

    peaks_mz = Vector{Float64}(undef,0)
    peaks_int = Vector{Float64}(undef,0)
    peaks_s = Vector{Float64}(undef,0)

    diff_prev = 0.0
    diff      = 0.0 

    # rough numerical differentiation of correlation vector to find its maximum
    for i =2:length(correlation)-2
        diff = (-correlation[i-1] +correlation[i+1]) / 2.0 
        if diff < 0.0
            if diff_prev > 0.0
                max_value = maximum( scan.int[i-3:i+3] )
                max_index = num2pnt(scan.int, max_value)
                push!(peaks_mz, scan.mz[max_index])
                push!(peaks_int, scan.int[max_index])
                if scan isa MSscans
                    push!(peaks_s, scan.s[max_index])
                end
            end
        end       
        diff_prev = diff
    end

    basePeakIntensity = maximum(peaks_int)
    basePeakMz = peaks_mz[ num2pnt(peaks_int, basePeakIntensity) ]
    
    if scan isa MSscans
        return MSscans(scan.num, scan.rt, sum(peaks_int), peaks_mz, peaks_int, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, peaks_s)
    elseif scan isa MSscan
        return MSscan(scan.num, scan.rt, sum(peaks_int), peaks_mz, peaks_int, scan.level, basePeakMz, basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy)
    end
end

"""
function snra(scan::MScontainer)                                        # Signal to Noise Ration Analysis
    error("SNR not implemented")
end

function cwt(scan::MScontainer)                                         # Continuous Wavelet Transform 
    error("Wavelets not implemented")
end
"""

