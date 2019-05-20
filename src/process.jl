  
# Mass spectra

"""
    msprocess(scan::MScontainer; resolution::Symbol = :medium, R::Real = 4500., shape::Symbol = :gauss, threshold::Real = 0.2)
Checks the file extension and calls the right function to load the mass spectra if it exists. Returns an array of individual mass spectra 
# Examples
```julia-repl
julia> reduced_data = msJ.find_peaks(filename)
msJ.MSscans
msJ.MSscans(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
```
"""
#function msprocess(scan::MScontainer; resolution::Symbol = :medium, R::Real = 4500., shape::Symbol = :gauss, threshold::Real = 0.2 )
function test(scan::MScontainer, method::MethodType...)
    for el in method
        if el isa PP
            println(el)
            if el.method == nothing
                println("default tbpd")
            end
            
        end
    end
end

function smooth(scan::MScontainer; method::MethodType=SG(7, 15))
    if method isa msJ.SG
        println(method.order, "  ", method.window)
        return savitzky_golay_filtering(scan, method.order, method.window)
    end
    
    
end

function savitzky_golay_filtering(scan::MScontainer, order::Int, window::Int, deriv::Int=0)
    println(order, " ",window )
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
        
    pad = vcat(scan.int[1]*ones(half_window), scan.int, scan.int[end]*ones(half_window))
    y = conv(coefs[end:-1:1], pad)[2 * half_window + 1 : end - 2 * half_window]
    
    if scan isa msJ.MSscan
        return MSscan(scan.num, scan.rt, scan.tic, scan.mz, y, scan.level, scan.basePeakMz, scan.basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy)
    elseif scan isa MSscans
        return msJ.MSscans(scan.num, scan.rt, scan.tic, scan.mz, y, scan.level, scan.basePeakMz, scan.basePeakIntensity, scan.precursor, scan.polarity, scan.activationMethod, scan.collisionEnergy, scan.s)
    end
end
 

function centroid(scan::MScontainer; resolution::Symbol = :medium, R::Real = 4500., shape::Symbol = :gauss, threshold::Real = 0.2 )
    if R == 4500.
        if resolution == :low
            R = 1000.
        elseif resolution == :medium
            R = 4500.
        elseif resolution == :high
            R = 50000.
        elseif resolution == :veryhigh
            R = 300000.
        else
            error("Unsupported resolution")
        end
    else
        resolution = :custom
    end
    return tbpd(scan, R, shape, threshold)
end


function snra(scan::MScontainer)                                        # Signal to Noise Ration Analysis
    error("SNR not implemented")
end

function cwt(scan::MScontainer)                                         # Continuous Wavelet Transform 
    error("Wavelets not implemented")
end

function tbpd(scan::MScontainer, R::Real, shape::Symbol, thres::Real)   #template based peak detection
    if shape == :gauss
        # Gaussian shape function
        # width            = p[1]
        # x0               = p[2]
        # height           = p[3]
        # background level = p[4]
        @. model(x, p) = p[4] + p[3] * exp(- ( (x-p[2])/p[1] )^2)
        ∆mz = 500.0 / R                  # according to ∆mz / mz  = R, we take the value @ m/z 500
    elseif shape == :lorentz
        # Lorentzian shape function
        # width            = p[1]
        # x0               = p[2]
        # height           = p[3]
        # background level = p[4]
        #@. model(x, p) = p[4] + (p[3] / ( p[1] * (x-p[2])^2) )
    elseif shape == :voigt
        # Voight shape function
    else
        error("Unsupported function")
    end
    box = num2pnt(scan.mz, scan.mz[1]+0.4) - 1        # taking a box of 0.5 width m/z
    correlation = zeros(length(scan.mz))
    maxi = maximum(scan.int)
    val = 0.0
    for i = 1:1:length(scan.mz)-box
        level = scan.int[i]
        if level >=  maxi * thres / 100. 
            bkg = 0.0
            p0 = [∆mz, scan.mz[i], level, bkg]
            ydata = model(scan.mz[i:i+box], p0)
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

    # numerical differentiation of correlation vector to find its maximum
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
