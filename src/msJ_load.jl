
"""
    load(filename::String)
Checks the file extension and calls the right function to load the mass spectra if it exists. Returns an array of individual mass spectra 
# Examples
```julia-repl
julia> scans = msJ.load("test.mzXML")
51-element Array{msJ.MSscan,1}:
 msJ.MSscan(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
```
"""

function load(filename::String)
    extension = split(filename, ".")[end]

    if Unicode.normalize(extension, casefold=true) == "mzxml"
        #println("loading mzXML file...")
        return load_mzxml!(filename)

#    elseif Unicode.normalize(extension, casefold=true) == "ascii"
#        #println("loading ascii file...")
#        return load_ascii!(filename, scans)

    else
        return ErrorException("File format not supported.")
    end
end



"""
    retention_time(filename::String)
Returns an array composed of the retention times of the individual mass spectra. 
# Examples
```julia-repl
julia> msJ.retention_time("test.mzXML")
51-element Array{Float64,1}:
  0.1384
  0.7307
  2.1379
....
```
"""
function retention_time(filename::String)
    extension = split(filename, ".")[end]
    rt  = Vector{Float64}(undef,0)

    if Unicode.normalize(extension, casefold=true) == "mzxml"
        xdoc = parse_file(filename)
        xroot = root(xdoc)
        if name(xroot) != "mzXML"
            ErrorException("Not an mzXML file.")
        end
        
        msRun = find_element(xroot, "msRun")    
        scanCount = attribute(msRun, "scanCount")
        
        rt = retention_time(msRun)
        free(xdoc)

#    elseif Unicode.normalize(extension, casefold=true) == "ascii"
#        rt = rt_ascii!(filename, rt)
    else
        ErrorException("File format not supported.")
    end
    return rt
end


function retention_time(scans::Vector{MSscan})
    rt  = Vector{Float64}(undef,0)
    for elem in scans
        push!(rt, elem.rt)
    end
    return rt
end



"""
    chromatogram(filename::String, filters::FilterType...; method::MethodType=TIC())
Returns the retention time and the total ion current by default for all the mass spectra within the file. Alternatively, other options may be supplied such as method = msJ.BasePeak, which returs the base peak intensity, method = msJ.∆MZ([500,5]), which returns the ion current for the range mz = 500 ± 5, or method = msJ.MZ([200,1000]) which return the ion current in the range from m/z 200 to m/z 1000.  The data may be filtered by ms level, precursor mass, activation methods, etc, using the arguments msJ.Level(N), msJ.Precursor(mz), msJ.Activation_Method("method")...
# Examples
```julia-repl
julia> rt, ic = msJ.chromatogram("test.mzxml")
([0.1384  …  60.4793], [4.74795e6  …  17.4918])
julia> rt, ic = msJ.chromatogram("test.mzxml", method = msJ.BasePeak() )
([0.1384  …  60.4793], [102558.0  …  1.23181])
julia> rt, ic = msJ.chromatogram("test.mzxml", method = msJ.∆MZ([500,5]) )
([0.1384  …  60.4793], [46036.6  …  14.2529])
julia> rt, ic = msJ.chromatogram("test.mzxml", method = msJ.MZ([200,1000]))
([0.1384  …  60.4793], [4.74795e6  …  17.4918])
```
"""
function chromatogram(filename::String, filters::FilterType...; method::MethodType=TIC())
    xrt = Vector{Float64}(undef,0)
    xic = Vector{Float64}(undef,0)
    index = Set{Int}()
    extension = split(filename, ".")[end]
    if Unicode.normalize(extension, casefold=true) == "mzxml"
        xdoc = parse_file(filename)
        xroot = root(xdoc)
        if name(xroot) != "mzXML"
            ErrorException("Not an mzXML file.")
        end
        
        msRun = find_element(xroot, "msRun")    
        scanCount = parse(Int, attribute(msRun, "scanCount"))

        index = Set( i for i in 1:scanCount )
        
        for el in filters
            subindex = filter(msRun, el)
            index = intersect(index, subindex) 
        end
        
#    elseif Unicode.normalize(extension, casefold=true) == "ascii"
#        xrt, xic = xtic_ascii!(filename, xrt, xtic, arg::FilterType )       
    else
        ErrorException("File format not supported.")
    end
    free(xdoc)

    indices = sort([ i for i in index])
    if length(indices) != 0
        return extracted_chromatogram(filename, indices, method)
    else
        ErrorException("No matching spectra.")
    end

end

"""
    msfilter(filename::String, filters::FilterType...; stats::Bool=true)
Returns the average mass spectrum container (MSscans) along with the sample standard deviation of the intensities with stats=true (default) for all the mass spectra within file. The data may be filtered by level, precursor mass, activation methods, etc, using the arguments msJ.Level(N), msJ.Precursor(mz), msJ.Activation_Method("method"), or any combination of these arguments.
# Examples
```julia-repl
julia> spectrum = msJ.msJ.msfilter("test.mzxml")
msJ.MSscans([1, 2, 3 ....
julia> spectrum = msJ.msJ.msfilter("test.mzxml", msJ.Level(1) )
msJ.MSscans([1, 4, 7, 10,
julia> spectrum = msJ.msJ.msfilter("test.mzxml", msJ.Precursor(1255.5) )
msJ.MSscans([2, 5, 8, 11, ...
julia> spectrum = msJ.msJ.msfilter("test.mzxml", msJ.Activation_Method("PQD") )
msJ.MSscans([3, 6, 9, 12, 15,
julia> spectrum = msJ.msJ.msfilter("test.mzxml", msJ.Activation_Method("PQD"), msJ.Polarity("+"), msJ.RT([10,20]))
msJ.MSscans([9, 12, 15, 18], ...
```
"""
function msfilter(filename::String, arguments::FilterType...; stats::Bool=true)
    index = Set{Int}()
    extension = split(filename, ".")[end]
    
    if Unicode.normalize(extension, casefold=true) == "mzxml"
        # MZ
        xdoc = parse_file(filename)
        xroot = root(xdoc)
        if name(xroot) != "mzXML"
            ErrorException("Not an mzXML file.")
        end        
        msRun = find_element(xroot, "msRun")    
        scanCount = parse(Int, attribute(msRun, "scanCount"))

        index = Set( i for i in 1:scanCount )
        
        for el in arguments
            subindex = filter(msRun, el)
            index = intersect(index, subindex)
        end

#    elseif Unicode.normalize(extension, casefold=true) == "ascii"
#        error("msfilter not supported for Bruker ascii")
    else
        ErrorException("File format not supported.")
    end

    free(xdoc)
    indices = sort([ i for i in index])
    if length(indices) >= 2
        return composite_spectra(filename, indices, stats)
    elseif length(indices) == 1
        return load_mzxml(filename, indices[1])
    else
        ErrorException("No matching spectra.")
    end
end

"""
    chromatogram(scans::Vector{MSscan}, filters::FilterType...; method::MethodType=TIC())
Returns the retention time and the total ion current by default for all the mass spectra within the Array of mass spectrum container MSscan. Alternatively, other options may be supplied such as method = msJ.BasePeak, which returs the base peak intensity, method = msJ.∆MZ([500,5]), which returns the ion current for the range mz = 500 ± 5, or method = msJ.MZ([200,1000]) which return the ion current in the range from m/z 200 to m/z 1000.  The data may be filtered by ms level, precursor mass, activation methods, etc, using the arguments msJ.Level(N), msJ.Precursor(mz), msJ.Activation_Method("method")...
# Examples
```julia-repl
julia> rt, ic = msJ.chromatogram("test.mzxml")
([0.1384  …  60.4793], [4.74795e6  …  17.4918])
julia> rt, ic = msJ.chromatogram("test.mzxml", method = msJ.BasePeak() )
([0.1384  …  60.4793], [102558.0  …  1.23181])
julia> rt, ic = msJ.chromatogram("test.mzxml", method = msJ.∆MZ([500,5]) )
([0.1384  …  60.4793], [46036.6  …  14.2529])
julia> rt, ic = msJ.chromatogram("test.mzxml", method = msJ.MZ([200,1000]))
([0.1384  …  60.4793], [4.74795e6  …  17.4918])
```
"""
function chromatogram(scans::Vector{MSscan}, filters::FilterType...; method::MethodType=TIC())
    # Ranges of mz value used to compute the tic from
    xrt = Vector{Float64}(undef,0)
    xic = Vector{Float64}(undef,0)
    index = Set( i for i in 1:length(scans) )

    for el in filters
        subindex = filter(scans, el)
        index = intersect(index, subindex)
    end

    indices = sort([ i for i in index])
    if length(indices) != 0
        return extracted_chromatogram(scans, indices, method)
    else
        ErrorException("No matching spectra.")
    end    
end


"""
    msfilter(filename::String, filters::FilterType...; stats::Bool=true)
Returns the average mass spectrum container (MSscans) along with the sample standard deviation of the intensities with stats=true (default) for all the mass spectra within the Array of mass spectrum container MSscan.. The data may be filtered by level, precursor mass, activation methods, etc, using the arguments msJ.Level(N), msJ.Precursor(mz), msJ.Activation_Method("method"), or any combination of these arguments.
# Examples
```julia-repl
julia> spectrum = msJ.msJ.msfilter("test.mzxml")
msJ.MSscans([1, 2, 3 ....
julia> spectrum = msJ.msJ.msfilter("test.mzxml", msJ.Level(1) )
msJ.MSscans([1, 4, 7, 10,
julia> spectrum = msJ.msJ.msfilter("test.mzxml", msJ.Precursor(1255.5) )
msJ.MSscans([2, 5, 8, 11, ...
julia> spectrum = msJ.msJ.msfilter("test.mzxml", msJ.Activation_Method("PQD") )
msJ.MSscans([3, 6, 9, 12, 15,
julia> spectrum = msJ.msJ.msfilter("test.mzxml", msJ.Activation_Method("PQD"), msJ.Polarity("+"), msJ.RT([10,20]))
msJ.MSscans([9, 12, 15, 18], ...
```
"""
function msfilter(scans::Vector{MSscan}, arguments::FilterType...; stats::Bool=true)
    index = Set( i for i in 1:length(scans) )
    
    for el in arguments
        subindex = filter(scans, el)
        index = intersect(index, subindex)
    end

    indices = sort([ i for i in index])
    if length(indices) >= 2
        return composite_spectra(scans, indices, stats)
    elseif length(indices) == 1
        return scans[indices[1] ]
    else
        ErrorException("No matching spectra.")
    end
end

