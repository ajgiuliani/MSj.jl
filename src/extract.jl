"""
Module for extracting subsets from a Vector{MSscan} according to specific conditions
"""

# User Interface.
# ---------------

export extract



"""
    extract(filename::String, arguments::FilterType...)
Search for scans matching the argument MS level and returns an array of matching MSscans otherwise returns an ErrorException: "No matching spectra found."
# Examples
```julia-repl
julia> sub_set = extract("test.mzxml")
6-element Array{MSj.MSscan,1}:
 MSj.MSscan(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
julia> sub_set = extract("test.mzxml", MSj.Level(2) )      # extract MS/MS spectra
MSj.MSscan(2, 0.7307, 9727.2, [345.083, 345.167, 345.25, 345.333, 345.417, 345.5, 345.583, 345.667, 345.75, 345.833  …  1999.25, 1999.33, 1999.42, 1999.5, 1999.58 ....
 MSj.MSscan(5, 4.3442, 12203.5, [345.083, 345.167, 345.25, 345.333, 345.417, 345.5, 345.583, 345.667, 345.75, 345.833  …  1999.25, 1999.33, 1999.42, 1999.5, 1999.58, ....
```
"""
function extract(filename::String, arguments::FilterType...)
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
    if length(indices) >= 1
        return build_subset(filename, indices)
    else
        ErrorException("No matching spectra.")
    end
end



"""
    build_subset(filename::String, indices::Vector{Int})
Returns a Vector of MSscan from the input file according to the scan num (indices).
"""
function build_subset(filename::String, indices::Vector{Int})
    sub_set = Vector{MSscan}(undef,0)   
    for i = 1:length(indices)
        push!(sub_set, load_mzxml(filename, indices[i]))
    end
    return sub_set
end


"""
    extract(scans::Vector{MSscan}, arguments::FilterType...)
Search for scans matching the argument MS level and returns an array of matching MSscans otherwise returns an ErrorException: "No matching spectra found."
# Examples
```julia-repl
julia> scans = load("test.mzxml")                          # load mass spectra
6-element Array{MSj.MSscan,1}:
 MSj.MSscan(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ...
julia> sub_set = extract(scans)                            # extract a sub_set without conditions returns the original data
6-element Array{MSj.MSscan,1}:
 MSj.MSscan(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....
julia> sub_set = extract(scans, MSj.Level(2) )      # extract MS/MS spectra
MSj.MSscan(2, 0.7307, 9727.2, [345.083, 345.167, 345.25, 345.333, 345.417, 345.5, 345.583, 345.667, 345.75, 345.833  …  1999.25, 1999.33, 1999.42, 1999.5, 1999.58 ....
 MSj.MSscan(5, 4.3442, 12203.5, [345.083, 345.167, 345.25, 345.333, 345.417, 345.5, 345.583, 345.667, 345.75, 345.833  …  1999.25, 1999.33, 1999.42, 1999.5, 1999.58, ....
```

"""
function extract(scans::Vector{MSscan}, arguments::FilterType...)
    index = Set( i for i in 1:length(scans) )
    
    for el in arguments
        subindex = filter(scans, el)
        index = intersect(index, subindex)
    end

    indices = sort([ i for i in index])
    
    if length(indices) >= 1
        return build_subset(scans, indices)
    else
        ErrorException("No matching spectra found.")
    end

end

"""
    build_subset(scans::Vector{MSscan}, indices::Vector{Int})
Returns a Vector of MSscan according to the input scan num.
"""
function build_subset(scans::Vector{MSscan}, indices::Vector{Int})
    sub_set = Vector{MSscan}(undef,0)
    for i = 1:length(indices)
        push!(sub_set, scans[indices[i]])
    end
    return sub_set
end


