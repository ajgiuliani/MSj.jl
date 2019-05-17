"""
     info(filename::String, verbose::Bool = false)
The function looks into an file and returns in an Array{String} containing the number of scans and the different scans described by their MS level, polarity and eventually the precursor m/z followed by the activation method and collision energy. Each entry is unique, which gives a summary of the input file. With verbose = true, the functions also returns the parentFile, msManufacturer, msModel, msIonisation, msMassAnalyzer, msDetector, software and dataProcessing if existing.
# Example
```julia-repl
julia> msJ.info("test1.mzXML")
4-element Array{String,1}:
 "51 scans"               
 "MS1+"                   
 "MS2+ 1255.5  CID(CE=18)"
 "MS3+ 902.33  PQD(CE=35)"
julia> msJ.info("test1.mzXML", verbose = true)
12-element Array{String,1}:
 "parentFile: test1_msJ_1.raw"
 "msManufacturer: Thermo Finnigan"
 "msModel: LTQ XL"      
 "msIonisation: ESI"
 "msMassAnalyzer: ITMS"
 "msDetector: unknown"
 "software: Xcalibur, 2.6.0 SP3"
 "dataProcessing: conversion, ReAdW 4.3.1(build Sep  9 2009 12:30:29)"
 "51 scans"
 "MS1+"
 "MS2+ 1255.5  CID(CE=18)"
 "MS3+ 902.33  PQD(CE=35)"
```
"""
function info(filename::String; verbose::Bool = false)
    info = Vector{String}(undef,0)
    extension = split(filename, ".")[end]

    if Unicode.normalize(extension, casefold=true) == "mzxml"
        return info_mzxml!(filename, info, verbose)

    elseif Unicode.normalize(extension, casefold=true) == "ascii"
        return info_ascii!(filename, info)

    else 
        error("File format not supported.") 
    end 
end


