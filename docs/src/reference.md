References
==========

```@meta
CurrentModule = msJ
DocTestSetup  = quote
    using LightXML
end
```

This page lists all the documented elements of the `msJ.jl` package covering all modules and submodules.

## Contents

```@contents
Pages = ["reference.md"]
```

# Main module
```@docs
msJ
```

## Types
--------
Submodule with types and structures used to stored the data and dispatch to the right methods.

### Data types
```@docs
msJ.MScontainer
msJ.MSscan
msJ.MSscans
msJ.Chromatogram
```

### Methods Types
```@docs
msJ.MethodType
msJ.BasePeak
msJ.TIC
msJ.∆MZ
msJ.MZ
msJ.SG
msJ.TBPD
```

### Filters
```@docs
msJ.FilterType
msJ.RT
msJ.IC
msJ.Level
msJ.Scan
msJ.Polarity
msJ.Activation_Method
msJ.Activation_Energy
msJ.Precursor
```


## I/O
------

Module for importing and exporting data. Dispatch to specific methods according to the file extension

```@docs
msJ.info(filename::String; verbose::Bool = false)
msJ.load(filename::String)
msJ.retention_time(filename::String)
msJ.chromatogram(filename::String, filters::FilterType...; method::MethodType=TIC())
msJ.msfilter(filename::String, arguments::FilterType...; stats::Bool=true)
```

### mzXML
---------
Interface to the mzxml file format

```@docs
msJ.info_mzxml
msJ.load_mzxml_all
msJ.load_mzxml
msJ.load_mzxml_spectrum
msJ.retention_time(msRun::XMLElement)
msJ.filter(msRun::XMLElement, argument::Level{<:Int})
msJ.filter(msRun::XMLElement, argument::Level{<:AbstractVector})
msJ.filter(msRun::XMLElement, argument::Scan{<:Int})
msJ.filter(msRun::XMLElement, argument::Scan{<:AbstractVector})
msJ.filter(msRun::XMLElement, argument::Polarity{<:String})
msJ.filter(msRun::XMLElement, argument::Polarity{<:AbstractVector})
msJ.filter(msRun::XMLElement, argument::RT{<:Real})
msJ.filter(msRun::XMLElement, argument::RT{<:AbstractVector})
msJ.filter(msRun::XMLElement, argument::RT{<:AbstractVector{<:AbstractVector} } )
msJ.filter(msRun::XMLElement, argument::IC{<:AbstractVector})
msJ.filter(msRun::XMLElement, argument::Precursor{<:Real})
msJ.filter(msRun::XMLElement, argument::Precursor{<:AbstractVector})
msJ.filter(msRun::XMLElement, argument::Activation_Energy{<:Real})
msJ.filter(msRun::XMLElement, argument::Activation_Energy{<:AbstractVector})
msJ.filter(msRun::XMLElement, argument::Activation_Method{<:String})
msJ.filter(msRun::XMLElement, argument::Activation_Method{<:AbstractVector})
msJ.extracted_chromatogram(filename::String, indices::Vector{Int},method::MethodType)
msJ.composite_spectra(filename::String, indices::Vector{Int}, stats::Bool)
```


## Filtering
```@docs
msJ.msfilter(scans::Vector{MSscan}, arguments::FilterType...; stats::Bool=true)
msJ.chromatogram(scans::Vector{MSscan}, filters::FilterType...; method::MethodType=TIC())
msJ.retention_time(scans::Vector{MSscan})
msJ.filter(scans::Vector{MSscan}, argument::Scan{<:Int})
msJ.filter(scans::Vector{MSscan}, argument::Scan{<:AbstractVector})
msJ.filter(scans::Vector{MSscan}, argument::Level{<:Int})
msJ.filter(scans::Vector{MSscan}, argument::Level{<:AbstractVector})
msJ.filter(scans::Vector{MSscan}, argument::Precursor{<:Real})
msJ.filter(scans::Vector{MSscan}, argument::Precursor{<:AbstractVector})
msJ.filter(scans::Vector{MSscan}, argument::Activation_Energy{<:Real})
msJ.filter(scans::Vector{MSscan}, argument::Activation_Energy{<:AbstractVector})
msJ.filter(scans::Vector{MSscan}, argument::Activation_Method{<:String})
msJ.filter(scans::Vector{MSscan}, argument::Activation_Method{<:AbstractVector})
msJ.filter(scans::Vector{MSscan}, argument::Polarity{<:String})
msJ.filter(scans::Vector{MSscan}, argument::Polarity{<:AbstractVector})
msJ.filter(scans::Vector{MSscan}, argument::RT{<:Real}) 
msJ.filter(scans::Vector{MSscan}, argument::RT{<:AbstractVector})
msJ.filter(scans::Vector{MSscan}, argument::RT{<:AbstractVector{<:AbstractVector} } )
msJ.filter(scans::Vector{MSscan}, argument::IC{<:AbstractVector})
msJ.extracted_chromatogram(scans::Vector{MSscan}, indices::Vector{Int},method::MethodType)
msJ.composite_spectra(scans::Vector{MSscan}, indices::Vector{Int}, stats::Bool)
```


## Process
----------

### Mass spectrum
-----------------
```@docs
msJ.smooth(scan::MScontainer; method::MethodType=SG(5, 9, 0))
msJ.savitzky_golay_filtering(scan::MScontainer, order::Int, window::Int, deriv::Int)
msJ.centroid(scan::MScontainer; method::MethodType=TBPD(:gauss, 4500., 0.2) )
msJ.tbpd(scan::MScontainer, shape::Symbol,  R::Real, thres::Real)
```


### Chromatogram
----------------
No functions yet. To be added.



## Plots
--------
```@autodocs
Modules = [msJ.plots]

```


## Utilities
------------

### Base overloaded
```@docs
+(a::MScontainer, b::MScontainer)
-(a::MScontainer, b::MScontainer)
/(a::MSscan, N::Real)
/(a::MSscans, N::Real)
*(a::MSscan, N::Real)
*(a::MSscans, N::Real)
*(N::Real, a::MScontainer)
*(a::MScontainer, b::MScontainer)
```

### Utility function
```@docs
msJ.avg(a::MScontainer, b::MScontainer)
msJ.add_ion_current(x::Vector{Float64}, y::Vector{Float64}, a::Float64, b::Float64)
msJ.num2pnt(x::Vector{Float64}, val::Real)
```
