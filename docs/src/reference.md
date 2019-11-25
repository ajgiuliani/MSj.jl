References
==========
MSJ.builoto
```@meta
CurrentModule = MSJ
DocTestSetup  = quote
    using LightXML
end
```

This page lists all the documented elements of the `MSJ.jl` package covering all modules and submodules.

## Contents

```@contents
Pages = ["reference.md"]
```

# Main module
```@docs
MSJ
```

## Types
--------
Submodule with types and structures used to stored the data and dispatch to the right methods.

### Data types
```@docs
MSJ.MScontainer
MSJ.MSscan
MSJ.MSscans
MSJ.Chromatogram
```

### Methods Types
```@docs
MSJ.MethodType
MSJ.BasePeak
MSJ.TIC
MSJ.∆MZ
MSJ.MZ
MSJ.SG
MSJ.TBPD
MSJ.SNRA
MSJ.TopHat
MSJ.LOESS
MSJ.IPSA
```

### Filters
```@docs
MSJ.FilterType
MSJ.RT
MSJ.IC
MSJ.Level
MSJ.Scan
MSJ.Polarity
MSJ.Activation_Method
MSJ.Activation_Energy
MSJ.Precursor
MSJ.Isotopes
```

## I/O
------

Module for importing and exporting data. Dispatch to specific methods according to the file extension

```@docs
MSJ.info(filename::String; verbose::Bool = false)
MSJ.load(filename::String)
MSJ.retention_time(filename::String)
MSJ.chromatogram(filename::String, filters::FilterType...; method::MethodType=TIC())
MSJ.average(filename::String, arguments::FilterType...; stats::Bool=true)
```

### mzXML
---------
Interface to the mzxml file format

```@docs
MSJ.info_mzxml
MSJ.load_mzxml_all
MSJ.load_mzxml
MSJ.load_mzxml_spectrum
MSJ.retention_time(msRun::XMLElement)
MSJ.average(msRun::XMLElement, argument::Level{<:Int})
MSJ.average(msRun::XMLElement, argument::Level{<:AbstractVector})
MSJ.average(msRun::XMLElement, argument::Scan{<:Int})
MSJ.average(msRun::XMLElement, argument::Scan{<:AbstractVector})
MSJ.average(msRun::XMLElement, argument::Polarity{<:String})
MSJ.average(msRun::XMLElement, argument::Polarity{<:AbstractVector})
MSJ.average(msRun::XMLElement, argument::RT{<:Real})
MSJ.average(msRun::XMLElement, argument::RT{<:AbstractVector})
MSJ.average(msRun::XMLElement, argument::RT{<:AbstractVector{<:AbstractVector} } )
MSJ.average(msRun::XMLElement, argument::IC{<:AbstractVector})
MSJ.average(msRun::XMLElement, argument::Precursor{<:Real})
MSJ.average(msRun::XMLElement, argument::Precursor{<:AbstractVector})
MSJ.average(msRun::XMLElement, argument::Activation_Energy{<:Real})
MSJ.average(msRun::XMLElement, argument::Activation_Energy{<:AbstractVector})
MSJ.average(msRun::XMLElement, argument::Activation_Method{<:String})
MSJ.average(msRun::XMLElement, argument::Activation_Method{<:AbstractVector})
MSJ.extracted_chromatogram(filename::String, indices::Vector{Int},method::MethodType)
MSJ.composite_spectra(filename::String, indices::Vector{Int}, stats::Bool)
```


## Filtering
```@docs
MSJ.average(scans::Vector{MSscan}, arguments::FilterType...; stats::Bool=true)
MSJ.chromatogram(scans::Vector{MSscan}, filters::FilterType...; method::MethodType=TIC())
MSJ.retention_time(scans::Vector{MSscan})
MSJ.average(scans::Vector{MSscan}, argument::Scan{<:Int})
MSJ.average(scans::Vector{MSscan}, argument::Scan{<:AbstractVector})
MSJ.average(scans::Vector{MSscan}, argument::Level{<:Int})
MSJ.average(scans::Vector{MSscan}, argument::Level{<:AbstractVector})
MSJ.average(scans::Vector{MSscan}, argument::Precursor{<:Real})
MSJ.average(scans::Vector{MSscan}, argument::Precursor{<:AbstractVector})
MSJ.average(scans::Vector{MSscan}, argument::Activation_Energy{<:Real})
MSJ.average(scans::Vector{MSscan}, argument::Activation_Energy{<:AbstractVector})
MSJ.average(scans::Vector{MSscan}, argument::Activation_Method{<:String})
MSJ.average(scans::Vector{MSscan}, argument::Activation_Method{<:AbstractVector})
MSJ.average(scans::Vector{MSscan}, argument::Polarity{<:String})
MSJ.average(scans::Vector{MSscan}, argument::Polarity{<:AbstractVector})
MSJ.average(scans::Vector{MSscan}, argument::RT{<:Real}) 
MSJ.average(scans::Vector{MSscan}, argument::RT{<:AbstractVector})
MSJ.average(scans::Vector{MSscan}, argument::RT{<:AbstractVector{<:AbstractVector} } )
MSJ.average(scans::Vector{MSscan}, argument::IC{<:AbstractVector})
MSJ.extracted_chromatogram(scans::Vector{MSscan}, indices::Vector{Int},method::MethodType)
MSJ.composite_spectra(scans::Vector{MSscan}, indices::Vector{Int}, stats::Bool)
```

## Extracting subsets
```@docs
MSJ.extract(filename::String, arguments::FilterType...)
MSJ.extract(scans::Vector{MSscan}, arguments::FilterType...)
MSJ.build_subset(filename::String, indices::Vector{Int})
MSJ.build_subset(scans::Vector{MSscan}, indices::Vector{Int})
```

## Process
----------

### Mass spectrum
-----------------
```@docs
MSJ.smooth(scan::MScontainer; method::MethodType=SG(5, 9, 0))
MSJ.smooth(scans::Vector{MSscan}; method::MethodType=SG(5, 9, 0))
MSJ.savitzky_golay_filtering(scan::MScontainer, order::Int, window::Int, deriv::Int)
MSJ.smooth(scans::Vector{MSscan}; method::MethodType=SG(5, 9, 0))
MSJ.centroid(scan::MScontainer; method::MethodType=SNRA(1., 100) )
MSJ.centroid(scans::Vector{MSscan}; method::MethodType=SNRA(1., 100) 
MSJ.snra(scan::MScontainer, thres::Real, region::Int)
MSJ.tbpd(scan::MScontainer, model::Function,  ∆mz::Real, thres::Real)
MSJ.gauss(x::Float64, p::Vector{Float64})
MSJ.lorentz(x::Float64, p::Vector{Float64})
MSJ.voigt(x::Float64, p::Vector{Float64})
MSJ.tbpd(scan::MScontainer, shape::Symbol,  R::Real, thres::Real)
MSJ.baseline_correction(scan::MScontainer; method::MethodType=TopHat(100) )
MSJ.baseline_correction(scans::Vector{MSscan}; method::MethodType=TopHat(100) )
MSJ.tophat_filter(scan::MScontainer, region::Int )
MSJ.tophat_filter(scans::Vector{MSscan}, region::Int )
MSJ.loess(scans::Vector{MSscan}, iter::Int )
MSJ.loess(scan::MScontainer, iter::Int )
MSJ.ipsa(scan::MScontainer, width::Real, maxiter::Int)
MSJ.ipsa(scans::Vector{MSscan}, width::Real, maxiter::Int)
```


### Chromatogram
----------------
No functions yet. To be added.


## Simulation
-------------
```@docs
MSJ.formula
MSJ.masses
MSJ.isotopic_distribution
MSJ.simulate
```


## Plots
--------
```@autodocs
Modules = [MSJ.plots]

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
MSJ.avg(a::MScontainer, b::MScontainer)
MSJ.add_ion_current(x::Vector{Float64}, y::Vector{Float64}, a::Float64, b::Float64)
MSJ.num2pnt(x::Vector{Float64}, val::Real)
MSJ.savitzky_golay(int::AbstractArray, order::Int, window::Int, deriv::Int)
MSJ.extremefilt(input::AbstractArray, minmax::Function, region::Int)
MSJ.morpholaplace(input::AbstractArray, region::Int)
MSJ.morphogradient(input::AbstractArray, region::Int)
MSJ.tophat(input::AbstractArray, region::Int)
MSJ.bottomhat(input::AbstractArray, region::Int) 
MSJ.opening(input::AbstractArray, region::Int)
MSJ.closing(input::AbstractArray, region::Int)
MSJ.erosion(input::AbstractArray, region::Int)
MSJ.dilatation(input::AbstractArray, region::Int)
MSJ.convolve(a::AbstractArray, b::AbstractArray)
```
