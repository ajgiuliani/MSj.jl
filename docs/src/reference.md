References
==========
msJ.builoto
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
msJ.SNRA
msJ.TopHat
msJ.LOESS
msJ.IPSA
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
msJ.Isotopes
```

## I/O
------

Module for importing and exporting data. Dispatch to specific methods according to the file extension

```@docs
msJ.info(filename::String; verbose::Bool = false)
msJ.load(filename::String)
msJ.retention_time(filename::String)
msJ.chromatogram(filename::String, filters::FilterType...; method::MethodType=TIC())
msJ.average(filename::String, arguments::FilterType...; stats::Bool=true)
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
msJ.average(msRun::XMLElement, argument::Level{<:Int})
msJ.average(msRun::XMLElement, argument::Level{<:AbstractVector})
msJ.average(msRun::XMLElement, argument::Scan{<:Int})
msJ.average(msRun::XMLElement, argument::Scan{<:AbstractVector})
msJ.average(msRun::XMLElement, argument::Polarity{<:String})
msJ.average(msRun::XMLElement, argument::Polarity{<:AbstractVector})
msJ.average(msRun::XMLElement, argument::RT{<:Real})
msJ.average(msRun::XMLElement, argument::RT{<:AbstractVector})
msJ.average(msRun::XMLElement, argument::RT{<:AbstractVector{<:AbstractVector} } )
msJ.average(msRun::XMLElement, argument::IC{<:AbstractVector})
msJ.average(msRun::XMLElement, argument::Precursor{<:Real})
msJ.average(msRun::XMLElement, argument::Precursor{<:AbstractVector})
msJ.average(msRun::XMLElement, argument::Activation_Energy{<:Real})
msJ.average(msRun::XMLElement, argument::Activation_Energy{<:AbstractVector})
msJ.average(msRun::XMLElement, argument::Activation_Method{<:String})
msJ.average(msRun::XMLElement, argument::Activation_Method{<:AbstractVector})
msJ.extracted_chromatogram(filename::String, indices::Vector{Int},method::MethodType)
msJ.composite_spectra(filename::String, indices::Vector{Int}, stats::Bool)
```


## Filtering
```@docs
msJ.average(scans::Vector{MSscan}, arguments::FilterType...; stats::Bool=true)
msJ.chromatogram(scans::Vector{MSscan}, filters::FilterType...; method::MethodType=TIC())
msJ.retention_time(scans::Vector{MSscan})
msJ.average(scans::Vector{MSscan}, argument::Scan{<:Int})
msJ.average(scans::Vector{MSscan}, argument::Scan{<:AbstractVector})
msJ.average(scans::Vector{MSscan}, argument::Level{<:Int})
msJ.average(scans::Vector{MSscan}, argument::Level{<:AbstractVector})
msJ.average(scans::Vector{MSscan}, argument::Precursor{<:Real})
msJ.average(scans::Vector{MSscan}, argument::Precursor{<:AbstractVector})
msJ.average(scans::Vector{MSscan}, argument::Activation_Energy{<:Real})
msJ.average(scans::Vector{MSscan}, argument::Activation_Energy{<:AbstractVector})
msJ.average(scans::Vector{MSscan}, argument::Activation_Method{<:String})
msJ.average(scans::Vector{MSscan}, argument::Activation_Method{<:AbstractVector})
msJ.average(scans::Vector{MSscan}, argument::Polarity{<:String})
msJ.average(scans::Vector{MSscan}, argument::Polarity{<:AbstractVector})
msJ.average(scans::Vector{MSscan}, argument::RT{<:Real}) 
msJ.average(scans::Vector{MSscan}, argument::RT{<:AbstractVector})
msJ.average(scans::Vector{MSscan}, argument::RT{<:AbstractVector{<:AbstractVector} } )
msJ.average(scans::Vector{MSscan}, argument::IC{<:AbstractVector})
msJ.extracted_chromatogram(scans::Vector{MSscan}, indices::Vector{Int},method::MethodType)
msJ.composite_spectra(scans::Vector{MSscan}, indices::Vector{Int}, stats::Bool)
```

## Extracting subsets
```@docs
msJ.extract(filename::String, arguments::FilterType...)
msJ.extract(scans::Vector{MSscan}, arguments::FilterType...)
msJ.build_subset(filename::String, indices::Vector{Int})
msJ.build_subset(scans::Vector{MSscan}, indices::Vector{Int})
```

## Process
----------

### Mass spectrum
-----------------
```@docs
msJ.smooth(scan::MScontainer; method::MethodType=SG(5, 9, 0))
msJ.smooth(scans::Vector{MSscan}; method::MethodType=SG(5, 9, 0))
msJ.savitzky_golay_filtering(scan::MScontainer, order::Int, window::Int, deriv::Int)
msJ.smooth(scans::Vector{MSscan}; method::MethodType=SG(5, 9, 0))
msJ.centroid(scan::MScontainer; method::MethodType=SNRA(1., 100) )
msJ.centroid(scans::Vector{MSscan}; method::MethodType=SNRA(1., 100) 
msJ.snra(scan::MScontainer, thres::Real, region::Int)
msJ.tbpd(scan::MScontainer, model::Function,  ∆mz::Real, thres::Real)
msJ.gauss(x::Float64, p::Vector{Float64})
msJ.lorentz(x::Float64, p::Vector{Float64})
msJ.voigt(x::Float64, p::Vector{Float64})
msJ.tbpd(scan::MScontainer, shape::Symbol,  R::Real, thres::Real)
msJ.baseline_correction(scan::MScontainer; method::MethodType=TopHat(100) )
msJ.baseline_correction(scans::Vector{MSscan}; method::MethodType=TopHat(100) )
msJ.tophat_filter(scan::MScontainer, region::Int )
msJ.tophat_filter(scans::Vector{MSscan}, region::Int )
msJ.loess(scans::Vector{MSscan}, iter::Int )
msJ.loess(scan::MScontainer, iter::Int )
msJ.ipsa(scan::MScontainer, width::Real, maxiter::Int)
msJ.ipsa(scans::Vector{MSscan}, width::Real, maxiter::Int)
```


### Chromatogram
----------------
No functions yet. To be added.


## Simulation
-------------
```@docs
msJ.formula
msJ.masses
msJ.isotopic_distribution
msJ.simulate
```


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
msJ.savitzky_golay(int::AbstractArray, order::Int, window::Int, deriv::Int)
msJ.extremefilt(input::AbstractArray, minmax::Function, region::Int)
msJ.morpholaplace(input::AbstractArray, region::Int)
msJ.morphogradient(input::AbstractArray, region::Int)
msJ.tophat(input::AbstractArray, region::Int)
msJ.bottomhat(input::AbstractArray, region::Int) 
msJ.opening(input::AbstractArray, region::Int)
msJ.closing(input::AbstractArray, region::Int)
msJ.erosion(input::AbstractArray, region::Int)
msJ.dilatation(input::AbstractArray, region::Int)
msJ.convolve(a::AbstractArray, b::AbstractArray)
```
