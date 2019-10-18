Manual
======

## Introduction
---------------
The `msJ` package aims at providing an API to the most common file format in mass spectrometry. The following file formats are currently supported:
- mzXML 

## Public elements
-----------------
The functions below are exported:
- [`info`](@ref) 
- [`load`](@ref) 
- [`chromatogram`](@ref)
- [`msfilter`](@ref)
- [`extract`](@ref)
- [`centroid`](@ref)
- [`smooth`](@ref)
- [`baseline_correction`](@ref)


## Data types
-------------
The main data type of the package is the abstract type [`msJ.MScontainer`](@ref).

Mass spectrometry scans are stored in the following structures, inspired from the mzXML format, which is a subtype of [`msJ.MScontainer`](@ref). 
```julia
struct MSscan <: MScontainer
    num::Int                          # num
    rt::Float64                       # retentionTime
    tic::Float64                      # totIonCurrent
    mz::Vector{Float64}               # m/z
    int::Vector{Float64}              # intensity
    level::Int                        # msLevel
    basePeakMz::Float64               # basePeakMz
    basePeakIntensity::Float64        # basePeakIntensity
    precursor::Float64                # precursorMz
    polarity::String                  # polarity
    activationMethod::String          # activationMethod
    collisionEnergy::Float64          # collisionEnergy
end
```

Another subtype, [`msJ.Chromatogram`](@ref), is used to store the retention time, the ionic current and the maximum value of the ion current.

```julia
struct Chromatogram  <: MScontainer
    rt::Vector{Float64}               # araay of retention times
    ic::Vector{Float64}               # array of ion current
    maxic::Float64                    # maximum ion current (used in plotting normalization)
end
```

Combination of mass spectra requires another subtype of [`msJ.MScontainer`](@ref) called [`msJ.MSscans`](@ref) (notice the ending s).

```julia
struct MSscans  <: MScontainer
    num::Vector{Int}                  # num
    rt::Vector{Float64}               # retentionTime
    tic::Float64                      # totIonCurrent
    mz::Vector{Float64}               # m/z
    int::Vector{Float64}              # intensity
    level::Vector{Int}                # msLevel
    basePeakMz::Float64               # basePeakMz
    basePeakIntensity::Float64        # basePeakIntensity
    precursor::Vector{Float64}        # precursorMz
    polarity::Vector{String}          # polarity
    activationMethod::Vector{String}  # activationMethod
    collisionEnergy::Vector{Float64}  # collisionEnergy
    s::Vector{Float64}                # variance
end
```
The [`msJ.MSscans`](@ref) structure is very similar to the [`msJ.MSscan`](@ref) one, except that the fields `num`, `rt`, `precursor`, `polarity`, `activationMethod` and `collisionEnergy` are vectors. The idea is to keep track of the *history* of the operations that have led to this result. For example, if a `MSscans` element is the result of the addition of two individual scans such as scans[1] + scans[2], then the `num`field of resulting `MSscans` is an array [1, 2]. The same applies to the retention time, precursor m/z, polarity, activation method and collision energy fields.

## Information
--------------

The [`info`](@ref) public function reads the content of a file, but without loading the mass spectrometry data, and returns a `Vector{String}`containing the number of scans, scans level and for MS/MS data, the precursor m/z, the activation method and energy. Additional information may be gained by setting `verbose = true`.
```julia-repl
info(filename)
4-element Array{String,1}:
 "51 scans"
 "MS1+"
 "MS2+ 1255.5  CID(CE=18)"
 "MS3+ 902.33  PQD(CE=35)"
```

## Importing data
-----------------

When loading a file containing more than a single acquisition, the individual mass spectrometry scans are pushed into an array of [`msJ.MSscan`](@ref).  The individual scans may be retrieve from the array the usual way:

```julia-repl
julia> scans = load("filename")
51-element Array{msJ.MSscan,1}:
 msJ.MSscan(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....)
...

julia> scans[1]
msJ.MSscan(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....)
```

As mentioned above, chromatograms may be retrieved from a file and imported in [`msJ.Chromatogram`](@ref) :
```julia-repl
julia> chromatogram("filename")
msJ.Chromatogram([0.1384, 0.7307, 2.1379, 3.7578, 4.3442, 5.7689], [5.08195e6, 9727.2, 11.3032, 4.8084e6, 12203.5, 4.84455], 5.08195e6)
```

The function [`msJ.retention_time`](@ref) reads the retention time of an input file and returns a `Vector{Float64}`containing the time in seconds.
```julia-repl
julia> msJ.retention_time("filename")
51-element Array{Float64,1}:
  0.1384
  0.7307
  2.1379
....msJ.FilterType

```


## Exporting data
-----------------

!!! note

    This feature is currently under development. 
 

## Combining and filtering data
---------------------------------------

### Average
-------------

The [`msfilter`](@ref) returns the average of the mass spectra directly from a `Vector{MSscan}` after [Importing data](@ref) data or directly from the `filename`.

```julia-repl
julia> msfilter("filename")
msJ.MSscans([1, 2, 3 ....

julia> scans = load("filename")
51-element Array{msJ.MSscan,1}:
 msJ.MSscan(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....)
...

julia> msfilter(scans)
msJ.MSscans([1, 2, 3 ....

```

Operating on files takes more time than working on `Vector{MSscan}` but may be useful to reduce the memory load.

Without any argument the [`msfilter`](@ref) function averages the entire content of the data and the [`chromatogram`](@ref) function operates on also on the entire data.


### Filtering
------------

The [`msfilter`](@ref) and [`chromatogram`](@ref) functions may takes arguments to select specific fields of interest within the data and operate on them. The argument belongs to the [`msJ.FilterType`](@ref). Their properties are listed below:

| FilterType            | Description       | Arguments                                | Specificity            |
|-----------------------|-------------------|------------------------------------------|------------------------|
| msJ.Scan              | Scan num          | Int, Vector{Int}                         | msfilter, chromatogram |
| msJ.Level             | MS level          | Int, Vector{Int}                         | msfilter, chromatogram |
| msJ.Polarity          | Polarity          | String, Vector{String}                   | msfilter, chromatogram |
| msJ.Activation_Method | Activation method | String, Vector{String}                   | msfilter, chromatogram |
| msJ.Activation_Energy | Activation energy | Real, Vector{Real}                       | msfilter, chromatogram |
| msJ.Precursor         | Precursor _m/z_   | Real, Vector{Real}                       | msfilter, chromatogram |
| msJ.RT                | Retention time    | Real, Vector{Real}, Vector{Vector{Real}} | msfilter               |
| msJ.IC                | Ion current       | Vector{Real}                             | msfilter               |



!!! note

    The filtering function goes first through the argument and setup an array of scan num that matches the conditions. Then it uses this array to calculate the average mass spectrum.  So this procedure needs two passes through the data, which is not very efficient. This is a point to make better in the future.
 

When the argument is restricted to a single value, such as `msJ.Scan(1)`, filtering is performed on that specific value. If the argument is a vector then filtering involves all the values within the range.  Filtering on `msJ.scan([1,10])` means that the result will be obtained for scans ranging from 1 to 10.  The same applies for all `FilterType`with the exception of `msJ.∆MZ`, for which the first value of the vector represents the *mz* and the second value represents the spread ∆mz, so that filtering is operated for all *mz* value in the range [m/z - ∆mz , m/z + ∆mz].  The `msJ.RT`type may take a vector or vectors as argument, such `msJ.RT([ [1,10], [20, 30] ]).  In that case, mass spectra will be averaged in [1,10] and [20,30] range.


These filters may be combined together if necessary. For example, the input below returns the average mass spectrum for:
- the MS2 scans (level = 2), 
- precursor m/z 1255.5, 
- upon CID activation conditions
- with an activation energy of 18 
- and for retention times in the 1 to 60 s range.

```julia
msfilter("filename", msJ.Precursor(1255.5),
                     msJ.Activation_Energy(18),
                     msJ.Activation_Method("CID"),
                     msJ.Level(2),
                     msJ.RT( [1, 60] ),
                     )
```

Several filter types may also be combined for `chromatograms`:
```julia
chromatogram("filename", msJ.Precursor(1255.5),
                         msJ.Activation_Energy(18),
                         msJ.Activation_Method("CID"),
                         msJ.Level(2),
                         )
```

If the condition does not match any existing data, then an `ErrorException` is returned with the `"No matching spectra."` message.


The `chromatogram` function has some methods using `msJ.MethodType` arguments:

| MethodType   | Description         | Arguments    | Remark  |
|--------------|---------------------|--------------|---------|
| msJ.TIC      | Total ion current   | None         | Default |
| msJ.BasePeak | Base peak intensity | None         |         |
| msJ.MZ       | *m/z* range         | Vector{Real} |         |
| msJ.∆MZ      | *m/z* ± ∆mz         | Vector{Real} |         |


These types control the way chromatograms are calculated: either using the total ionic current, the base peak intensity or using a *m/z* range.  The `method` argument of the `msJ.chromatogram`function is set to msJ.TIC() by default. This setting may be overruled by setting the method to desired value:

```julia
chromatogram("filename", method = msJ.BasePeak())
chromatogram("filename", method = msJ.MZ( [257, 259] ) ) 
chromatogram("filename", method = msJ.∆MZ( [258, 1] ) ) 
```

### Extracting subsets
----------------------

The [`extract`](@ref) returns a Vector of `MSscan`from either a file of from a Vector{MSscan} following a ['load'](@ref) command, which corresponds to the filter conditions. See the [filtering](Filtering) part above.

```julia
sub_set = extract("filename")                     # extracting without any conditions returns a vector identical to the output 
sub_set = extract("filename", msJ.Level(2) )      # extract MS/MS spectra
scans = load("test.mzxml")                          # load mass spectra
sub_set = extract(scans)                            # extract a sub_set without conditions returns the original data
```


## Processing
------------
### Smooth
----------
The [`smooth`](@ref) function is public and applies on `MSscan`or `MSscans` objects, with an optional `method` argument set to `msJ.SG(5, 9, 0)`.  Smoothing is performed on the `int` field using the [Savinsky and Golay](https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter). The first argument is the order (5 by default), the second is the number of points (default 9)  and the last, is the derivative level (0).
The function returns an `MScontainer` type identical to the input. 

Other smoothing algorithms will be implemented in the future.


### Base line correction
------------------------
Base line correction is performed using the [`baseline_correction`](@ref) function. This function as two methods and operates either on [`MScontainer`](@ref) or on Array of [`MSscan`](@ref) such as obtained after [importing data](Importing data).
```julia
baseline_correction(scans)
baseline_correction(scans, method = msJ.IPSA(51, 100))
```
The `method` argument allows choosing the algorithm. 

#### Top Hat
This filter is based Top Hat transform used in image processing ([wikipedia](https://en.wikipedia.org/wiki/Top-hat_transform), [Sauve et al. (2004)](https://pdfs.semanticscholar.org/c04c/afc9b2670edd1ea38f0f724cadbe2ec321e9.pdf). The region onto which the operation is performed is set using the `region`field of the [`msJ.TopHat`](@ref). This filter removes every structure from the input which are smaller in size than the structuring element. Usually a region of 100 points is enough.This filter is fast and works quite well on large and complex backgrounds.

#### Iterative polynomial smoothing algorithm (IPSA)
The default algorithm is the IPSA for iterative polynomial smoothing algorithm ([Wang et al. (2017)](https://doi.org/10.1177/0003702816670915). This iterative algorithm use a zero ordre Savinsly and Golay smoothing to estimate a background. Then a new input, constructed by taking the minimum of either the original spectrum or the background, is smooth again. The process is repeated until the maximum iteration is reached or when the background does not change much. The termination criteria has been changed from the original paper.

##### Locally weighted error sum of squares regression (LOESS)
The LOESS family of algorithm is based on non-parametric linear [local regression](https://en.wikipedia.org/wiki/Local_regression) where the regression is weighted to reduced the influence of more distant data. We use here the iterative robust estimation procedure where the weights are updated with a bisquare function of the median of the residuals.
This algorithm takes the number of iteration to be performed. Usually 3 iteration is enough. This algorithm is slow and is not recommended. The implementation will be improved in future versions.


### Peak picking
----------------
Pick-picking is performed using the public [`centroid`](@ref) function. It operates on `MSscan`or `MSscans`type of data and return a similar type. It takes a method argument, set by default to the Signal to Noise Analysis method: `msJ.SNRA`.
```julia
centroid(scan)
```

#### Signal to Noise Ratio Analysis (SNRA)
Signal to noise ratio analysis is a very general approach, which relies on the definition of noise. Here, we use TopHat filter to define the noise. Then the signal to noise ratio is calculated. Peaks are found by searching for a local maximum for which the signal to noise ratio is above the threshold. By defaults the `msJ.SNRA` uses a threshold = 1.0 and a structuring element of 10 points.
```julia
centroid(scan, method = SNRA(1., 10)
```

#### Threshold base peak detection algorithm (TBPD)
The TBPD method identifies features based on their similarity (as described by the Pearson correlation coefficient) with a [template peak](https://doi.org/10.1007/978-1-60761-987-1_22). By default the `msJ.TBPD` method type uses a Gaussian function, with 1000 mass resolving power and a threshold level set to 0.2% as :
```julia
centroid(scan, method = TBPD(:gauss, 1000, 0.2)
```
Two other shape functions are available:
- `:loretz` which uses a Cauchy-Lorentz function and
- `:voigt` which implements a pseudo-voigt profile ([Ida et al., J. Appl. Cryst. (2000)](https://doi.org/10.1107%2Fs0021889800010219), [Wikipedia](https://en.wikipedia.org/wiki/Voigt_profile#Pseudo-Voigt_approximation))

The `:lorentz` profile fits better Fourrier Transform mass spectra. The `:voigt` shape is the result of the convolution of gaussi and Cauchy-Lorentz shape.

```julia
centroid(scan, method = TBPD(:lorentz, 1000., 0.1)
centroid(scan, method = TBPD(:voight,  1000., 0.1)
```



## Plotting
-----------

Plotting facilities are available as a submodule to the `msJ` package.  The [`msJ.plots`](@ref) module relies on the [RecipesBase package](https://github.com/JuliaPlots/RecipesBase.jl), which allows writing recipes to plot users' data types. Hence, recipes have been created for `MSscan`, `Msscans` and `Chromatogram`:

```julia-repl
julia> plot(scans[1], method = :relative))

```

By default plotting is made in relative intensities, which may be changed by setting method to `:absolute`.
