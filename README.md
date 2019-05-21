# msJ.jl

[![Build Status](https://travis-ci.org/ajgiuliani/msJ.jl.svg?branch=master)](https://travis-ci.org/ajgiuliani/msJ.jl)
[![Coverage Status](https://coveralls.io/repos/github/ajgiuliani/msJ.jl/badge.svg?branch=master)](https://coveralls.io/github/ajgiuliani/msJ.jl?branch=master)
[![codecov](https://codecov.io/gh/ajgiuliani/msJ.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ajgiuliani/msJ.jl)



This package provides an API to the mzXML file format for mass spectrometry data.

## Installation
This package is unregistered, so it should be installed by:
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/ajgiuliani/msJ.jl"))
```
or by:
```julia
]

(v1.1) pkg> add https://github.com/ajgiuliani/msJ.jl
```

## Usage

### getting information on an mzXML file

```julia
msJ.info(filename)
```

### Loading mzXML data

```julia
scans = msJ.load(filename)
```

Individual mass spectra are stored in structures called MSscan reflecting the organization of the mzXML files.

| Field             | Description                                           | Type            |
|-------------------|-------------------------------------------------------|-----------------|
| num               | scan number                                           | Int             |
| rt                | retention time (seconds)                              | Float64         |
| tic               | total ion current                                     | Float64         |
| mz                | m/z data                                              | Vector{Float64} |
| int               | intensity data                                        | Vector{Float64  |
| level             | MS level                                              | Int             |
| basePeakMz        | m/z ratio of the base peak                            | Float64         |
| basePeakIntensity | intensity of the base peak                            | Float64         |
| precursor         | m/z ratio of the precursor ion for MS/MS acquisitions | Float64         |
| polarity          | Ion mode of the scan ("+" or "-")                     | String          |
| activationMethod  | Activation method used in MS/MS acquisitions          | String          |
| collisionEnergy   | Collision energy used in MS/MS acquisitions           | Float64         |


The spectra contained in the mzXML files are loaded in an array of MSscan.


### Chromatograms
Chromatograms may be obtained from a file or from an array of MSscans previously loaded. The function returns an struct containing an array for both retention time and ion current and the maximum value of the ion current vector.

```julia
run1 = msJ.chromatogram(filename)
run2 = msJ.chromatogram(scans)
```
Both `run1.rt, run2.rt` and `run1.ic, run2.ic` sets are identical if `scans` has been loaded from `filename`.

### Average mass spectra
Average mass spectra may be obtained from a file or from an array of MSscans:

```julia
ms1 = msJ.msfilter(filename)
ms2 = msJ.msfilter(scans)
```
The `msfilter` functions return an `MSscans` structure (! notice the s at the end) which is similar to the previous `MSscan` structure, but keeps tracks of the individual scans used to get the average. Some fields, such as `num`, now become vectors.

### Filters
Chromatogram and average mass spectra may be filtered to the structure fields, such as:

```julia
msJ.msfilter(filename, msJ.Precursor(1255.5))             # where 1255.5 is the m/z ratio of the precursor
msJ.msfilter(filename, msJ.Activation_Energy(18))         # where 18 is the collision energy of the MS/MS spectra
msJ.msfilter(filename, msJ.Activation_Method("CID"))      # where CID is the activation method
msJ.msfilter(filename, msJ.Level(2))                      # where 2 corresponds to the MS^n level
msJ.msfilter(filename, msJ.RT([1,5))                      # where [1,5] stands for retention time from 1 to 5s.
...
```
```julia
msJ.chromatogram(filename, method = msJ.MZ([0, 500]))     # where [0,500] specifies the m/z range to get the total ion current
msJ.chromatogram(filename, method = msJ.MZ([500, 2000]))  # same as above with m/z in the [500,2000] range
msJ.chromatogram(filename, msJ.Level(1))                  # chromaogram for all the MS spectra
...
```

or combined together. For example the script below may be used to obtain the average mass spectrum for the MS2 scans (level = 2), for which the precursor is m/z 1255.5, obtained upon CID with an activation energy of 18 and for retention times in the 1 to 60 s range.

```julia
msJ.msfilter(filename, msJ.Precursor(1255.5),
                       msJ.Activation_Energy(18),
	   	       msJ.Activation_Method("CID"),
		       msJ.Level(2),
		       msJ.RT( [1, 60] ))
```

### Basic operations
Mass spectral data may be added, subtracted, multiplied together or multiplied or divided by a number.

```julia
ms450 = msJ.msfilter(data, msJ.RT([430, 470]))
bk1 = msJ.msfilter(data, msJ.RT([500,550]))
ms450_bkg = ms450 - 0.95 * bk1
```

## Other packages
* [mzXML](https://github.com/timholy/mzXML.jl): Load mass spectrometry mzXML files.
