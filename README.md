# msJ.jl

[![Build Status](https://travis-ci.org/ajgiuliani/msJ.jl.svg?branch=master)](https://travis-ci.org/ajgiuliani/msJ.jl)

[![Coverage Status](https://coveralls.io/repos/ajgiuliani/msJ.j/badge.svg?branch=master&service=github)](https://coveralls.io/github//ajgiuliani/msJ.j?branch=master)

[![codecov.io](http://codecov.io/github//ajgiuliani/msJ.j/coverage.svg?branch=master)](http://codecov.io/github//ajgiuliani/msJ.j?branch=master)

This package provides an API to the mzXML file format for mass spectrometry data.

## Installation
This package is unregistered, so it should be installed by:
```julia
]
```
```pkg
 add https://github.com/ajgiuliani/msJ.jl
```

## Usage

### getting information on an mzXML file

```julia
msJ.info(filename)
```

### Loading mzXML data

```julia
scans = msJ.load("/Users/alex/hubiC/development/msJ/test.mzXML")
```

Individual mass spectra a stored in a structure called MSscan reflecting that of the mzXML files, containing:

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
|-------------------|-------------------------------------------------------|-----------------|

The spectra contained in the mzXML file are loaded in an array of MSscan.

### Chromatograms
Chromatograms may be obtained from a file or from an array of MSscans and return an array for both retention time and ion current

```julia
rt1, ic1 = msJ.load(filename)
rt2, ic2 = msJ.load(scans)
```
Both sets of rt and ic are identical if scans has been obtained from loading the 'filename' file.

### Average mass spectra
Average mass spectra may be obtained from a file or from an array of MSscans:

```julia
ms1 = msJ.msfilter(filename)
ms2 = msJ.msfilter(scans)
```
The msfilter functions return an MSscans (! notice the s) which is based on the previous MSscan structure, but keeps tracks of the individual scans used to get the average. Some fields, such as num, are  now vectors.

### Filters
Chromatogram and average mass spectra may be filtered to the structure fields, such as

```julia
msJ.msfilter(filename, msJ.Precursor(1255.5))
msJ.msfilter(filename, msJ.Activation_Energy(18))
msJ.msfilter(filename, msJ.Activation_Method("CID"))
msJ.msfilter(filename, msJ.Level(2))
...
```

```julia
msJ.chromatogram(filename, method = msJ.MZ([0, 500]))
msJ.chromatogram(filename, method = msJ.MZ([500, 2000]))
msJ.chromatogram(filename, msJ.Level(1))
...
```

### Basic operations
Mass spectral data may added, subtracted, multiplied or divided by a number.

```julia
ms450 = msJ.msfilter(data, msJ.RT([430, 470]))
bk1 = msJ.msfilter(data, msJ.RT([500,550]))
ms450_bkg = ms450 - 0.95 * bk1
```

## Other packages
* [mzXML][https://github.com/timholy/mzXML.jl]
