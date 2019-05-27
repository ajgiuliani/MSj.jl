Tutorial
========

The msJ package intends to provide an access to the common open source mass spectrometry file format.

# Loading mass spectrometry data
For this example, the mzXML file from the test folder of the package will be used. First, the [`info`](@ref) function may be used to see what's inside the file.

```@example 1
using msJ
info("../../test/test.mzXML")

```
From the output, we learn that the file contains MS scans in the positive ion mode, MS/MS scans also in positive mode, where the precursor m/z 1255.5 is activated under CID conditions with 18 collision energy and the m/Z 902.33 precursor is activated using the PQD method with the activation energy set to 35.

Then, the data can be loaded into a variable that we call `scans`.

```@example 1
scans = load("../../test/test.mzXML")
```

As expected, the `scans` variable contains an array of 6 `MSscans`. The fields of the individual scans may be accessed by:

```@example 1
println("scan num      : ", scans[1].num)
```
```@example 1
println("Polarity      : ", scans[1].polarity)
```
```@example 1
println("retention time: ", scans[1].rt)
```
```@example 1
println("MS level      : ", scans[1].level)
```
```@example 1
println("Base peak m/z : ", scans[1].basePeakMz)
```
```@example 1
println("Base peak int.: ", scans[1].basePeakIntensity)
```

The mass spectrum is stored in two arrays:
```julia
mz_2  =  scans[2].mz
int_2 =  scans[2].int ; nothing
```

Getting chromatograms is straightforward using the [`chromatogram`](@ref) method:
```@example 1
full_TIC   = chromatogram(scans)
```
```@example 1
MS1_TIC    = chromatogram(scans, msJ.Level(1))
```
```@example 1
CID_TIC    = chromatogram(scans, msJ.Activation_Method("CID"))
```
```@example 1
mz902_TIC  = chromatogram(scans, msJ.Precursor(902.33))
```

The individual extracted chromatogram may be plotted. First we need to import the `Julia` package for plotting:
```@example 1
Using Plots
```
Then, we can make a figure with the chromatograms:
```@example 1
using Plots
pyplot()
p1 = plot(full_TIC, label = "full tic")
p2 = plot(MS1_TIC, label = " MS1_TIC")
p3 = plot(CID_TIC, label = "CID_TIC")
p4 = plot(mz902_TIC, label = "mz902_TIC")

plot(p1, p2, p3, p4, layout = (4,1))
#savefig("chromatograms.png"); nothing # hide
savefig("chromatograms.png"); nothing # hide
```
![](chromatograms.png)

Average mass spectra may be obtained using the proper [`msfilter`](@ref) functions:
```@example 1
ms1 = msfilter(scans, msJ.Level(1))
ms2_CID = msfilter(scans, msJ.Activation_Method("CID"))
ms2_PQD = msfilter(scans, msJ.Activation_Method("PQD"))
ms2_1255 = msfilter(scans, msJ.Precursor(1255.5))
```

and then plotted similarly:
```@example 1
p5 = plot(ms1, label = "MS", color = :blue)
p6 = plot(ms2_CID, label = "CID", color = :green)
p7 = plot(ms2_PQD, label = "PQD", color = :purple)
p8 = plot(ms2_1255, label = "mz 1255", color = :orange)
plot(p5, p6, p7, p8, layout = (4,1), size = (800,600))
savefig("ms.png"); nothing # hide
```
![](ms.png)

The average mass spectrum for CID (green) and for the m/z 1255.5 precursor ion (orange) are identical, which is not surprising.
