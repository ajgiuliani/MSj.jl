# Information
The [`info`](@ref) public function reads the content of a file, but without loading the mass spectrometry data, and returns a `Vector{String}`containing the number of scans, scans level and for MS/MS data, the precursor m/z, the activation method and energy. Additional information may be gained by setting `verbose = true`.
```julia-repl
info(filename)
4-element Array{String,1}:
 "51 scans"
 "MS1+"
 "MS2+ 1255.5  CID(CE=18)"
 "MS3+ 902.33  PQD(CE=35)"
```
