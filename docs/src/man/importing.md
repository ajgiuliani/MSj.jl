# Importing data
When loading a file containing more than a single acquisition, the individual mass spectrometry scans are pushed into an array of [`MSj.MSscan`](@ref).  The individual scans may be retrieve from the array the usual way:

```julia-repl
julia> scans = load("filename")
51-element Array{MSj.MSscan,1}:
 MSj.MSscan(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....)
...

julia> scans[1]
MSj.MSscan(1, 0.1384, 5.08195e6, [140.083, 140.167, 140.25, 140.333, 140.417, 140.5, 140.583, 140.667, 140.75, 140.833  …  1999.25, 1999.33, 1999.42, ....)
```

As mentioned above, chromatograms may be retrieved from a file and imported in [`MSj.Chromatogram`](@ref) :
```julia-repl
julia> chromatogram("filename")
MSj.Chromatogram([0.1384, 0.7307, 2.1379, 3.7578, 4.3442, 5.7689], [5.08195e6, 9727.2, 11.3032, 4.8084e6, 12203.5, 4.84455], 5.08195e6)
```

The function [`MSj.retention_time`](@ref) reads the retention time of an input file and returns a `Vector{Float64}`containing the time in seconds.
```julia-repl
julia> MSj.retention_time("filename")
51-element Array{Float64,1}:
  0.1384
  0.7307
  2.1379
....MSj.FilterType

```
