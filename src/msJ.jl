"""
Main module for `msJ.jl`-- A Julia package to load and process mass spectrometry data

"""
module msJ

using DataFrames
using DataFramesMeta


include("types.jl")
include("Io.jl")
include("msscans.jl")
include("mzxml.jl")
include("process.jl")
include("extract.jl")
include("utilities.jl")


# Submodules
# ----------

include("plots.jl")



end # module
