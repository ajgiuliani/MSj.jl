"""
Main module for `msJ.jl`-- A Julia package to load and process mass spectrometry data

This module exports :

$(EXPORTS)

"""
module msJ

using DataFrames
using DataFramesMeta

# User Interface.
# ---------------

#export msJ.load, msfilter, chromatogram, smooth, centroid
#export FilterType, MScontainer, MSscan, MSscan, SG

#greets() = print("Hello World!")


# Submodules
# ----------

include("types.jl")
include("msscans.jl")
include("mzxml.jl")
include("process.jl")
include("plots.jl")
include("io.jl")
include("utilities.jl")


end # module
