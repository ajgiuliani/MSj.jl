"""
Main module for `msJ.jl`-- A Julia package to load and process mass spectrometry data.

"""
module msJ


using Statistics           # used for Perasons correlation calculation
using LsqFit               # used for curve fitting
using DSP                  # used for convolution
using Interpolations       # used for interpolation
using LinearAlgebra        # used for matrix operation 
using LightXML, Codecs     # used for mzXML file import
using Unicode              # used for file io


import Base: +, -, *, /


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
