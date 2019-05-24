"""
Submodule with types and structures used to stored the data and dispatch to the right methods.
"""
module types


### Containers

"""
Abstract type containing any imported data belongs to the MScontainer type.
"""
abstract type MScontainer  end


"""
Data structure used to store individual mass spectrometry scans organized following the structure of mzXML files.
"""
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

"""
Data structure designed to store mass spectra obtained after filtering operation along with the history of these operation.
"""
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


"""
Data structure used to retrieve chromatography data.
"""
struct Chromatogram  <: MScontainer
    rt::Vector{Float64}
    ic::Vector{Float64}
    maxic::Float64
end



### Methods
"""
This type contains all the methods used for filtering the data.
"""
abstract type MethodType  end


# chromatogram

"""
Structure for multiple dispatching to retrieve base peak chromatogram.
"""
struct BasePeak <: MethodType
   #field = "base peak"
   BasePeak() = new()
end


"""
Structure for multiple dispatching to retrieve total ion current chromatogram.
"""
struct TIC <: MethodType
   #field = "TIC"
   TIC() = new()
end


"""
Structure for multiple dispatching to retrieve extracted ion current chromatogram around an m/z ± ∆mz value given by arg = [mz, ∆mz]
"""
struct ∆MZ{argT <: Union{Real, AbstractVector{<:Real} }} <: MethodType
   arg::argT
   #field = "∆mz range"
   ∆MZ(arg::argT) where{argT} = new{argT}(arg)
end


"""
Structure for multiple dispatching to retrieve extracted ion current chromatogram around for m/z in the range arg = [mz1, mz2].
"""
struct MZ{argT <: Union{Real, AbstractVector{<:Real} }} <: MethodType
   arg::argT
   #field = "mz range"
   MZ(arg::argT) where{argT} = new{argT}(arg)
end


"""
Structure for multiple dispatching to Savinsky & Golay filtering, providing the order, window size and derivative to be performed.  Defaults values are provided in functions calls.
"""
struct SG{argT <: Int} <: MethodType   #Savinsky & Golay filtering
    order::argT
    window::argT
    derivative::argT
    SG(order::argT, window::argT, derivative::argT) where{argT} = new{argT}(order, window, derivative)
end


"""
Structure for multiple dispatching to Template Base Peak Detection centroiding, providing the shape of the template function, the resolution and threshold.  Defaults values are provided in functions calls.
"""
struct TBPD{argT1 <: Symbol, argT2 <: Real}  <: MethodType
    shape::argT1
    resolution::argT2
    threshold::argT2
    TBPD(shape::argT1, resolution::argT2, threshold::argT2) where{argT1, argT2} = new{argT1, argT2}(shape, resolution, threshold)
end


"""
struct SNRA{argT <: Real}  <: MethodType
    threshold::argT
    SNRA(threshold::argT) where{argT} = new{argT}(threshold)
end

struct CWT{argT <: Real}  <: MethodType
    threshold::argT
    CWT(threshold::argT) where{argT} = new{argT}(threshold)
end
"""

### Filters

"""
    FilterType 
This type contains  the structures for filtering the data.
"""
abstract type FilterType end


"""
    RT{argT <: Union{Real, AbstractVector{<:Real},  AbstractVector{<:AbstractVector{<:Real}} }}
This type contains  the structures for filtering the data.
"""
struct RT{argT <: Union{Real, AbstractVector{<:Real},  AbstractVector{<:AbstractVector{<:Real}} }} <: FilterType
   arg::argT
   RT(arg::argT) where{argT} = new{argT}(arg)
end

"""
    IC{argT <: Union{Real, AbstractVector{<:Real} }}
Structure for multiple dispatching to Template Base Peak Detection centroiding, providing the shape of the template function, the resolution and threshold.  Defaults values are provided in functions calls.
"""
struct IC{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
   arg::argT
   IC(arg::argT) where{argT} = new{argT}(arg)
end

struct Level{argT <: Union{Int, AbstractVector{<:Int} }} <: FilterType
   arg::argT
   #field = "level"
   Level(arg::argT) where{argT} = new{argT}(arg)
end

struct Scan{argT <: Union{Int, AbstractVector{<:Int} }} <: FilterType
   arg::argT
   #field = "num"
   Scan(arg::argT) where{argT} = new{argT}(arg)
end

struct Polarity{argT <: Union{String, AbstractVector{<:String} }} <: FilterType
   arg::argT
   #field = "polarity"
   Polarity(arg::argT) where{argT} = new{argT}(arg)
end

struct Activation_Method{argT <: Union{String, AbstractVector{<:String} }} <: FilterType
   arg::argT
   #field = "activationMethod"
   Activation_Method(arg::argT) where{argT} = new{argT}(arg)
end

struct Activation_Energy{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
   arg::argT
   #field = "collisionEnergy"
   Activation_Energy(arg::argT) where{argT} = new{argT}(arg)
end

struct Precursor{argT <: Union{Real, AbstractVector{<:Real} }} <: FilterType
   arg::argT
   #field = "precursorMz"
   Precursor(arg::argT) where{argT} = new{argT}(arg)
end



end # submodule
