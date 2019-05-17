module msJ

using LightXML, Codecs
using Interpolations
using Unicode
using LsqFit

using Statistics
import Base: +, -, *, /

### Containers
    
abstract type MScontainer  end
    
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


### Methods
# chromatogram
abstract type MethodType  end

struct BasePeak <: MethodType
   #field = "base peak"
   BasePeak() = new()
end

struct TIC <: MethodType
   #field = "TIC"
   TIC() = new()
end

struct ∆MZ{argT <: Union{Real, AbstractVector{<:Real} }} <: MethodType
   arg::argT
   #field = "mz range"
   ∆MZ(arg::argT) where{argT} = new{argT}(arg)
end

struct MZ{argT <: Union{Real, AbstractVector{<:Real} }} <: MethodType
   arg::argT
   #field = "mz range"
   MZ(arg::argT) where{argT} = new{argT}(arg)
end


### Filters

abstract type FilterType  end

struct RT{argT <: Union{Real, AbstractVector{<:Real},  AbstractVector{<:AbstractVector{<:Real}} }} <: FilterType
   arg::argT
   RT(arg::argT) where{argT} = new{argT}(arg)
end

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



greets() = print("Hello World!")


### Overloaded Base function

function /(a::MSscan, N::Real)
    return MSscan(a.num, a.rt, a.tic / N, a.mz, a.int ./ N, a.level, a.basePeakMz, a.basePeakIntensity / N, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy)
end

function /(a::MSscans, N::Real)
    return MSscans(a.num, a.rt, a.tic / N, a.mz, a.int ./ N, a.level, a.basePeakMz, a.basePeakIntensity / N, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy, a.s / N)
end

function *(a::MSscan, N::Real)
    return MSscan(a.num, a.rt, a.tic * N, a.mz, a.int .* N, a.level, a.basePeakMz, a.basePeakIntensity * N, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy)
end

function *(a::MSscans, N::Real)
    return MSscans(a.num, a.rt, a.tic * N, a.mz, a.int .* N, a.level, a.basePeakMz, a.basePeakIntensity * N, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy, a.s * N)
end

function *(N::Real, a::MScontainer)
    return a * N
end

function *(a::MScontainer, b::MScontainer)
    num = vcat(a.num, b.num)
    rt  = vcat(a.rt, b.rt)
    tic = a.tic * b.tic
    eps = 1e-5
    m   = Vector{Float64}(undef,0)
    s   = Vector{Float64}(undef,0)
    int = Vector{Float64}(undef,0)
    mz  = Vector{Float64}(undef,0)
    
    if length(a.mz) > length(b.mz)
        extrap = LinearInterpolation(b.mz, b.int, extrapolation_bc = Line())
        int = a.int .* [extrap(x) for x in a.mz ]
        mz = a.mz
      
    elseif length(a.mz) < length(b.mz)
        extrap = LinearInterpolation(a.mz, a.int, extrapolation_bc = Line())
        int = [extrap(x) for x in b.mz ] .* b.int
        mz = b.mz
            
    elseif length(a.mz) == length(b.mz)
        mz = a.mzq
        int = a.int .* b.int
    end

    basePeakIntensity = maximum(int)
    basePeakMz = mz[num2pnt(int, basePeakIntensity)]

    level = vcat(a.level, b.level)
    precursor = vcat(a.precursor, b.precursor)
    polarity = vcat(a.polarity, b.polarity)
    activationMethod = vcat( a.activationMethod,  b.activationMethod)
    collisionEnergy = vcat(a.collisionEnergy, b.collisionEnergy)
    
    return  MSscans(num, rt, tic, mz, int, level, basePeakMz, basePeakIntensity, precursor, polarity, activationMethod, collisionEnergy, s)
end

function -(a::MScontainer, b::MScontainer)
    #num = vcat(a.num, -b.num)
    num = vcat(a.num)
    rt  = vcat(a.rt, -b.rt)
    tic = a.tic - b.tic

    s    = Vector{Float64}(undef,0)
    int  = Vector{Float64}(undef,0)
    mz   = Vector{Float64}(undef,0)
    
    if length(a.mz) > length(b.mz)
        extrap = LinearInterpolation(b.mz, b.int, extrapolation_bc = Line())
        int = a.int - [extrap(x) for x in a.mz ]
        mz = a.mz
        if a isa MSscans
            if b isa MSscans
                extrap2 = LinearInterpolation(b.mz, b.s, extrapolation_bc = Line())
                s = a.s + [extrap2(x) for x in a.mz]
            else
                s = a.s
            end
        elseif a isa MSscan
            if b isa MSscans
                extrap2 = LinearInterpolation(b.mz, b.s, extrapolation_bc = Line())
                s = [extrap2(x) for x in a.mz]
            end            
        end

    elseif length(a.mz) < length(b.mz)
        extrap = LinearInterpolation(a.mz, a.int, extrapolation_bc = Line())
        int = [extrap(x) for x in b.mz ] - b.int
        mz = b.mz
        if a isa MSscans
            extrap2 = LinearInterpolation(a.mz, a.s, extrapolation_bc = Line())
            if b isa MSscans
                s = [extrap2(x) for x in b.mz] + b.s
            else
                s = [extrap2(x) for x in b.mz]
            end
        elseif a isa MSscan
            if b isa MSscans
                s = b.s
            end            
        end

    elseif length(a.mz) == length(b.mz)
        mz = a.mz
        int = a.int - b.int
        if a isa MSscans
            if b isa MSscans
                s = a.s + b.s
            else
                s = a.s
            end
        elseif a isa MSscan
            if b isa MSscans
                s = b.s
            end
        end         

    end

    basePeakIntensity = maximum(int)
    basePeakMz = mz[num2pnt(int, basePeakIntensity)]
    
    level = vcat(a.level, b.level)
    precursor = vcat(a.precursor, b.precursor)
    polarity = vcat(a.polarity, b.polarity)
    activationMethod = vcat( a.activationMethod,  b.activationMethod)
    collisionEnergy = vcat(a.collisionEnergy, b.collisionEnergy)
    
    return  MSscans(num, rt, tic, mz, int, level, basePeakMz, basePeakIntensity, precursor, polarity, activationMethod, collisionEnergy, s)
end

function +(a::MScontainer, b::MScontainer)
    num = vcat(a.num, b.num)
    rt  = vcat(a.rt, b.rt)
    tic = a.tic + b.tic

    m   = Vector{Float64}(undef,0)
    s   = Vector{Float64}(undef,0)
    int = Vector{Float64}(undef,0)
    mz  = Vector{Float64}(undef,0)
    
    if length(a.mz) > length(b.mz)
        extrap = LinearInterpolation(b.mz, b.int, extrapolation_bc = Line())
        int = a.int + [extrap(x) for x in a.mz ]
        mz = a.mz
        if a isa MSscans
            if b isa MSscans
                extrap2 = LinearInterpolation(b.mz, b.s, extrapolation_bc = Line())
                s = a.s + [extrap2(x) for x in a.mz]
            else
                s = a.s
            end
        elseif a isa MSscan
            if b isa MSscans
                extrap2 = LinearInterpolation(b.mz, b.s, extrapolation_bc = Line())
                s = [extrap2(x) for x in a.mz]
            end            
        end
      
    elseif length(a.mz) < length(b.mz)
        extrap = LinearInterpolation(a.mz, a.int, extrapolation_bc = Line())
        int = [extrap(x) for x in b.mz ] + b.int
        mz = b.mz
        if a isa MSscans
            extrap2 = LinearInterpolation(a.mz, a.s, extrapolation_bc = Line())
            if b isa MSscans
                s = [extrap2(x) for x in b.mz] + b.s
            else
                s = [extrap2(x) for x in b.mz]
            end
        elseif a isa MSscan
            if b isa MSscans
                s = b.s
            end            
        end
            
    elseif length(a.mz) == length(b.mz)
        mz = a.mz
        int = a.int + b.int
        if a isa MSscans
            if b isa MSscans
                s = a.s + b.s
            else
                s = a.s
            end
        elseif a isa MSscan
            if b isa MSscans
                s = b.s
            end
        end         
    end

    basePeakIntensity = maximum(int)
    basePeakMz = mz[num2pnt(int, basePeakIntensity)]

    level = vcat(a.level, b.level)
    precursor = vcat(a.precursor, b.precursor)
    polarity = vcat(a.polarity, b.polarity)
    activationMethod = vcat( a.activationMethod,  b.activationMethod)
    collisionEnergy = vcat(a.collisionEnergy, b.collisionEnergy)
    
    return  MSscans(num, rt, tic, mz, int, level, basePeakMz, basePeakIntensity, precursor, polarity, activationMethod, collisionEnergy, s)
end



### Average & Variance calculation using an incremental Welford algorithm

function avg(a::MScontainer, b::MScontainer)
    num = vcat(a.num, b.num)
    rt  = vcat(a.rt, b.rt)
    tic = a.tic + b.tic

    m   = Vector{Float64}(undef,0)
    s   = Vector{Float64}(undef,0)
#    int = Vector{Float64}(undef,0)
    mz  = Vector{Float64}(undef,0)
    
    if length(a.mz) > length(b.mz)
        # interpolating b.int to the a.mz values and adding the result to a.int
        extrap = LinearInterpolation(b.mz, b.int, extrapolation_bc = Line())
#        int = a.int + [extrap(x) for x in a.mz]
        mz = a.mz
        if a isa MSscans
            m = a.int + ([extrap(x) for x in a.mz] - a.int) / (length(a.num)+1 )
            s = a.s + ([extrap(x) for x in a.mz] - m) .* ([extrap(x) for x in a.mz] - a.int)
#            s = (a.s .* a.s) + ([extrap(x) for x in a.mz] - m) .* ([extrap(x) for x in a.mz] - a.int)
        elseif a isa MSscan
            m = a.int + ([extrap(x) for x in a.mz] - a.int) / 2
            s = ([extrap(x) for x in a.mz] - m) .* ([extrap(x) for x in a.mz] - a.int)
        end
    elseif length(a.mz) < length(b.mz)
        # interpoling a.int to the b.mz values and adding the result to b.int
        extrap = LinearInterpolation(a.mz, a.int, extrapolation_bc = Line())
#        int = [extrap(x) for x in b.mz] + b.int
        mz = b.mz
        if a isa MSscans
            extrap2 = LinearInterpolation(a.mz, a.s, extrapolation_bc = Line())
            m = [extrap(x) for x in b.mz]  + (b.int -[extrap(x) for x in b.mz] ) / (length(a.num)+1 )
            s = [extrap2(x) for x in b.mz] + (b.int - m) .* (b.int - [extrap(x) for x in b.mz ])
#            s = ([extrap2(x) for x in b.mz] .* [extrap2(x) for x in b.mz]) + (b.int - m) .* (b.int - [extrap(x) for x in b.mz ])
        elseif a isa MSscan
            m = [extrap(x) for x in b.mz] + (b.int - [extrap(x) for x in b.mz] ) / 2
            s = (b.int - m) .* (b.int - [extrap(x) for x in b.mz ])
        end
    elseif length(a.mz) == length(b.mz)
        #    else
        mz = a.mz
#        int = a.int + b.int
        if a isa MSscans
            m = a.int + (b.int - a.int) / (length(a.num)+1 )
            s = a.s  + ( (b.int - m) .* (b.int - a.int) )
#            s = (a.s .* a.s)  + ( (b.int - m) .* (b.int - a.int) )
        elseif a isa MSscan
            m = a.int + (b.int - a.int) / 2
            s = (b.int - m) .* (b.int - a.int)
        end
    end

    basePeakIntensity = maximum(m)
    basePeakMz = mz[num2pnt(m, basePeakIntensity)]


    level = vcat(a.level, b.level)
    precursor = vcat(a.precursor, b.precursor)
    polarity = vcat(a.polarity, b.polarity)
    activationMethod = vcat(a.activationMethod,  b.activationMethod)
    collisionEnergy = vcat(a.collisionEnergy, b.collisionEnergy)
    
    return  MSscans(num, rt, tic, mz, m, level, basePeakMz, basePeakIntensity, precursor, polarity, activationMethod, collisionEnergy, s)
end

function standard_deviation(a::MSscans, N::Int)
    return MSscans(a.num, a.rt, a.tic, a.mz, a.int, a.level, a.basePeakMz, a.basePeakIntensity, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy, map(sqrt, (a.s / (N -1)) ) )
end



### Utility functions

function add_ion_current(x::Vector{Float64}, y::Vector{Float64}, a::Float64, b::Float64)
    ia = num2pnt(x, a)
    ib = num2pnt(x, b)
    return sum( y[ia:ib] )
end


function num2pnt(x::Vector{Float64}, val::Real)
    ibest = 1
    dxbest = abs(x[ibest] - val)
    for i in eachindex(x)
        dx = abs(x[i] - val)
        if dx < dxbest
            dxbest = dx
            ibest = i
        end
    end
    return ibest
end 


include("msJ_info.jl")
include("msJ_load.jl")
include("msscans.jl")
include("mzxml.jl")
include("analysis.jl")

end # module
