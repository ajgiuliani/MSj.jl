"""
Utility functions to work on MScontainer data types.

"""
module utilities

using Interpolations
using LinearAlgebra

import Base: +, -, *, /

# User Interface.
# --------------

export avg    

# Overloaded Base function
# ------------------------

"""
    /(a::MSscan, N::Real)
Divide the intensity and the tic data of a MSscan by a number
```julia-repl
julia> scans[1] / 1.0e2
msJ.MSscan(1, 0.1384, 50819.5, [140. ....
```
"""
function /(a::MSscan, N::Real)
    return MSscan(a.num, a.rt, a.tic / N, a.mz, a.int ./ N, a.level, a.basePeakMz, a.basePeakIntensity / N, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy)
end


"""
    /(a::MSscans, N::Real)
Divide in the intenisty, tic and variance of a MSscans by a number
```julia-repl
julia> a / 1.0e2
msJ.MSscans(1, 0.1384, 50819.5, [140. ....
```
"""
function /(a::MSscans, N::Real)
    return MSscans(a.num, a.rt, a.tic / N, a.mz, a.int ./ N, a.level, a.basePeakMz, a.basePeakIntensity / N, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy, a.s / N)
end

"""
    *(a::MSscan, N::Real)
Multiply the intensity and the tic data of a MSscan by a number
```julia-repl
julia> scans[1] * 1.0e2
msJ.MSscan(1, 0.1384, 50819.5, [140. ....
```
"""
function *(a::MSscan, N::Real)
    return MSscan(a.num, a.rt, a.tic * N, a.mz, a.int .* N, a.level, a.basePeakMz, a.basePeakIntensity * N, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy)
end

"""
    *(a::MSscans, N::Real)
Multiply in the intenisty, tic and variance of a MSscans by a number
```julia-repl
julia> a * 1.0e2
msJ.MSscans(1, 0.1384, 50819.5, [140. ....
```
"""
function *(a::MSscans, N::Real)
    return MSscans(a.num, a.rt, a.tic * N, a.mz, a.int .* N, a.level, a.basePeakMz, a.basePeakIntensity * N, a.precursor, a.polarity, a.activationMethod, a.collisionEnergy, a.s * N)
end

"""
    *(N::Real, a::MScontainer)
Commutation of multiplication of number with MSscontainer
"""
function *(N::Real, a::MScontainer)
    return a * N
end


"""
    *(a::MScontainer, b::MScontainer)
Multiplication of mass spectra elementwise.
```julia-repl
julia> a * b
msJ.MSscans([2, 5], [0.7307, 4.344
```
"""
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
        mz = a.mz
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


"""
    -(a::MScontainer, b::MScontainer)
Substraction of mass spectra elementwise. Negative scan num refers the 'b' MScontainer
```julia-repl
julia> a - b
msJ.MSscans([1, 4], [0.1384, 3.7578, -0.1384, -3.7578]...
```
"""
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


"""
    +(a::MScontainer, b::MScontainer)
Addition of mass spectra elementwise.
```julia-repl
julia> scans[1] - scans[2]
msJ.MSscans([1, 2], [0.1384, 0.7307]
```
"""
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

"""
    avg(a::MScontainer, b::MScontainer)
Returns the average of the input mass spectra and compute the variance using an incremental Welford algorithm.
```julia-repl
julia> msJ.avg(scans[1], scans[4])
msJ.MSscans([1, 4], [0.1384, 3.7578], ....
```
"""
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

"""
    add_ion_current(x::Vector{Float64}, y::Vector{Float64}, a::Float64, b::Float64)
Returns sum the ion current (int) within the m/z range defined by the a and b input values.
"""
function add_ion_current(x::Vector{Float64}, y::Vector{Float64}, a::Float64, b::Float64)
    ia = num2pnt(x, a)
    ib = num2pnt(x, b)
    return sum( y[ia:ib] )
end


"""
    num2pnt(x::Vector{Float64}, val::Real)
General purpose utility function used to retrieve the index of an array for which the value is closest to the input.
"""
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



end
