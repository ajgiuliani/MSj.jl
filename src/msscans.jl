"""
Interface to the MScontainer data for filtering scans
"""
module msscans


"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.Level{<:Int})
Search for scans matching the argument MS level and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.Level{<:Int})
    subindex = Set{Int}()
    for elem in scans
        if elem.level == argument.arg
            push!(subindex, elem.num)
        end
    end
    return subindex
end

"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.Level{<:AbstractVector})
Search for scans matching the argument MS levels and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.Level{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg
        for elem in scans
            if elem.level == i
                push!(subindex, elem.num)
            end
        end
    end
    return subindex
end

"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.Precursor{<:Real})
Search for scans matching the argument precursor mz and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.Precursor{<:Real})
    subindex = Set{Int}()
    for elem in scans
        if elem.precursor == argument.arg
            push!(subindex, elem.num)
        end
    end
    return subindex
end


"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.Precursor{<:AbstractVector})
Search for scans matching the argument precursors mz and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.Precursor{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg
        for elem in scans
            if elem.precursor == i
                push!(subindex, elem.num)
            end
        end        
    end
    return subindex
end

"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.Activation_Energy{<:Real})
Search for scans matching the argument activation energy and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.Activation_Energy{<:Real})
    subindex = Set{Int}()
    for elem in scans
        if elem.collisionEnergy == argument.arg
            push!(subindex, elem.num)
        end
    end
    return subindex
end

"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.Activation_Energy{<:AbstractVector})
Search for scans matching the argument activation energies and returns a list of the corresponding indexes

"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.Activation_Energy{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg       
        for elem in scans
            if elem.collisionEnergy == i
                push!(subindex, elem.num)
            end
        end
    end
    return subindex
end

"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.Activation_Method{<:String})
Search for scans matching the argument activation method and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.Activation_Method{<:String})
    subindex = Set{Int}()
    for elem in scans
        if elem.activationMethod == argument.arg
            push!(subindex, elem.num)
        end
    end
    return subindex
end

"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.Activation_Method{<:AbstractVector})
Search for scans matching the argument activation methods and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.Activation_Method{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg
        for elem in scans
            if elem.activationMethod == i
                push!(subindex, elem.num)
            end
        end
    end
    return subindex
end

"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.Polarity{<:AbstractVector})
Search for scans matching the argument polarities and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.Polarity{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg
        for elem in scans
            if elem.polarity == i
                push!(subindex, elem.num)
            end
        end
    end
    return subindex
end


"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.Polarity{<:String})
Search for scans matching the argument polarity and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.Polarity{<:String})
    subindex = Set{Int}()
    for elem in scans
        if elem.polarity == argument.arg
            push!(subindex, elem.num)
        end
    end
    return subindex
end


"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.Scan{<:Int})
Search for scans matching the argument scan num and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.Scan{<:Int})
    subindex = Set{Int}()
    for elem in scans
        if elem.num == argument.arg
            push!(subindex, elem.num)
        end
    end
    return subindex
end


"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.Scan{<:AbstractVector})
Search for scans matching the argument scan nums and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.Scan{<:AbstractVector})
    subindex = Set{Int}()
    for i in argument.arg       
        for elem in scans
            if elem.num == i
                push!(subindex, elem.num)
            end
        end
    end
    return subindex
end


"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.RT{<:Real}) 
Search for scans matching the argument retention time and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.RT{<:Real}) 
    subindex = Set{Int}()
    rt = retention_time(scans)
    index = num2pnt(rt, argument.arg)
    for elem in scans
        if elem.num == index
            push!(subindex, elem.num)
        end
    end
    return subindex
end


"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.RT{<:AbstractVector})
Search for scans matching the argument retention time in the specified range and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.RT{<:AbstractVector})
    subindex = Set{Int}()
    rt = retention_time(scans)
    for i =1:2:length(argument.arg)
        index_low  = num2pnt(rt, argument.arg[i])
        index_high = num2pnt(rt, argument.arg[i+1])
        for elem in scans
            if  index_low <= elem.num <= index_high
                push!(subindex, elem.num)
            end
        end
    end
    return subindex
end


"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.RT{<:AbstractVector{<:AbstractVector} } )
Search for scans matching the argument retention time in the specified ranges and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.RT{<:AbstractVector{<:AbstractVector} } )
    subindex = Set{Int}()
    rt = retention_time(scans)
    for el  in argument.arg
        index_low  = num2pnt(rt, el[1])
        index_high = num2pnt(rt, el[2])
        for elem in scans
            if  index_low <= elem.num <= index_high
                push!(subindex, elem.num)
            end
        end
    end
    return subindex
end


"""
    filter(scans::Vector{msJ.MSscan}, argument::msJ.IC{<:AbstractVector})
Search for scans matching the argument total ion current within the specified ranges and returns a list of the corresponding indexes
"""
function filter(scans::Vector{msJ.MSscan}, argument::msJ.IC{<:AbstractVector})
    subindex = Set{Int}()
    for elem in scans
        if  argument.arg[1] <= elem.tic <= argument.arg[2]
            push!(subindex, elem.num)
        end
    end
    return subindex
end



### Extraction of the ion current according to the selected filters and method

"""
    extracted_chromatogram(scans::Vector{msJ.MSscan}, indices::Vector{Int},method::msJ.MethodType)
Returns the extracted chromatogram for input Array of MSscan according to the selected method and for set of scan num as input
"""
function extracted_chromatogram(scans::Vector{msJ.MSscan}, indices::Vector{Int},method::msJ.MethodType)
    xrt = Vector{Float64}(undef,0)
    xic = Vector{Float64}(undef,0)
    if method isa BasePeak
        for i = 1:length(indices)
            push!(xrt, scans[indices[i]].rt)
            push!(xic, scans[indices[i]].basePeakIntensity)
        end
    elseif method isa ∆MZ
        mz1 = convert(Float64, method.arg[1] - method.arg[2] )  # mz - ∆mz
        if(mz1 < 0.0)
            return ErrorException("Bad mz ± ∆mz values.")
            #error("Bad mz ± ∆mz values")
        end
        mz2 = convert(Float64, method.arg[1] + method.arg[2] ) # mz + ∆mz
        for i = 1:length(indices)
            value = add_ion_current(scans[indices[i]].mz, scans[indices[i]].int, mz1, mz2)
            push!(xrt, scans[indices[i]].rt)
            push!(xic, value)
        end
        
    elseif method isa MZ
        mz1 = convert(Float64,method.arg[1])
        mz2 = convert(Float64,method.arg[2])
        for i = 1:length(indices)
            value = add_ion_current(scans[indices[i]].mz, scans[indices[i]].int, mz1, mz2)
            push!(xrt, scans[indices[i]].rt)
            push!(xic, value)
        end
    else
        for i = 1:length(indices)
            push!(xrt, scans[indices[i]].rt)
            push!(xic, scans[indices[i]].tic)
        end

    end
    return Chromatogram(xrt, xic, maximum(xic)) 
end




### Calculation of the composite MSscans composite mass spectrum according to the selected filters and method


"""
    composite_spectra(scans::Vector{msJ.MSscan}, indices::Vector{Int}, stats::Bool)
Returns the average MSscans for input Array of MSscan and according to the input scan num. Calculation of variance is controlled by the stats Boolean variable.
"""
function composite_spectra(scans::Vector{msJ.MSscan}, indices::Vector{Int}, stats::Bool)
    if stats == false
        result = scans[indices[1]]
        for i = 2:length(indices)
            result += scans[indices[i]]
        end
        return result / length(indices)
    elseif stats == true
        result = scans[indices[1]]
        for i = 2:length(indices)
            result = avg(result, scans[indices[i]])
        end
        return standard_deviation(result, length(indices))
    end
end





end
