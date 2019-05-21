
### Filter functions for MSscan objetcs

function filter(scans::Vector{MSscan}, argument::Level{<:Int})
    subindex = Set{Int}()
    for elem in scans
        if elem.level == argument.arg
            push!(subindex, elem.num)
        end
    end
    return subindex
end

function filter(scans::Vector{MSscan}, argument::Level{<:AbstractVector})
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

function filter(scans::Vector{MSscan}, argument::Precursor{<:Real})
    subindex = Set{Int}()
    for elem in scans
        if elem.precursor == argument.arg
            push!(subindex, elem.num)
        end
    end
    return subindex
end
function filter(scans::Vector{MSscan}, argument::Precursor{<:AbstractVector})
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

function filter(scans::Vector{MSscan}, argument::Activation_Energy{<:Real})
    subindex = Set{Int}()
    for elem in scans
        if elem.collisionEnergy == argument.arg
            push!(subindex, elem.num)
        end
    end
    return subindex
end

function filter(scans::Vector{MSscan}, argument::Activation_Energy{<:AbstractVector})
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

function filter(scans::Vector{MSscan}, argument::Activation_Method{<:String})
    subindex = Set{Int}()
    for elem in scans
        if elem.activationMethod == argument.arg
            push!(subindex, elem.num)
        end
    end
    return subindex
end

function filter(scans::Vector{MSscan}, argument::Activation_Method{<:AbstractVector})
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

function filter(scans::Vector{MSscan}, argument::Polarity{<:AbstractVector})
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

function filter(scans::Vector{MSscan}, argument::Polarity{<:String})
    subindex = Set{Int}()
    for elem in scans
        if elem.polarity == argument.arg
            push!(subindex, elem.num)
        end
    end
    return subindex
end

function filter(scans::Vector{MSscan}, argument::Scan{<:Int})
    subindex = Set{Int}()
    for elem in scans
        if elem.num == argument.arg
            push!(subindex, elem.num)
        end
    end
    return subindex
end

function filter(scans::Vector{MSscan}, argument::Scan{<:AbstractVector})
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

function filter(scans::Vector{MSscan}, argument::RT{<:Real}) 
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

function filter(scans::Vector{MSscan}, argument::RT{<:AbstractVector})
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

function filter(scans::Vector{MSscan}, argument::RT{<:AbstractVector{<:AbstractVector} } )
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

function filter(scans::Vector{MSscan}, argument::IC{<:AbstractVector})
    subindex = Set{Int}()
    for elem in scans
        if  argument.arg[1] <= elem.tic <= argument.arg[2]
            push!(subindex, elem.num)
        end
    end
    return subindex
end

### Extraction of the ion current according to the selected filters and method

function extracted_chromatogram(scans::Vector{MSscan}, indices::Vector{Int},method::MethodType)
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


function composite_spectra(scans::Vector{MSscan}, indices::Vector{Int}, stats::Bool)
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
