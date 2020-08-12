#=
cpa_filter.jl

Filter functions for CpaAndTime datatypes

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#
include("./cpa_structs.jl")

################################################################################
# Filter by CpaFilter Type
################################################################################

function filterCpasAndTimeVec(cpasandtimevec::Vector{CpasAndTime}, filter::CpaFilter)

    filtered_cpasandtimevec = CpasAndTime[]

    for cpasandtime in cpasandtimevec
        filteredCpas = filterCpasAndTime(cpasandtime, filter)

        if !isnothing(filteredCpas)
            push!(filtered_cpasandtimevec, filteredCpas)
        end
    end

    return filtered_cpasandtimevec
end

function filterCpasAndTime(cpasandtime::CpasAndTime, filter::CpaFilter)

    cpas = filterCpas(cpasandtime.cpas, filter)

    if !isempty(cpas)
        return CpasAndTime(cpasandtime.time, cpas)
    else
        return nothing
    end
end

function filterCpas(cpas::Vector{Cpa}, filter::CpaFilter)

    filtervec = falses(length(cpas))

    for (i,cpa) in enumerate(cpas)
        filtervec[i] = filterCpa(cpa, filter)
    end

    return cpas[filtervec]
end

function filterCpa(cpa::Cpa, filter::CpaFilter)::Bool

    if cpa.c_count > filter.c_count_min && cpa.fmax > filter.fmax_min
        return true
    else
        return false
    end
end


################################################################################
# Filter by boundary names
################################################################################

function is_bname_in_cpa(cpa::Cpa, bname::String)::Bool

    for str1 in cpa.bnames
        if str1 == bname
            return true
        end
    end

    return false
end


function are_bnames_in_cpa(cpa::Cpa, bnames::Vector{String})::Bool

    for bname in bnames
        if is_bname_in_cpa(cpa, bname)
            return true
        end
    end

    return false
end


function get_cpas_with_bname(cpasandtime::CpasAndTime, bname::String)

    hascpawithbname = false
    cpaswithbname = Vector{Cpa}(undef,0)

    for cpa in cpasandtime.cpas
        if is_bname_in_cpa(cpa, bname)
            hascpawithbname = true
            push!(cpaswithbname, cpa)
        end
    end

    return hascpawithbname, CpasAndTime(cpasandtime.time, cpaswithbname)
end


function get_cpasandtimevec_with_bname(cpasandtimevec::Vector{CpasAndTime}, bname::String)::Vector{CpasAndTime}

    cpasandtimevecwithbname = Vector{CpasAndTime}(undef,0)

    for cpasandtime in cpasandtimevec
        hascpawithbname, cpasandtimewithbname = get_cpas_with_boundary(cpasandtime, bname)

        if hascpawithbname
            push!(cpasandtimevecwithbname, cpasandtimewithbname)
        end
    end

    return cpasandtimevecwithbname
end


function get_cpas_without_boundaries(cpasandtime::CpasAndTime)

    hasnoboundcpas = false
    cpaswithoutbound = Vector{Cpa}(undef,0)

    for cpa in cpasandtime.cpas
        if cpa.bcount == 0
            hasnoboundcpas = true
            push!(cpaswithoutbound, cpa)
        end
    end

    return hasnoboundcpas, CpasAndTime(cpasandtime.time, cpaswithoutbound)
end


function get_cpasandtimevec_no_boundaries(cpasandtimevec::Vector{CpasAndTime})::Vector{CpasAndTime}

        cpasandtimevecwithoutbound = Vector{CpasAndTime}(undef,0)

        for cpasandtime in cpasandtimevec
            hasnoboundcpas, cpasandtimewithoutbound = get_cpas_without_boundaries(cpasandtime)

            if hasnoboundcpas
                push!(cpasandtimevecwithoutbound, cpasandtimewithoutbound)
            end
        end

        return cpasandtimevecwithoutbound
end


function get_cpas_without_bnames(cpasandtime::CpasAndTime, bnames::Vector{String})

    iscpawithoutbname = false
    cpaswithoutbnames = Vector{Cpa}(undef,0)

    for cpa in cpasandtime.cpas
        if !are_bnames_in_cpa(cpa, bnames)
            iscpawithoutbname = true
            push!(cpaswithoutbnames, cpa)
        end
    end

    return iscpawithoutbname, CpasAndTime(cpasandtime.time, cpaswithoutbnames)
end


function get_cpasandtimevec_without_bnames(cpasandtimevec::Vector{CpasAndTime}, bnames::Vector{String})::Vector{CpasAndTime}

    cpasandtimevecwithoutbnames = Vector{CpasAndTime}(undef,0)

    for cpasandtime in cpasandtimevec
        iscpawithoutbname, cpasandtimewithbname = get_cpas_without_bnames(cpasandtime, boundaryname)

        if iscpawithoutbname
            push!(cpasandtimevecwithoutbnames, cpasandtimewithbname)
        end
    end

    return cpasandtimevecwithoutbnames
end
