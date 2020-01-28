#=
cpa_eval_help_functions.jl


License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
=#

include("./cpa_eval_structs.jl")

function orthogonal_dist_cpa_axis(cpa::CpaInfo, origin::Vector{Float64}, dirvec::Vector{Float64})

    return orthogonal_dist_point_axis(
                                        cpa.x, cpa.y, cpa.z,
                                        origin[1], origin[2], origin[3], 
                                        dirvec[1], dirvec[2], dirvec[3]
                                     )

end



function orthogonal_dist_point_axis(
                                    xp, yp, zp, 
                                    originx, originy, originz, 
                                    dirvecx, dirvecy, dirvecz
                                    )
    od = norm(
                cross([originx,originy,originz] .-
                [xp,yp,zp],[dirvecx,dirvecy,dirvecz])
            ) ./ norm([dirvecx,dirvecy,dirvecz])

    return od
end



function getboundarypidsandcells(adjboundstr::Array{String,1}, myid::Int64)
# function which returns 
# 1. boundary processor ids (pids) for processes boundaries adjacent for the 
# processor with id mypid
    bstr_matches = map(adjboundstr) do mystr
        match(r"procBoundary[0-9]+to[0-9]+", mystr)
    end

    bstr_matches = bstr_matches[bstr_matches .!= nothing]

    if length(bstr_matches) > 0
        pidsfrom = zeros(Int64, length(bstr_matches))
        pidsto = zeros(Int64, length(bstr_matches))

        for istrmatches = 1:1:length(bstr_matches)
            pfromtostr = bstr_matches[istrmatches].match

            fromstr = match(r"y[0-9]+t", pfromtostr)
            tostr = match(r"o[0-9]+$", pfromtostr)

            fromstr = fromstr.match[2:end-1]
            tostr = tostr.match[2:end]

            fromid = parse(Int64, fromstr) + 1
            toid = parse(Int64, tostr) + 1

            if fromid != myid 
                error("Error getboundarypidsandcells(): fromid not equal ",
                	"to processor id, this should not happen!")
            end

            pidsto[istrmatches] = toid
        end
    else
        error("Error getboundarypidsandcells(): Here should be something ",
        "to gahter ...")
    end

    return pidsto
end



function prefilter_np_rawinfo(rawinfo::CpaInfoRaw, no_cells_min::Int64, alpha_max_min::Float64)
    # pre filters not parallel  data to reduce memory requirements...
    
    if no_cells_min > rawinfo.no_cells || alpha_max_min > rawinfo.alpha_max
        return nothing 
    else
        return rawinfo
    end    
end



function valofvec1_occursin_vec2(vec1, vec2)

    occurs_valofvec1_in_vec2 = false

    for i = 1:1:length(vec1)
        a = findfirst(isequal(vec1[i]),vec2)

        if !isnothing(a)
            occurs_valofvec1_in_vec2 = true
            break
        end
    end

    return occurs_valofvec1_in_vec2
end



function has_cpa_adj_boundary(cpainfo::CpaInfo, boundaryname::Union{String,Nothing})::Bool

    if isnothing(cpainfo.adj_boundaries) && isnothing(boundaryname)
        return true
    elseif isnothing(cpainfo.adj_boundaries) || isnothing(boundaryname)
        return false
    else
        for str1 in cpainfo.adj_boundaries
            if str1 == boundaryname
                return true
                break
            end
        end

        return false
    end
end



function get_cpas_with_boundary(cpas::Cpas,  boundaryname::Union{String,Nothing})

    hascpa = false
    cpas_with_bname = Array{CpaInfo,1}(undef,0)

    for cpainfo in cpas.cpas
        if has_cpa_adj_boundary(cpainfo, boundaryname)
            hascpa = true
            push!(cpas_with_bname, cpainfo)
        end
    end

    return hascpa, cpas_with_bname
end



function get_cpasarray_with_boundary(cpasarr::Array{Cpas,1}, boundaryname::Union{String,Nothing})
    cpasarr_with_bname = Array{Cpas,1}(undef,0)

    for cpas in cpasarr
        hascpas, cpaswithbname = get_cpas_with_boundary(cpas, boundaryname)

        if hascpas
            push!(cpasarr_with_bname, Cpas(cpas.time, cpaswithbname))
        end
    end

    return cpasarr_with_bname
end



function get_cpasarray_without_boundaries(cpasarr::Array{Cpas,1}, boundarynames_arr::Array{String,1})

    cpasarr_without_bnames = Array{Cpas,1}(undef,0)

    for cpas in cpasarr

        curret_cpas_have_bname = false
        cpainfo_withoutbnames = Array{CpaInfo,1}(undef,0) 
        for cpainfo in cpas.cpas
            cpahasbname = false
            for bname in boundarynames_arr
                if has_cpa_adj_boundary(cpainfo, bname)
                    cpahasbname =true
                    break
                end
            end
            if !cpahasbname
                curret_cpas_have_bname = true
                push!(cpainfo_withoutbnames, cpainfo)
            end
        end 

        if curret_cpas_have_bname
            push!(cpasarr_without_bnames, Cpas(cpas.time, cpainfo_withoutbnames))
        end
    end

    return cpasarr_without_bnames
end


function printarrayarray(arr)
    println("[")
    if !isnothing(arr)
        for i=1:length(arr)
            print("\t[")

            if !isnothing(arr[i])
                print(arr[i][1])

                for j=2:length(arr[i])
                    print(", " , arr[i][j])
                end
            else
                print(arr[i])
            end

            print("]\n")
        end
    else
        println(arr)
    end
    println("]")
end