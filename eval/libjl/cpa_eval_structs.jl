#=
cpa_eval_structs.jl


License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
=#

struct CpaParInfo
    no_pcells::Int64
    pcellids::Array{Int64,1} 
end


struct CpaInfoRaw
    no_cells::Int64 
    x::Float64
    y::Float64
    z::Float64
    alpha_mean::Float64
    alpha_max::Float64
    mass::Float64
    vol::Float64
    adj_boundaries::Union{Array{String,1}, Nothing}
    pinfo::Union{CpaParInfo, Nothing}
end


struct CpasRawPrefiltered
   time::Float64 
   cpas::Union{Array{CpaInfoRaw,1}, Nothing}
end


struct CpaInfo
    no_cells::Int64 
    x::Float64
    y::Float64
    z::Float64
    alpha_mean::Float64
    alpha_max::Float64
    mass::Float64
    vol::Float64
    adj_boundaries::Union{Array{String,1}, Nothing}
end


struct Cpas
   time::Float64 
   cpas::Array{CpaInfo,1}
end


function get_cpa_from_rawcpa(cparaw::CpaInfoRaw)
    cpa = CpaInfo(
                    cparaw.no_cells, 
                    cparaw.x, 
                    cparaw.y, 
                    cparaw.z, 
                    cparaw.alpha_mean, 
                    cparaw.alpha_max, 
                    cparaw.mass, 
                    cparaw.vol,
                    cparaw.adj_boundaries
                    )
    return cpa
end


function Base.:+(a::CpaInfo, b::CpaInfo)

    no_cells = a.no_cells + b.no_cells
    mass = a.mass + b.mass
    vol = a.vol + b.vol
    x = (a.x * a.mass + b.x * b.mass) / mass
    y = (a.y * a.mass + b.y * b.mass) / mass
    z = (a.z * a.mass + b.z * b.mass) / mass
    alpha_max = maximum([a.alpha_max b.alpha_max])
    alpha_mean = (a.alpha_mean * a.vol + b.alpha_mean * b.vol) / vol

    adj_boundaries = add_uniques(a.adj_boundaries, b.adj_boundaries)

    c = CpaInfo(no_cells,x,y,z,alpha_mean,alpha_max, mass, vol, adj_boundaries)
    return c
end


function add_uniques(a::Union{Array{String,1}, Nothing}, b::Union{Array{String,1}, Nothing})

    res = a

    if isnothing(a)
        res = b
    else
        if !isnothing(b)
            res = add_unique_str_to_str_list(a,b)
        end
    end

    return res
end



function add_unique_str_to_str_list(a::Array{String,1}, b::Array{String,1})
    res = a

    for str1 in b
        match = false

        for str2 in a
            if str1 == str2
                match = true
                break
            end
        end

        if !match
            push!(res, str1)
        end
    end

    return res
end



function convert_raw_to_final(raw_cpas_over_t::Array{CpasRawPrefiltered,1})
    # mainly for reading of resulting  files

    cpas_over_t = Array{Cpas,1}(undef, length(raw_cpas_over_t))

    for (i,raw_cpas) in enumerate(raw_cpas_over_t)

        time = raw_cpas.time
        cpas =  Array{CpaInfo,1}(undef, length(raw_cpas.cpas))
        for (ii,raw_cpa) in enumerate(raw_cpas.cpas)

            cpas[ii] = CpaInfo(
                                raw_cpa.no_cells, 
                                raw_cpa.x, 
                                raw_cpa.y, 
                                raw_cpa.z, 
                                raw_cpa.alpha_mean, 
                                raw_cpa.alpha_max, 
                                raw_cpa.mass, 
                                raw_cpa.vol,
                                raw_cpa.adj_boundaries
                                )
            if !isnothing(raw_cpa.pinfo)
                println("Warining convert_raw_to_final(): Raw Cpa contains ",
                "parallel information, this function is probably not what you ", "want to use right now ...")
            end
        end

        cpas_over_t[i]=Cpas(time, cpas)
    end

    return cpas_over_t
end