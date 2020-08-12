#=
cpa_import.jl

Import rawcaps from cpa files written by OpenFOAM or AnsysFluent Conitnuos Phase Area Detection Libraries.

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

include("./cpa_structs.jl")

function read_raw_par_cpas_from_files(dirstr, cpa_idfier::String)

    if isdir(dirstr)
        dir_str = readdir(dirstr)

        rawcpasandtime = Vector{Vector{RawCpasAndTime}}(undef,0)

        cpa_idces = findall(x-> occursin(cpa_idfier, x), dir_str)
        cpa_files = dir_str[cpa_idces]
        cpa_file_idces_str = split.(cpa_files,"_")

        cpa_file_idces = map(cpa_file_idces_str) do fidx
         parse(Int64, fidx[1])
        end

        cpa_file_order = sortperm(cpa_file_idces)

        for i in cpa_file_order
            f = joinpath(dirstr,cpa_files[i])
            rawcpas_of_f = read_raw_cpas_from_file(f)
            push!(rawcpasandtime, rawcpas_of_f)
        end

        return rawcpasandtime
    else
        error("In get_cpa_files_ordered(): ", _dir, " is not a valid directory string")
    end
end


function parse_cpalinestr_to_cpa(cell_count::Int64, cpa_str::SubString{String})::RawCpa

    cpa_str_no_brackets = split(cpa_str, r"\[(.*?)\]", keepempty=:false)
    cpa_str_brackets = eachmatch( r"\[(.*?)\]", cpa_str)
    cpa_str_brackets = collect(cpa_str_brackets)
    bnames_str = cpa_str_brackets[1].captures
    pinfo_str = cpa_str_brackets[2].captures
    cpa_info_str = split(cpa_str_no_brackets[1], ",", keepempty=:false)

    x = parse(Float64, cpa_info_str[1])
    y = parse(Float64, cpa_info_str[2])
    z = parse(Float64, cpa_info_str[3])

    mass = parse(Float64, cpa_info_str[4])
    vol = parse(Float64, cpa_info_str[5])
    fmean = parse(Float64, cpa_info_str[6])
    fmax = parse(Float64, cpa_info_str[7])

    bnames_list = Vector{String}(undef,0)

    if bnames_str[1] != nothing && bnames_str[1] !=""
        bnames_list = split(bnames_str[1], ",", keepempty=:false)
        bnames_list = string.(bnames_list)
            bno = length(bnames_list)
    else
        bno = 0
    end

    no_par_cell_info_str = split(cpa_str_no_brackets[2], ",", keepempty=:false)
    if length(no_par_cell_info_str) == 1
        no_par_cells = parse(Int64, no_par_cell_info_str[1])
    else
        no_par_cells = 0
    end

    if no_par_cells > 1
        pinfo = string.(pinfo_str)
        pinfo_str_list = split(pinfo[1], ",", keepempty=:false)
    elseif no_par_cells > 0
        pinfo = string.(pinfo_str)
        pinfo_str_list = [pinfo[1]]
    else
        pinfo_str_list =[""]
    end

    inforaw = RawCpa(cell_count, x, y, z, fmean,
                      fmax, mass, vol, bno, bnames_list, no_par_cells, pinfo_str_list)

    return inforaw
end


function parse_cpa_entry(cpa_line_str::SubString{String})

    strbuff = IOBuffer(UInt8[], read=true, write=true, append=true)
    istart_info_str = 0

    for c = 1:1:length(cpa_line_str)
        if cpa_line_str[c] == '['
            istart_info_str = c
            break
        else
            write(strbuff, cpa_line_str[c])
        end
    end

    cpa_cell_count_str =  read(strbuff, String)
    cpa_cell_count = parse(Int64, cpa_cell_count_str)

    rawcpa = parse_cpalinestr_to_cpa(cpa_cell_count, cpa_line_str[istart_info_str+1:end-1])
    return rawcpa
end


function read_raw_cpas_from_file(filename::String)::Array{RawCpasAndTime,1}

    print("Reading ", filename, " ... ")
    io = open(filename, "r")

    cpas_at_t_str = String[]
    iobuff = IOBuffer(UInt8[], read=true, write=true, append=true)

    for line in eachline(io)
        if line[1] =='{'
            flush(iobuff)
        elseif line[1] =='}'
            push!(cpas_at_t_str, read(iobuff,String))
        else
            write(iobuff, line * "\n")
        end
    end

   rawcpasandtime = RawCpasAndTime[]

    for i=1:1:length(cpas_at_t_str)
        cpa_str_line = split(cpas_at_t_str[i], "\n")
        cpa_str_line = cpa_str_line[1:end-1]

        t_cpas = parse(Float64, cpa_str_line[1])

        rawcpas_at_t = RawCpa[]

        for j=3:1:length(cpa_str_line)
            rawcpa = parse_cpa_entry(cpa_str_line[j])

            if !isnothing(rawcpa)
                push!(rawcpas_at_t, rawcpa)
            end
        end

        if length(rawcpas_at_t) > 0
            push!(rawcpasandtime, RawCpasAndTime(t_cpas,rawcpas_at_t))
        else
            push!(rawcpasandtime, RawCpasAndTime(t_cpas,nothing))
        end
    end

    print("Finished\n")
    return rawcpasandtime
end


function Base.convert(::Type{Cpa},rawcpa::RawCpa)::Cpa

    return Cpa(rawcpa.c_count,
               rawcpa.x,
               rawcpa.y,
               rawcpa.z,
               rawcpa.fmean,
               rawcpa.fmax,
               rawcpa.mass,
               rawcpa.vol,
               rawcpa.bcount,
               rawcpa.bnames)
end


function Base.convert(::Type{CpasAndTime},rawcpasandtime::RawCpasAndTime)::CpasAndTime

#     println("Warning Convert(): this function should only be used for "*
#             "allready reconstructed or single core evaluated cpa files!")

    return CpasAndTime(rawcpasandtime.time, convert.(Cpa, rawcpasandtime.cpas))
end
