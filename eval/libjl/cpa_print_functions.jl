#=
cpa_print_functions.jl


License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
=#

include("./cpa_eval_structs.jl")

function printrawcpa(r::CpaInfoRaw)
    println("No Cells: \t" * string(r.no_cells))
    println("COM:\t \t" * string(r.x) * ", "
                        * string(r.y) * ", "
                        * string(r.z))
    println("Adj Bound.: \t", r.adj_boundaries)

    if isnothing(r.pinfo)
        
    else
        println("With Parallel Boundaries!");
    end
    println("--")

end


function printrawcpas(rpfs::CpasRawPrefiltered)
    println("->Time: " * string(rpfs.time) * "\n")

    if !isnothing(rpfs.s)
        println("Containing " * string(length(rpfs.s)) * " Cpas:\n--")
        for i=1:1:length(rpfs.s)
            printrawcpa(rpfs.s[i])
        end
    else
        println("Empty")
    end
    println("------\n")
end


function printrawcpasarrays(rpfsarr::Array{CpasRawPrefiltered,1})
    println("Printing S Array, containing " * string(length(rpfsarr))
    * " time points: \n")

    for i=1:1:length(rpfsarr)
        printrawcpas(rpfsarr[i])
    end
end



function printcpainfo(r::CpaInfo)
    println("No Cells: \t" * string(r.no_cells))
    println("COM:\t \t" * string(r.x) * ", "
                        * string(r.y) * ", "
                        * string(r.z))
    println("Mass:\t\t" * string(r.mass) * ", ")
    println("Adj Bound.: \t", r.adj_boundaries)

    println("--")

end


function printcpas(cpa::Cpas)
    println("->Time: " * string(cpa.time) * "\n")

    if !isnothing(cpa.s)
        println("Containing " * string(length(cpa.s)) * " Cpas:\n--")
        for i=1:1:length(cpa.s)
            printcpainfo(cpa.s[i])
        end
    else
        println("Empty")
    end
    println("------\n")
end


function printcpasarrays(cpas_over_t::Array{Cpas,1})
    println("Printing S Array, containing " * string(length(cpas_over_t))
    * " time points: \n")

    for i=1:1:length(cpas_over_t)
        printcpas(cpas_over_t[i])
    end
end


function printcpastofile(cpas_over_t::Array{Cpas,1}, file::String)
#todo
    io = open(file, "w")

    for cpas in cpas_over_t
        write(io, "{\n")
        write(io, string(cpas.time)*"\n")
        write(io, "cell_count[com.(x)_in_m,com.(y),com.(z),mass_in_kg,volume_in_m³,cellAveragedAlpha,maxAlpha,[list_of_boundary_names],[list_of_processor_cells]]\n")
        for cpa in cpas.s
            write(io, string(cpa.no_cells))
            write(io, "[")
            write(io, string(cpa.x)*",")
            write(io, string(cpa.y)*",")
            write(io, string(cpa.z)*",")
            write(io, string(cpa.mass)*",")
            write(io, string(cpa.vol)*",")
            write(io, string(cpa.alpha_mean)*",")
            write(io, string(cpa.alpha_max)*",")
            write(io, "[")
            if !isnothing(cpa.adj_boundaries)
                for (i,b) in enumerate(cpa.adj_boundaries)
                    if i > 1
                        write(io, ",")
                    end
                    write(io, string(b))
                end
            end
            write(io, "],")
            write(io,"[]")
            write(io, "]\n")
        end
        write(io, "}\n")
    end

    close(io)
end


