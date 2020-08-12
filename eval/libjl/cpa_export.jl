#=
cpa_filter.jl

Various export functions for CpaAndTime datatypes

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#
include("./cpa_structs.jl")

function cpavec_to_csv(cpas::Vector{Cpa}, filename::String)

    io = open(filename, "w")
    write(io, "x;y;z;mass;vol\n")
    for cpa in cpas
        write(io, string(cpa.x) * ";")
        write(io, string(cpa.y) * ";")
        write(io, string(cpa.z) * ";")
        write(io, string(cpa.mass) * ";")
        write(io, string(cpa.vol) * "\n")
    end

    close(io)
end


function cpasandtimevec_to_file(cpasandtimevec::Vector{CpasAndTime}, file::String)
    print("Writing CpasAndTime vector to file ... ")
    io = open(file, "w")

    for cpasandtime in cpasandtimevec
        write(io, "{\n")
        write(io, string(cpasandtime.time)*"\n")
        write(io, "cell_count[com.(x)_in_m,com.(y),com.(z),mass_in_kg,volume_in_m³,cellAveragedAlpha,maxAlpha,[list_of_boundary_names],[list_of_processor_cells]]\n")
        for cpa in cpasandtime.cpas
            write(io, string(cpa.c_count))
            write(io, "[")
            write(io, string(cpa.x)*",")
            write(io, string(cpa.y)*",")
            write(io, string(cpa.z)*",")
            write(io, string(cpa.mass)*",")
            write(io, string(cpa.vol)*",")
            write(io, string(cpa.fmean)*",")
            write(io, string(cpa.fmax)*",")
            write(io, "[")
            if !isnothing(cpa.bnames)
                for (i,b) in enumerate(cpa.bnames)
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
    print("Finished\n")
end


function cpasandtimevec_to_paraview(cpasandtimevec::Vector{CpasAndTime}, folder::String ,filenameprefix::String)

    print("Starting PV-export ... ")
    if !isdir(folder)
        mkdir(folder)
    end

    k = 1

    for (i, cpas) in enumerate(cpasandtimevec)
            file = joinpath(folder, filenameprefix*"-"*replace(string(cpas.time),"." => "_")*".csv")
            cpavec_to_csv(cpas.cpas, file::String)

            k += 1
    end

    print("Finished ", string(k-1)," cpa file(s) written\n")
end
