#=
eval_par_cpas.jl


License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
=#
include("./libjl/cpa_eval_structs.jl")
include("./libjl/cpa_eval_functions.jl")
include("./libjl/cpa_print_functions.jl")
include("./libjl/cpa_eval_help_functions.jl")

_dir = joinpath(joinpath(@__DIR__, "tests"),"test1")
# _dir = joinpath(joinpath(@__DIR__, "tests"),"test3") # another test

cpa_parallel_identifier = "_cpa"

prefilter_alpha_max_min = 0.0
finfilter_alpha_max_min = 0.99
prefilter_num_cells_min = 1
finfilter_num_cells_min = 8


dir_str = readdir(_dir)
cpa_idces = findall(x-> occursin(cpa_parallel_identifier, x), dir_str)
cpa_files = dir_str[cpa_idces]

raw_cpas_over_t = Array{CpasRawPrefiltered,1}[]

global pid

for pid=1:1:length(cpa_files)
    println("File: " * string(pid))
    raw_cpas = read_filter_cpasrawfile(
                                                _dir * '/' * cpa_files[pid], 
                                                prefilter_num_cells_min, 
                                                prefilter_alpha_max_min
                                             )

    push!(raw_cpas_over_t, raw_cpas)
    println("EOF ----------------")
end

# various prints:
printrawcpasarrays(raw_cpas_over_t[1]) # per processor [pid]

rec_cpas_over_t = reconstruct_filter_parcpas(
                                                raw_cpas_over_t, 
                                                finfilter_num_cells_min, 
                                                finfilter_alpha_max_min
                                            )

# print reconstructed cpas
# printcpasarrays(rec_cpas_over_t)

#save reconstructed cpas to file
printcpastofile(rec_cpas_over_t, "testcpa_saves.txt")

GC.gc()


