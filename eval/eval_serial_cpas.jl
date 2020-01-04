#=
eval_serial_cpas.jl


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

cpa_file = joinpath(joinpath(joinpath(@__DIR__, "tests"),"test2"),"cpa.txt")

finfilter_alpha_max_min = 0.99
finfilter_num_cells_min = 8

rcpas_over_t = read_filter_cpasrawfile(
                                        cpa_file, 
                                        finfilter_num_cells_min, 
                                        finfilter_alpha_max_min
                                      )
println("EOF ----------------")

# various prints:
# printrawcpasarrays(rcpas_over_t)

cpas_over_t = convert_raw_to_final(rcpas_over_t)
# printcpasarrays(cpas_over_t)

#save to file
printcpastofile(cpas_over_t, "testcpa_saves.txt")

GC.gc()


