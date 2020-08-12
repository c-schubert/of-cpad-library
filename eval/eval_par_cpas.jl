#=
eval_par_cpas.jl

example for reconstruction of parallel written CPAs
=#

include("./libjl/cpalib.jl")

_dir = joinpath(joinpath(@__DIR__, "tests"),"test5") # fluent data

cpa_parallel_identifier = "_cpa"

rawcpasandtime = read_raw_par_cpas_from_files(_dir, cpa_parallel_identifier)

rawcpasforparrecandtime_arr = prepare_parrec(rawcpasandtime)

rawcpasandtime = []
GC.gc()

cpasandtime = convert.(CpasAndTime,rawcpasforparrecandtime_arr)

rawcpasforparrecandtime_arr = []
GC.gc()

cpasandtimevec_to_file(cpasandtime,"testfile_big.txt")
cpasandtimevec_to_paraview(cpasandtime, "pv_export","cpa_to_pv")
