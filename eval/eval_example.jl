#=
eval_example.jl

Example for CPA postprocessing
=#

include("./libjl/cpalib.jl")

# Read (reconstructed) CPA file
cpasfile = joinpath(@__DIR__, "./tests/test4/testcpa_saves.txt")

rawcpasandtime = read_raw_cpas_from_file(cpasfile)

cpasandtime = convert(Vector{CpasAndTime}, rawcpasandtime)

# CpaFilter(Minimal Number Cells in CPA, Minimal Volume Fraction of Cell in CPA)
cpasandtime_fil = filterCpasAndTimeVec(cpasandtime, CpaFilter(24,0.99))

# Get only CPAs not adjecent to boundaries
dropletcpasandtime = get_cpasandtimevec_no_boundaries(cpasandtime_fil)

rotaxes_origin = [0.,0.,0.]
rotaxes_dirvec = [0.,0.,1.]

dropletpostcpasandtime = followDroplets(dropletcpasandtime,0.02)
trackeddropletpostcpasandtime= exractTrackedCpas(dropletpostcpasandtime, 3)

trackeddropletstats = get_tracked_cpa_stats(trackeddropletpostcpasandtime,
 												rotaxes_origin, rotaxes_dirvec)

meantime = [td.mean_occurence_time for td in trackeddropletstats]
mass_final = [td.mass_final for td in trackeddropletstats]
restime = [td.residencetime for td in trackeddropletstats]
orthdis_final = [td.orthdist_final for td in trackeddropletstats]
phi_final = [td.phi_final for td in trackeddropletstats]
tlength = [td.trajectorylength for td in trackeddropletstats]

# plot or do something else with it ...
