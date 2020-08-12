
include("./cpa_structs.jl")
include("./cpa_postprocess_sub.jl")

function get_tracked_cpa_stats(tcpavec::Vector{TrackedCpa}, rot_orig::Vector{Float64},
							rot_dir::Vector{Float64})


	len_tcpavec = length(tcpavec)
	tcpa_statsvec = Vector{TrackedCpaStats}(undef, len_tcpavec)

	for (j,tcpa) in enumerate(tcpavec)
		orthdis_tmp = [orthogonal_dist_cpa_axis(i,rot_orig,rot_dir ) for i in tcpa.cpa]
		phi_tmp = [phi_polar_cpa(i,orthogonal_dist_cpa_axis(i,rot_orig,rot_dir )) for i in tcpa.cpa]
		mean_occurence_time = (tcpa.time[end]+tcpa.time[1])/2
		orthdis_final = orthdis_tmp[end]
		phi_final = phi_tmp[end]
		mass_final = tcpa.cpa[end].mass
		restime = tcpa.residencetime
		trajectlen = tcpa.trajectorylength
		tcpa_statsvec[j] = TrackedCpaStats(orthdis_final,phi_final,mass_final,mean_occurence_time,trajectlen,restime)
	end

	println("Fin TrackedCpa stats!")
	return tcpa_statsvec
end



function get_tracked_cpa_stats(cpa::Union{Cpa,Nothing}, rot_orig::Vector{Float64},
							rot_dir::Vector{Float64})
# simplfication to get orthis and phi for single cpa ...

		orthdis = orthogonal_dist_cpa_axis(cpa,rot_orig,rot_dir)
		phi = phi_polar_cpa(cpa,orthogonal_dist_cpa_axis(cpa,rot_orig,rot_dir))
		mean_occurence_time = 0
		mass = cpa.mass
		trajectlen = 0
        restime = 0
		return TrackedCpaStats(orthdis,phi,mass,mean_occurence_time,trajectlen,restime)
end


function get_tracked_cpa_stats(cpavec::Vector{Union{Cpa,Nothing}}, rot_orig::Vector{Float64},
							rot_dir::Vector{Float64})
# simplfication to get orthis and phi for single cpa ...
    len_cpavec = length(cpavec)
	tcpa_statsvec = Vector{Union{TrackedCpaStats,Nothing}}(nothing, len_cpavec)
    
    for (i,cpa) in enumerate(cpavec)
        if !isnothing(cpa)
            tcpa_statsvec[i] = get_tracked_cpa_stats(cpa, rot_orig, rot_dir)
        end
    end
    
    return tcpa_statsvec
end




