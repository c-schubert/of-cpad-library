#=
cpa_postprocess_sub.jl

Functions used for postprocessing of cpaandtime datatypes

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#
include("./cpa_structs.jl")
include("./cpa_postprocess_sub.jl")


function followDroplets(dropletsandtimevec::Vector{CpasAndTime}, Δt_max::Float64)::Vector{CpasAndTimePost}
# filter unqiue droplets over time
# detect droplet belonging to droplet of timestep before

    Δm_tol = 0.05 
    Δt_tol = 1E-6
    id_counter = 0
    cpaspost = Vector{CpasAndTimePost}(undef,0)
    droplet_ids = []
    v_stokes_fac = 1.2

    for (i,dropletsandtime) in enumerate(dropletsandtimevec)
        # println("i: ", i)
        nodroplets  = length(dropletsandtime.cpas)

        if id_counter == 0 && nodroplets > 0
            droplet_ids_temp = [j for j = 1 : 1 : nodroplets]
            id_counter = nodroplets

        elseif id_counter > 0 && i > 1
            droplet_ids_temp = zeros(Int64, nodroplets)
            cpas_before = dropletsandtimevec[i-1].cpas
            cpas_now = dropletsandtime.cpas
            Δt =  dropletsandtime.time - dropletsandtimevec[i-1].time

            for j = 1 : 1 : nodroplets
                rdist = cpa_dist(cpas_now[j] , cpas_before)

                i_minrdist =argmin(rdist)
                minrdist = rdist[i_minrdist]

#                 vdropletmax = v_stokes_fac * abs(v_stokes_droplet(cpas_now[j]))
                vdropletmax = 1.5
                # println("vdroplet ", vdropletmax)
                dist = (vdropletmax*Δt)

                if dist > 2E-2
                    dist = 2E-2
                end

                if ( minrdist < dist &&
                    cpa_is_mdiff_small(cpas_before[i_minrdist], cpas_now[j], Δm_tol) &&
                     Δt <= (Δt_max + Δt_tol) )
                    droplet_ids_temp[j] = droplet_ids[i_minrdist]

#                     println("minrdist: ", minrdist)
#                     println("vdropletmax: ", vdropletmax)
#                     println("Δt: ", Δt)
#                     println("vdropletmax*Δt: ", vdropletmax*Δt)
                else
                    id_counter += 1
                    droplet_ids_temp[j] = id_counter
                end
            end
        end
        droplet_ids = droplet_ids_temp

        push!(cpaspost, CpasAndTimePost(dropletsandtime.time,
                						 dropletsandtime.cpas,
                						 droplet_ids))

    end
    println("Fin give ids")
    return cpaspost
end


function get_droplet_id_max(droplets::Vector{CpasAndTimePost})::Int64
    max = 1

    for i=1:1:length(droplets)
        max = maximum(droplets[i].id)
    end

    return max
end


function get_droplet_id_idces(droplets::Vector{CpasAndTimePost}, idx::Int64)::Array{Int64,1}

    b = zeros(length(droplets))
    for i=1:1:length(droplets)
        tmp = findfirst(isequal(idx),droplets[i].id)
        if !isnothing(tmp)
            b[i] = tmp
        end
    end

    return b
end

function exractTrackedCpas(droplets::Vector{CpasAndTimePost}, minlength::Int64)::Array{TrackedCpa,1}
# minlenght -> minimal occurence in droplet
# function track trajectory of these droplets
# determine redidence time
    tracked_droplets = Array{TrackedCpa,1}(undef,0)

    for id = 1:get_droplet_id_max(droplets)

        tracked_droplet_tmp = Array{Cpa,1}(undef,0)
        tracked_time_tmp = Array{Float64,1}(undef,0)

        id_arr = get_droplet_id_idces(droplets, id)
        lines_arr = findall(!isequal(0), id_arr)
        id_arr_short = id_arr[id_arr .!= 0]

        if length(id_arr_short) > minlength
            for i = 1:length(id_arr_short)
                push!(tracked_droplet_tmp,droplets[lines_arr[i]].cpas[id_arr_short[i]])
                push!(tracked_time_tmp,droplets[lines_arr[i]].time)
            end

            #calc residence time
            residence_time = tracked_time_tmp[end] - tracked_time_tmp[1]
            trajectorylength = 0

            #calc trajectory
            for i = 2:1:length(tracked_droplet_tmp)
                trajectorylength += cpa_dist(tracked_droplet_tmp[i-1], tracked_droplet_tmp[i])
            end

            push!(tracked_droplets, TrackedCpa(id,tracked_time_tmp,tracked_droplet_tmp,trajectorylength,residence_time))
        end
    end

    println("Fin tracing!")
    return tracked_droplets
end
