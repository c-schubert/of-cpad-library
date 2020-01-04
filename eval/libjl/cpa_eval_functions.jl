#=
cpa_eval_functions.jl


License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
=#


include("./cpa_eval_structs.jl")
include("./cpa_eval_help_functions.jl")



function reconstruct_filter_parcpas(
                                    cpas::Array{Array{CpasRawPrefiltered,1},1}, 
                                    filter_no_cells_min, 
                                    filter_alpha_max_min
                                   )
    # nur "upstream" suche
    # z.B. von processor0 -> processor10

    no_procs = length(cpas)
    println("No processors: ", no_procs)
    no_t_steps = 0

    for pid=1:1:no_procs
        if pid == 1
        no_t_steps = length(cpas[pid])
        else
            if length(cpas[pid]) != no_t_steps
                error("Error reconstructParCpas (): All Cpas should have ", 
                " same amount of time points, cannot proceed!")
            end
        end
    end

    rec_cpas_over_t = Array{Cpas, 1}(undef, no_t_steps)

    # loop over  times 
    # all cpas should have equal amount of times (some may contain nothing
    # for some times)
    for tid=1:no_t_steps
        println("Time id - tid: ", tid, "\n")

        haspidcpas = falses(no_procs)
        pid_iscpapar = Array{Union{Array{Bool,1},Nothing},1}(undef, no_procs)
        pid_iscpareconstructed = Array{Union{Array{Bool,1},Nothing},1}(undef, no_procs)
    
        pid_parcpa_ptoids = Array{
                                    Union{
                                        Nothing,  
                                        Array{
                                                Union{Nothing, Array{Union{Nothing, Int64}, 1}}
                                                , 1
                                                }
                                         }
                                    , 1
                                 }(undef,no_procs)
    
        for pid=1:1:no_procs
            if cpas[pid][tid].s != nothing
                haspidcpas[pid] = true
    
                nocpas = length(cpas[pid][tid].s)
                iscpareconstructed = falses(nocpas)
                iscpapar = falses(nocpas)
                cpas_toids_arr = Array{Union{Nothing, Array{Int64, 1}}, 1}(undef,nocpas)
    
                for i=1:nocpas
                    if !isnothing(cpas[pid][tid].s[i].pinfo)
                        # which processor boundaries to look for?
                        iscpapar[i] = true
                        bstr_arr = cpas[pid][tid].s[i].adj_boundaries
                        topids = getboundarypidsandcells(bstr_arr, pid)
                        cpas_toids_arr[i] = topids 
                    else
                        cpas_toids_arr[i] = nothing
                    end
                end
                
                pid_iscpareconstructed[pid] = iscpareconstructed
                pid_iscpapar[pid] = iscpapar
                pid_parcpa_ptoids[pid] = cpas_toids_arr
            else
                pid_iscpareconstructed[pid] = nothing
                pid_iscpapar[pid] = nothing
                pid_parcpa_ptoids[pid] = nothing
            end
        end
    
        # println("-------Info---------")
        # println(length(pid_iscpapar))
        # for pid=1:1:no_procs
        #     println("My pid:", pid)
        #     println("My par bound.  ids:")
        #     println(pid_iscpapar[pid])
        #     println("Parallel to processors:")
        #     println(pid_parcpa_ptoids[pid])
        #     println("--- End pid")
        # end
        # println("------End Info--------\n\n")
    
    
        pid_parcpa_ptocpaids = Array{
                                        Union{
                                            Nothing,  
                                            Array{
                                                    Union{Nothing, Array{Union{Nothing, Int64}, 1}}
                                                    , 1
                                                    }
                                            }
                                        , 1
                                        }(undef,no_procs)
    
        for pfromid=1:1:no_procs
            println("pfromid: ", pfromid)
            ispfromcpapar = pid_iscpapar[pfromid]
            if !isnothing(ispfromcpapar)
                println(length(ispfromcpapar)) 
                ptocpaids =  Array{Union{Nothing, Array{Union{Nothing, Int64}, 1}}, 1}(undef,length(ispfromcpapar)) 
    
                for cpafromid in eachindex(ispfromcpapar)
                    arrid_toids = Array{Union{Nothing, Int64}, 1}(undef,0)     
    
                    if ispfromcpapar[cpafromid]
                        println("cpafromid: ",  cpafromid)
    
                        bcellsfrom = cpas[pfromid][tid].s[cpafromid].pinfo.pcellids
                        ptoids = pid_parcpa_ptoids[pfromid][cpafromid]
                        # possible parallel boundary cpa's ppbcpas
                        for ptoid in ptoids
                            if ptoid > pfromid
                                isptocpapar = pid_iscpapar[ptoid]
                                # going the other way should not be necessary
                                for cpatoid in eachindex(pid_iscpapar[ptoid])
                                    if isptocpapar[cpatoid]
                                        println("ptoid: ", ptoid)
                                        println("cpatoid: ", cpatoid)
                            
                                        if  (
                                            Bool(sum(pid_parcpa_ptoids[ptoid][cpatoid] .== pfromid)) 
                                            && !isnothing(cpas[ptoid][tid].s[cpatoid].pinfo)
                                            )
                                            println("possible match?")
                                            bcellsto = cpas[ptoid][tid].s[cpatoid].pinfo.pcellids
    
                                            if valofvec1_occursin_vec2(abs.(bcellsfrom), abs.(bcellsto))
                                                println("match!")
                                                # append
                                                push!(arrid_toids, cpatoid)
                                            end
                                        end
                                    end
                                end     
                                if !isassigned(arrid_toids,1)
                                    arrid_toids = nothing
                                    println("Warning no assignment of arrid_toids!")
                                end  
                            else
                                push!(arrid_toids, nothing)
                            end
                        end
                    else
                        arrid_toids = nothing
                    end
    
                    ptocpaids[cpafromid] = arrid_toids
                end
    
    
                pid_parcpa_ptocpaids[pfromid] = ptocpaids
            else
                pid_parcpa_ptocpaids[pfromid] = nothing
            end
            println("-------------------")
        end
    
        println("pid_parcpa_ptoids: ", pid_parcpa_ptoids)
            # Set Recontructed S
        println("----------------------------------------------")
        println("----------------------------------------------")

        #######################    
        #RECONSTRUCT 
        ####################### 
        # s_at_t = CpaInfo[]

        rec_cpas = CpaInfo[]
        time = cpas[1][tid].time
        for pid=1:1:no_procs
            if haspidcpas[pid]
                nocpas = length(cpas[pid][tid].s)
                for icpa=1:nocpas
                    if !pid_iscpareconstructed[pid][icpa]
                        if pid_iscpapar[pid][icpa] 
                            #reconstruct
                            rec_cpa = reconstruct_cpa(
                                                        cpas, 
                                                        tid, 
                                                        pid, 
                                                        icpa, 
                                                        pid_parcpa_ptoids, 
                                                        pid_parcpa_ptocpaids, 
                                                        pid_iscpareconstructed
                                                    )
                        else
                            rec_cpa = get_cpa_from_rawcpa(cpas[pid][tid].s[icpa])
                        end

                        if (filter_add_cpa(rec_cpa, filter_alpha_max_min, filter_no_cells_min))
                            push!(rec_cpas, rec_cpa)
                        end
                    end
                end
            end

            if time != cpas[pid][tid].time
                println("Inconsitencys in times of cpas ...")
            end
        end

        rec_cpas_over_t[tid] = Cpas(cpas[1][tid].time, rec_cpas)
    end

    return rec_cpas_over_t
end



function filter_add_cpa(cpa::CpaInfo, alpha_max_min::Float64, no_cells_min::Int64)

    if (cpa.alpha_max >= alpha_max_min && cpa.no_cells >= no_cells_min)
        return true
    else
        return false
    end
end



function reconstruct_cpa(
                        cpas::Array{Array{CpasRawPrefiltered,1},1},
                        tid::Int64, 
                        pid::Int64,
                        i_cpa::Int64,
                        pid_parcpa_ptoids::Array{Union{Nothing, Array{ Union{Nothing, Array{Union{Nothing, Int64}, 1}}, 1}}, 1},
                        pid_parcpa_ptocpaids::Array{Union{Nothing, Array{Union{Nothing, Array{Union{Nothing, Int64}, 1}}, 1}}, 1},
                        pid_iscpareconstructed::Array{Union{Array{Bool,1},Nothing},1}
                        )

    println("Reconstruction of cpa ", i_cpa, " of processor ", pid)
    cpa_from = get_cpa_from_rawcpa(cpas[pid][tid].s[i_cpa])
    cpa_res = cpa_from

    if !(pid_iscpareconstructed[pid][i_cpa])
        ptoids = pid_parcpa_ptoids[pid][i_cpa]
        i_cpato_idces = pid_parcpa_ptocpaids[pid][i_cpa]

        pid_iscpareconstructed[pid][i_cpa] = true

        for i in eachindex(ptoids)
            ptoid = ptoids[i]
            if ptoid > pid
                for i_cpato in i_cpato_idces[i]
                    if (isnothing(ptoid) || isnothing(i_cpato))
                        println("Warning reconstruct_cpa(): Something is not right here...")
                    else

                        if (!pid_iscpareconstructed[ptoid][i_cpato]
                            && !isnothing(pid_parcpa_ptoids[ptoid][i_cpato]) 
                            && (length(pid_parcpa_ptoids[ptoid][i_cpato]) > 1)
                            )
                            #recursion
                            println("recursion ...")
                            cpa_to = reconstruct_cpa(   
                                cpas, 
                                tid, 
                                ptoid, 
                                i_cpato, 
                                pid_parcpa_ptoids, 
                                pid_parcpa_ptocpaids, 
                                pid_iscpareconstructed
                            )
                        else
                            cpa_to = get_cpa_from_rawcpa(cpas[ptoid][tid].s[i_cpato])

                            pid_iscpareconstructed[ptoid][i_cpato] = true
                        end

                        cpa_res = cpa_res + cpa_to
                    end
                end
            end
        end

    end

    return cpa_res
end



function parse_rawinfo(cell_no::Int64, info_str)
    info_str_no_brackets = split(info_str, r"\[(.*?)\]", keepempty=:false)
    
    info_str_brackets = eachmatch( r"\[(.*?)\]", info_str)
    info_str_brackets = collect(info_str_brackets)

    adj_boundaries = info_str_brackets[1].captures
    pinfo = info_str_brackets[2].captures
    
        
    cpa_info_str = split(info_str_no_brackets[1], ",", keepempty=:false)
    
    x = parse(Float64, cpa_info_str[1])
    y = parse(Float64, cpa_info_str[2])
    z = parse(Float64, cpa_info_str[3])
    
    mass = parse(Float64, cpa_info_str[4])
    vol = parse(Float64, cpa_info_str[5])
    alpha_mean = parse(Float64, cpa_info_str[6])
    alpha_max = parse(Float64, cpa_info_str[7])
    
    if adj_boundaries[1] != nothing && adj_boundaries[1] !=""
        adj_boundaries_list = split(adj_boundaries[1], ",", keepempty=:false)
        adj_boundaries_list = string.(adj_boundaries_list)
    else
        adj_boundaries_list = nothing
    end
    
    no_par_cell_info_str = split(info_str_no_brackets[2], ",", keepempty=:false)
    if length(no_par_cell_info_str) == 1
        no_par_cells = parse(Int64, no_par_cell_info_str[1])
    else
        no_par_cells = 0
    end
    
    if no_par_cells > 0 
        pinfo = string.(pinfo)
        pinfo_str_list = split(pinfo[1], ",", keepempty=:false)

        pinfo_vec = parse.(Int64, pinfo_str_list)
        
        parinfo = CpaParInfo(no_par_cells, pinfo_vec)
    else
        parinfo = nothing
    end

    inforaw = CpaInfoRaw(cell_no, x, y, z, alpha_mean,
                      alpha_max, mass, vol, adj_boundaries_list, parinfo)
    
    return inforaw
end

function parse_entry(line_str)
    
    strbuff = IOBuffer(UInt8[], read=true, write=true, append=true)
    istrat_info_str = 0
    
    for c = 1:1:length(line_str)
        
        if line_str[c] == '['
            istrat_info_str = c
            break
        else
            write(strbuff, line_str[c])
        end
    end
    
    cell_no_str =  read(strbuff, String)
    cell_no = parse(Int64, cell_no_str)

    inforaw = parse_rawinfo(cell_no, line_str[istrat_info_str+1:end-1])
    return inforaw
end


function read_filter_cpasrawfile(filename::String, 
                        no_cells_min::Int64, alpha_max_min::Float64)
    io = open(filename, "r")
    
    cpas_at_t = String[]
    iobuff = IOBuffer(UInt8[], read=true, write=true, append=true)

    for line in eachline(io)
        if line[1] =='{'
            flush(iobuff)
        elseif line[1] =='}'
            push!(cpas_at_t, read(iobuff,String))
        else
            write(iobuff, line * "\n")
        end
    end
      
   rawprefiltereds_over_t = CpasRawPrefiltered[]

    for i=1:1:length(cpas_at_t)
        cpa_str_line = split(cpas_at_t[i], "\n") 
        cpa_str_line = cpa_str_line[1:end-1]        
        
        time = parse(Float64, cpa_str_line[1])

        rawprefiltered_at_t = CpaInfoRaw[]
        
        for j=3:1:length(cpa_str_line)
            
            inforaw = parse_entry(cpa_str_line[j])
            
            #prefilter non-par-boundary  entrys
            if inforaw.pinfo == nothing
                # apply prefiltering
                inforaw = prefilter_np_rawinfo(inforaw, no_cells_min, alpha_max_min)
            else
                # leave everything till reconstruction
            end

            if isnothing(inforaw)
                # nothing
            else
                push!(rawprefiltered_at_t, inforaw)
            end
        end
        
        #println("<------- End read S for Timestep -------->")

        if length(rawprefiltered_at_t) > 0
            push!(rawprefiltereds_over_t, CpasRawPrefiltered(time,rawprefiltered_at_t))
        else
            push!(rawprefiltereds_over_t, CpasRawPrefiltered(time,nothing))
        end
    end
    
    return rawprefiltereds_over_t
end


