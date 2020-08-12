
#=
cpa_import.jl

Reconstruction functions for RawCpasAndTime to CpasAndTime

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#
include("./cpa_structs.jl")
include("./cpa_parallel_reconstruct_prep.jl")
include("./cpa_tuple_funcs.jl")


function recon_lookup_id(nproc_icpa::Vector{Tuple{Int64,Int64}},A::Tuple{Int64,Int64})
    isequal_arr = tupleArrEntryIsEqualTo(nproc_icpa,A)

    if sum(isequal_arr) == 1
        return findfirst(isequal_arr.== true)
    elseif sum(isequal_arr) > 1
        error("In recon_lookup_id(): Found multiple indices")
    else
        error("In recon_lookup_id(): Index not found")
    end
end


function recon_getConnectedRawCpas(cpaA::RawCpaForParRec,cpaArrB::Vector{RawCpaForParRec})

    con_cpa_arr = Vector{Int64}(undef,0)
    i = 1
    for cpaB in cpaArrB
        if recon_areCpasConnected(cpaA,cpaB)
            push!(con_cpa_arr,i)
        end
        i += 1
    end

    return con_cpa_arr
end


function recon_getDirectConnectionsOfRawCpa(cpaA::RawCpaForParRec,rcpaforparrec_per_proc::Vector{Union{Nothing,Vector{RawCpaForParRec}}})

conn_procs = cpaA.toprocids

direct_conections = Vector{Union{Nothing,Tuple{Int64,Int64}}}(undef, 0)

if !isnothing(conn_procs)
    for proc in conn_procs
        sub_cpa_ids = recon_getConnectedRawCpas(cpaA, rcpaforparrec_per_proc[proc])

        for i in sub_cpa_ids
            push!(direct_conections, (proc, i))
        end
    end
else
    direct_conections = nothing
end

return direct_conections
end


function recon_areCpasConnected(cpaA::RawCpaForParRec, cpaB::RawCpaForParRec)

    for paridferA in cpaA.paridfer
        if !isnothing(paridferA)
            if !isnothing(cpaB.paridfer) && sum(abs(paridferA) .== abs.(cpaB.paridfer)) > 0
                return true
            end
        else
            println("This should not happen ...")
        end
    end

    return false
end


function recon_chainDirectConnections(cpa_direct_connections::Vector{Vector{Tuple{Int64,Int64}}}, nproc_icpa::Vector{Tuple{Int64,Int64}})

    if length(cpa_direct_connections) != length(nproc_icpa)
        error("In recon_chainDirectConnections(): Should have same size!")
    end

    n_conn_cpas = length(nproc_icpa)
    cpa_chained_connections = Vector{Union{Nothing,Vector{Tuple{Int64,Int64}}}}(undef,n_conn_cpas)
    isconnected = falses(n_conn_cpas)

    for j = 1 : 1 : n_conn_cpas
        if ( !isconnected[j] &&
             !( isempty(cpa_direct_connections[j]) || isnothing(cpa_direct_connections[j]) ) )

            pot_followed_conns = cpa_direct_connections[j]
            no_pot_followed_conns = length(pot_followed_conns)
            followed_conns = [nproc_icpa[j]]
            isconnected[j] = true

            i = 1
            while i <= no_pot_followed_conns
                # verfolge verbindung und hänge verbindungen an
                i_new_conns = recon_lookup_id(nproc_icpa,pot_followed_conns[i])
                if !isconnected[i_new_conns]
                    pot_followed_conns = [pot_followed_conns; cpa_direct_connections[i_new_conns]]
                    pot_followed_conns = unique(pot_followed_conns)
                    no_pot_followed_conns = length(pot_followed_conns)
                    push!(followed_conns, nproc_icpa[i_new_conns])
                    isconnected[i_new_conns] = true
                end
                i += 1
            end

            cpa_chained_connections[j] = followed_conns
        else
            cpa_chained_connections[j] = nothing
        end
    end

    return cpa_chained_connections
end


function Base.convert(::Type{Cpa}, rawcpaforparrec::RawCpaForParRec)

    cpa = Cpa(
                rawcpaforparrec.c_count,
                rawcpaforparrec.x,
                rawcpaforparrec.y,
                rawcpaforparrec.z,
                rawcpaforparrec.fmean,
                rawcpaforparrec.fmax,
                rawcpaforparrec.mass,
                rawcpaforparrec.vol,
                rawcpaforparrec.bcount,
                rawcpaforparrec.bnames,
                )

    return cpa
end


function Base.convert(::Type{CpasAndTime}, rawcpaforparrecandtime::RawCpaForParRecAndTime)

    print("Reconstructing cpas for time ", rawcpaforparrecandtime.time, " s ...")
    cpasandtime = CpasAndTime(rawcpaforparrecandtime.time, convert(Vector{Cpa},rawcpaforparrecandtime.cpas_per_proc))

    print("Finished\n")
    return cpasandtime
end


function Base.convert(::Type{Vector{Cpa}}, rcpaforparrec_per_proc::Vector{Union{Nothing,Vector{RawCpaForParRec}}})
    # reconstruct par ...

    nprocs = length(rcpaforparrec_per_proc)
    nproc_icpa = Vector{Tuple{Int64,Int64}}(undef,0)

    for n = 1 : nprocs
        if !isnothing(rcpaforparrec_per_proc[n])
            for i = 1 : length(rcpaforparrec_per_proc[n])
                push!(nproc_icpa, (n, i))
            end
        end
    end

    # filter out non parallel cpas -> toproccount = 0 / procid = 0 etc.
    # and empty cpas (no cpa for time and core - should orrcur sometimes?)
    cpas = Cpa[]

    del_idx = []
    idx = 1
    for n = 1 : nprocs
        proc_cpa_arr = rcpaforparrec_per_proc[n]
        if !isnothing(proc_cpa_arr)
            for  i = 1 : length(proc_cpa_arr)
                cpa_of_proc = proc_cpa_arr[i]

                if ( cpa_of_proc.toproccount == 0 || isnothing(cpa_of_proc.toprocids) )
                    push!(cpas, convert(Cpa, cpa_of_proc))
                    push!(del_idx,idx)
                end
                idx += 1
            end
        end
    end
    deleteat!(nproc_icpa, del_idx)

    #
    cpas_direct_connections = Vector{Vector{Tuple{Int64,Int64}}}(undef,0)

    for (p,c) in nproc_icpa
        if !(isnothing(rcpaforparrec_per_proc[p][c].toprocids) || isnothing(rcpaforparrec_per_proc[p][c].paridfer))
            cpa_direct_connection = recon_getDirectConnectionsOfRawCpa(rcpaforparrec_per_proc[p][c],rcpaforparrec_per_proc)
        else
            error("This should not happen!")
        end
        push!(cpas_direct_connections, cpa_direct_connection)
    end

    # println("Direkt Connections")
    # # debug pring
    # for (i,(p,c)) in enumerate(nproc_icpa)
    #     println("proc: ", p, " cpa: ", c, "\t", cpas_direct_connections[i])
    # end

    # concat connections
    cpas_chained_connections = recon_chainDirectConnections(cpas_direct_connections, nproc_icpa)
    #
    # println("Chained Connections")
    # println(cpas_chained_connections)

    # cpas = []
    tempcpa = Cpa(0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,nothing)
    no_chained = 0
    for chained_con in cpas_chained_connections
        if !isnothing(chained_con)
            resettempcpa = true
            no_chained += 1
            for (proc,cpa_id) in chained_con
                if resettempcpa
                    tempcpa = convert(Cpa,rcpaforparrec_per_proc[proc][cpa_id])
                    resettempcpa = false
                else
                    tempcpa = tempcpa +  convert(Cpa,rcpaforparrec_per_proc[proc][cpa_id])
                end
            end
            push!(cpas, tempcpa)
        end
    end

    return cpas
end
