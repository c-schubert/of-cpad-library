#=
cpa_import.jl

Preparation functions for RawCpasAndTime reconstruction
Converting of RawCpasAndTime to RawCpaForParRec

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

include("./cpa_structs.jl")

function Base.convert(::Type{RawCpaForParRec}, rcpa::RawCpa)

    proc_idsfrom, proc_idsto = getProcIdsFromStringArr(rcpa.bnames)
    if isnothing(proc_idsfrom)
        procid = nothing
    else
        procid = proc_idsfrom[1]
    end

    if isnothing(proc_idsto)
        toproccount = 0
    else
        toproccount = length(proc_idsto)
    end

    bnames = removeProcStrings(rcpa.bnames)
    bcount = length(bnames)

    if rcpa.paridfercount > 0
        paridfer = parse.(Int64,rcpa.paridferstr)
    else
        paridfer = nothing
    end

    rawcpaforparrec =
        RawCpaForParRec(
                            rcpa.c_count,
                            rcpa.x,
                            rcpa.y,
                            rcpa.z,
                            rcpa.fmean,
                            rcpa.fmax,
                            rcpa.mass,
                            rcpa.vol,
                            bcount,
                            bnames,
                            procid,
                            toproccount,
                            proc_idsto,
                            rcpa.paridfercount,
                            paridfer
                        )

    return rawcpaforparrec
end


function isProcBoundaryString(str::String)
    return match(r"procBoundary[0-9]+to[0-9]+", str)
end


function removeProcStrings(strarr::Array{String,1})
    isProcString = isProcBoundaryString.(strarr)
    cleanarr = strarr[isProcString .== nothing]

    return cleanarr
end


function getProcIdsFromStringArr(bnames::Array{String,1})
# function which returns
# 1. boundary processor ids (pids) for processes boundaries adjacent for the
# processor with id mypid
    bstr_matches = isProcBoundaryString.(bnames)
    bstr_matches = bstr_matches[bstr_matches .!= nothing]

    if length(bstr_matches) > 0
        pidsfrom = zeros(Int64, length(bstr_matches))
        pidsto = zeros(Int64, length(bstr_matches))

        for istrmatches = 1:1:length(bstr_matches)
            pfromtostr = bstr_matches[istrmatches].match

            fromstr = match(r"y[0-9]+t", pfromtostr)
            tostr = match(r"o[0-9]+$", pfromtostr)

            fromstr = fromstr.match[2:end-1]
            tostr = tostr.match[2:end]

            fromid = parse(Int64, fromstr) + 1
            toid = parse(Int64, tostr) + 1

            pidsto[istrmatches] = toid
            pidsfrom[istrmatches] = fromid
        end
    else
        pidsfrom = nothing
        pidsto = nothing
    end

    if !isnothing(pidsfrom) && sum(pidsfrom .== pidsfrom[1]) !=  length(pidsfrom)
        error("More than one pidsfrom id. This should not happen!")
    end

    return pidsfrom, pidsto
end


function prepare_parrec( rcpa_at_t_procarr::Vector{Vector{RawCpasAndTime}})::Vector{RawCpaForParRecAndTime}

    print("Preparing Parallel Reconstruction ... ")
    n = length(rcpa_at_t_procarr)

    #check all procs have same time stamps
    tsteps = 0
    tsteps = length(rcpa_at_t_procarr[1])
    for (i,rcpa_at_t) in enumerate(rcpa_at_t_procarr)
        if tsteps != length(rcpa_at_t_procarr[i])
            error("Error prepare_parrec(): Unequal amount of time steps in processors. Maybe cpa file is corrupted?")
        end
    end

    procs = length(rcpa_at_t_procarr)

    rawcpaforparrecandtimearr = Array{RawCpaForParRecAndTime,1}(undef, tsteps)
    t = rcpa_at_t_procarr[1][1].time

    for it = 1 : tsteps
        cpas_per_proc=Array{Union{Array{RawCpaForParRec,1}, Nothing},1}(undef,procs)

        for p = 1 : procs
            if p == 1
                t = rcpa_at_t_procarr[p][it].time
            else
                if t != rcpa_at_t_procarr[p][it].time
                    error("Error prepare_parrec(): Timestamps are not equal. Maybe cpa file is corrupted?")
                end
            end

            if !isnothing(rcpa_at_t_procarr[p][it].cpas)
                cpas_per_proc[p] = convert.(RawCpaForParRec,rcpa_at_t_procarr[p][it].cpas)
            else
                cpas_per_proc[p] = nothing
            end
        end

     rawcpaforparrecandtimearr[it] = RawCpaForParRecAndTime(t, cpas_per_proc)
    end

    print("Finished \n")
    return rawcpaforparrecandtimearr
end
