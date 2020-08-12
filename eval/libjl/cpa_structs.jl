#=
cpa_structs.jl

Core CPA structs

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

struct RawCpa
    c_count::Int64
    x::Float64
    y::Float64
    z::Float64
    fmean::Float64
    fmax::Float64
    mass::Float64
    vol::Float64
    bcount::Int64
    bnames::Union{Vector{String}, Nothing}
    paridfercount::Int64
    paridferstr::Vector{String}
end

struct RawCpasAndTime
    time::Float64
    cpas::Union{Vector{RawCpa}, Nothing}
end

struct RawCpaForParRec
    c_count::Int64
    x::Float64
    y::Float64
    z::Float64
    fmean::Float64
    fmax::Float64
    mass::Float64
    vol::Float64
    bcount::Int64
    bnames::Union{Vector{String}, Nothing}
    procid::Union{Int64,Nothing}
    toproccount::Union{Int64,Nothing}
    toprocids::Union{Array{Int64,1}, Nothing}
    paridfer_count::Int64
    paridfer::Union{Array{Int64,1},Nothing}
end


struct RawCpaForParRecAndTime
    time::Float64
    cpas_per_proc::Vector{Union{Vector{RawCpaForParRec}, Nothing}}
end


struct Cpa
    c_count::Int64
    x::Float64
    y::Float64
    z::Float64
    fmean::Float64
    fmax::Float64
    mass::Float64
    vol::Float64
    bcount::Int64
    bnames::Union{Vector{String}, Nothing}
end


struct CpaFilter
    c_count_min::Int64
    fmax_min::Float64
end


struct CpasAndTime
   time::Float64
   cpas::Vector{Cpa}
end


struct CpasAndTimePost
    time::Float64
    cpas::Array{Cpa,1}
    id::Array{Int64,1}
end


struct TrackedCpa
    id::Int64
    time::Array{Float64,1}
    cpa::Array{Cpa,1}
    trajectorylength::Float64
    residencetime::Float64
end

struct TrackedCpaStats
	orthdist_final::Float64
	phi_final::Float64
	mass_final::Float64
	mean_occurence_time::Float64
	trajectorylength::Float64
	residencetime::Float64
end


include("./cpa_stringvec_funcs.jl")

function Base.:+(a::Cpa, b::Cpa)::Cpa

    c_count = a.c_count + b.c_count
    mass = a.mass + b.mass
    vol = a.vol + b.vol
    x = (a.x * a.mass + b.x * b.mass) / mass
    y = (a.y * a.mass + b.y * b.mass) / mass
    z = (a.z * a.mass + b.z * b.mass) / mass
    fmax = maximum([a.fmax b.fmax])
    fmean = (a.fmean * a.vol + b.fmean * b.vol) / vol

    bnames = add_uniques(a.bnames, b.bnames)

    if isnothing(bnames)
        bcount = 0
    else
        bcount = length(bnames)
    end

    c = Cpa(c_count,x,y,z,fmean,fmax, mass, vol,bcount, bnames)
    return c
end
