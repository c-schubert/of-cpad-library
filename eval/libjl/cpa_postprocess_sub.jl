#=
cpa_postprocess_sub.jl

Additional functions used for postprocessing of cpaandtime datatypes

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

include("./cpa_structs.jl")
using LinearAlgebra
using Plots
pyplot()

function phi_polar_cpa(cpa::Cpa, r::Float64)
# geht davon aus, dass die Hauptachse so liegt wie sie liegt ...

x = cpa.x
y = cpa.y

if x > 0
    phi = atan(y/x)
elseif x < 0 && y >= 0
    phi = atan(y/x) + pi
elseif x < 0 && y < 0
    phi = atan(y/x) - pi
elseif x == 0 && y > 0
    phi = pi/2
elseif x == 0 && y < 0
    phi = -pi/2
end

return phi

end


function orthogonal_dist_cpa_axis(cpa::Cpa, origin::Vector{Float64}, dirvec::Vector{Float64})

    return orthogonal_dist_point_axis(
                                        cpa.x, cpa.y, cpa.z,
                                        origin[1], origin[2], origin[3],
                                        dirvec[1], dirvec[2], dirvec[3]
                                     )
end


function orthogonal_dist_point_axis(
                                    xp, yp, zp,
                                    originx, originy, originz,
                                    dirvecx, dirvecy, dirvecz
                                    )
    od = norm(
                cross([originx,originy,originz] .-
                [xp,yp,zp],[dirvecx,dirvecy,dirvecz])
            ) ./ norm([dirvecx,dirvecy,dirvecz])

    return od
end


function cpa_dist(a::Cpa, b::Vector{Cpa})

    dist = zeros(length(b))

    for i = 1 : 1 : length(b)
        dist[i] = sqrt( (b[i].x-a.x)^2 + (b[i].y-a.y)^2 + (b[i].z-a.z)^2 )
    end

    return dist
end


function cpa_dist(a::Cpa, b::Cpa)

    dist = 0
    dist = sqrt( (b.x-a.x)^2 + (b.y-a.y)^2 + (b.z-a.z)^2 )

    return dist
end


function cpa_is_mdiff_small(a::Cpa, b::Cpa, mtol_rel::Float64)

    mtol = mtol_rel * max(a.mass,b.mass)
    mdiff = abs( a.mass - b.mass)

    if mdiff < mtol
        return true
    else
        return false
    end

end

function cpa_requiv(a::Cpa)

    requiv = ( (3.0*a.vol)/(4.0*pi) )^(1.0 / 3.0)

    return requiv
end


function v_stokes_droplet(droplet::Cpa)::Float64

    d = 2.0*cpa_requiv(droplet)
    rhof = 2600.0
    rhod = 6950.0
    nu_slag = 0.025
    g = 9.81

    v = (g * (rhof-rhod)* d^2)/(18.0 * nu_slag)

    return v
end



function plot_cpasandtime(dropletcpasandtime::Vector{CpasAndTime})

    n_droplet_times = length(dropletcpasandtime)
    droplet_t_arr = zeros(n_droplet_times)
    droplet_count = zeros(n_droplet_times)

    for i = 1:1:n_droplet_times
      droplet_t_arr[i] =  dropletcpasandtime[i].time
      droplet_count[i] = length(dropletcpasandtime[i].cpas)
    end

    pl1=scatter(droplet_t_arr,droplet_count, xlabel = "Zeit in s", ylabel="Anzahl Tropfen")
    
    return pl1
end
