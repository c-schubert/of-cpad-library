#=
cpa_postprocess_sub.jl

Helper functions for Vector{String} operations used in cpa_structs.jl

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

function add_unique_str_to_str_list(a::Vector{String}, b::Vector{String})
    res = a

    for str1 in b
        match = false

        for str2 in a
            if str1 == str2
                match = true
                break
            end
        end

        if !match
            push!(res, str1)
        end
    end

    return res
end


function add_uniques(a::Union{Vector{String}, Nothing}, b::Union{Vector{String}, Nothing})

    res = a

    if isnothing(a)
        res = b
    else
        if !isnothing(b)
            res = add_unique_str_to_str_list(a,b)
        end
    end

    return res
end
