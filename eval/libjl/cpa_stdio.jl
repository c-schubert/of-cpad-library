#=
cpa_stdio.jl

Various print functions ...

License (MIT):

Copyright (c) 2019 Christian Schubert

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

include("./cpa_structs.jl")

function printcpa(cpa::Cpa)
    println("No Cells: \t" * string(cpa.c_count))
    println("COM:\t \t" * string(cpa.x) * ", "
                        * string(cpa.y) * ", "
                        * string(cpa.z))
    println("Mass:\t\t" * string(cpa.mass) * ", ")
    println("Adj Bound.: \t", cpa.bnames)

    println("--")
end


function printcpasandtime(cpasandtime::CpasAndTime)
    println("->Time: " * string(cpasandtime.time) * "\n")

    if !isnothing(cpasandtime.cpas)
        println("Containing " * string(length(cpasandtime.cpas)) * " Cpas:\n--")
        for cpa in cpasandtime.cpas
            printcpa(cpa)
        end
    else
        println("Empty")
    end
    println("------\n")
end


function printcpasandtimevec(cpasandtimevec::Vector{CpasAndTime})
    println("Printing cpasandtimevec, containing " * string(length(cpasandtimevec))
    * " time points: \n")

    for cpasandtime in cpasandtimevec
        printcpasandtime(cpasandtime)
    end
end



function printarrayarray(arr)
    println("[")
    if !isnothing(arr)
        for i=1:length(arr)
            print("\t[")

            if !isnothing(arr[i])
                print(arr[i][1])

                for j=2:length(arr[i])
                    print(", " , arr[i][j])
                end
            else
                print(arr[i])
            end

            print("]\n")
        end
    else
        println(arr)
    end
    println("]")
end


function printarray(arr)
    if !isnothing(arr)
    print("[")
    for i=1:length(arr)-1
        print(arr[i],", ")
    end
    print(arr[length(arr)],"]\n")
    else
        print(arr)
    end
end
