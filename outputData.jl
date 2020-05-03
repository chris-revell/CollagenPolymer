#
#  outputData.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module OutputData

using Printf
using DelimitedFiles
using StaticArrays

@inline function outputData(pos::MMatrix,outfile::IOStream,t::Float64,tmax::Float64)

    writedlm(outfile,pos,", ")
    flush(outfile)
    run(`clear`)
    println("Simulating: $t/$tmax")

end

export outputData

end
