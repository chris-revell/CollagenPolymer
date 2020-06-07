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

    # Measure trimer lengths
    for i=0:Ntrimers-1
        l = 0
        for j=1:Ndomains-1
            dx = pos[i*Ndomains+j+1,:]-pos[i*Ndomains+j,:]
            l += sqrt(dot(dx,dx))
        end
        println("$(i+1) => $l, $(l/σ), re=$re, Ndomains=$Ndomains")
    end

end

export outputData

end
