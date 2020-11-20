#
#  OutputData.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module OutputData

using DelimitedFiles
using LinearAlgebra

@inline function outputData(pos,outfile,t,tMax,nTrimers,nDomains,σ)

    writedlm(outfile,pos,", ")
    flush(outfile)
    run(`clear`)
    println("Simulating: $t/$tMax")

    # Measure trimer lengths
    # for i=0:nTrimers-1
    #     l = 0
    #     for j=1:nDomains-1
    #         dx = pos[i*nDomains+j+1,:]-pos[i*nDomains+j,:]
    #         l += sqrt(dot(dx,dx))
    #     end
    #     println("$(i+1) => $(l+σ), nDomains=$nDomains")
    # end

    return nothing

end

export outputData

end
