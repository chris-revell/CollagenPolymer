#
#  outputData.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module OutputData

using DelimitedFiles
using LinearAlgebra

@inline function outputData(pos,outfile,t,tmax,Ntrimers,Ndomains,σ)

    writedlm(outfile,pos,", ")
    flush(outfile)
    run(`clear`)
    println("Simulating: $t/$tmax")

    # Measure trimer lengths
    # for i=0:Ntrimers-1
        # l = 0
        # for j=1:Ndomains-1
            # dx = pos[i*Ndomains+j+1,:]-pos[i*Ndomains+j,:]
            # l += sqrt(dot(dx,dx))
        # end
        # println("$(i+1) => $(l+σ), Ndomains=$Ndomains")
    # end

    return nothing

end

export outputData

end
