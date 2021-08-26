#
#  OutputData.jl
#  CollagenPolymer
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module OutputData

# Import Julia packages
using DelimitedFiles
using LinearAlgebra

# Import local program modules


@inline @views function outputData(pos,outfile)#,t,tMax)#,nMonomers,nDomains,σ)

    writedlm(outfile,pos,", ")
    flush(outfile)

    # Measure monomer lengths
    # for i=0:nMonomers-1
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
