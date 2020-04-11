#
#  vanderWaalsForces.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module VanderWaalsForces

include("lennardJones.jl")
using LinearAlgebra

function vanderWaalsForces!(pos,v,Nmonomers,epsilon,re)
    for ii=1:(2*Nmonomers)
        for jj=1:(2*Nmonomers)
            if ii!=jj
                F = pos[jj,:] - pos[ii,:]
                LennardJones.lennardJones!(F,epsilon,re)
                v[ii,:] = v[ii,:] - F
                v[jj,:] = v[jj,:] + F
            end
        end
    end
end

export vanderWaalsForces!

end
