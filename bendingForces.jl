#
#  bendingForces.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module BendingForces

using LinearAlgebra
using .Threads

@inline function bendingForces!(pos,F,Ntrimers,Ndomains,Ebend,AC,AB,D)

    @threads for jj=1:Ndomains*Ntrimers

        if (jj+2)%Ndomains == 0
            #skip
        else
            AC[:,threadid()] = pos[jj+2,:] - pos[jj,:]
            AB[:,threadid()] = pos[jj+1,:] - pos[jj,:]
            D[:,threadid()]  = Ebend*(AC[:,threadid()]/2.0 - AB[:,threadid()])

            F[ii*Ndomains+jj,:,threadid()]   -= D[:,threadid()]./2.0
            F[ii*Ndomains+jj+1,:,threadid()] += D[:,threadid()]
            F[ii*Ndomains+jj+2,:,threadid()] -= D[:,threadid()]./2.0
        end
    end
end

export bendingForces!

end
