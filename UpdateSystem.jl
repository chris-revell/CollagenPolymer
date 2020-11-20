#
#  UpdateSystem.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module UpdateSystem

using LinearAlgebra
using StaticArrays

@inline @views function updateSystem!(pos,F,W,t,Δt,D,kT,nParticles)

    @threads for i=1:nParticles
        pos[i] = SVector{3}(F[:,:,1].*(Δt*D/kT) .+ W[i].*sqrt(2.0*D*Δt))
    end
    F .= 0

    return t += Δt

end

export updateSystem!

end
