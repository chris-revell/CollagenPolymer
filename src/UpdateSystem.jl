#
#  UpdateSystem.jl
#  CollagenPolymer
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module UpdateSystem

# Import Julia packages
using LinearAlgebra
using StaticArrays

# Import local program modules


@views function updateSystem!(pos,F,W,t,Δt,D,kT,nParticles)

    pos .+= F[:,1]*(Δt*D/kT) .+ W.*sqrt(2.0*D*Δt)
    fill!(F,@SVector zeros(3))

    return t += Δt

end

export updateSystem!

end
