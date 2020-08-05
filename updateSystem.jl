#
#  updateSystem.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module UpdateSystem

using LinearAlgebra

@inline function updateSystem!(pos,F,W,t,dt,D,kT)

    pos .+= F[:,:,1].*(dt*D/kT) .+ W[:,:,1].*sqrt(2.0*D)
    F .= 0
    W .= 0

    return t += dt

end

export updateSystem!

end
