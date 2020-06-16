#
#  updateSystem.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module UpdateSystem

using LinearAlgebra
using StaticArrays

@inline function updateSystem!(pos,F,W,t,dt,D,kT)

    pos .+= F.*(dt*D/kT) .+ W.*sqrt(2.0*D)
    F .= 0
    W .= 0

    return t += dt

end

export updateSystem!

end
