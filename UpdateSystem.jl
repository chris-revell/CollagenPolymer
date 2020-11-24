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
using Base.Threads

@inline @views function updateSystem!(pos,F,W,t,Δt,D,kT,nParticles)

    pos .+= F[:,:,1].*(Δt*D/kT) .+ W.*sqrt(2.0*D*Δt)
    fill!(F,0.0)

    return t += Δt

end

export updateSystem!

end
