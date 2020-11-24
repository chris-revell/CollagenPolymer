#
#  AdaptTimestep.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 06/07/2020.
#
#

module AdaptTimestep

using LinearAlgebra
using Base.Threads

@inline @views function adaptTimestep!(F,ξ,nParticles,σ,D,kT)

    F[:,:,1] = sum(F,dims=3)

    Fmax_sq = maximum(sum(F[:,:,1].*F[:,:,1],dims=2))
    ξmax_sq = maximum(sum(ξ.*ξ,dims=2))

    Δt = min(σ^2/(4.0*ξmax_sq),σ/(2.0*sqrt(Fmax_sq)))

    return Δt

end

export adaptTimestep!

end
