#
#  AdaptTimestep.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 06/07/2020.
#
#

module AdaptTimestep

using LinearAlgebra
using Base.Threads

@inline function adaptTimestep!(F,Fmags,nTrimers,nDomains,σ,D,kT)

    F[:,:,1] = sum(F,dims=3)
    @threads for i=1:nParticles
        Fmags[i] = sum(F[i,:,1].*F[i,:,1])
    end
    Fmax_sq = maximum(Fmags)
    Δt = min(σ^2/(32*D),kT*σ/(2.0*D*sqrt(Fmax_sq)))

    return Δt

end

export adaptTimestep!

end
