#
#  calculateNoise.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module CalculateNoise

using Random
using Distributions
using LinearAlgebra
using .Threads

@inline function calculateNoise!(W,Ntrimers,Ndomains,dt,RNG)

    # Loop over all trimers
    @threads for ii=1:Ntrimers*Ndomains
        ranTheta = 2.0*pi*rand(RNG[threadid()],Uniform(0.0,1.0))
        ranPhi = pi*rand(RNG[threadid()],Uniform(0.0,1.0))
        ranR = abs(rand(RNG[threadid()],Normal(0.0,sqrt(dt))))

        W[ii,1,threadid()] += ranR*cos(ranTheta)*sin(ranPhi)
        W[ii,2,threadid()] += ranR*sin(ranTheta)*sin(ranPhi)
        W[ii,3,threadid()] += ranR*cos(ranPhi)
    end
    W[:,:,1] = sum(W,dims=3)
    return nothing
end

export calculateNoise!

end
