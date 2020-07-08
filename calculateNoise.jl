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

@inline function calculateNoise!(W,Ntrimers,Ndomains,dt,ThreadRNG)

    # Loop over all trimers
    @threads for ii=1:Ntrimers*Ndomains
        ranTheta = 2.0*pi*rand(ThreadRNG[threadid()],Uniform(0.0,1.0))
        ranPhi = pi*rand(ThreadRNG[threadid()],Uniform(0.0,1.0))
        ranR = abs(rand(ThreadRNG[threadid()],Normal(0.0,sqrt(dt))))

        W[ii,1] += ranR*cos(ranTheta)*sin(ranPhi)
        W[ii,2] += ranR*sin(ranTheta)*sin(ranPhi)
        W[ii,3] += ranR*cos(ranPhi)
    end
    return nothing
end

export calculateNoise!

end
