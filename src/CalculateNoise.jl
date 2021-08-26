#
#  CalculateNoise.jl
#  CollagenPolymer
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module CalculateNoise

# Import Julia packages
using Random
using Distributions
using LinearAlgebra
using StaticArrays

# Import local program modules


function calculateNoise!(W,nParticles,threadRNG)

    # Loop over all monomers
    for ii=1:nParticles
        W[ii] = SVector{3}(rand(threadRNG,Uniform(-1.0,1.0),3))
        magW = norm(W[ii])
        W[ii] = W[ii]*rand(threadRNG,Normal(0.0,1.0))/magW
    end

    return nothing
end

export calculateNoise!

end
