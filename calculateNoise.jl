#
#  CalculateNoise.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module CalculateNoise

using Random
using Distributions
using LinearAlgebra
using StaticArrays
using Base.Threads

@inline function calculateNoise!(W,N,Δt,RNG)

    # Loop over all trimers
    @threads for ii=1:N

        W[ii] = SVector{3}(rand(RNG[threadid()],Normal(0.0,1.0)).* normalize!(rand(RNG[threadid()],Uniform(-1.0,1.0),3)))

    end

    return nothing
end

export calculateNoise!

end
