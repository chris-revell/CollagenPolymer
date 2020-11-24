#
#  CalculateNoise.jl
#  collagen-brownian-polymer
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

@inline @views function calculateNoise!(W,N,threadRNG)

    # Loop over all trimers
    @threads for ii=1:N
        W[ii,:] .= rand(threadRNG[threadid()],Normal(0.0,1.0)).* normalize!(rand(threadRNG[threadid()],Uniform(-1.0,1.0),3))
    end

    return nothing
end

export calculateNoise!

end
