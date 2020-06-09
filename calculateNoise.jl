#
#  calculateNoise.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module CalculateNoise

using Random
using Distributions
using LinearAlgebra
using StaticArrays

@inline function calculateNoise!(W::MMatrix,Ntrimers::Int64,Ndomains::Int64,dt::Float64)

    # Loop over all trimers
    for ii=1:Ntrimers*Ndomains
        ranTheta = 2.0*pi*rand(Uniform(0.0,1.0))
        ranPhi = pi*rand(Uniform(0.0,1.0))
        ranR = abs(rand(Normal(0.0,sqrt(dt))))

        W[ii,1] += ranR*cos(ranTheta)*sin(ranPhi)
        W[ii,2] += ranR*sin(ranTheta)*sin(ranPhi)
        W[ii,3] += ranR*cos(ranPhi)
    end
end

export calculateNoise!

end
