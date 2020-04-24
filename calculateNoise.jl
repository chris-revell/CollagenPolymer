#
#  calculateNoise.jl
#  collagen-model
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

function calculateNoise!(v::MMatrix,Nmonomers,Ndomains,zetaMag)

    Random.seed!()
    # Loop over all monomers
    @inbounds for ii=1:Nmonomers*Ndomains
        ranTheta = 2.0*pi*rand(Uniform(0.0,1.0))
        ranPhi = pi*rand(Uniform(0.0,1.0))
        ranR = abs(rand(Normal(0.0,zetaMag)))

        v[ii,1] = v[ii,1] + ranR*cos(ranTheta)*sin(ranPhi)
        v[ii,2] = v[ii,2] + ranR*sin(ranTheta)*sin(ranPhi)
        v[ii,3] = v[ii,3] + ranR*cos(ranPhi)
    end
end

export calculateNoise!

end
