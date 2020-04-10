#
#  calculateNoise.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module CalculateNoise

using Random, Distributions

function calculateNoise!(v,zetaMag,Nmonomers)

    Random.seed!()
    # Loop over all monomers
    for ii=1:Nmonomers
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
