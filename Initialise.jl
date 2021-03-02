#
#  Initialise.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 01/05/2020.
#
#

module Initialise

using Random
using Distributions
using LinearAlgebra

@inline @views function initialise!(pos,nMonomers,nDomains,domainLength,boxSize)

    # Arrays to prevent reallocation
    initialx = zeros(3)
    dx       = zeros(3)

    # Create each monomer
    for ii=0:nMonomers-1
        notFound = true
        # Loop to find random monomer position and orientation such that it will fit within system box
        while notFound
            initialx .= rand(Uniform(-boxSize/2.0,boxSize/2.0),3)
            dx       .= rand(Float64,3).-0.5
            normalize!(dx)
            if false in (-boxSize/2.0 .< (initialx .+ dx.*(nDomains-1)*domainLength) .< boxSize/2.0)
                # Repeat
            else
                notFound = false
            end
        end
        # Once position and orientation is found, initialise all particles within monomer
        for jj=1:nDomains
            pos[ii*nDomains+jj,:] .= initialx .+ (jj-1)*domainLength.*dx
        end
    end
    return nothing
end

export initialise!

end
