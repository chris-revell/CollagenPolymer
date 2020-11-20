#
#  Initialise.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 01/05/2020.
#
#

module Initialise

using Random
using Distributions

@inline function initialise!(pos,nTrimers,nDomains,domainLength,boxSize)

    # Arrays to prevent reallocation
    initialx = zeros(3)
    initialAngle = zeros(2)
    dx = zeros(3)

    # Create each trimer
    for ii=0:nTrimers-1

        notFound = true
        # Loop to find random trimer position and orientation such that it will fit within system box
        while notFound
            initialx     .= rand(Uniform(-boxSize/2.0,boxSize/2.0),3)
            initialAngle .= rand(Uniform(0.0,π),2) .* [2.0,1.0]
            dx .= [cos(initialAngle[1])*sin(initialAngle[2]), sin(initialAngle[1])*sin(initialAngle[2]), cos(initialAngle[2])]
            if false in (-boxSize/2.0 .< (initialx .+ dx.*nDomains*domainLength) .< boxSize/2.0)
                # Repeat
            else
                notFound = false
            end
        end
        # Once position and orientation is found, initialise all particles within trimer
        for jj=1:nDomains
            pos[ii*nDomains+jj] = initialx .+ (jj-1)*domainLength.*dx
        end
    end
    return nothing
end

export initialise!

end
