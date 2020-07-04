#
#  initialise.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 01/05/2020.
#
#

module Initialise

using Random
using Distributions

@inline function initialise(pos,Ntrimers,Ndomains,domainLength,boxSize)
    initialx = zeros(3)
    initialAngle = zeros(2)
    dx = zeros(3)

    for ii=0:Ntrimers-1

        notFound = true
        # initialx     .= rand(Uniform(-boxSize/2.0,boxSize/2.0),3)
        while notFound
            initialx     .= rand(Uniform(-boxSize/2.0,boxSize/2.0),3)
            initialAngle .= rand(Uniform(0.0,π),2) .* [2.0,1.0]
            dx .= [cos(initialAngle[1])*sin(initialAngle[2]), sin(initialAngle[1])*sin(initialAngle[2]), cos(initialAngle[2])]
            if false in (-boxSize/2.0 .< (initialx .+ dx.*Ndomains*domainLength) .< boxSize/2.0)
                # Repeat
            else
                notFound = false
            end
        end

        for jj=1:Ndomains
            pos[ii*Ndomains+jj,:] = initialx .+ (jj-1)*domainLength.*dx
        end
    end
end

export initialise

end
