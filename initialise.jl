#
#  initialise.jl
#  collagen-model
#
#  Created by Christopher Revell on 01/05/2020.
#
#
__precompile__()
module Initialise

using StaticArrays
using Random
using Distributions

@inline function initialise(pos::MMatrix,Ntrimers::Int64,Ndomains::Int64,domainLength::Float64,boxsize::Float64)
    initialx = zeros(3)
    initialAngle = zeros(2)

    for ii=0:Ntrimers-1
        initialx     .= rand(Uniform(0,boxsize),3)
        initialAngle .= rand(Uniform(0.0,π),2) .* [2.0,1.0]
        for jj=1:Ndomains
            pos[ii*Ndomains+jj,1] = initialx[1]+(jj-1)*domainLength*cos(initialAngle[1])*sin(initialAngle[2]);
            pos[ii*Ndomains+jj,2] = initialx[2]+(jj-1)*domainLength*sin(initialAngle[1])*sin(initialAngle[2]);
            pos[ii*Ndomains+jj,3] = initialx[3]+(jj-1)*domainLength*cos(initialAngle[2]);
        end
    end
end

export initialise

end
