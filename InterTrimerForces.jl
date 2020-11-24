#
#  InterTrimerForces.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module InterTrimerForces

using LinearAlgebra
using LennardJones
using StaticArrays
using Base.Threads

@inline @views function interTrimerForces!(pairsList,pos,F,nDomains,ϵ,σ,dx,WCAthreshSq,intrctnThrshld)

    @threads for (ii,jj) in pairsList
        if floor(Int8,(ii-1)/nDomains)==floor(Int8,(jj-1)/nDomains) && abs(ii-jj)<=1
            # Skip adjacent particles in same trimer
        else
            dx[:,threadid()] .= pos[jj,:] .- pos[ii,:]
            dxmag_sq = dot(dx,dx)
        	if (ii+3)%nDomains == (jj-1)%nDomains && floor(Int8,(ii-1)/nDomains)!=floor(Int8,(jj-1)/nDomains)
            	# Apply adhesive van der waals force in stepped fashion between trimers
                lennardJones!(dx,ϵ,σ)
                F[ii,:,threadid()] .+= dx[:,threadid()]
                F[jj,:,threadid()] .-= dx[:,threadid()]
            elseif dxmag_sq < WCAthreshSq
                # For all other particles, apply WCA potential (truncated repulsive Lennard-Jones)
                lennardJones!(dx,ϵ,σ)
                F[ii,:,threadid()] .+= dx[:,threadid()]
                F[jj,:,threadid()] .-= dx[:,threadid()]
            else
                # Skip any pairs within interaction range, beyond WCA range, and without specified adhesive rule
            end
        end
    end
end

export interTrimerForces!

end
