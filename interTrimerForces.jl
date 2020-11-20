#
#  interTrimerForces.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module InterTrimerForces

include("lennardJones.jl")
using LinearAlgebra
using LennardJones
using Base.Threads

@inline function interTrimerForces!(pairs_list,pos,F,Ndomains,ϵ,σ,dx,WCAthresh_sq,intrctnThrshld)

    @threads for (ii,jj) in pairs_list
        if floor(Int8,(ii-1)/Ndomains)==floor(Int8,(jj-1)/Ndomains) && abs(ii-jj)<=1
            # Skip adjacent particles in same trimer
        else
            dx[:,threadid()] .= pos[jj,:] - pos[ii,:]
            dxmag_sq = dot(dx,dx)
        	if (ii+3)%Ndomains == (jj-1)%Ndomains && floor(Int8,(ii-1)/Ndomains)!=floor(Int8,(jj-1)/Ndomains)
            	# Apply adhesive van der waals force in stepped fashion between trimers
                lennardJones!(dx,ϵ,σ)
                F[ii,:,threadid()] .+= dx[:,threadid()]
                F[jj,:,threadid()] .-= dx[:,threadid()]
            elseif dxmag_sq < WCAthresh_sq
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
