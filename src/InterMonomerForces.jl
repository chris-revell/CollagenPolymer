#
#  InterMonomerForces.jl
#  CollagenPolymer
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module InterMonomerForces

# Import Julia packages
using LinearAlgebra
using StaticArrays

# Import local program modules
using LennardJones

@inline @views function interMonomerForces!(pairsList,pos,F,nDomains,ϵ,σ,dx,WCAthreshSq,intrctnThrshld)

    for (ii,jj) in pairsList
        if floor(Int8,(ii-1)/nDomains)==floor(Int8,(jj-1)/nDomains) && abs(ii-jj)<=1
            # Skip adjacent particles in same monomer
        else
            dx = pos[jj] - pos[ii]
            dxMag = norm(dx)
            if dxMag < intrctnThrshld
                dx = dx*lennardJones(dxMag,ϵ,σ)/dxMag
                F[ii] += dx
                F[jj] -= dx
            end


        	# if (ii+3)%nDomains == (jj-1)%nDomains && floor(Int8,(ii-1)/nDomains)!=floor(Int8,(jj-1)/nDomains)
            # 	# Apply adhesive van der waals force in stepped fashion between monomers
            #     dx = dx*lennardJones(dx,ϵ,σ)
            #     F[ii] += dx
            #     F[jj] -= dx
            # elseif dxmag_sq < WCAthreshSq
            #     # For all other particles, apply WCA potential (truncated repulsive Lennard-Jones)
            #     dx = dx*lennardJones(dx,ϵ,σ)
            #     F[ii] += dx
            #     F[jj] -= dx
            # else
            #     # Skip any pairs within interaction range, beyond WCA range, and without specified adhesive rule
            # end
        end
    end

    return nothing

end

export interMonomerForces!

end
