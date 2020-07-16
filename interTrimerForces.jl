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
using DataStructures
using .LennardJones

@inline function interTrimerForces!(nonEmptyGridPoints,pos,F,Ndomains,ϵ,σ,dx,Ng,WCAthresh_sq,intrctnThrshld,boxSize,dxMatrix,r_m)

	nextindex = zeros(Int64,3)

    for (index,list) in nonEmptyGridPoints

        for ii in list
            for xx=-1:1
                for yy=-1:1
                    for zz=-1:1
						#if haskey(nonEmptyGridPoints,index.+[xx,yy,zz])
							nextindex = index+[xx,yy,zz]
	                        for jj in nonEmptyGridPoints[nextindex]
	                            if floor(Int8,(ii-1)/Ndomains)==floor(Int8,(jj-1)/Ndomains) && abs(ii-jj)<=1
	                                # Skip adjacent particles in same trimer
	                            else
	                                dx .= pos[jj,:] - pos[ii,:]
	                                dxmag_sq = dot(dx,dx)
	                                if dxmag_sq > intrctnThrshld^2
	                                    # Skip pairs with separation beyond threshold (technically some may exist despite cell list)
	                                else
	                                	if (ii+3)%Ndomains == (jj-1)%Ndomains && floor(Int8,(ii-1)/Ndomains)!=floor(Int8,(jj-1)/Ndomains)
	                                    	# Apply adhesive van der waals force in stepped fashion between trimers
	                                        lennardJones!(dx,ϵ,σ)
	                                        F[ii,:] .+= dx
	                                        F[jj,:] .-= dx
	                                    elseif dxmag_sq < WCAthresh_sq
	                                        # For all other particles, apply WCA potential (truncated repulsive Lennard-Jones)
	                                        lennardJones!(dx,ϵ,σ)
	                                        F[ii,:] .+= dx
	                                        F[jj,:] .-= dx
	                                    else
	                                        # Skip any pairs within interaction range, beyond WCA range, and without specified adhesive rule
	                                    end
	                                end
	                            end
	                        end
						#else
							# skip
						#end
                    end
                end
            end
        end


		# Boundary forces
        for jj=1:3
			if index[jj]==1
				for kk in list
					dxmag = (boxSize/2.0 + pos[kk,jj])
					if dxmag < r_m
						# Use morse potential approximation of Lennard Jones because it is valid over values less than zero. Approximation from http://www.znaturforsch.com/aa/v58a/s58a0615.pdf
						Fmag = (12.0*ϵ/r_m)*(exp(6.0*(1.0-dxmag/r_m))-exp(12.0*(1.0-dxmag/r_m)))
						F[kk,:] .-= (Fmag/dxmag).*dxMatrix[jj,:]
					end
				end
			end
			if index[jj]==Ng
				for kk in list
					dxmag = (boxSize/2.0 - pos[kk,jj])
					if dxmag < r_m
						# Use morse potential approximation of Lennard Jones because it is valid over values less than zero. Approximation from http://www.znaturforsch.com/aa/v58a/s58a0615.pdf
						Fmag = (12.0*ϵ/r_m)*(exp(6.0*(1.0-dxmag/r_m))-exp(12.0*(1.0-dxmag/r_m)))
						F[kk,:] .+= (Fmag/dxmag).*dxMatrix[jj,:]
					end
				end
			end
		end



    end
end

export interTrimerForces!

end
