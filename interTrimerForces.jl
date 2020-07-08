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
using .LennardJones
using .Threads

@inline function interTrimerForces!(pos,F,Ntrimers,Ndomains,ϵ,σ,dx,cellLists,Ng,WCAthresh_sq,intrctnThrshld,nonZeroGrids,Nfilled,boxSize,dxMatrix,r_m)

    @threads for nn=1:Nfilled

		# Inter-trimer forces
        kk,ll,mm = nonZeroGrids[nn]
        for ii=1:cellLists[kk,ll,mm,1]
            celllabel1 = cellLists[kk,ll,mm,ii+1]
            for xx=-1:1
                if (kk+xx)==0 || (kk+xx)==Ng+1
                    #skip
                else
                    for yy=-1:1
                        if (ll+yy)==0 || (ll+yy)==Ng+1
                            #skip
                        else
                            for zz=-1:1
                                if (mm+zz)==0 || (mm+zz)==Ng+1
                                    #skip
                                else
                                    for jj=1:cellLists[kk+xx,ll+yy,mm+zz,1]
                                        celllabel2 = cellLists[kk+xx,ll+yy,mm+zz,jj+1]
                                        if floor(Int8,(celllabel1-1)/Ndomains)==floor(Int8,(celllabel2-1)/Ndomains) && abs(celllabel1-celllabel2)<=1
                                            # Skip adjacent particles in same trimer
                                        else
                                            dx[:,threadid()] .= pos[celllabel2,:] - pos[celllabel1,:]
                                            dxmag_sq = dot(dx[:,threadid()],dx[:,threadid()])
                                            if dxmag_sq > intrctnThrshld^2
                                                # Skip pairs with separation beyond threshold (technically some may exist despite cell list)
                                            else
                                                if (celllabel1+3)%Ndomains == (celllabel2-1)%Ndomains && floor(Int8,(celllabel1-1)/Ndomains)!=floor(Int8,(celllabel2-1)/Ndomains)
                                                    # Apply adhesive van der waals force in stepped fashion between trimers
                                                    lennardJones!(dx[:,threadid()],ϵ,σ)
                                                    F[celllabel1,:,threadid()] .+= dx[:,threadid()]
                                                    F[celllabel2,:,threadid()] .-= dx[:,threadid()]
                                                elseif dxmag_sq < WCAthresh_sq
                                                    # For all other particles, apply WCA potential (truncated repulsive Lennard-Jones)
                                                    lennardJones!(dx[:,threadid()],ϵ,σ)
                                                    F[celllabel1,:,threadid()] .+= dx[:,threadid()]
                                                    F[celllabel2,:,threadid()] .-= dx[:,threadid()]
                                                else
                                                    # Skip any pairs within interaction range, beyond WCA range, and without specified adhesive rule
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end


		# Boundary forces
        for jj=1:3
			if nonZeroGrids[nn][jj]==1
				for kk in 1:cellLists[nonZeroGrids[nn]...,1]
					dxmag = (boxSize/2.0 + pos[cellLists[nonZeroGrids[nn]...,1+kk],jj])
					if dxmag < r_m
						# Use morse potential approximation of Lennard Jones because it is valid over values less than zero. Approximation from http://www.znaturforsch.com/aa/v58a/s58a0615.pdf
						Fmag = (12.0*ϵ/r_m)*(exp(6.0*(1.0-dxmag/r_m))-exp(12.0*(1.0-dxmag/r_m)))
						F[cellLists[nonZeroGrids[nn]...,1+kk],:,threadid()] .-= (Fmag/dxmag).*dxMatrix[jj,:]
					end
				end
			end
			if nonZeroGrids[nn][jj]==Ng
				for kk in 1:cellLists[nonZeroGrids[nn]...,1]
					dxmag = (boxSize/2.0 - pos[cellLists[nonZeroGrids[nn]...,1+kk],jj])
					if dxmag < r_m
						# Use morse potential approximation of Lennard Jones because it is valid over values less than zero. Approximation from http://www.znaturforsch.com/aa/v58a/s58a0615.pdf
						Fmag = (12.0*ϵ/r_m)*(exp(6.0*(1.0-dxmag/r_m))-exp(12.0*(1.0-dxmag/r_m)))
						F[cellLists[nonZeroGrids[nn]...,1+kk],:,threadid()] .+= (Fmag/dxmag).*dxMatrix[jj,:]
					end
				end
			end
		end



    end
end

export interTrimerForces!

end
