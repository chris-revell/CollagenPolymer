#
#  interTrimerForces.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module InterTrimerForces

include("lennardJones.jl")
using LinearAlgebra
using StaticArrays
using .LennardJones

@inline function interTrimerForces!(pos::MMatrix,F::MMatrix,Ntrimers::Int64,Ndomains::Int64,ϵ::Float64,σ::Float64,D::MArray{Tuple{3},Float64,1,3},cellLists::Array{Int64},Ng::Int64,WCAthresh_sq::Float64,NgThreshold::Float64)

    for kk=2:Ng-1
        for ll=2:Ng-1
            for mm=2:Ng-1
                for ii=1:cellLists[kk,ll,mm,1]
                    celllabel1 = cellLists[kk,ll,mm,ii+1]
                    for xx=-1:1
                        for yy=-1:1
                            for zz=-1:1
                                for jj=1:cellLists[kk+xx,ll+yy,mm+zz,1]

                                    celllabel2 = cellLists[kk+xx,ll+yy,mm+zz,jj+1]

                                    if floor(Int8,(celllabel1-1)/Ndomains)==floor(Int8,(celllabel2-1)/Ndomains) && abs(celllabel1-celllabel2)<=1
                                        # Skip adjacent particles in same trimer
                                    else
                                        D = pos[celllabel2,:] .- pos[celllabel1,:]
                                        Dmag_sq = dot(D,D)
                                        if Dmag_sq > NgThreshold^2
                                            # Skip pairs with separation beyond threshold (technically some may exist despite cell list)
                                        else
                                            if (celllabel1+3)%Ndomains == (celllabel2-1)%Ndomains && floor(Int8,(celllabel1-1)/Ndomains)!=floor(Int8,(celllabel2-1)/Ndomains)
                                                # Apply adhesive van der waals force in stepped fashion between trimers
                                                lennardJones!(D,ϵ,σ)
                                                F[celllabel1,:] .+= D
                                                F[celllabel2,:] .-= D
                                            elseif Dmag_sq < WCAthresh_sq
                                                # For all other particles, apply WCA potential (truncated repulsive Lennard-Jones)
                                                lennardJones!(D,ϵ,σ)
                                                F[celllabel1,:] .+= D
                                                F[celllabel2,:] .-= D
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

end

export interTrimerForces!

end
