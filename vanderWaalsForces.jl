#
#  vanderWaalsForces.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module VanderWaalsforces

include("lennardJones.jl")
using LinearAlgebra
using StaticArrays
using .LennardJones

@inline function vanderWaalsForces!(pos::MMatrix,F::MMatrix,Nmonomers::Int64,Ndomains::Int64,ϵ::Float64,σ::Float64,D::MArray{Tuple{3},Float64,1,3})

    # Loop over all particles
    @inbounds for ii=1:Ndomains*Nmonomers
        # For each particle, loop over all other particles
        @inbounds for jj=1:Ndomains*Nmonomers
            # Exclude particles in the same monomer
        #    if floor(Int8,(ii-1)/Ndomains)!=floor(Int8,(jj-1)/Ndomains)
        #
        #        mod(ii-1,Ndomains)+1) $(mod(jj-1,Ndomains)+1)
            D = pos[jj,:] .- pos[ii,:]
            Dmag_sq = dot(D,D)
            if (ii-1)%Ndomains+1 == (jj-1)%Ndomains+4
                lennardJones!(D,ϵ,σ)
                F[ii,:] .+= D
                F[jj,:] .-= D
            elseif Dmag_sq < σ^2 && ii!=jj
                lennardJones!(D,ϵ,σ)
                F[ii,:] .+= D
                F[jj,:] .-= D
            else
            end

        end
    end
end

export vanderWaalsForces!

end
