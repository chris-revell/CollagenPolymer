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

@inline function vanderWaalsForces!(pos::MMatrix,F::MMatrix,Ntrimers::Int64,Ndomains::Int64,ϵ::Float64,σ::Float64,D::MArray{Tuple{3},Float64,1,3})

    WCAthresh_sq = (2.0^(1.0/6.0)*σ)^2 # Cutoff threshold for WCA potential - separation at which force is zero - squared

    # Loop over all particles
    @inbounds for ii=1:Ndomains*Ntrimers
        # For each particle, loop over all other particles
        @inbounds for jj=1:Ndomains*Ntrimers
            if floor(Int8,(ii-1)/Ndomains)==floor(Int8,(jj-1)/Ndomains)
                # Skip particles in same trimer
            else
                D = pos[jj,:] .- pos[ii,:]
                Dmag_sq = dot(D,D)
                if (ii+3)%Ndomains == (jj-1)%Ndomains
                    # Apply adhesive van der waals force in stepped fashion between trimers
                    lennardJones!(D,ϵ,σ)
                    F[ii,:] .+= D
                    F[jj,:] .-= D
                elseif Dmag_sq < WCAthresh_sq
                    # For all other particles, apply WCA potential (truncated repulsive Lennard-Jones)
                    lennardJones!(D,ϵ/100.0,σ)
                    F[ii,:] .+= D
                    F[jj,:] .-= D
                else
                    # Outside range of WCA so do nothing
                end
            end

        end
    end
end

export vanderWaalsForces!

end
