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

    @inbounds for ii=1:Ndomains*Nmonomers
        @inbounds for jj=1:Ndomains*Nmonomers
            if floor(Int8,(ii-1)/Ndomains)!=floor(Int8,(jj-1)/Ndomains)
                D = pos[jj,:] .- pos[ii,:]
                lennardJones!(D,ϵ,σ)
                F[ii,:] .+= D
                F[jj,:] .-= D
            end
        end
    end
end

export vanderWaalsForces!

end
