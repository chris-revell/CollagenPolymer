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
using StaticArrays

@inline function vanderWaalsForces!(pos::MMatrix,v::MMatrix,Nmonomers::Int64,Ndomains::Int64,epsilon::Float64,re::Float64,F::MArray{Tuple{3},Float64,1,3})

    @inbounds for ii=1:Ndomains*Nmonomers
        @inbounds for jj=1:Ndomains*Nmonomers
            if ii!=jj
                F = pos[jj,:] .- pos[ii,:]
                LennardJones.lennardJones!(F,epsilon,re)
                v[ii,:] .-= F
                v[jj,:] .+= F
            end
        end
    end
end

export vanderWaalsForces!

end
