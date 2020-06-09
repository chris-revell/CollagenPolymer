#
#  bendingForces.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#

__precompile__()
module BendingForces

using LinearAlgebra
using StaticArrays

@inline function bendingForces!(pos::MMatrix,F::MMatrix,Ntrimers::Int64,Ndomains::Int64,Ebend::Float64,AC::Array{Float64,1},AB::Array{Float64,1},D::Array{Float64,1})

    # Loop over all trimers
    for ii=0:Ntrimers-1
        # Loop over all sets of 3 in each trimer chain
        for jj=1:Ndomains-2

            AC = pos[ii*Ndomains+jj+2,:] .- pos[ii*Ndomains+jj,:]
            AB = pos[ii*Ndomains+jj+1,:] .- pos[ii*Ndomains+jj,:]
            D  = Ebend*(AC/2.0 .- AB)

            F[ii*Ndomains+jj,:]   .-= D/2.0
            F[ii*Ndomains+jj+1,:] .+= D
            F[ii*Ndomains+jj+2,:] .-= D/2.0
        end
    end
end

export bendingForces!

end
