#
#  bendingModulus.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

__precompile__()
module BendingModulus

using LinearAlgebra
using StaticArrays

@inline function bendingModulus!(pos::MMatrix,F::MMatrix,Nmonomers::Int64,Ndomains::Int64,Ebend::Float64,AC::Array{Float64,1},AB::Array{Float64,1},D::Array{Float64,1})

    # Loop over all monomers
    @inbounds for ii=1:Nmonomers
        # Loop over all sets of 3 in each monomer chain
        @inbounds for jj=1:Ndomains-2

            AC = pos[(ii-1)*Ndomains+jj+2,:] .- pos[(ii-1)*Ndomains+jj,:]
            AB = pos[(ii-1)*Ndomains+jj+1,:] .- pos[(ii-1)*Ndomains+jj,:]
            D  = Ebend*(AC/2.0 .- AB)

            F[(ii-1)*Ndomains+jj,:]   .-= D/2.0
            F[(ii-1)*Ndomains+jj+1,:] .+= D
            F[(ii-1)*Ndomains+jj+2,:] .-= D/2.0
        end
    end
end

export bendingModulus!

end
