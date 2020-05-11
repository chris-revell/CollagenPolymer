#
#  updateSystem.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module UpdateSystem

using LinearAlgebra
using StaticArrays

@inline function updateSystem!(pos::MMatrix,F::MMatrix,W::MMatrix,Ntrimers::Int64,Ndomains::Int64,t::Float64,dt::Float64,D::Float64,kT::Float64)

    pos .+= F.*(dt*D/kT) .+ W.*sqrt(2.0*D)
    F .= 0
    W .= 0

    return t += dt

end

export updateSystem!

end
