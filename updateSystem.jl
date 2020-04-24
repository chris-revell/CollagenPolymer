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

function updateSystem!(pos::MMatrix,v::MMatrix,Nmonomers::Int64,Ndomains::Int64,t::Float64,dt::Float64)

    pos .+= v.*dt
    v .= 0

    return t += dt

end

export updateSystem!

end
