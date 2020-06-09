#
#  lennardJones.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module LennardJones

using LinearAlgebra
using StaticArrays

@inline function lennardJones!(dx::MArray{Tuple{3},Float64,1,3},epsilon::Float64,σ::Float64)

  dx_mag = sqrt(dot(dx,dx))
  F_mag  = 24.0*epsilon*((σ^6.0)/(dx_mag^7.0) - 2.0*(σ^12.0)/(dx_mag^13.0))
  dx .= F_mag.*dx./dx_mag

end

export lennardJones!

end
