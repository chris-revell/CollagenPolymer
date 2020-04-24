#
#  lennardJones.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module LennardJones

using LinearAlgebra
using StaticArrays

function lennardJones!(dx::MArray{Tuple{3},Float64,1,3},epsilon::Float64,re::Float64)

  dx_mag = sqrt(dot(dx,dx))
  F_mag  = -12.0*epsilon*((re^6.0)/(dx_mag^7.0) - (re^12.0)/(dx_mag^13.0))
  dx .= F_mag.*dx./dx_mag

end

export lennardJones!

end
