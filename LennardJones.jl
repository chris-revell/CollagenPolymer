#
#  LennardJones.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module LennardJones

using LinearAlgebra

@inline function lennardJones!(dx,ϵ,σ)

  dx_mag = sqrt(dot(dx,dx))
  F_mag  = 24.0*ϵ*((σ^6.0)/(dx_mag^7.0) - 2.0*(σ^12.0)/(dx_mag^13.0))
  dx .= (F_mag/dx_mag).*dx

  return nothing
end

export lennardJones!

end
