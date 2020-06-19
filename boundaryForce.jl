#
#  boundaryForce.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 18/06/2020.
#
#
__precompile__()
module BoundaryForce

using LinearAlgebra
using StaticArrays

@inline function boundaryForce!(pos,F,Ntrimers,Ndomains,ϵ,dx)


	# Loop over edges of cell lists grid


end

export boundaryForce!

end
