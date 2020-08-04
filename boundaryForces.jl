#
#  boundarForces.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 04/08/2020.
#
#

module BoundaryForces

using LinearAlgebra

@inline function boundaryForces!(boundary_list,pos,F,ϵ,Ng,boxSize,dxMatrix,r_m)
	# Boundary forces
	for (particle,dimension,edge) in boundary_list
		if edge == 1
			dxmag = boxSize/2.0 + pos[particle,dimension]
			if dxmag < r_m
				Fmag = (12.0*ϵ/r_m)*(exp(6.0*(1.0-dxmag/r_m))-exp(12.0*(1.0-dxmag/r_m)))
				F[particle,:] .-= (Fmag/dxmag).*dxMatrix[dimension,:]
			end
		end
		if edge == Ng
			dxmag = boxSize/2.0 - pos[particle,dimension]
			if dxmag < r_m
				Fmag = (12.0*ϵ/r_m)*(exp(6.0*(1.0-dxmag/r_m))-exp(12.0*(1.0-dxmag/r_m)))
				F[particle,:] .+= (Fmag/dxmag).*dxMatrix[dimension,:]
			end
		end
	end
end

export boundaryForces!

end
