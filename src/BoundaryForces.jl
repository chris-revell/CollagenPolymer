#
#  BoundarForces.jl
#  CollagenPolymer
#
#  Created by Christopher Revell on 04/08/2020.
#
#

module BoundaryForces

# Import Julia packages
using LinearAlgebra
using StaticArrays

# Import local program modules


@inline @views function boundaryForces!(boundaryList,pos,F,ϵ,nGrid,boxSize,dxMatrix,rₘ)

	# Boundary forces
	for (particle,dimension,edge) in boundaryList
		if edge == 1
			dxMag = boxSize/2.0 + pos[particle][dimension]
			if dxMag < rₘ
				Fmag = (12.0*ϵ/rₘ)*(exp(6.0*(1.0-dxMag/rₘ))-exp(12.0*(1.0-dxMag/rₘ)))
				F[particle] -= (Fmag/dxMag)*dxMatrix[dimension,:]
			end
		end
		if edge == nGrid
			dxMag = boxSize/2.0 - pos[particle][dimension]
			if dxMag < rₘ
				Fmag = (12.0*ϵ/rₘ)*(exp(6.0*(1.0-dxMag/rₘ))-exp(12.0*(1.0-dxMag/rₘ)))
				F[particle] += (Fmag/dxMag)*dxMatrix[dimension,:]
			end
		end
	end

	return nothing

end

export boundaryForces!

end
