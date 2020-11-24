#
#  BoundarForces.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 04/08/2020.
#
#

module BoundaryForces

using LinearAlgebra

@inline @views function boundaryForces!(boundaryList,pos,F,ϵ,nGrid,boxSize,dxMatrix,rₘ)
	# Boundary forces
	for (particle,dimension,edge) in boundaryList
		if edge == 1
			dxmag = boxSize/2.0 + pos[particle,dimension]
			if dxmag < rₘ
				Fmag = (12.0*ϵ/rₘ)*(exp(6.0*(1.0-dxmag/rₘ))-exp(12.0*(1.0-dxmag/rₘ)))
				F[particle,:,1] .-= (Fmag/dxmag).*dxMatrix[dimension,:]
			end
		end
		if edge == nGrid
			dxmag = boxSize/2.0 - pos[particle,dimension]
			if dxmag < rₘ
				Fmag = (12.0*ϵ/rₘ)*(exp(6.0*(1.0-dxmag/rₘ))-exp(12.0*(1.0-dxmag/rₘ)))
				F[particle,:,1] .+= (Fmag/dxmag).*dxMatrix[dimension,:]
			end
		end
	end
end

export boundaryForces!

end
