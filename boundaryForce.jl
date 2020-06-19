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

@inline function boundaryForce!(pos,F,cellLists,nonZeroGrids,Nfilled,Ng)


	# Loop over edges of cell lists grid
	for i in 1:Nfilled
		if 1 in nonZeroGrids[i] || Ng in nonZeroGrids[i]
			kk,ll,mm = nonZeroGrids[i]
			#println(nonZeroGrids[i])
			for j in cellLists[kk,ll,mm,1]

########################################################################################
				# Set direction and magnitude of force properly
				F[j,:] = -10000.0*pos[j,:]
########################################################################################
			end
		end
	end

end

export boundaryForce!

end
