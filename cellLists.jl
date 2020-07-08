#
#  cellLists.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 23/06/2020.
#
#

module CellLists

using LinearAlgebra

@inline function cellLists!(pos,allDomains,cellLists,nonZeroGrids,cellIndex,boxSize,intrctnThrshld)

	# Create cell list to identify trimer pairs within interaction range
	# Refresh numbers of trimers in each grid point to zero
	cellLists[:,:,:,1] .= 0
	Nfilled = 0
	# Loop over all particles
	for jj=1:allDomains
		# Find grid point that particle sits within
		cellIndex = ceil.(Int64,(pos[jj,:] .+ boxSize/2.0)/intrctnThrshld)
		# Update grid point with this particle
		cellLists[cellIndex...,1] += 1
		cellLists[cellIndex...,cellLists[cellIndex...,1]+1] = jj
		# Store this grid point in nonZeroGrids array if it isn't already there
		if cellIndex in nonZeroGrids[1:Nfilled]
			# Skip
		else
			Nfilled+=1
			nonZeroGrids[Nfilled] = cellIndex
		end
	end
	# Return the number of grid points that contain particles
	return Nfilled
end

export cellLists!

end
