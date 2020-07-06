#
#  cellLists.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 23/06/2020.
#
#

module CellLists

using LinearAlgebra

@inline function cellLists!(pos,allDomains,cellLists,nonZeroGrids,iₓ,boxSize,intrctnThrshld)

	# Create cell list to identify trimer pairs within interaction range
	# Refresh numbers of trimers in each grid point to zero
	cellLists[:,:,:,1] .= 0
	Nfilled = 0
	# Loop over all particles
	for jj=1:allDomains
		# Find grid point that particle sits within
		iₓ = ceil.(Int64,(pos[jj,:] .+ boxSize/2.0)/intrctnThrshld)
		# Update grid point with this particle
		cellLists[iₓ...,1] += 1
		cellLists[iₓ...,cellLists[iₓ...,1]+1] = jj
		# Store this grid point in nonZeroGrids array if it isn't already there
		if iₓ in nonZeroGrids[1:Nfilled]
			# Skip
		else
			Nfilled+=1
			nonZeroGrids[Nfilled] = iₓ
		end
	end
	# Return the number of grid points that contain particles 
	return Nfilled
end

export cellLists!

end
