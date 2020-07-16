#
#  cellLists.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 23/06/2020.
#
#

module CellLists

using LinearAlgebra
using DataStructures

@inline function cellLists(nonEmptyGridPoints,pos,allDomains,halfBoxSize,intrctnThrshld)

	empty!(nonEmptyGridPoints)
	# Loop over all particles
	for jj=1:allDomains
		# Find grid point that particle sits within
		cellIndex = ceil.(Int64,(pos[jj,:] .+ halfBoxSize)/intrctnThrshld)
		# Update cell list with this particle
		push!(nonEmptyGridPoints[cellIndex],jj)
	end
	# Return the number of grid points that contain particles
	return nonEmptyGridPoints

end

export cellLists

end
