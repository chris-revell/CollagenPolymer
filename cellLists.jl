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
	cellLists[:,:,:,1] .= 0
	Nfilled = 0
	for jj=1:allDomains
		iₓ = ceil.(Int64,(pos[jj,:] .+ boxSize/2.0)/intrctnThrshld)
		cellLists[iₓ...,1] += 1
		cellLists[iₓ...,cellLists[iₓ...,1]+1] = jj
		if iₓ in nonZeroGrids[1:Nfilled]
			# Skip
		else
			Nfilled+=1
			nonZeroGrids[Nfilled] = iₓ
		end
	end
	return Nfilled
end

export cellLists!

end
