#
#  boundaryForce.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 18/06/2020.
#
#
__precompile__()
module BoundaryForce

include("./lennardJones.jl")
using LinearAlgebra
using StaticArrays
using .LennardJones

@inline function boundaryForce!(pos,F,cellLists,nonZeroGrids,Nfilled,Ng,boxSize,σ,ϵ)

	dxMatrix = Matrix(1I, 3, 3)
	r_m = σ*2.0^(1.0/6.0)

	# Loop over edges of cell lists grid
	for ii in 1:Nfilled
		for jj=1:3
			if nonZeroGrids[ii][jj]==1
				for kk in 1:cellLists[nonZeroGrids[ii]...,1]
					dxmag = (boxSize/2.0 + pos[cellLists[nonZeroGrids[ii]...,1+kk],jj])
					if dxmag < r_m
						# Use morse potential approximation of Lennard Jones because it is valid over values less than zero. Approximation from http://www.znaturforsch.com/aa/v58a/s58a0615.pdf
						Fmag = (12.0*ϵ/r_m)*(exp(6.0*(1.0-dxmag/r_m))-exp(12.0*(1.0-dxmag/r_m)))
						F[cellLists[nonZeroGrids[ii]...,1+kk],:] .-= (Fmag/dxmag).*dxMatrix[jj,:]
					end
				end
			elseif nonZeroGrids[ii][jj]==Ng
				for kk in 1:cellLists[nonZeroGrids[ii]...,1]
					dxmag = (boxSize/2.0 - pos[cellLists[nonZeroGrids[ii]...,1+kk],jj])
					if dxmag < r_m
						# Use morse potential approximation of Lennard Jones because it is valid over values less than zero. Approximation from http://www.znaturforsch.com/aa/v58a/s58a0615.pdf
						Fmag = (12.0*ϵ/r_m)*(exp(6.0*(1.0-dxmag/r_m))-exp(12.0*(1.0-dxmag/r_m)))
						F[cellLists[nonZeroGrids[ii]...,1+kk],:] .+= (Fmag/dxmag).*dxMatrix[jj,:]
					end
				end
			else
				#skip
			end
		end
	end

end

export boundaryForce!

end
