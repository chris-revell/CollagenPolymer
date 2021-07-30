#
#  CellListFunctions.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 03/08/2020.
#
#

# Functions from https://discourse.julialang.org/t/cell-list-algorithm-is-slower-than-double-for-loop/36621/2

module CellListFunctions

    using LinearAlgebra
    using Dictionaries
    using StaticArrays

    function find_pairs(nParticles, pos, s, nGrid, neighbourCells)

        pairsList = Tuple{Int64, Int64}[]           # Array of tuples storing neighbour pairs
        boundaryList = Tuple{Int64, Int64, Int64}[]

        # Allocate all particles in matrix pos to grid points
        cellLists = gridAllocate(pos,nParticles,s)

        for key in keys(cellLists)
            # Find neighbour pair list from particles in the same grid cell and neighbouring grid cells
            sameCellPairs!(cellLists[key], pairsList, pos, s)
            loopNeighbour!(key, neighbourCells, cellLists, pairsList, pos, s)

            # Check edge grid cells to find list of particles interacting with boundary
            for jj=1:3
        		if key[jj]==1
        			for kk in cellLists[key]
        				push!(boundaryList,(kk,jj,1))
        			end
        		end
        		if key[jj]==nGrid
        			for kk in cellLists[key]
        				push!(boundaryList,(kk,jj,nGrid))
        			end
        		end
        	end

        end

        return pairsList,boundaryList

    end

    # Allocate all particles in matrix pos to grid points. Return Dictionary cellLists mapping (x,y,z) indices to list of particles.
    function gridAllocate(pos, N, interactionThresh)

        cellLists = Dictionary{SVector{3,Int64},Vector{Int64}}()

        for i = 1:N
            indexSVector = ceil.(Int64,pos[i]/interactionThresh)
            if indexSVector ∉ keys(cellLists)
                insert!(cellLists,indexSVector,[i])
            else
                push!(cellLists[indexSVector],i)
            end
        end

        return cellLists
    end

    # Find list of neighbouring cells around given cell index
    function getForwardCell!(neighbourCells, index)
        neighbourCells[1]  = (index[1]+1, index[2]  , index[3]  )
        neighbourCells[2]  = (index[1]+1, index[2]  , index[3]-1)
        neighbourCells[3]  = (index[1]+1, index[2]  , index[3]+1)
        neighbourCells[4]  = (index[1]+1, index[2]+1, index[3]  )
        neighbourCells[5]  = (index[1]+1, index[2]+1, index[3]-1)
        neighbourCells[6]  = (index[1]+1, index[2]+1, index[3]+1)
        neighbourCells[7]  = (index[1]+1, index[2]-1, index[3]  )
        neighbourCells[8]  = (index[1]+1, index[2]-1, index[3]-1)
        neighbourCells[9]  = (index[1]+1, index[2]-1, index[3]+1)
        neighbourCells[10] = (index[1]  , index[2]  , index[3]+1)
        neighbourCells[11] = (index[1]  , index[2]+1, index[3]  )
        neighbourCells[12] = (index[1]  , index[2]+1, index[3]+1)
        neighbourCells[13] = (index[1]  , index[2]+1, index[3]-1)

        return neighbourCells
    end

    function sameCellPairs!(cellList, pairsList, pos, interactionThresh)
        len = length(cellList)
        if len >= 2
            for i=1:len-1
                for j=i+1:len
                    addOrNot!(cellList[i], cellList[j], pos, interactionThresh, pairsList)
                end
            end
        end
        return pairsList
    end

    #
    function loopNeighbour!(key, neighbourCells, cellLists, pairsList, pos,interactionThresh)
        getForwardCell!(neighbourCells, key)
        for neighbourIndex in neighbourCells
            if neighbourIndex ∈ keys(cellLists)
                findPair!(cellLists[key], cellLists[neighbourIndex], pairsList, pos, interactionThresh)
            end
        end
        return pairsList
    end

    # Function to loop over all particles in two given grid points and then call function  to test whether each particle pair is in interaction range
    function findPair!(cellList1, cellList2, pairsList, pos, interactionThresh)
        for P1 in cellList1
            for P2 in cellList2
                addOrNot!(P1, P2, pos, interactionThresh, pairsList)
            end
        end
        return pairsList
    end

    # Test whether two points P1 and P2 are within separation limit interactionThresh; if so, add tuple (P1,P2) to vector of pairs pairsList
    function addOrNot!(P1, P2, pos, interactionThresh, pairsList)
        dx = pos[P1] - pos[P2]
        dist = norm(dx)
        if dist <= interactionThresh
            push!(pairsList, (P1,P2)) # Add this pair to array of neighbour pairs if within range
        end
        return pairsList
    end

export find_pairs

end
