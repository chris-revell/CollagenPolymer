#
#  Initialise.jl
#  CollagenPolymer
#
#  Created by Christopher Revell on 01/05/2020.
#
#

module Initialise

# Import Julia packages
using Random
using Distributions
using LinearAlgebra
using StaticArrays

# Import local program modules


function initialise(nMonomers,nDomains,nParticles,domainLength,boxSize)

    # Allocate system arrays
    pos            = Array{SVector{3,Float64}}(undef,nParticles)            # xyz positions of all particles
    F              = Array{SVector{3,Float64}}(undef,nParticles)            # xyz dimensions of all forces applied to particles
    fill!(F,@SVector zeros(3))
    magsF          = zeros(Float64,nParticles)                              # Vector of force magnitudes for all particules
    W              = Array{SVector{3,Float64}}(undef,nParticles)            # xyz values of stochastic Wiener process for all particles
    magsW          = zeros(Float64,nParticles)                              # Vector of force magnitudes for all particules
    pairsList      = Tuple{Int64, Int64}[]                                  # Array storing tuple of particle interaction pairs eg pairsList[2]=(1,5) => 2nd element of array shows that particles 1 and 5 are in interaction range
    boundaryList   = Tuple{Int64,Int64,Int64}[]                             # Array storing list of particles in boundary cells, with 2nd and 3rd components of tuples storing which part of boundary cell is at
    dxMatrix       = SMatrix{3,3}(Matrix(1I, 3, 3))                         # Identity matrix for later calculations
    neighbourCells = MVector{13}(Vector{Tuple{Int64,Int64,Int64}}(undef, 13))     # Vector storing 13 neighbouring cells for a given cell

    # Arrays to prevent reallocation
    initialx = zeros(3)
    dx       = zeros(3)

    # Create each monomer
    for ii=0:nMonomers-1
        notFound = true
        # Loop to find random monomer position and orientation such that it will fit within system box
        while notFound
            initialx .= rand(Uniform(-boxSize/2.0,boxSize/2.0),3)
            dx       .= rand(Float64,3).-0.5
            normalize!(dx)
            if false in (-boxSize/2.0 .< (initialx .+ dx.*(nDomains-1)*domainLength) .< boxSize/2.0)
                # Repeat
            else
                notFound = false
            end
        end
        # Once position and orientation is found, initialise all particles within monomer
        for jj=1:nDomains
            pos[ii*nDomains+jj] = SVector{3}(initialx .+ (jj-1)*domainLength.*dx)
        end
    end

    @info "Initialised $nMonomers monomers with $(nDomains*nMonomers) particles"

    # Create random number generators for each thread
    # threadRNG = Vector{Random.MersenneTwister}(undef, nthreads())
    # for i in 1:nthreads()
    #     threadRNG[i] = MersenneTwister()
    # end

    threadRNG = MersenneTwister()

    # Allocate variables to reuse in calculations and prevent memory reallocations
    AA     = @SVector zeros(3)
    BB     = @SVector zeros(3)
    CC     = @SVector zeros(3)
    DD     = @SVector zeros(3)
    EE     = @SVector zeros(3)


    return pos,F,W,magsF,magsW,pairsList,boundaryList,dxMatrix,neighbourCells,threadRNG,AA,BB,CC,DD,EE

end

export initialise

end
