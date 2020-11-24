#
#  main.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module Simulate

# Import Julia packages
using Distributions
using LinearAlgebra
using Random
using StaticArrays
using Base.Threads

# Import local program modules
using InternalForces
using BoundaryForces
using InterTrimerForces
using CalculateNoise
using UpdateSystem
using OutputData
using Initialise
using CreateRunDirectory
using AdaptTimestep
using CellListFunctions
using Visualise

# Define function for bringing together modules to run simulation
@inline @views function simulate(nTrimers::Int64,L::Float64,σ::Float64,ϵLJ_in::Float64,k_in::Float64,Ebend_in::Float64,boxSize::Float64,tMax::Float64,outputFlag::Int64,renderFlag::Int64)

    # Thermodynamic parameters
    μ              = 1.0                                # Fluid viscosity
    kT             = 1.0                                # Boltzmann constant*Temperature

    # Force parameters
    ϵLJ            = ϵLJ_in*kT                          # External Lennard-Jones energy
    k              = k_in*kT                            # Internal spring stiffness for forces between adjacent particles within a trimer
    Ebend          = Ebend_in*kT                        # Internal bending modulus of trimer

    # Derived parameters
    outputInterval = tMax/100.0                         # Time interval for writing position data to file
    intrctnThrshld = 2.0*σ                              # Threshold separation for van der Waals interactions
    nDomains       = 1+ceil(Int64,(5.0*L)/(4.0*σ)+1.0)  # Number of particles per trimer, ensuring rₑ<=0.4σ and thus preventing corrugation issues
    rₑ             = (L-σ)/(nDomains-1)                 # Equilibrium separation of springs connecting adjacent particles within a trimer
    nParticles     = nTrimers*nDomains                  # Total number of particles in system
    rₘ             = σ*2.0^(1.0/6.0)                    # Separation at which LJ potential is minimised
    WCAthreshSq    = (rₘ)^2                             # Cutoff threshold for WCA potential - separation at which force is zero - squared
    D              = kT/(6.0*π*μ*σ)                     # Diffusion constant
    nGrid          = ceil(Int64,boxSize/intrctnThrshld) # Dimensions of cellLists array

    # Data arrays
    pos            = SizedVector{nParticles}(fill(SVector{3}(zeros(Float64,3)),nParticles))             # xyz positions of all particles
    F              = SizedArray{Tuple{nParticles,3,nthreads()}}(zeros(Float64,nParticles,3,nthreads())) # xyz dimensions of all forces applied to particles
    Fmags          = SizedVector{nParticles}(zeros(Float64,nParticles))                                                          # Vector of force magnitudes for all particules
    W              = SizedVector{nParticles}(fill(SVector{3}(zeros(Float64,3)),nParticles))                               # xyz values of stochastic Wiener process for all particles
    pairsList      = Tuple{Int64, Int64}[]                                                               # Array storing tuple of particle interaction pairs eg pairsList[2]=(1,5) => 2nd element of array shows that particles 1 and 5 are in interaction range
    boundaryList   = Tuple{Int64,Int64,Int64}[]                                                          # Array storing list of particles in boundary cells, with 2nd and 3rd components of tuples storing which part of boundary cell is at
    dxMatrix       = SMatrix{3,3}(Matrix(1I, 3, 3))                   # Identity matrix for later calculations
    neighbourCells = MVector{13}(Vector{Tuple{Int64,Int64,Int64}}(undef, 13))     # Vector storing 13 neighbouring cells for a given cell

    # Initialise system time
    t  = 0.0
    Δt = 0.0

    # Allocate variables to reuse in calculations and prevent memory reallocations
    AA = zeros(Float64,3,nthreads())
    AA_bar = zeros(Float64,3,nthreads())
    BB = zeros(Float64,3,nthreads())
    BB_bar = zeros(Float64,3,nthreads())
    CC = zeros(Float64,3,nthreads())
    DD = zeros(Float64,3,nthreads())
    DD_bar = zeros(Float64,3,nthreads())
    EE = zeros(Float64,3,nthreads())
    EE_bar = zeros(Float64,3,nthreads())

    # Create random number generators for each thread
    threadRNG = Vector{Random.MersenneTwister}(undef, nthreads())
    for i in 1:nthreads()
        threadRNG[i] = MersenneTwister()
    end

    # Initialise trimers within boxSize space
    initialise!(pos,nTrimers,nDomains,rₑ,boxSize)

    # Setup data folder and output files; output initial state
    if outputFlag == 1
        foldername = createRunDirectory(nTrimers,L,nDomains,μ,kT,ϵLJ,σ,k,rₑ,Ebend,D,tMax,outputInterval,nGrid,boxSize)
        outfile = open("output/$(foldername)/output.txt","w")
        outputData(pos,outfile,t,tMax,nTrimers,nDomains,σ)
    end

    #Iterate over time until max system time is reached
    while t<tMax

        # Create list of particle interaction pairs based on cell lists algorithm
        pairsList,boundaryList = find_pairs(nParticles,pos,intrctnThrshld,nGrid,neighbourCells)

        # Calculate tension and bending forces within each trimer
        internalForces!(pos,F,nTrimers,nDomains,nParticles,k,rₑ,Ebend,AA,AA_bar,BB,BB_bar,CC,DD,DD_bar,EE,EE_bar)

        # Calculate van der Waals/electrostatic interactions between nearby trimer domains
        interTrimerForces!(pairsList,pos,F,nDomains,ϵLJ,σ,AA,WCAthreshSq,intrctnThrshld)

        # Calculate forces on particles from system boundary
        boundaryForces!(boundaryList,pos,F,ϵLJ,nGrid,boxSize,dxMatrix,rₘ)

        # Find stochastic term (Wiener process) for all monomers
        calculateNoise!(W,nParticles,threadRNG)

        # Adapt timestep to maximum force value
        Δt = adaptTimestep!(F,Fmags,W,nParticles,σ,D,kT)

        # Integrate system with forward euler
        t = updateSystem!(pos,F,W,t,Δt,D,kT,nParticles)

        if (t%outputInterval)<Δt && outputFlag == 1
            outputData(pos,outfile,t,tMax,nTrimers,nDomains,σ)
        end

    end

    if outputFlag == 1
        close(outfile)
        renderFlag == 1 ? visualise("output/"*foldername,nParticles,σ) : nothing
    end

end

export simulate

end
