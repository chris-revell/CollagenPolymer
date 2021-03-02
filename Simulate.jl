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
using InterMonomerForces
using CalculateNoise
using UpdateSystem
using OutputData
using Initialise
using CreateRunDirectory
using AdaptTimestep
using CellListFunctions
using Visualise

# Define function for bringing together modules to run simulation
@inline @views function simulate()

    # Define run parameters
    nMonomers      = 8        # Number of collagen monomers
    L              = 0.5      # Length of one monomer
    σ              = 0.005    # Diameter of monomer = Diameter of particle within monomer/External LJ length scale (separation at which V=0) = 2*particle radius
    ϵLJ_in         = 00.0     # External Lennard-Jones energy
    k_in           = 1000.0 # Internal spring stiffness for Fraenkel spring forces between adjacent particles within a monomer
    Ebend_in       = 0.0      # Internal bending modulus of monomer
    boxSize        = 1.0      # Dimensions of cube in which particles are initialised
    tMax           = 0.01      # Total simulation time
    outputFlag     = 1        # Controls whether or not data is printed to file
    renderFlag     = 1        # Controls whether or not system is visualised with povRay automatically


    # Thermodynamic parameters
    μ              = 1.0                                # Fluid viscosity
    kT             = 1.0                                # Boltzmann constant*Temperature

    # Force parameters
    ϵLJ            = ϵLJ_in*kT                          # External Lennard-Jones energy
    k              = k_in*kT                            # Internal spring stiffness for forces between adjacent particles within a monomer
    Ebend          = Ebend_in*kT                        # Internal bending modulus of monomer

    # Derived parameters
    outputInterval = tMax/100.0                         # Time interval for writing position data to file
    intrctnThrshld = 2.0*σ                              # Threshold separation for van der Waals interactions
    nDomains       = 1+ceil(Int64,(5.0*L)/(4.0*σ)+1.0)  # Number of particles per monomer, ensuring rₑ<=0.4σ and thus preventing corrugation issues
    rₑ             = (L-σ)/(nDomains-1)                 # Equilibrium separation of springs connecting adjacent particles within a monomer
    nParticles     = nMonomers*nDomains                 # Total number of particles in system
    rₘ             = σ*2.0^(1.0/6.0)                    # Separation at which LJ potential is minimised
    WCAthreshSq    = (rₘ)^2                             # Cutoff threshold for WCA potential - separation at which force is zero - squared
    D              = kT/(6.0*π*μ*σ)                     # Diffusion constant
    nGrid          = ceil(Int64,boxSize/intrctnThrshld) # Dimensions of cellLists array, for boundary interaction.

    # Data arrays
    pos            = SizedArray{Tuple{nParticles,3}}(zeros(nParticles,3))             # xyz positions of all particles
    F              = SizedArray{Tuple{nParticles,3,nthreads()}}(zeros(nParticles,3,nthreads())) # xyz dimensions of all forces applied to particles
    #Fmags          = SizedVector{nParticles}(zeros(Float64,nParticles))                                                          # Vector of force magnitudes for all particules
    W              = SizedArray{Tuple{nParticles,3}}(zeros(nParticles,3))                               # xyz values of stochastic Wiener process for all particles
    pairsList      = Tuple{Int64, Int64}[]                                                               # Array storing tuple of particle interaction pairs eg pairsList[2]=(1,5) => 2nd element of array shows that particles 1 and 5 are in interaction range
    boundaryList   = Tuple{Int64,Int64,Int64}[]                                                          # Array storing list of particles in boundary cells, with 2nd and 3rd components of tuples storing which part of boundary cell is at
    dxMatrix       = SMatrix{3,3}(Matrix(1I, 3, 3))                   # Identity matrix for later calculations
    neighbourCells = MVector{13}(Vector{Tuple{Int64,Int64,Int64}}(undef, 13))     # Vector storing 13 neighbouring cells for a given cell

    # Initialise system time
    t  = 1E-10
    Δt = 0.0

    # Allocate variables to reuse in calculations and prevent memory reallocations
    AA     = SizedArray{Tuple{nthreads(),3}}(zeros(nthreads(),3))
    AA_bar = SizedArray{Tuple{nthreads(),3}}(zeros(nthreads(),3))
    BB     = SizedArray{Tuple{nthreads(),3}}(zeros(nthreads(),3))
    BB_bar = SizedArray{Tuple{nthreads(),3}}(zeros(nthreads(),3))
    CC     = SizedArray{Tuple{nthreads(),3}}(zeros(nthreads(),3))
    DD     = SizedArray{Tuple{nthreads(),3}}(zeros(nthreads(),3))
    DD_bar = SizedArray{Tuple{nthreads(),3}}(zeros(nthreads(),3))
    EE     = SizedArray{Tuple{nthreads(),3}}(zeros(nthreads(),3))
    EE_bar = SizedArray{Tuple{nthreads(),3}}(zeros(nthreads(),3))

    # Create random number generators for each thread
    threadRNG = Vector{Random.MersenneTwister}(undef, nthreads())
    for i in 1:nthreads()
        threadRNG[i] = MersenneTwister()
    end

    # Initialise monomers within boxSize space
    initialise!(pos,nMonomers,nDomains,rₑ,boxSize)

    # Setup data folder and output files; output initial state
    if outputFlag == 1
        foldername = createRunDirectory(nMonomers,L,nDomains,μ,kT,ϵLJ,σ,k,rₑ,Ebend,D,tMax,outputInterval,nGrid,boxSize,intrctnThrshld)
        outfile = open("output/$(foldername)/output.txt","w")
        outputData(pos,outfile,t,tMax,nMonomers,nDomains,σ)
    end

    #Iterate over time until max system time is reached
    while t<tMax

        # Create list of particle interaction pairs based on cell lists algorithm
        pairsList,boundaryList = find_pairs(nParticles,pos,intrctnThrshld,nGrid,neighbourCells)

        # Calculate tension and bending forces within each monomer
        internalForces!(pos,F,nMonomers,nDomains,nParticles,k,rₑ,Ebend,AA,AA_bar,BB,BB_bar,CC,DD,DD_bar,EE,EE_bar)

        # Calculate van der Waals/electrostatic interactions between nearby monomer domains
        interMonomerForces!(pairsList,pos,F,nDomains,ϵLJ,σ,AA,WCAthreshSq,intrctnThrshld)

        # Calculate forces on particles from system boundary
        boundaryForces!(boundaryList,pos,F,ϵLJ,nGrid,boxSize,dxMatrix,rₘ)

        # Find stochastic term (Wiener process) for all monomers
        calculateNoise!(W,nParticles,threadRNG)

        # Adapt timestep to maximum force value
        Δt = adaptTimestep!(F,W,nParticles,σ,D,kT)

        # Integrate system with forward euler
        t = updateSystem!(pos,F,W,t,Δt,D,kT,nParticles)

        if (t%outputInterval)<Δt && outputFlag == 1
            outputData(pos,outfile,t,tMax,nMonomers,nDomains,σ)
        end

    end

    if outputFlag == 1
        close(outfile)
        renderFlag == 1 ? run(`python3 visualise.py output/$foldername`) : nothing
    end

end

export simulate

end
