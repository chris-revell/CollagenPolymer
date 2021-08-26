#
#  CollagenPolymer.jl
#  CollagenPolymer
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module CollagenPolymer

# Import Julia packages
using Distributions
using LinearAlgebra
using Random
using StaticArrays

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
@views function collagenPolymer(nMonomers,L,σ,ϵLJ_in,k_in,Ebend_in,boxSize,tMax,outputFlag,renderFlag)

    # Run parameters
    # nMonomers      (eg. 2    )   Number of collagen monomers
    # L              (eg. 0.5  )   Length of one monomer
    # σ              (eg. 0.05 )   Diameter of monomer = Diameter of particle within monomer/External LJ length scale (separation at which V=0) = 2*particle radius
    # ϵLJ_in         (eg. 00.0 )   External Lennard-Jones energy
    # k_in           (eg. 100.0)   Internal spring stiffness for Fraenkel spring forces between adjacent particles within a monomer
    # Ebend_in       (eg. 0.0  )   Internal bending modulus of monomer
    # boxSize        (eg. 1.0  )   Dimensions of cube in which particles are initialised
    # tMax           (eg. 0.1  )   Total simulation time
    # outputFlag     (eg. 1    )   Controls whether or not data is printed to file
    # renderFlag     (eg. 1    )   Controls whether or not system is visualised with povRay automatically

    # Thermodynamic parameters
    μ              = 1.0                                # Fluid viscosity
    kT             = 1.0                                # Boltzmann constant*Temperature

    # Derived parameters
    ϵLJ            = ϵLJ_in*kT                          # External Lennard-Jones energy
    k              = k_in*kT                            # Internal spring stiffness for forces between adjacent particles within a monomer
    Ebend          = Ebend_in*kT                        # Internal bending modulus of monomer
    outputInterval = tMax/100.0                         # Time interval for writing position data to file
    intrctnThrshld = 2.0*σ                              # Threshold separation for van der Waals interactions
    nDomains       = 1+ceil(Int64,(5.0*L)/(4.0*σ)+1.0)  # Number of particles per monomer, ensuring rₑ<=0.4σ and thus preventing corrugation issues
    rₑ             = (L-σ)/(nDomains-1)                 # Equilibrium separation of springs connecting adjacent particles within a monomer
    nParticles     = nMonomers*nDomains                 # Total number of particles in system
    rₘ             = σ*2.0^(1.0/6.0)                    # Separation at which LJ potential is minimised
    WCAthreshSq    = (rₘ)^2                             # Cutoff threshold for WCA potential - separation at which force is zero - squared
    D              = kT/(6.0*π*μ*σ)                     # Diffusion constant
    nGrid          = ceil(Int64,boxSize/intrctnThrshld) # Dimensions of cellLists array, for boundary interaction.

    # Initialise system time
    t  = 1E-10

    # Allocate system arrays and initialise monomers locations and orientations
    pos,F,W,magsF,magsW,pairsList,boundaryList,dxMatrix,neighbourCells,threadRNG,AA,BB,CC,DD,EE = initialise(nMonomers,nDomains,nParticles,rₑ,boxSize)

    # Setup data folder and output files; output initial state
    if outputFlag == 1
        foldername = createRunDirectory(nMonomers,L,nDomains,μ,kT,ϵLJ,σ,k,rₑ,Ebend,D,tMax,outputInterval,nGrid,boxSize,intrctnThrshld)
        outfile = open("output/$(foldername)/output.txt","w")
        outputData(pos,outfile)#,t,tMax)#,nMonomers,nDomains,σ)
    end

    println("|                   |                   |                   |                   |                  |")
    #Iterate over time until max system time is reached
    while t<tMax

        # Create list of particle interaction pairs based on cell lists algorithm
        pairsList,boundaryList = find_pairs(nParticles,pos,intrctnThrshld,nGrid,neighbourCells)

        # Calculate tension and bending forces within each monomer
        internalForces!(pos,F,nMonomers,nDomains,nParticles,k,rₑ,Ebend,AA,BB,CC,DD,EE)

        # Calculate van der Waals/electrostatic interactions between nearby monomer domains
        interMonomerForces!(pairsList,pos,F,nDomains,ϵLJ,σ,AA,WCAthreshSq,intrctnThrshld)

        # Calculate forces on particles from system boundary
        boundaryForces!(boundaryList,pos,F,ϵLJ,nGrid,boxSize,dxMatrix,rₘ)

        # Find stochastic term (Wiener process) for all monomers
        calculateNoise!(W,nParticles,threadRNG)

        # Adapt timestep to maximum force value
        Δt = adaptTimestep(F,W,magsF,magsW,σ)

        # Integrate system with forward euler
        t = updateSystem!(pos,F,W,t,Δt,D,kT,nParticles)

        if (t%outputInterval)<Δt
            print(">")
            outputFlag == 1 ? outputData(pos,outfile) : nothing
        end

    end
    print("\n")

    if outputFlag == 1
        close(outfile)
        renderFlag == 1 ? run(`python3 src/visualise.py output/$foldername`) : nothing
    end

end

export collagenPolymer

end
