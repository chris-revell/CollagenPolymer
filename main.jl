#
#  main.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#



# Import Julia packages
using Distributions
using StaticArrays
using LinearAlgebra
# Import program modules
include("./outputData.jl")
include("./calculateNoise.jl")
include("./updateSystem.jl")
include("./tensionForces.jl")
include("./interTrimerForces.jl")
include("./bendingForces.jl")
include("./initialise.jl")
include("./createRunDirectory.jl")
include("./boundaryForce.jl")
include("./cellLists.jl")
using .TensionForces
using .BendingForces
using .InterTrimerForces
using .CalculateNoise
using .UpdateSystem
using .OutputData
using .Initialise
using .CreateRunDirectory
using .BoundaryForce
using .CellLists


# Define run parameters
const Ntrimers       = 5                  # Number of collagen trimers
const L              = 0.5                # Length of one trimer
const σ              = 0.05               # Diameter of trimer = Diameter of particle within trimer/External LJ length scale (separation at which V=0) = 2*particle radius
const ϵLJ_in         = 1.0                # External Lennard-Jones energy
const k_in           = 100000.0           # Internal spring stiffness for forces between adjacent particles within a trimer
const Ebend_in       = 0.0                # Internal bending modulus of trimer
const boxSize        = 1.0                # Dimensions of cube in which particles are initialised
const tmax           = 0.5                # Total simulation time
const renderFlag     = 1                  # Controls whether or not system is visualised with povRay automatically


#%%

# Define function for bringing together modules to run simulation
@inline function main(Ntrimers::Int64,L::Float64,σ::Float64,ϵLJ_in::Float64,k_in::Float64,Ebend_in::Float64,boxSize::Float64,tmax::Float64,renderFlag::Int64)

    # Thermodynamic parameters
    μ              = 1.0                # Fluid viscosity
    kT             = 1.0                # Boltzmann constant*Temperature

    # Force parameters
    ϵLJ            = ϵLJ_in*kT          # External Lennard-Jones energy
    k              = k_in*kT            # Internal spring stiffness for forces between adjacent particles within a trimer
    Ebend          = Ebend_in*kT        # Internal bending modulus of trimer

    # Derived parameters
    outputInterval = tmax/100.0                        # Time interval for writing position data to file
    intrctnThrshld = 2.0*σ                             # Threshold separation for van der Waals interactions
    Ndomains       = 1+ceil(Int64,(5.0*L)/(4.0*σ)+1.0) # Number of particles per trimer, ensuring re<=0.4σ and thus preventing corrugation issues
    re             = (L-σ)/(Ndomains-1)                # Equilibrium separation of springs connecting adjacent particles within a trimer
    allDomains     = Ntrimers*Ndomains                 # Total number of particles in system
    r_m            = σ*2.0^(1.0/6.0)                   # Separation at which LJ potential is minimised
    WCAthresh_sq   = (r_m)^2                           # Cutoff threshold for WCA potential - separation at which force is zero - squared
    D              = kT/(6.0*π*μ*σ)                    # Diffusion constant
    Ng             = ceil(Int64,boxSize/intrctnThrshld)# Dimensions of cellLists array

    # Data arrays
    pos            = zeros(Float64,Ntrimers*Ndomains,3)     # xyz positions of all particles
    F              = zeros(Float64,Ndomains*Ntrimers,3)     # xyz dimensions of all forces applied to particles
    W              = zeros(Float64,Ndomains*Ntrimers,3)     # xyz values of stochastic Wiener process for all particles
    cellLists      = zeros(Int64,Ng,Ng,Ng,10)               # Cell list grid
    nonZeroGrids   = fill(zeros(Int64,3), Ndomains*Ntrimers)# Stores locations of non-empty grid points in cellLists
    Nfilled        = 0                                      # Number of non-empty grid points in cellLists
    dxMatrix       = Matrix(1I, 3, 3)                       # Identity matrix

    # Initialise system time
    t = 0.0
    dt = 0.0

    # Allocate variables to reuse in calculations and prevent memory reallocations
    AA,BB,CC = fill(zeros(Float64,3),3)
    DD = zeros(Int64,3)

    # Setup data folder and output files
    foldername = createRunDirectory(Ntrimers,L,Ndomains,μ,kT,ϵLJ,σ,k,re,Ebend,D,tmax,outputInterval,Ng,boxSize)
    outfile = open("output/$(foldername)/output.txt","w")

    # Initialise trimers within boxSize space
    initialise(pos,Ntrimers,Ndomains,re,boxSize)

    # Output initial state
    outputData(pos,outfile,t,tmax,Ntrimers,Ndomains,σ)

    #Iterate over time until max system time is reached
    while t<tmax

        # Create cell lists array for interactions
        Nfilled = cellLists!(pos,allDomains,cellLists,nonZeroGrids,DD,boxSize,intrctnThrshld)

        # Spring forces between particles along trimer chain
        tensionForces!(pos,F,Ntrimers,Ndomains,k,re,AA)

        # Calculate forces from trimer bending stiffness
        bendingForces!(pos,F,Ntrimers,Ndomains,Ebend,AA,BB,CC)

        # Calculate van der Waals/electrostatic interactions between nearby trimer domains
        interTrimerForces!(pos,F,Ntrimers,Ndomains,ϵLJ,σ,AA,cellLists,Ng,WCAthresh_sq,intrctnThrshld,nonZeroGrids,Nfilled)

        boundaryForce!(pos,F,cellLists,nonZeroGrids,Nfilled,Ng,boxSize,σ,ϵLJ,dxMatrix,r_m)

        # Adapt timestep to maximum force value
        Fmax_sq = max(sum(F.*F,dims=2)...)
        dt = min(σ^2/(32*D),kT*σ/(2.0*D*sqrt(Fmax_sq)))

        # Find stochastic term (Wiener process) for all monomers
        calculateNoise!(W,Ntrimers,Ndomains,dt)

        # Integrate system with forward euler
        t = updateSystem!(pos,F,W,t,dt,D,kT)

        if (t%outputInterval)<dt
            outputData(pos,outfile,t,tmax,Ntrimers,Ndomains,σ)
        end

    end
    close(outfile)

    # Render images from output files
    if renderFlag == 1
        run(`python3 visualise.py output/$foldername`)
    end
end

#%%

main(Ntrimers,L,σ,ϵLJ_in,k_in,Ebend_in,boxSize,tmax,renderFlag)

# using Profile
#
# Profile.clear()
# @profile main(Ntrimers,L,σ,ϵLJ_in,k_in,Ebend_in,boxSize,tmax,renderFlag)
# Juno.profiler(; C=true)

# using BenchmarkTools
# @benchmark main(Ntrimers,L,σ,ϵLJ_in,k_in,Ebend_in,boxSize,tmax,renderFlag)
