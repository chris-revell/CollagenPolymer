#
#  main.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

# Import Julia packages
using Distributions
using StaticArrays
# Import program modules
include("./outputData.jl")
include("./calculateNoise.jl")
include("./updateSystem.jl")
include("./intraTrimerForces.jl")
include("./vanderWaalsForces.jl")
include("./bendingModulus.jl")
include("./initialise.jl")
include("./createRunDirectory.jl")
using .IntraTrimerForces
using .BendingModulus
using .VanderWaalsforces
using .CalculateNoise
using .UpdateSystem
using .OutputData
using .Initialise
using .CreateRunDirectory

# Define run parameters
const Ntrimers       = 20             # Number of collagen trimers
const L              = 0.5           # Length of one trimer
const a              = 0.05          # Diameter of trimer = Diameter of particle within trimer
const Ndomains       = 1+ceil(Int64,(5.0*L)/(4.0*a)+1.0) # Number of particles per trimer, ensuring re<=0.4σ
const boxSize        = 1.0           # Dimensions of cube in which particles are initialised

# Thermodynamic parameters
const μ              = 1.0           # Fluid viscosity
const kT             = 1.0           # Boltzmann constant*Temperature  3.76917046 × 10-21

# Force parameters
const ϵLJ            = 1.0*kT        # External Lennard-Jones energy
const σ              = a             # External LJ length scale (separation at which V=0) = 2*particle radius
const k              = 1000.0*kT     # Internal spring stiffness for forces between adjacent particles within a trimer
const re             = L/(Ndomains-1)# Equilibrium separation of springs connecting adjacent particles within a trimer
const Ebend          = 100000.0*kT   # Internal bending modulus of trimer

# Derived parameters
const D              = kT/(6.0*π*μ*a)# Diffusion constant
const trimerVolume   = L*π*a^2       # Volume of one trimer
#const ϕ              = trimerVolume/(2.0*boxSize)^3 # Volume fraction

# Simulation parameters
const tmax           = 5.00          # Total simulation time
const outputInterval = tmax/100.0    # Time interval for writing position data to file
const renderFlag     = 1             # Controls whether or not system is visualised with povRay automatically
const interactionThreshold = 2.0*σ   # Threshold for van der Waals interactions
const boxMultiples   = 4             # Multiple of boxsizes over which to define cell list grid to allow for system expansion
const Ng             = ceil(Int64,boxMultiples*boxSize/interactionThreshold)+1 #

# Data arrays
const pos            = MMatrix{Ntrimers*Ndomains,3}(zeros(Ntrimers*Ndomains,3)) # xyz positions of all particles
const F              = MMatrix{Ntrimers*Ndomains,3}(zeros(Ndomains*Ntrimers,3)) # xyz dimensions of all forces applied to particles
const W              = MMatrix{Ntrimers*Ndomains,3}(zeros(Ndomains*Ntrimers,3)) # xyz values of stochastic Wiener process for all particles
const cellLists      = zeros(Int64,Ng,Ng,Ng,50) # Cell list grid. Too many components for static array?


# Define function for bringing together modules to run simulation
@inline function runsim(Ntrimers::Int64,Ndomains::Int64,tmax::Float64,outputInterval::Float64,boxSize::Float64,σ::Float64,k::Float64,Ebend::Float64,ϵLJ::Float64,re::Float64,D::Float64,kT::Float64,pos::MMatrix,F::MMatrix,W::MMatrix,renderFlag::Int64,Ng::Int64,cellLists::Array{Int64},interactionThreshold::Float64,boxMultiples::Int64)

    # Initialise system time
    t = 0.0

    # Allocate variables to reuse in calculations and prevent memory reallocations
    AA = zeros(3)
    BB = zeros(3)
    CC  = zeros(3)
    DD = MArray{Tuple{3},Float64,1,3}(zeros(3))

    # Setup data folder and output files
    foldername = createRunDirectory(Ntrimers,L,a,Ndomains,μ,kT,ϵLJ,σ,k,re,Ebend,D,tmax,outputInterval)
    outfile = open("output/$(foldername)/output.txt","w")

    # Initialise trimers within boxSize space
    initialise(pos,Ntrimers,Ndomains,re,boxSize)

    # Output initial state
    outputData(pos,outfile,t,tmax)

    #Iterate over time until max system time is reached
    while t<tmax

        # Apply periodic boundary conditions
        #pos .= (pos .+ boxSize).%boxSize

        # Create cell list to identify trimer pairs within interaction range
        fill!(cellLists,0)
        for i=1:Ntrimers*Ndomains
            iₓ = ceil.(Int64,(pos[i,:] .+ boxMultiples*boxSize/2.0)./interactionThreshold)
            #println(iₓ)
            cellLists[iₓ...,1] += 1
            cellLists[iₓ...,cellLists[iₓ...,1]+1] = i

        end

        # Spring forces between particles along trimer chain
        intraTrimerforces!(pos,F,Ntrimers,Ndomains,k,re,AA)

        # Calculate forces from trimer bending stiffness
        bendingModulus!(pos,F,Ntrimers,Ndomains,Ebend,AA,BB,CC)

        # Calculate van der Waals/electrostatic interactions between nearby trimer domains
        vanderWaalsForces!(pos,F,Ntrimers,Ndomains,ϵLJ,σ,DD,cellLists,Ng)

        # Adapt timestep to maximum force value
        Fmax_sq = max(sum(F.*F,dims=2)...)
        dt = min(σ^2/(32*D),kT*σ/(2.0*D*sqrt(Fmax_sq)))

        # Find stochastic term (Wiener process) for all monomers
        calculateNoise!(W,Ntrimers,Ndomains,dt)

        # Integrate system with forward euler
        t = updateSystem!(pos,F,W,Ntrimers,Ndomains,t,dt,D,kT)

        if (t%outputInterval)<dt
            outputData(pos,outfile,t,tmax)            
            # Measure trimer lengths
        #    for i=1:Ntrimers
        #        for j=1:Ndomains
        end

    end
    close(outfile)

    # Render images from output files
    if renderFlag == 1
        run(`python3 visualise.py output/$foldername`)
    end
end

runsim(Ntrimers,Ndomains,tmax,outputInterval,boxSize,σ,k,Ebend,ϵLJ,re,D,kT,pos,F,W,renderFlag,Ng,cellLists,interactionThreshold,boxMultiples)
