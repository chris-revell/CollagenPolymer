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
using .TensionForces
using .BendingForces
using .InterTrimerForces
using .CalculateNoise
using .UpdateSystem
using .OutputData
using .Initialise
using .CreateRunDirectory
using .BoundaryForce


# Define run parameters
const Ntrimers       = 10                  # Number of collagen trimers
const L              = 0.5                # Length of one trimer
const σ              = 0.05               # Diameter of trimer = Diameter of particle within trimer/External LJ length scale (separation at which V=0) = 2*particle radius
const ϵLJ_in         = 10.0                # External Lennard-Jones energy
const k_in           = 100000.0           # Internal spring stiffness for forces between adjacent particles within a trimer
const Ebend_in       = 0.0                # Internal bending modulus of trimer
const boxSize        = 1.0                # Dimensions of cube in which particles are initialised
const tmax           = 10.0                # Total simulation time
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
    outputInterval = tmax/100.0         # Time interval for writing position data to file
    intrctnThrshld = 2.0*σ              # Threshold for van der Waals interactions
    Ndomains       = 1+ceil(Int64,(5.0*L)/(4.0*σ)+1.0) # Number of particles per trimer, ensuring re<=0.4σ
    WCAthresh_sq   = (2.0^(1.0/6.0)*σ)^2# Cutoff threshold for WCA potential - separation at which force is zero - squared
    re             = (L-σ)/(Ndomains-1) # Equilibrium separation of springs connecting adjacent particles within a trimer
    D              = kT/(6.0*π*μ*σ)     # Diffusion constant
    Ng             = ceil(Int64,boxSize/intrctnThrshld) #
    #trimerVolume   = L*π*σ^2            # Volume of one trimer
    #ϕ              = trimerVolume/(2.0*boxSize)^3 # Volume fraction
    # Data arrays
    pos            = MMatrix{Ntrimers*Ndomains,3}(zeros(Float64,Ntrimers*Ndomains,3)) # xyz positions of all particles
    F              = MMatrix{Ntrimers*Ndomains,3}(zeros(Float64,Ndomains*Ntrimers,3)) # xyz dimensions of all forces applied to particles
    W              = MMatrix{Ntrimers*Ndomains,3}(zeros(Float64,Ndomains*Ntrimers,3)) # xyz values of stochastic Wiener process for all particles
    cellLists      = zeros(Int16,Ng,Ng,Ng,20)                                 # Cell list grid. Too many components for static array?
    nonZeroGrids   = fill(zeros(Int16,3), Ndomains*Ntrimers)


    # Initialise system time
    t = 0.0

    # Allocate variables to reuse in calculations and prevent memory reallocations
    AA = zeros(Float64,3)
    BB = zeros(Float64,3)
    CC = zeros(Float64,3)
    DD = MArray{Tuple{3},Float64,1,3}(zeros(3))

    # Setup data folder and output files
    foldername = createRunDirectory(Ntrimers,L,Ndomains,μ,kT,ϵLJ,σ,k,re,Ebend,D,tmax,outputInterval,Ng,boxSize)
    outfile = open("output/$(foldername)/output.txt","w")

    # Initialise trimers within boxSize space
    initialise(pos,Ntrimers,Ndomains,re,boxSize)

    # Output initial state
    outputData(pos,outfile,t,tmax,Ntrimers,Ndomains,σ)

    #Iterate over time until max system time is reached
    while t<tmax

        # Create cell list to identify trimer pairs within interaction range
        cellLists[:,:,:,1] .= 0
        Nfilled = 0
        fill!(nonZeroGrids,zeros(Int16,3))
        for i::Int16=1:Ntrimers*Ndomains
            iₓ = ceil.(Int64,(pos[i,:] .+ boxSize/2.0)./intrctnThrshld)
            cellLists[iₓ...,1] += 1
            cellLists[iₓ...,cellLists[iₓ...,1]+1] = i
            if iₓ in nonZeroGrids
                # Skip
            else
                Nfilled+=1
                nonZeroGrids[Nfilled] = iₓ
            end
        end

        # Spring forces between particles along trimer chain
        tensionForces!(pos,F,Ntrimers,Ndomains,k,re,AA)

        # Calculate forces from trimer bending stiffness
        bendingForces!(pos,F,Ntrimers,Ndomains,Ebend,AA,BB,CC)

        # Calculate van der Waals/electrostatic interactions between nearby trimer domains
        interTrimerForces!(pos,F,Ntrimers,Ndomains,ϵLJ,σ,DD,cellLists,Ng,WCAthresh_sq,intrctnThrshld,nonZeroGrids,Nfilled)

        boundaryForce!(pos,F,cellLists,nonZeroGrids,Nfilled,Ng,boxSize,σ,ϵLJ)

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
# @profile main(Ntrimers,L,σ,ϵLJ,k,Ebend,boxSize,tmax,renderFlag)
# Juno.profiler(; C=true)

# using BenchmarkTools
# @benchmark main(Ntrimers,L,σ,ϵLJ,k,Ebend,boxSize,tmax,renderFlag)
