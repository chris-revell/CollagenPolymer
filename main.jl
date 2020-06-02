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
const Ntrimers       = 2             # Number of collagen trimers
const L              = 0.5           # Length of one trimer
const a              = 0.05          # Diameter of trimer = Diameter of particle within trimer
const Ndomains       = 1+ceil(Int64,(5.0*L)/(4.0*a)+1.0) # Number of particles per trimer, ensuring re<=0.4σ
const boxSize        = 1.0           # Dimensions of cube in which particles are initialised

# Thermodynamic parameters
const μ              = 1.0           # Fluid viscosity
const kT             = 1.0           # Boltzmann constant*Temperature  3.76917046 × 10-21

# Force parameters
const ϵLJ            = 100.0*kT      # External Lennard-Jones energy
#const ϵWCA           = ϵLJ/10.0     # External Weeks-Chandler-Anderson (purely repulsive) energy. WCA = LJ+ϵ for r<re (r<σ^(1/6)) and 0 otherwise.
const σ              = a#2.0*a       # External LJ length scale (separation at which V=0) = 2*particle radius
const k              = 1000.0*kT     # Internal spring stiffness for forces between adjacent particles within a trimer
const re             = L/(Ndomains-1)# Equilibrium separation of springs connecting adjacent particles within a trimer
const Ebend          = 10.0*kT       # Internal bending modulus of trimer

# Derived parameters
const D              = kT/(6.0*π*μ*a)# Diffusion constant
const trimerVolume   = L*π*a^2       # Volume of one trimer
#const ϕ              = trimerVolume/(2.0*boxSize)^3 # Volume fraction

# Simulation parameters
const tmax           = 10.00         # Total simulation time
dt             = 0.0001              # Time step between iterations
const outputInterval = tmax/100.0    # Time interval for writing position data to file
const renderFlag     = 1             # Controls whether or not system is visualised with povRay automatically

# Data arrays
const pos            = MMatrix{Ntrimers*Ndomains,3}(zeros(Ntrimers*Ndomains,3)) # xyz positions of all particles
const F              = MMatrix{Ntrimers*Ndomains,3}(zeros(Ndomains*Ntrimers,3)) # xyz dimensions of all forces applied to particles
const W              = MMatrix{Ntrimers*Ndomains,3}(zeros(Ndomains*Ntrimers,3)) # xyz values of stochastic Wiener process for all particles
const Ng             = 2*ceil(Int64,2.0*boxSize/(4.0*σ))+1
const cellLists      = zeros(Int64,Ng,Ng,Ng,50) # Cell list grid. Too many components for static array?


# Define function for bringing together modules to run simulation
@inline function runsim(Ntrimers::Int64,Ndomains::Int64,tmax::Float64,dt::Float64,outputInterval::Float64,boxSize::Float64,σ::Float64,k::Float64,Ebend::Float64,ϵLJ::Float64,re::Float64,D::Float64,kT::Float64,pos::MMatrix,F::MMatrix,W::MMatrix,renderFlag::Int64,Ng::Int64,cellLists::Array{Int64})

    # Allocate variables needed for calculations
    t = 0.0
    AA = zeros(3)
    BB = zeros(3)
    CC  = zeros(3)
    DD = MArray{Tuple{3},Float64,1,3}(zeros(3))
    foldername = createRunDirectory(Ntrimers,L,a,Ndomains,μ,kT,ϵLJ,σ,k,re,Ebend,D,tmax,dt,outputInterval)
    outfile = open("output/$(foldername)/output.txt","w")

    initialise(pos,Ntrimers,Ndomains,re,boxSize)

    while t<tmax

        if (t%outputInterval)<dt
            outputData(pos,outfile,t,tmax)

            for i=1:Ntrimers
                for j=1:Ndomains
                    

        end

        # Create cell list to identify trimer pairs within interaction range
        fill!(cellLists,0)
        # Apply periodic boundary conditions
        pos .= (pos .+ boxSize).%boxSize

        for i=1:Ntrimers*Ndomains
            iₓ = ceil.(Int64,(pos[i,:] .+ 2.0*boxSize)./(4.0*σ))
            cellLists[iₓ...,1] += 1
            cellLists[iₓ...,cellLists[iₓ...,1]+1] = i
        end

        # Spring forces along trimer chain
        intraTrimerforces!(pos,F,Ntrimers,Ndomains,k,re,AA)

        #bendingModulus!(pos,F,Ntrimers,Ndomains,Ebend,AA,BB,CC)

        # Calculate van der Waals/electrostatic interactions between nearby trimer domains
        vanderWaalsForces!(pos,F,Ntrimers,Ndomains,ϵLJ,σ,DD,cellLists,Ng)

        # Adapt timestep to maximum force value
        Fmax_sq = max(sum(F.*F,dims=2)...)
        dt = min(σ^2/(32*D),kT*σ/(2.0*D*sqrt(Fmax_sq)))
        #println(dt)
        calculateNoise!(W,Ntrimers,Ndomains,dt)

        t = updateSystem!(pos,F,W,Ntrimers,Ndomains,t,dt,D,kT)

    end
    close(outfile)

    if renderFlag == 1
        run(`python3 visualise.py output/$foldername`)
    end
end

runsim(Ntrimers,Ndomains,tmax,dt,outputInterval,boxSize,σ,k,Ebend,ϵLJ,re,D,kT,pos,F,W,renderFlag,Ng,cellLists)
