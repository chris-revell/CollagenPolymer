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
include("./intraMonomerForces.jl")
include("./vanderWaalsForces.jl")
include("./bendingModulus.jl")
include("./initialise.jl")
include("./createRunDirectory.jl")
using .IntraMonomerForces
using .BendingModulus
using .VanderWaalsforces
using .CalculateNoise
using .UpdateSystem
using .OutputData
using .Initialise
using .CreateRunDirectory

# Define run parameters
const Nmonomers      = 15            # Number of collagen monomers
const L              = 0.5           # Length of one monomer
const a              = 0.05          # Diameter of monomer = Diameter of particle within monomer
const Ndomains       = ceil(Int64,(5.0*L)/(4.0*a)+1.0) # Number of particles per monomer, ensuring re<=0.4σ
const boxSize        = 1.0           # Dimensions of cube in which particles are initialised

# Thermodynamic parameters
const μ              = 1.0           # Fluid viscosity
const kT             = 2.0           # Boltzmann constant*Temperature  3.76917046 × 10-21

# Force parameters
const ϵLJ            = kT/1.0       # External Lennard-Jones energy
#const ϵWCA           = ϵLJ/10.0     # External Weeks-Chandler-Anderson (purely repulsive) energy. WCA = LJ+ϵ for r<re (r<σ^(1/6)) and 0 otherwise.
const σ              = 2.0*a         # External LJ
const k              = 100.0*kT      # Internal spring stiffness
const re             = L/(Ndomains-1)# Equilibrium separation of springs connecting internal particles
const Ebend          = 100.0*kT      # Internal bending modulus

# Derived parameters
const D = kT/(6.0*π*μ*a)             # Diffusion constant

# Simulation parameters
const tmax           = 2000.00       # Total simulation time
const dt             = 0.0001        # Time step between iterations
const outputInterval = tmax/100.0           # Time interval for writing position data to file
const renderFlag     = 1             # Controls whether or not system is visualised with povRay automatically

# Data arrays
const pos            = MMatrix{Nmonomers*Ndomains,3}(zeros(Nmonomers*Ndomains,3)) # xyz positions of all particles
const F              = MMatrix{Nmonomers*Ndomains,3}(zeros(Ndomains*Nmonomers,3)) # xyz dimensions of all forces applied to particles
const W              = MMatrix{Nmonomers*Ndomains,3}(zeros(Ndomains*Nmonomers,3)) # xyz values of stochastic Wiener process for all particles

# Define function for bringing together modules to run simulation
@inline function runsim(Nmonomers::Int64,Ndomains::Int64,tmax::Float64,dt::Float64,outputInterval::Float64,boxSize::Float64,σ::Float64,k::Float64,Ebend::Float64,ϵLJ::Float64,re::Float64,D::Float64,kT::Float64,pos::MMatrix,F::MMatrix,W::MMatrix,renderFlag::Int64)

    # Allocate variables needed for calculations
    t = 0.0
    AA = zeros(3)
    BB = zeros(3)
    CC  = zeros(3)
    DD = MArray{Tuple{3},Float64,1,3}(zeros(3))
    foldername = createRunDirectory(Nmonomers,L,a,Ndomains,μ,kT,ϵLJ,σ,k,re,Ebend,D,tmax,dt,outputInterval)
    outfile = open("output/$(foldername)/output.txt","w")

    initialise(pos,Nmonomers,Ndomains,re,boxSize)

    while t<tmax

        if (t%outputInterval)<dt
            outputData(pos,outfile,t,tmax)
        end

        intraMonomerforces!(pos,F,Nmonomers,Ndomains,k,re,AA)

        bendingModulus!(pos,F,Nmonomers,Ndomains,Ebend,AA,BB,CC)

        vanderWaalsForces!(pos,F,Nmonomers,Ndomains,ϵLJ,σ,DD)

        calculateNoise!(W,Nmonomers,Ndomains,dt)

        t = updateSystem!(pos,F,W,Nmonomers,Ndomains,t,dt,D,kT)

    end
    close(outfile)

    if renderFlag == 1
        run(`python3 visualise.py output/$foldername`)
    end
end

runsim(Nmonomers,Ndomains,tmax,dt,outputInterval,boxSize,σ,k,Ebend,ϵLJ,re,D,kT,pos,F,W,renderFlag)

#using BenchmarkTools

#bm = @benchmark runsim(Nmonomers,Ndomains,tmax,dt,outputInterval,zetaMag,k,Ebend,epsilon,re,l,pos,v)
#println(minimum(bm))
#println(mean(bm))
#println(maximum(bm))
#println(allocs(bm))
