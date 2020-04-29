#
#  main.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

# Import Julia packages
using Distributions
using Dates
using StaticArrays
# Import program modules
include("./outputData.jl")
include("./calculateNoise.jl")
include("./updateSystem.jl")
include("./intraMonomerForces.jl")
include("./vanderWaalsForces.jl")
include("./bendingModulus.jl")
using .IntraMonomerForces
using .BendingModulus
using .VanderWaalsforces
using .CalculateNoise
using .UpdateSystem
using .OutputData


# Define run parameters
const Nmonomers      = 2            # Number of collagen monomers
const L              = 0.2          # Length of one monomer
const a              = 0.02        # Diameter of monomer = Diameter of particle within monomer
const Ndomains       = ceil(Int64,(5.0*L)/(4.0*a)+1.0) # Number of particles per monomer, ensuring re<=0.4σ

# Thermodynamic parameters
const μ              = 1.0          # Fluid viscosity
const kT             = 1.0          # Boltzmann constant*Temperature

# Force parameters
const ϵLJ            = 10*kT        # External Lennard-Jones energy
const ϵWCA           = ϵLJ/100.0    # External Weeks-Chandler-Anderson (purely repulsive) energy. WCA = LJ+ϵ for r<re (r<σ^(1/6)) and 0 otherwise.
const σ              = 2.0*a        # External LJ
const k              = 2.0          # Internal spring stiffness
const re             = L/(Ndomains-1)      # Equilibrium separation of springs connecting internal particles
const Ebend          = 5.0          # Internal bending modulus

# Derived parameters
const D = kT/(6.0*π*μ*a)                  # Diffusion constant

# Simulation parameters
const tmax           = 100.0      #
const dt             = 0.01         #
const outputInterval = 50.0         #

# Data arrays
pos                  = MMatrix{Nmonomers*Ndomains,3}(rand(Uniform(-2.0,2.0),Nmonomers*Ndomains,3))  #
F                    = MMatrix{Nmonomers*Ndomains,3}(zeros(Ndomains*Nmonomers,3))                   #
W                    = MMatrix{Nmonomers*Ndomains,3}(zeros(Ndomains*Nmonomers,3))                   #
#%%
# Define function for bringing together modules to run simulation
@inline function runsim(Nmonomers::Int64,Ndomains::Int64,tmax::Float64,dt::Float64,outputInterval::Float64,σ::Float64,k::Float64,Ebend::Float64,ϵLJ::Float64,re::Float64,D::Float64,kT::Float64,pos::MMatrix,F::MMatrix,W::MMatrix)

    foldername = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")
    mkdir("output/$(foldername)")
    open("output/$(foldername)/conditions.txt","w") do conditionsfile
        println(conditionsfile,"Nmonomers      ",Nmonomers      )
        println(conditionsfile,"Ndomains       ",Ndomains       )
        println(conditionsfile,"tmax           ",tmax           )
        println(conditionsfile,"dt             ",dt             )
        println(conditionsfile,"outputInterval ",outputInterval )
        #println(conditionsfile,"zetaMag        ",zetaMag        )
        println(conditionsfile,"k              ",k              )
        println(conditionsfile,"Ebend          ",Ebend          )
        println(conditionsfile,"ϵLJ            ",ϵLJ            )
        println(conditionsfile,"re             ",re             )
        #println(conditionsfile,"l              ",l              )
    end

    # Allocate variables needed for calculations
    t = 0.0
    AA = zeros(3)
    BB = zeros(3)
    CC  = zeros(3)
    DD = MArray{Tuple{3},Float64,1,3}(zeros(3))
    outfile = open("output/$(foldername)/output.txt","w")

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
end

runsim(Nmonomers,Ndomains,tmax,dt,outputInterval,σ,k,Ebend,ϵLJ,re,D,kT,pos,F,W)

#using BenchmarkTools

#bm = @benchmark runsim(Nmonomers,Ndomains,tmax,dt,outputInterval,zetaMag,k,Ebend,epsilon,re,l,pos,v)
#println(minimum(bm))
#println(mean(bm))
#println(maximum(bm))
#println(allocs(bm))
