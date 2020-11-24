#
#  main.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 20/11/2020.
#
#

# Import local program modules
push!(LOAD_PATH, "./")
using Simulate


# Define run parameters
nTrimers       = 4        # Number of collagen trimers
L              = 0.5      # Length of one trimer
σ              = 0.005     # Diameter of trimer = Diameter of particle within trimer/External LJ length scale (separation at which V=0) = 2*particle radius
ϵLJ_in         = 100.0     # External Lennard-Jones energy
k_in           = 1000.0     # Internal spring stiffness for Fraenkel spring forces between adjacent particles within a trimer
Ebend_in       = 0.0     # Internal bending modulus of trimer
boxSize        = 1.0      # Dimensions of cube in which particles are initialised
tMax           = 0.01   # Total simulation time
outputFlag     = 0        # Controls whether or not data is printed to file
renderFlag     = 0        # Controls whether or not system is visualised with povRay automatically


# Quick run to precompile
#simulate(1,0.5,0.05,1.0,1.0,1.0,1.0,0.00001,0,0)

#using BenchmarkTools
#println("Timing")
#@benchmark simulate(nTrimers,L,σ,ϵLJ_in,k_in,Ebend_in,boxSize,tMax,0,0) seconds=300

#using Profile
#Profile.clear()
#@profile simulate(nTrimers,L,σ,ϵLJ_in,k_in,Ebend_in,boxSize,tMax,0,0)
#Juno.profiler(; C=true)
