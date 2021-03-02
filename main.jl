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


simulate(nMonomers,L,σ,ϵLJ_in,k_in,Ebend_in,boxSize,tMax,outputFlag,renderFlag)

# Quick run to precompile
#simulate(1,0.5,0.05,1.0,1.0,1.0,1.0,0.00001,0,0)

#using BenchmarkTools
#println("Timing")
#@benchmark simulate(nMonomers,L,σ,ϵLJ_in,k_in,Ebend_in,boxSize,tMax,0,0) seconds=300

#using Profile
#Profile.clear()
#@profile simulate(nMonomers,L,σ,ϵLJ_in,k_in,Ebend_in,boxSize,tMax,0,0)
#Juno.profiler(; C=true)
