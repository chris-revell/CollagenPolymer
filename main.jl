#
#  main.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

#include("./monomer.jl")
include("./calculateNoise.jl")
include("./updateSystem.jl")
include("./intraMonomerForces.jl")
include("./vanderWaalsForces.jl")
include("./bendingModulus.jl")

using Distributions
using DelimitedFiles
using Printf
using Dates
using StaticArrays

using .IntraMonomerForces
using .BendingModulus
using .VanderWaalsForces
using .CalculateNoise
using .UpdateSystem

const Nmonomers      = 5
const Ndomains       = 5
const tmax           = 10000.0
const dt             = 0.01
const outputInterval = 50.0
const zetaMag        = 1.0
const k              = 2.0
const Ebend          = 5.0
const epsilon        = 0.0005
const re             = 0.1
const l              = 0.2
pos                  = MMatrix{Nmonomers*Ndomains,3}(rand(Uniform(-1.0,1.0),Nmonomers*Ndomains,3))
v                    = MMatrix{Nmonomers*Ndomains,3}(zeros(Ndomains*Nmonomers,3))

foldername = Dates.format(Dates.now(),"yyyy-mm-dd-HH-MM-SS")

mkdir("output/$(foldername)")
open("output/$(foldername)/conditions.txt","w") do conditionsfile
    println(conditionsfile,"Nmonomers      ",Nmonomers      )
    println(conditionsfile,"Ndomains       ",Ndomains       )
    println(conditionsfile,"tmax           ",tmax           )
    println(conditionsfile,"dt             ",dt             )
    println(conditionsfile,"outputInterval ",outputInterval )
    println(conditionsfile,"zetaMag        ",zetaMag        )
    println(conditionsfile,"k              ",k              )
    println(conditionsfile,"Ebend          ",Ebend          )
    println(conditionsfile,"epsilon        ",epsilon        )
    println(conditionsfile,"re             ",re             )
    println(conditionsfile,"l              ",l              )
end


@inline function runsim(Nmonomers::Int64,Ndomains::Int64,tmax::Float64,dt::Float64,outputInterval::Float64,zetaMag::Float64,k::Float64,Ebend::Float64,epsilon::Float64,re::Float64,l::Float64,pos::MMatrix,v::MMatrix,foldername::String)

    t=0.0
    outfile = open("output/$(foldername)/output.txt","w")

    while t<tmax

        if (t%outputInterval)<dt
            writedlm(outfile,pos,", ")
            Printf.@printf("Simulating: %f/%f\n",t,tmax)
        end

        intraMonomerForces!(pos,v,Nmonomers,Ndomains,k,l)

        bendingModulus!(pos,v,Nmonomers,Ndomains,Ebend)

        vanderWaalsForces!(pos,v,Nmonomers,Ndomains,epsilon,re)

        calculateNoise!(v,Nmonomers,Ndomains,zetaMag)

        t = updateSystem!(pos,v,Nmonomers,Ndomains,t,dt)

    end

    close(outfile)
end

runsim(Nmonomers,Ndomains,tmax,dt,outputInterval,zetaMag,k,Ebend,epsilon,re,l,pos,v,foldername)
