#
#  main.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

#include("./monomer.jl")
include("./calculateNoise.jl")
#include("outputData.jl")
#include("readParameters.jl")
include("./intraMonomerForces.jl")
include("./vanderWaalsForces.jl")

using Distributions
using DelimitedFiles
using Printf

Nmonomers      = 10
tmax           = 1000.0
outputInterval = 100.0
dt             = 0.01
zetaMag        = 1.0
k              = 1.0
epsilon        = 1.0
re             = 0.1
t              = 0.0
l              = 0.2
pos            = rand(Uniform(-1.0,1.0),2*Nmonomers,3)
v              = zeros(2*Nmonomers,3)

open("output/conditions.txt","w") do conditionsfile
    println(conditionsfile,"Nmonomers      ",Nmonomers      )
    println(conditionsfile,"tmax           ",tmax           )
    println(conditionsfile,"outputInterval ",outputInterval )
    println(conditionsfile,"dt             ",dt             )
    println(conditionsfile,"zetaMag        ",zetaMag        )
    println(conditionsfile,"epsilon        ",epsilon        )
    println(conditionsfile,"re             ",re             )
end

outfile = open("output/output.txt","w")

while t<tmax

    if (t%outputInterval)<dt
        writedlm(outfile,pos,", ")
        Printf.@printf("Simulating: %f/%f\n",t,tmax)
    end

    IntraMonomerForces.intraMonomerForces!(pos,v,Nmonomers,k,l)

    VanderWaalsForces.vanderWaalsForces!(pos,v,Nmonomers,epsilon,re)

    CalculateNoise.calculateNoise!(v,zetaMag,Nmonomers)

    global pos .= pos .+ v.*dt/100
    global t += dt
    global v .= 0

end

close(outfile)
