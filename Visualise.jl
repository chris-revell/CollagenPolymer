#
#  Visualise.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 22/11/2020.
#
#

module Visualise

using DelimitedFiles
using Makie
using LinearAlgebra
using AbstractPlotting
using GeometryBasics
using Colors
using Printf


function visualise(foldername,nParticles,σ)
    nImages = 100
    data = readdlm(foldername*"/output.txt",',',Float64)
    r = fill(Point3(0.0,0.0,0.0),nParticles)
    set_theme!(show_axis = false, scale_plot = false, resolution = (800, 800))
    for i in 0:nImages
        scene = Scene()
        meshscatter!(data[i*nParticles+1:(i+1)*nParticles,:],markersize=σ,color=:green)
        #r .= Point3.(eachslice(data[i*nParticles+1:(i+1)*nParticles,:],dims=1))
        #for j=1:nParticles
        #    mesh!(Sphere(r[j],σ),color=:green)
        #end
        save("$foldername/$(@sprintf("%03d",i)).png",scene)
    end
end

export visualise

end
