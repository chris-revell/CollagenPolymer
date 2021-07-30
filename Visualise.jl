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
using GeometryBasics
using Colors
using Printf


function visualise(foldername,nParticles,σ)
    data = readdlm(foldername*"/output.txt",',',Float64)
    nImages = floor(Int64,size(data)[1]/nParticles)
    r = fill(Point3(0.0,0.0,0.0),nParticles)
    set_theme!(show_axis = false, scale_plot = false, resolution = (800, 800))
    for i in 0:nImages-1
        scene = Scene()
        meshscatter!(data[i*nParticles+1:(i+1)*nParticles,:],markersize=σ/2.0,color=:green)
        #r .= Point3.(eachslice(data[i*nParticles+1:(i+1)*nParticles,:],dims=1))
        #for j=1:nParticles
        #    mesh!(Sphere(r[j],σ),color=:green)
        #end
        save("$foldername/$(@sprintf("%03d",i)).png",scene)
    end
end

export visualise

end
