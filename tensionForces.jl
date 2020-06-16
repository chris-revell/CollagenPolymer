#
#  intraTrimerForces.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module TensionForces

using LinearAlgebra
using StaticArrays

@inline function tensionForces!(pos,F,Ntrimers,Ndomains,k,re,dx)
	# Tension forces between trimer domains
	for ii=0:Ntrimers-1
		for jj=1:Ndomains-1 #Threads.@threads for jj=1:Ndomains-1
		    dx = pos[ii*Ndomains+jj+1,:] .- pos[ii*Ndomains+jj,:]
		    dx_mag = sqrt(dot(dx,dx))
		    dif = dx_mag - re
		    F[ii*Ndomains+jj,:]   .+= dif*k.*dx./dx_mag
		    F[ii*Ndomains+jj+1,:] .-= dif*k.*dx./dx_mag
		end
	end
end

export tensionForces!

end
