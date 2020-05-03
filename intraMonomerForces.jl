#
#  intraMonomerForces.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#
__precompile__()
module IntraMonomerForces

using LinearAlgebra
using StaticArrays

@inline function intraMonomerforces!(pos::MMatrix,F::MMatrix,Nmonomers::Int64,Ndomains::Int64,k::Float64,re::Float64,dx::Array{Float64,1})
	# Tension forces between monomer domains
	@inbounds for ii=0:Nmonomers-1
		@inbounds for jj=1:Ndomains-1
		    dx = pos[ii*Ndomains+jj+1,:] .- pos[ii*Ndomains+jj,:]
		    dx_mag = sqrt(dot(dx,dx))
		    dif = dx_mag - re
		    F[ii*Ndomains+jj,:]   .+= dif*k.*dx./dx_mag
		    F[ii*Ndomains+jj+1,:] .-= dif*k.*dx./dx_mag
		end
	end
end

export intraMonomerforces!

end
