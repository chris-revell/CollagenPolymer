#
#  intraMonomerForces.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module IntraMonomerForces

using LinearAlgebra

function intraMonomerForces!(pos,v,Nmonomers,k,l)
	# Tension forces between monomer domains
	for ii=1:Nmonomers
	    dx = pos[2*ii+1] - pos[2*ii]
	    dx_mag = sqrt(dot(dx,dx))
	    dif = dx_mag - l
	    v[2*ii]   = v[2*ii] + dif*k*dx/dx_mag
	    v[2*ii+1] = v[2*ii+1] - dif*k*dx/dx_mag
	end
end

export intraMonomerForces!

end
