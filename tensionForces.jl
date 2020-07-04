#
#  intraTrimerForces.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module TensionForces

using LinearAlgebra
using .Threads

@inline function tensionForces!(pos,F,Ntrimers,Ndomains,k,re,dx)
	# Tension forces between trimer domains

	@threads for jj=1:Ndomains*Ntrimers
		if (jj+1)%Ndomains == 0
			#skip
		else
		    dx[:,threadid()] = pos[jj+1,:] - pos[jj,:]
		    dx_mag = sqrt(dot(dx[:,threadid()],dx[:,threadid()]))
		    dif = dx_mag - re
		    F[jj,:,threadid()]   .+= dif*k.*dx[:,threadid()]./dx_mag
		    F[jj+1,:,threadid()] .-= dif*k.*dx[:,threadid()]./dx_mag
		end
	end
end

export tensionForces!

end
