#
#  internalForces.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 04/07/2020.
#
#

module InternalForces

using LinearAlgebra
using .Threads

@inline function internalForces!(pos,F,Ntrimers,Ndomains,k,re,Ebend,AA,BB,CC)

	@threads for jj=1:Ndomains*Ntrimers

		# Tension forces between trimer domains
		if jj%Ndomains == 0
			#skip
		else
		    AA[:,threadid()] = pos[jj+1,:] .- pos[jj,:]
		    AA_mag = sqrt(dot(AA[:,threadid()],AA[:,threadid()]))
		    dif = AA_mag - re
		    F[jj,:,threadid()]   += dif*k.*AA[:,threadid()]./AA_mag
		    F[jj+1,:,threadid()] -= dif*k.*AA[:,threadid()]./AA_mag

			# Bending forces
			if (jj+1)%Ndomains == 0
				#skip
			else
				AA[:,threadid()] = pos[jj+2,:] .- pos[jj,:]
				BB[:,threadid()] = pos[jj+1,:] .- pos[jj,:]
				CC[:,threadid()] = Ebend*(AA[:,threadid()]/2.0 .- BB[:,threadid()])

				F[jj,:,threadid()]   -= CC[:,threadid()]./2.0
				F[jj+1,:,threadid()] += CC[:,threadid()]
				F[jj+2,:,threadid()] -= CC[:,threadid()]./2.0
			end
		end
	end
end

export internalForces!

end
