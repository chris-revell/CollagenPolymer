#
#  internalForces.jl
#  collagen-model-jl
#
#  Created by Christopher Revell on 04/07/2020.
#
#

module InternalForces

using LinearAlgebra

@inline function internalForces!(pos,F,Ntrimers,Ndomains,k,re,Ebend,AA,BB,CC)

	# Loop over all particles
	@threads for jj=1:Ndomains*Ntrimers

		# Tension forces between trimer domains
		if jj%Ndomains == 0
			# Skip at the end of each trimer
		else
			# Calculate tension forces between adjacent particles with Hookean spring potential
			AA[:,threadid()] .= pos[jj+1,:] .- pos[jj,:]
		    AA_mag = sqrt(dot(AA[:,threadid()],AA[:,threadid()]))
		    dif = AA_mag - re
		    F[jj,:,threadid()]   .+= dif*k.*AA[:,threadid()]./AA_mag
		    F[jj+1,:,threadid()] .-= dif*k.*AA[:,threadid()]./AA_mag

			# Bending forces
			if (jj+1)%Ndomains == 0
				# Skip at the end of each trimer
			else
				# Bending force from small displacement potential
				BB[:,threadid()] .= pos[jj+2,:] .- pos[jj+1,:]
				CC[:,threadid()] .= Ebend.*(BB[:,threadid()]./2.0 .- AA[:,threadid()])
				F[jj,:,threadid()]   .-= CC[:,threadid()]./2.0
				F[jj+1,:,threadid()] .+= CC[:,threadid()]
				F[jj+2,:,threadid()] .-= CC[:,threadid()]./2.0
			end
		end
	end
	return nothing
end

export internalForces!

end
