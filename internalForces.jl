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

	for jj=1:Ndomains*Ntrimers

		# Tension forces between trimer domains
		if jj%Ndomains == 0
			#skip
		else
			AA = pos[jj+1,:] .- pos[jj,:]
		    AA_mag = sqrt(dot(AA,AA))
		    dif = AA_mag - re
		    F[jj,:]   .+= (dif*k/AA_mag).*AA
		    F[jj+1,:] .-= (dif*k/AA_mag).*AA

			# Bending forces
			if (jj+1)%Ndomains == 0
				#skip
			else
				BB = pos[jj+2,:] .- pos[jj+1,:]
				CC = Ebend.*(BB./2.0 .- AA)
				F[jj,:]   .-= CC./2.0
				F[jj+1,:] .+= CC
				F[jj+2,:] .-= CC./2.0
			end
		end
	end
end

export internalForces!

end
