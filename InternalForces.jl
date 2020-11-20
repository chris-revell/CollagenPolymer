#
#  InternalForces.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 04/07/2020.
#
#

module InternalForces

using LinearAlgebra
using Base.Threads

function internalForces!(pos,F,nTrimers,nDomains,nParticles,k,rₑ,Ebend,AA,AA_bar,BB,BB_bar,CC,DD,DD_bar,EE,EE_bar)

	# Loop over all particles
	@threads for jj=1:nParticles

		# Tension forces between trimer domains
		if jj%nDomains == 0
			# Skip at the end of each trimer
		else
			# Calculate tension forces between adjacent particles with Hookean spring potential
			AA[:,threadid()] .= pos[jj+1] .- pos[jj]
		    AA_mag = sqrt(dot(AA[:,threadid()],AA[:,threadid()]))
			AA_bar[:,threadid()] .= AA[:,threadid()]./AA_mag

			dif = AA_mag - rₑ
		    F[jj,:,threadid()]   .+= dif*k.*AA[:,threadid()]./AA_mag
		    F[jj+1,:,threadid()] .-= dif*k.*AA[:,threadid()]./AA_mag

			# Bending forces
			if (jj+1)%nDomains == 0
				# Skip at the end of each trimer
			else
				# Vector from second particle to third
				BB[:,threadid()] = pos[jj+2 .- pos[jj+1]
				BB_mag = sqrt(dot(BB[:,threadid()],BB[:,threadid()]))
				BB_bar[:,threadid()] .= BB[:,threadid()]./BB_mag

				if abs.(BB_bar[:,threadid()].-AA_bar[:,threadid()]) < [0.000001,0.000001,0.000001]
					# skip if colinear
				else

					# Vector perpendicular to BB and in plane defined by AA and BB
					CC[:,threadid()] = cross(AA[:,threadid()],BB[:,threadid()])

					# Vector perpendicular to AA and in plane defined by AA and BB
					DD[:,threadid()] = cross(AA[:,threadid()],CC[:,threadid()])
					DD_mag = sqrt(dot(DD[:,threadid()],DD[:,threadid()]))
					DD_bar[:,threadid()] .= DD[:,threadid()]./DD_mag

					# Vector opposing the sum of CC and DD
					EE[:,threadid()] = cross(BB[:,threadid()],CC[:,threadid()])
					EE_mag = sqrt(dot(EE[:,threadid()],EE[:,threadid()]))
					EE_bar[:,threadid()] .= EE[:,threadid()]./EE_mag

					θ=acos(dot(-AA_bar[:,threadid()],BB_bar[:,threadid()]))
					Fmag = Ebend*(π-θ)

					F[jj,:,threadid()]   .+= Fmag.*DD_bar[:,threadid()]
					#Fbend[jj,:]   .+= Fmag.*DD_bar
					F[jj+1,:,threadid()] .-= Fmag.*(DD_bar[:,threadid()].+EE_bar[:,threadid()])
					#Fbend[jj+1,:] .-= Fmag.*(DD_bar.+EE_bar)
					F[jj+2,:,threadid()] .+= Fmag.*EE_bar[:,threadid()]
					#Fbend[jj+2,:] .+= Fmag.*EE_bar
				end
			end
		end
	end
	return nothing
end

export internalForces!

end
