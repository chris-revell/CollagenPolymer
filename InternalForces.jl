#
#  InternalForces.jl
#  collagen-brownian-polymer
#
#  Created by Christopher Revell on 04/07/2020.
#
#

module InternalForces

using LinearAlgebra
using StaticArrays

function internalForces!(pos,F,nMonomers,nDomains,nParticles,k,rₑ,Ebend,AA,BB,CC,DD,EE)

	# Loop over all particles
	for jj=1:nParticles

		# Tension forces between monomer domains
		if jj%nDomains == 0
			# Skip at the end of each monomer
		else
			# Calculate tension forces between adjacent particles with Hookean spring potential

			# Find vector AA between 1st and 2nd particle
			AA = pos[jj+1] - pos[jj]
		    AA_mag = norm(AA)
			AA = AA/AA_mag
			#println("AA_mag = $AA_mag")

			dif = AA_mag - rₑ
			#println("dif = $dif")
		    F[jj]   += dif*k*AA
			#println("F[jj]   = $(F[jj]  )")
		    F[jj+1] -= dif*k*AA
			#println("F[jj+1] = $(F[jj+1])")

			# Bending forces
			# if (jj+1)%nDomains == 0
			# 	# Skip at the end of each monomer
			# else
			# 	# Vector BB between 2nd and 3rd particle
			# 	BB = pos[jj+2] - pos[jj+1]
			# 	BB_mag = norm(BB)
			# 	BB = BB/BB_mag
			#
			# 	if dot(BB,AA) > 0.9999999999
			# 		# Skip if colinear, ie cosθ=BB⋅AA≈1
			# 	else
			#
			# 		# Cross produce of AA and BB gives vector CC, perpendicular to plane defined by AA and BB
			# 		CC = AA×BB
			#
			# 		# Find vector DD perpendicular to AA and lying in plane of AA and BB with cross product of AA and CC
			# 		DD = AA×CC
			# 		DD_mag = norm(DD)
			# 		DD = DD/DD_mag
			#
			# 		# Vector opposing the sum of CC and DD
			# 		EE = BB×CC
			# 		EE_mag = norm(EE)
			# 		EE = EE/EE_mag
			#
			# 		θ    = acos(-AA⋅BB)
			# 		Fmag = Ebend*θ
			# 		F[jj]   += Fmag*DD
			# 		F[jj+1] -= Fmag*(DD+EE)
			# 		F[jj+2] += Fmag*EE
			#
			# 	end
			# end
		end

	end
	return nothing
end

export internalForces!

end
