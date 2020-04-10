#
#  intraMonomerForces.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

function intraMonomerForces!(pos,Nmonomers,k){


  # Tension forces between monomer domains
for ii=1:Nmonomers
for (int jj=0 jj<(Monomers[ii].Ndomains-1) jj++){
            dx = Monomers[ii].pos.row(3.0*(jj+1)+kk) - Monomers[ii].pos.row(3.0*jj+kk)
            dx_mag = sqrt(dot(dx,dx))
            dif = dx_mag - Monomers[ii].domainLength
            Monomers[ii].v.row(3.0*jj+kk) = Monomers[ii].v.row(3.0*jj+kk) + dif*k*dx/dx_mag
            Monomers[ii].v.row(3.0*(jj+1)+kk) = Monomers[ii].v.row(3.0*(jj+1)+kk) - dif*k*dx/dx_mag
        end
    end


end
