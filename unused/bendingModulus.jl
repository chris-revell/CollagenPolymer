#
#  bendingModulus.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

#include "bendingModulus.hpp"
#include "monomer.hpp"
#include <armadillo>
#include <omp.h>

using namespace std
using namespace arma

void bendingModulus(vector<monomer>& Monomers,const int& Nmonomers,const float& Ebend){

  rowvec AC = rowvec(3,fill::zeros)
  rowvec AB = rowvec(3,fill::zeros)
  rowvec F = rowvec(3,fill::zeros)

  #pragma omp parallel for private(AC,AB,F) shared(Monomers,Nmonomers,Ebend)
  {
  # Loop over all monomers
  for (int ii=0 ii<Nmonomersii++){
    # Loop over all sets of 3 in each monomer chain
    for (int jj=0 jj<(Monomers[ii].Ndomains-2) jj++){
      for (int kk=0 kk<3 kk++){
        AC = Monomers[ii].pos.row(3.0*(jj+2)+kk)-Monomers[ii].pos.row(3.0*jj+kk)
        AB = Monomers[ii].pos.row(3.0*(jj+1)+kk)-Monomers[ii].pos.row(3.0*jj+kk)
        F = Ebend*(AC/2.0 - AB)

        Monomers[ii].v.row(3.0*jj+kk) = Monomers[ii].v.row(3.0*jj+kk) - F/2.0
        Monomers[ii].v.row(3.0*(jj+1)+kk) = Monomers[ii].v.row(3.0*(jj+1)+kk) + F
        Monomers[ii].v.row(3.0*(jj+2)+kk) = Monomers[ii].v.row(3.0*(jj+2)+kk) - F/2.0
      }
    }
  }
  }
}
