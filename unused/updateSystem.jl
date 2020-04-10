#
#  updateSystem.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

#include "updateSystem.hpp"
#include <armadillo>
#include "monomer.hpp"
#include <vector>
#include <omp.h>

using namespace std
using namespace arma

void updateSystem(vector<monomer>& Monomers,const int& Nmonomers,float& t,const float& dt,const float& gamma){

  #pragma omp parallel for shared(Monomers,Nmonomers,dt,gamma)
  {
  for (int ii=0 ii<Nmonomers ii++){
    Monomers[ii].pos = Monomers[ii].pos + Monomers[ii].v*dt/gamma
    Monomers[ii].v.zeros()
  }
  }
  t = t+dt
}
