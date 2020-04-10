#
#  rotateVector.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

#include "rotateVector.hpp"
#include <armadillo>
#include <cmath>

using namespace std
using namespace arma

arma::rowvec rotateVector(subview_row<double> x,const rowvec& u,const float& psi)  {

  mat R = mat(3,3,fill::zeros)
  rowvec v
  R(0,0) = cos(psi) + u(0)*u(0)*(1.0-cos(psi))
  R(0,1) = u(0)*u(1)*(1.0-cos(psi)) - u(2)*sin(psi)
  R(0,2) = u(0)*u(2)*(1.0-cos(psi)) + u(1)*sin(psi)
  R(1,0) = u(0)*u(1)*(1.0-cos(psi)) + u(2)*sin(psi)
  R(1,1) = cos(psi) + u(1)*u(1)*(1.0-cos(psi))
  R(1,2) = u(1)*u(2)*(1-cos(psi)) - u(0)*sin(psi)
  R(2,0) = u(0)*u(2)*(1.0-cos(psi)) - u(1)*sin(psi)
  R(2,1) = u(1)*u(2)*(1-cos(psi)) + u(0)*sin(psi)
  R(2,2) = cos(psi) + u(2)*u(2)*(1.0-cos(psi))

  v = x*R
  return v
}
