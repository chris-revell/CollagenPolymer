#
#  outputData.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

#include "outputData.hpp"
#include "monomer.hpp"
#include <armadillo>
#include <vector>
#include <iostream>
#include <fstream>

function outputData!(const vector<monomer>& Monomers,ofstream& povrayFile,const int& Nmonomers,const float& t,const float& tmax){
  for (int jj=0 jj<Nmonomers jj++){
    for (int ii=0 ii<Monomers[jj].Ndomains*3 ii++){
      povrayFile << Monomers[jj].pos(ii,0) << ", " << Monomers[jj].pos(ii,1) << ", " << Monomers[jj].pos(ii,2) << ", " << Monomers[jj].domainLength/2.0 << ", " << t << endl
    }
  }
  system("clear")
  cout << "Simulating: " << t << "/" << tmax << endl
}
