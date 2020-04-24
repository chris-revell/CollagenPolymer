#
#  outputData.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module OutputData

function outputData!(pos,povrayFile,Nmonomers,t,tmax)
  for (int jj=0 jj<Nmonomers jj++){
    for (int ii=0 ii<Monomers[jj].Ndomains*3 ii++){
      povrayFile << Monomers[jj].pos(ii,0) << ", " << Monomers[jj].pos(ii,1) << ", " << Monomers[jj].pos(ii,2) << ", " << Monomers[jj].domainLength/2.0 << ", " << t << endl
    }
  }
  cout << "Simulating: " << t << "/" << tmax << endl
  if (t%outputInterval)<dt
      writedlm(outfile,pos,", ")
      Printf.@printf("Simulating: %f/%f\n",t,tmax)
  end
}
