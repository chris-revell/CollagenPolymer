#
#  readParameters.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

using namespace std

# Read simulation parameters from the command line.
void readParameters(ofstream& positionFile,char* buffer,int argc,char *argv[],int& Nmonomers,int& Ndomains,float& tmax,float& dt,float& outputInterval,float& gamma,float& zetaMag,float& k,float& Ebend,float& epsilon,int& visualiseFlag){

	ofstream conditionsFile

  std::time_t result = std::time(nullptr)
  struct tm * timeinfo
  timeinfo = localtime (&result)
  # Create directory to store simulation results
  strftime (buffer,27, "output/%F-%H-%M-%S\0",timeinfo)
  system(("mkdir "+string(buffer)).c_str())
  positionFile.open((string(buffer)+"/positionData.txt").c_str(),fstream::out)
  conditionsFile.open((string(buffer)+"/conditionsData.txt").c_str(),fstream::out)

	# Read from command line
	Nmonomers							= atoi(argv[1])
  Ndomains 							= atoi(argv[2])
  tmax                  = atof(argv[3])
  dt                    = atof(argv[4])
	outputInterval				= atof(argv[5])
  gamma                 = atof(argv[6])
  zetaMag               = atof(argv[7])
  k                     = atof(argv[8])
	Ebend									= atof(argv[9])
  epsilon               = atof(argv[10])
  visualiseFlag 				= atoi(argv[11])

	conditionsFile << Nmonomers << endl
	conditionsFile << Ndomains << endl
	conditionsFile << tmax << endl
	conditionsFile << dt << endl
	conditionsFile << outputInterval << endl
  conditionsFile << gamma << endl
  conditionsFile << zetaMag << endl
  conditionsFile << k << endl
	conditionsFile << Ebend << endl
  conditionsFile << epsilon << endl
  conditionsFile << visualiseFlag << endl
	conditionsFile.close()

}
