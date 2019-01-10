#include<fstream>
#include"model.h"
#include<string>

using namespace std;

int main(int argc, char **argv){
  int Size, Nstep;
  Size = atoi(argv[1]);
  Nstep = atoi(argv[2]);
  double T;
  string filename;

  fstream file;
  filename = "critical.txt";
  file.open(filename, ios::out);
  file << "Temperature\tField\tEnergy\tMagnetization\tSpecificHeat\tSusceptibility\tM_x"<<endl;
  file.close();

  /* Run a MC simulation at each temperature. */
  /* To derive critical exponents by fitting the data points,
  more points close to the critical are recommended, by means
  of modifying the varying step of temperature. */

  for(T=0.1;T<=5.1;T+=0.005){
    XYModel xy(Size, T);
    xy.Filename = filename;
    xy.MCsteps(Nstep);
  }
}
