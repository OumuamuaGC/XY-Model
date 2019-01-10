#include<iostream>
#include<fstream>
#include<string>
#include"model.h"

using namespace std;

/* Please input Size and MC steps when running the program. */
int main(int argc, char **argv){
  int Size, Nstep;
  Size = atoi(argv[1]);
  Nstep = atoi(argv[2]);
  double T[11] = {0.1, 0.3, 0.5, 0.7, 1.1, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0};
  /* Low T: 0.1~0.7   Critical T: 1.1   High T: 1.5~10.0 */

  fstream file;
  string filename;
  filename = "field.txt";
  file.open(filename, ios::out);
  file << "Temperature\tField\tEnergy\tMagnetization\tSpecificHeat\tSusceptibility\tM_x"<<endl;
  file.close();

  double h;
  double Hmin = -0.1;
  double Hmax = 0.1;

  /* Magnetic field ascends from minimun to maximum, then descends back to minimum. */
  /* New object is created only before the field-varying circle to realize the evolutional process. */
  for(int i=0;i<11;i++){
    XYModel xy(Size, T[i]);
    xy.Filename = filename;
    for(h=Hmin;h<=Hmax;h+=0.005){
      xy.H = h;
      xy.MCsteps(Nstep);
    }
    for(h=Hmax;h>=Hmin;h-=0.005){
      xy.H = h;
      xy.MCsteps(Nstep);
    }
  }
}
