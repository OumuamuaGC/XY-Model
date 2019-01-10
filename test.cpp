#include"model.h"

using namespace std;

/* Test one MC simulation, where Size, temperature and MC steps need to be input. */
int main(int argc, char **argv){
  int Size, Nstep;
  double T;

  Size = atoi(argv[1]);
  T = atof(argv[2]);
  Nstep = atoi(argv[3]);

  XYModel xy(Size, T);
  xy.MCsteps(Nstep);
}
