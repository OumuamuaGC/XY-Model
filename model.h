#include<stdlib.h>
#include<cmath>
#include<ctime>
#include<iostream>
#include<fstream>
#include<string>

using namespace std;

#define PI 3.1416

/* Size is not allowed to exceed MAXLen. */
/* The angles of spin vector are set as 72 dicrete values
to reduce degrees of freedom, i.e. {0, 5, ..., 355}, so that
the system can reach its quilibrium in fewer MC steps. */
#define MAXLen 129
#define Angles 72

/* J is the exchange interaction parameter, positive for ferromagnetism,
negative for anti-ferromagnetism. */
/* cosls[], sinls[] make use of arrays to store the cos/sin values,
with the angles as its index, instead of calculating them every time. */
class XYModel{
public:
  int Size;
  double Temp;
  double H;
  double J;
  double cosls[Angles];
  double sinls[Angles];
  int grid[MAXLen][MAXLen];
  double Energy;
  double Mfrac[2];
  string Filename;

  XYModel(int size, double t);
  void MCsteps(int N);

private:
  double Esite(int theta, int x, int y);
  void Flipper(int x, int y);
  void Collector(int i, double *Es, double *Ms, double *Esqrs, double *Msqrs, double *Mxs);
  void Calculator(int Nsam, double Es[], double Ms[], double Esqrs[], double Msqrs[], double Mxs[], double *output);

};
