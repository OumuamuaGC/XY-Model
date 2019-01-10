#include"model.h"

using namespace std;

/* Construction function to initialize the object. */
XYModel::XYModel(int size, double temp){
  Size = size;
  Temp = temp;

  /* Default values: ferromagnetism, no external field. */
  H = 0.0;
  J = 1.0;

  /* Default name of output file, allowed to rename in "int main()". */
  Filename = "output.txt";

  Energy = 0.0;
  Mfrac[0] = 0.0;
  Mfrac[1] = 0.0;

  srand((unsigned)time(NULL));

  /* All the initial directions are "+x". */
  /* Random directions can be realized easily with "rand()%Angles". */
  for(int i=0;i<MAXLen;i++){
    for(int j=0;j<MAXLen;j++){
      grid[i][j] = 0;
      //grid[i][j] = rand() % Angles;
    }
  }

  for(int i=0;i<Angles;i++){
    cosls[i] = cos(i * (360/Angles) * PI / 180.0);
    sinls[i] = sin(i * (360/Angles) * PI / 180.0);
  }

  /* Esite includes both interaction energy and magnetic energy, while
  the former will contribute twice to the total energy. */
  for(int i=0;i<Size;i++){
    for(int j=0;j<Size;j++){
      Energy += (Esite(grid[i][j], i, j) - H * cosls[grid[i][j]])/2.0;
      Mfrac[0] += cosls[grid[i][j]];
      Mfrac[1] += sinls[grid[i][j]];
    }
  }
}

/* Calculate the interaction energy and magnetic energy for each site. */
double XYModel::Esite(int theta, int x, int y){
  double e = 0;
  // Four adjacent spin interactions.
  e += - J * cosls[abs(theta - grid[x][(y + 1 + Size) % Size])];
  e += - J * cosls[abs(theta - grid[x][(y - 1 + Size) % Size])];
  e += - J * cosls[abs(theta - grid[(x + 1 + Size) % Size][y])];
  e += - J * cosls[abs(theta - grid[(x - 1 + Size) % Size][y])];
  // Magnetic energy.
  e += - H * cosls[theta];
  return e;
}

/* Collect the instant thermodynamics quantities at the ith MC step. */
/* All the quantities are averaged to each sites. */
void XYModel::Collector(int i, double *Es, double *Ms, double *Esqrs, double *Msqrs, double *Mxs){
  double E1, M1, M1x;
  E1 = Energy / double(Size * Size);
  M1 = sqrt(Mfrac[0] * Mfrac[0] + Mfrac[1] * Mfrac[1]) / double(Size * Size);
  M1x = Mfrac[0] / double(Size * Size);
  Es[i] = E1;
  Ms[i] = M1;
  Esqrs[i] = E1 * E1;
  Msqrs[i] = M1 * M1;
  Mxs[i] = M1x;
}

/* Calculate the mean values of various thermodynamics quantities after the whole MC simulation. */
void XYModel::Calculator(int Nsam, double Es[], double Ms[], double Esqrs[], double Msqrs[], double Mxs[], double *output){
  double Esig, Msig, Esqrsig, Msqrsig, Mxsig;
  Esig = 0.0, Msig = 0.0, Esqrsig = 0.0, Msqrsig = 0.0, Mxsig = 0.0;

  for(int i=0;i<Nsam;i++){
    Esig += Es[i];
    Msig += Ms[i];
    Esqrsig += Esqrs[i];
    Msqrsig += Msqrs[i];
    Mxsig += Mxs[i];
  }

  output[0] = Esig / double(Nsam);
  output[1] = Msig / double(Nsam);
  output[2] = ((Esqrsig / double(Nsam)) - (Esig / double(Nsam)) * (Esig / double(Nsam))) / (Temp * Temp);
  output[3] = ((Msqrsig / double(Nsam)) - (Msig / double(Nsam)) * (Msig / double(Nsam))) / Temp;
  output[4] = Mxsig / double(Nsam);
}

/* Decide to flip or not. Update the angle, the total energy
and the total magnetization if flip. */
/* Flipping is based on Metropolis Algorithm, please modify this
function if different algorithm needs to be applied. */
void XYModel::Flipper(int x, int y){
  int flipper, flip;
  flipper = rand() % Angles;
  flip  = (grid[x][y] + flipper) % Angles;

  double dE, prob;
  dE = Esite(flip, x, y) - Esite(grid[x][y], x, y);
  if(dE < 0){
    Mfrac[0] += cosls[flip] - cosls[grid[x][y]];
    Mfrac[1] += sinls[flip] - sinls[grid[x][y]];
    grid[x][y] = flip;
    Energy += dE;
  }
  else{
    prob = rand() / double(RAND_MAX);
    if(prob <= exp(- dE/Temp)){
      Mfrac[0] += cosls[flip] - cosls[grid[x][y]];
      Mfrac[1] += sinls[flip] - sinls[grid[x][y]];
      grid[x][y] = flip;
      Energy += dE;
    }
  }
}

/* Run N MC steps. */
void XYModel::MCsteps(int N){
  int x, y;
  int Nsample;
  Nsample = N / 10;
  /* Calculate the averages with the data from the last N/10 steps. */

  double Es[Nsample]; double Ms[Nsample];
  double Esqrs[Nsample]; double Msqrs[Nsample];
  double Mxs[Nsample];
  /* Collected data are stored into arrays. */

  double output[5]; /* Results of 5 thermodynamics quantities. */

  for(int i=0;i<N;i++){
    x = rand()%Size;
    y = rand()%Size;
    Flipper(x, y);
    if(i >= (N - Nsample)){
      Collector(i - N + Nsample, Es, Ms, Esqrs, Msqrs, Mxs);
    }
  }
  Calculator(Nsample, Es, Ms, Esqrs, Msqrs, Mxs, output);
  cout << "Temperature\tField\tEnergy\tMagnetization\tSpecificHeat\tSusceptibility\tM_x"<<endl;
  cout << Temp <<"\t"<< H <<"\t"<< output[0] <<"\t"<< output[1] <<"\t"<< output[2] <<"\t"<< output[3] <<"\t"<< output[4] <<endl;
  fstream file;
  file.open(Filename, ios::app);
  file << Temp <<"\t"<< H <<"\t"<< output[0] <<"\t"<< output[1] <<"\t"<< output[2] <<"\t"<< output[3] <<"\t"<< output[4] <<endl;
  file.close();
}
