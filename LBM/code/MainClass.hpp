#ifndef MAINCLASS
#define MAINCLASS

//import what we need
#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <string>
#include <iomanip>
#include <fstream>

//make it easier
using namespace std;
using namespace arma;

class MainClass
{
//the public variables
public:

  //writing to file
  string filename;
  int data_lines;

  int Nx;
  int Ny; 
  double alpha;
  double beta; 
  double gamma;
  double tol;
  double tau;

 Cube<double> f;
 Cube<double> f_prev;
 Cube<double> f_equil;
 Cube<double> S;

  Col<double> F;
  Col<double> omega;

  Mat<double> c;
  Mat<double> velocity;
  Mat<double> density;

  MainClass();
  MainClass(int NX, int NY, double TAU, double FX, double FY, double tolerence, string filename, int amount_of_data);

  //running and equilibrating
  void initialize(double rho);
};

#endif
