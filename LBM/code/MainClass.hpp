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

  int x_next;
  int x_prev;
  int y_next;
  int y_prev;
  
  double alpha;
  double beta; 
  double gamma;
  double tol;
  double tau;

  double u_squared;

  Cube<double> f;
  Cube<double> f_star;
  Cube<double> f_eq;
  Cube<double> S;
  Cube<double> u;

  Mat<double> rho;

  Col<double> F;


  MainClass();
  MainClass(int NX, int NY, double TAU, double FX, double FY, double tolerence, string filename, int amount_of_data);
  void initialize(double rho);
  void run();
};

#endif
