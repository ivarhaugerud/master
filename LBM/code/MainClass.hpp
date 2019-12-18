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
#include <tuple>
#include <vector>

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
  int counter;

  int x_next;
  int x_prev;
  int y_next;
  int y_prev;

  double alpha;
  double beta; 
  double gamma;
  double delta;
  double tol;
  double tau;
  double FU;

  double u_squared;
  double current_max_u;
  double prev_max_u;

  Cube<double> f;
  Cube<double> f_prev;
  Cube<double> f_star;
  Cube<double> f_eq;
  Cube<double> S;
  Cube<double> u;
  Cube<double> prev_u;

  Mat<double> rho;
  Col<double> F;

  int x;
  int y;
  vector<tuple<int, int>> boundary;

  MainClass();
  MainClass(int NX, int NY, double TAU, double FX, double FY, double tolerence, string filename, int amount_of_data);
  void initialize(double rho);
  void run();
  void write_u();
  void set_boundary();
  void test_mass_cons();
  void initialize_other(int x, int y, int i, double rho);
  void boundary_disc(int x, int y, double R)
};

#endif
