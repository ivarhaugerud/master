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
  double eta;
  double zeta;
  double tol;
  double tau;
  double tau_g;
  double FU;

  double three_u_squared;
  double current_max_u;
  double prev_max_u;

  //advection
  Cube<double> f;
  Cube<double> f_prev;
  Cube<double> f_star;
  Cube<double> f_eq;
  Cube<double> S;
  Cube<double> u;
  Cube<double> prev_u;

  //diffusion
  Cube<double> g;
  Cube<double> g_eq;
  Cube<double> g_star;
  Cube<double> g_prev;

  Mat<double> rho;
  Mat<double> C;
  Col<double> F;

  int x;
  int y;
  vector<tuple<int, int>> boundary;
  vector<tuple<int, int>> rest;

  MainClass();
  MainClass(int NX, int NY, double TAU, double TAU_G, double FX, double FY, double tolerence, string filename, int amount_of_data);
  void initialize(double rho);
  void run();
  void ADE(int t);
  void write_u();
  void write_C();
  void set_boundary();
  void open();  
  void test_mass_cons();
  void test_mass_diffusion();
  void initialize_other(int x, int y, int i, double rho);
  void initialize_C(int x, int y, int i, double rho);
  void boundary_disc(int x, int y, double R);
};

#endif
