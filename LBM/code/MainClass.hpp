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
  Cube<double> u;
  Cube<double> prev_u;

  //diffusion
  Cube<double> g;
  Cube<double> g_star;

  Mat<double> rho;
  Mat<double> C;
  Col<double> F;

  int x;
  int y;
  vector<tuple<int, int>> boundary;
  vector<tuple<int, int>> rest;
  vector<tuple<int, int>> source;

  MainClass();
  MainClass(int NX, int NY, double TAU, double TAU_G, double FX, double FY, double tolerence, string FILENAME, int amount_of_data);
  void initialize(double rho);
  void run();
  mat ADE(int t);
  void ADE_back(int t, mat C_in);
  mat ADE_heat(int T, double wall_T, string name);
  void write_u(string name);
  void write_C(int T, string filename2);
  void change_D(double D_factor);
  void change_F(double FX, double FY);

  void set_boundary();
  void heat_fluid(double heat_fluid);
  void define_sources(int x, int y);
  void open();  
  void clear_g();
  void test_mass_cons();
  void test_mass_diffusion();
  void initialize_other(int x, int y, int i, double rho);
  void initialize_C(int x, int y, int i, double rho);
  void boundary_disc(int x, int y, double R);
  void update_g();
  void update_g_reversed();
  void bounce_back(int x, int y);
  void anti_bounce_back(int x, int y, double wall_T);
  void propegate(int x, int y);
};

#endif
