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

  //Number of lattice sites and equilibration counter
  int Nx;
  int Ny;
  int counter;

  //indicies for propegating
  int x_next;
  int x_prev;
  int y_next;
  int y_prev;

  //simulation parameters derived from tau and tau_g, and the external force
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

  //to simplify local equilibrium calculations
  double three_u_squared;

  //advection
  Cube<double> f;
  Cube<double> f_prev;
  Cube<double> f_star;
  Cube<double> u;
  Cube<double> prev_u;

  //diffusion
  Cube<double> g;
  Cube<double> g_star;


  Mat<double> rho; //fluid density spanning domain
  Mat<double> C;   //concentration field
  Col<double> F;   //external force

  int x;
  int y;

  //define the domain into boundary and fluid (rest), in addition to defining the sources/sinks of the system
  vector<tuple<int, int>> boundary;
  vector<tuple<int, int>> rest;
  vector<tuple<int, int>> source;

  MainClass();
  MainClass(int NX, int NY, double TAU, double TAU_G, double FX, double FY, double tolerence, string FILENAME, int amount_of_data);

  void initialize(double rho); //set density in fluid
  

  //NSE and ADE solvers
  void run();                  //find stationary velocity field
  mat ADE(int t);              //simulate advection-diffusion equation
  void ADE_back(int t, mat C_in, string name, int injection_T);  //ade with exactly reversed velocity field
  void ADE_back_no_source(int t, string name);                   //ade with exactly reversed, without sources
  mat ADE_heat(int T, double wall_T, string name);               //ade with diffusion of heat BCs

  //writing data to file
  void write_u(string name);                                     //save velocity field to file
  void write_C(int T, string filename2);                         //save concentration to file

  //change parameters
  void change_D(double D_factor);                                //change the diffusion coefficient
  void change_F(double FX, double FY);                           //change the external force

  //set boundaries
  void set_boundary();                                           //create a boundary around two sides of the domain, e.g. to replicate Poiseuille flow
  void boundary_disc(int x, int y, double R);                    //place disk boundaries in fluid
  void open();                                                   //get indices of boundary points

  //benchmarking
  void test_mass_cons();                                         //verify conservation of fluid mass
  void test_mass_diffusion();                                    //verify conservation of solute mass
  void heat_fluid(double heat_fluid);                            //set a temperature of the fluid
  void define_sources(int x, int y);                             //define indices for sources in the fluid
  void clear_g();                                                //set the concentration to zero for a new run

  //initialize fluid and concentration
  void initialize_other(int x, int y, int i, double rho);        //initialize fluid outside of regular procedure
  void initialize_C(int x, int y, double rho);                   //initialize concentration

  //BCs
  void bounce_back(int x, int y);                                //bounce back BCs
  void anti_bounce_back(int x, int y, double wall_T);            //anti-bounce back BCs

  //streaming of fluid
  void propegate(int x, int y);                                  //propagate fluid

  //streaming of solute
  void update_g();                                               //stream concentration
  void update_g_reversed();                                      //stream concentration exactly reversed
  void update_g_oscillate(double t_rel);                         //stream concentration oscillating velocity field
  void update_g_reversed_oscillate(double t_rel);                //stream reversed oscillating velocity field
};

#endif
