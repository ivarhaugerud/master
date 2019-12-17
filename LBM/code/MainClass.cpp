#include "MainClass.hpp"

MainClass::MainClass()
{}

//all variabels we neeed for the class, in addition to the two functions for the energy and wave function
MainClass::MainClass(int NX, int NY, double TAU, double FX, double FY, double tolerence, string filename, int amount_of_data)

{
  //the file name for saving data, and number of lines when saving data
  filename   = filename;
  data_lines = amount_of_data;

  Nx = NX;
  Ny = NY;
  tau = TAU;

  f       = Cube<double>(Nx, Ny, 9);
  f_prev  = Cube<double>(Nx, Ny, 9);
  f_equil = Cube<double>(Nx, Ny, 9);
  S       = Cube<double>(Nx, Ny, 9);

  velocity = Mat<double>(Nx, Ny);
  density  = Mat<double>(Nx, Ny);
  c        = Mat<double>(2, 9);

  F        = Col<double>(2);
  omega    = Col<double>(9);

  alpha = 1 - 1/tau;
  beta  = 1/tau;
  gamma = 3*(1 - 1/(2*tau));
  tol = tolerence;

  omega(0) = 4/9;
  omega(1) = 1/9;
  omega(2) = 1/9;
  omega(3) = 1/9;
  omega(4) = 1/9;
  omega(5) = 1/36;
  omega(6) = 1/36;
  omega(7) = 1/36;
  omega(8) = 1/36;

  c(1, 1) =  1;
  c(0, 2) =  1;
  c(1, 3) = -1;
  c(0, 4) = -1;
  c(0, 5) =  1;
  c(1, 5) =  1;
  c(1, 6) = -1;
  c(0, 6) =  1;
  c(0, 7) = -1;
  c(1, 7) = -1;
  c(1, 8) =  1;
  c(0, 8) = -1;
}

//the metropolis algorithm
void MainClass::initialize(double rho)
{
  for (int i = 0; i < Nx; i ++)
  {
    for (int j = 0; i < Ny; j ++)
    {
      f(i, j, 0) = 4*rho/9;
      f(i, j, 1) = rho/9;
      f(i, j, 2) = rho/9;
      f(i, j, 3) = rho/9;
      f(i, j, 4) = rho/9;
      f(i, j, 5) = rho/36;
      f(i, j, 6) = rho/36;
      f(i, j, 7) = rho/36;
      f(i, j, 8) = rho/36;
    }
  }
}

