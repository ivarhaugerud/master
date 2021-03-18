#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdio>
#include <random>
#include <ctime>
#include <string>
#include <typeinfo>

#include "MainClass.hpp" //the class

using namespace std;
using namespace arma;

//simple run function for general input
int main(int argc, char const *argv[])
  {
  	int Nx = 128;  //atoi(argv[1]);
  	int Ny = 128;   //atoi(argv[2]);
    //                 Nx  Ny tau  tau_g               Fx           Fy    tolerance  file_name           datafiles
    MainClass instance(Nx, Ny, 2, 0.50 + 2*pow(10,-4), 5*pow(10,-8), 0, pow(10, -7), "benchmark_testing", 900);

    //instance.set_boundary();
    instance.boundary_disc(64, 64, 16);
    
    instance.open();
    instance.initialize(1);
    instance.run();

    cout << "equilibrated" << endl;
    instance.write_u("stokes");
    instance.test_mass_cons();
    instance.test_mass_diffusion();
    //instance.test_pousielle();
    //instance.initialize_C(10, 25 , 0, 0.001);
    //mat C_in = instance.ADE(T);

    return 0;
  }
