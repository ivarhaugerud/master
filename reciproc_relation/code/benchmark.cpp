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
  	int Nx = 10;  //atoi(argv[1]);
  	int Ny = atoi(argv[1]);
    //                 Nx  Ny tau  tau_g               Fx           Fy    tolerance  file_name           datafiles
    MainClass instance(Nx, Ny, 2, 0.50 + 2*pow(10,-4), pow(10,-8), 0, pow(10, -9), "benchmark_testing", 900);

    instance.set_boundary();
    //instance.boundary_disc(48, 128, 8);
    //instance.boundary_disc(16, 128, 8);

    instance.open();
    instance.initialize(1);
    instance.run();

    cout << "equilibrated" << endl;
    instance.write_u("pousielle" + to_string(Ny));
    //instance.test_mass_cons();
    //instance.test_mass_diffusion();
    //instance.initialize_C(10, 25 , 0, 0.001);
    //mat C_in = instance.ADE(T);

    return 0;
  }
