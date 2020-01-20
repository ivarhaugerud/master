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
  	int Nx = 512;//atoi(argv[1]);
  	int Ny = 128;//atoi(argv[2]);

    MainClass instance(Nx, Ny, 2, 0.53 + 3*pow(10,-5), pow(10,-7), 0, 5*pow(10, -3), "filename", 5);
    //instance.set_boundary();
    instance.boundary_disc(40, 70, 30);
    instance.boundary_disc(405, 55, 55);
    instance.boundary_disc(213, 25, 25);

    instance.open();
    instance.initialize(1);
    //instance.initialize_other(2, 2, 5, 2);
    instance.run();
    instance.initialize_C(256, 64 , 0, 100);

    //instance.ADE(5000);
    //instance.test_mass_cons();
    //instance.write_u();
    instance.test_mass_diffusion();
    //instance.write_C();
    return 0;
  }
