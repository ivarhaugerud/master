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
  	int Nx = atoi(argv[1]);
  	int Ny = atoi(argv[2]);

    MainClass instance(Nx, Ny, 2, 5*pow(10,-8), 0, 5*pow(10, -3), "filename", 5);
    //MainClass instance(Nx, Ny, 2, 0, 0, pow(10, -9), "filename", 5);
    //           (file name, matrix_size)

    instance.set_boundary();
    instance.boundary_disc(20, 45, 10);
    instance.boundary_disc(40, 90, 20);
    instance.boundary_disc(65, 135, 14);
    instance.boundary_disc(80, 180, 10);
    instance.boundary_disc(100, 25, 14);
    instance.boundary_disc(120, 65, 10);
    instance.boundary_disc(135, 115, 14);
    instance.boundary_disc(160, 140, 10);

    instance.open();
    instance.initialize(1);
    //instance.initialize_other(2, 2, 5, 2);
    instance.run();
    //instance.test_mass_cons();
    instance.write_u();
    return 0;
  }
