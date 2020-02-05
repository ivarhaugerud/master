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
    int Nx = 256;  //atoi(argv[1]);
    int Ny = 64;   //atoi(argv[2]);

    MainClass instance(Nx, Ny, 2, 0.50 + 6*pow(10,-5), 5*pow(10,-7), 0, 5*pow(10, -7), "050220", 100);

    
    instance.boundary_disc(134,  14, 5);
    instance.boundary_disc(240,  37, 6);
    instance.boundary_disc(33,   44, 7);
    instance.boundary_disc(202,  41, 8);
    instance.boundary_disc(115,  27, 9);
    instance.boundary_disc(201,  28, 4);
    instance.boundary_disc(16,   19, 3);
    instance.boundary_disc(58,   16, 4);
    instance.boundary_disc(101,  13, 5);
    instance.boundary_disc(42,   37, 6);
    instance.boundary_disc(130,  26, 7);
    instance.boundary_disc(187,  19, 8);
    instance.boundary_disc(192,  48, 9);
    instance.boundary_disc(175,  52, 4);
    
    instance.open();
    instance.initialize(1);
    instance.test_mass_cons();
    instance.run();
    instance.test_mass_diffusion();
    return 0;
  }