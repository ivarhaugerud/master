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

    MainClass instance(Nx, Ny, 2, 0.50 + 3*pow(10,-3), pow(10,-7), 0, 5*pow(10, -3), "filename", 5);
    //instance.set_boundary();
    instance.boundary_disc(20, 20, 10);
    //instance.boundary_disc(60, 40, 30);
    //instance.boundary_disc(37, 60, 20);

    instance.open();
    instance.initialize(1);
    //instance.initialize_other(2, 2, 5, 2);
    instance.run();
    cout << "equilibrated" << endl;
    instance.initialize_C(1, 1 , 5, 5);

    //instance.ADE(12000);
    //instance.test_mass_cons();
    //instance.write_u();
    instance.test_mass_diffusion();
    instance.write_C();
    return 0;
  }
