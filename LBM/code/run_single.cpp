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

    MainClass instance(Nx, Ny, 2, 0.7, pow(10,-7), 0, 5*pow(10, -3), "filename", 5);
    //MainClass instance(Nx, Ny, 2, 0, 0, pow(10, -9), "filename", 5);
    //           (file name, matrix_size)

    instance.set_boundary();
    //instance.boundary_disc(20, 20, 10);

    instance.open();
    instance.initialize(1);
    //instance.initialize_other(2, 2, 5, 2);
    instance.run();
    cout << "equilibrated" << endl;
    instance.initialize_C(40, 40 , 1, 0);
    instance.initialize_C(40, 40 , 1, 1);
    instance.initialize_C(40, 40 , 1, 2);
    instance.initialize_C(40, 40 , 1, 3);
    instance.initialize_C(40, 40 , 1, 4);
    instance.initialize_C(40, 40 , 1, 5);
    instance.initialize_C(40, 40 , 1, 6);
    instance.initialize_C(40, 40 , 1, 7);
    instance.initialize_C(40, 40 , 1, 8);

    //instance.ADE(12000);
    //instance.test_mass_cons();
    //instance.write_u();
    //instance.write_C();
    instance.test_mass_diffusion();
    return 0;
  }
