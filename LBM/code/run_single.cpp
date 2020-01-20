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

    MainClass instance(Nx, Ny, 2, 0.51 + 3*pow(10,-3), pow(10,-7), 0, 5*pow(10, -5), "filename", 10);
    //instance.set_boundary();
    instance.boundary_disc(12,  13, 10);
    instance.boundary_disc(15,  44, 10);
    instance.boundary_disc(60,  33, 10);
    instance.boundary_disc(83,  53, 10);
    instance.boundary_disc(72,  10, 10);
    instance.boundary_disc(91,  22, 10);
    instance.boundary_disc(122, 25, 10);
    instance.boundary_disc(180, 13, 10);
    instance.boundary_disc(213, 42, 10);
    instance.boundary_disc(234, 17, 10);
    instance.boundary_disc(199, 20, 10);
    instance.boundary_disc(212, 31, 10);

    instance.open();

    instance.initialize(1);
    instance.run();

    instance.initialize_C(32, 0 , 0, 100);

    instance.ADE(50000);
    instance.write_u();
    return 0;
  }
