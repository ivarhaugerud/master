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

    MainClass instance(Nx, Ny, 2, 5*pow(10, -8), 0, pow(10, -9), "filename", 5);
    //           (file name, matrix_size)

    instance.initialize(1);
    instance.set_boundary();
    instance.run();
    instance.write_u();
    return 0;
  }
