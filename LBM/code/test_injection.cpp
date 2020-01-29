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
  	int Nx = 64;  //atoi(argv[1]);
  	int Ny = 64;   //atoi(argv[2]);

    MainClass instance(Nx, Ny, 2, 0.50 + 6*pow(10,-3), 5*pow(10,-6), 0, 5*pow(10, -5), "filename", 100);

    instance.boundary_disc(20,  33, 8);
    instance.boundary_disc(43,  50, 12);
    instance.boundary_disc(32,  10, 9);
    instance.boundary_disc(51,  22, 11);

    instance.open();
    instance.initialize(1);
    instance.run();
    cout << "equilibrated" << endl;

    instance.initialize_C(35, 32 , 0, 100);

    for (int i = 0; i < Ny; i++)
    {
    instance.define_sources(0, i);
    }

    mat C_in = instance.ADE(100);
    instance.clear_g();

    instance.ADE_back(100, C_in);
    instance.write_u();
    return 0;
  }
