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
    int Nx = 140;  //atoi(argv[1]);
    int Ny = 64;   //atoi(argv[2]);
    int T  = 600000;
    MainClass instance(Nx, Ny, 2, 0.50 + 6*pow(10,-5), 5*pow(10,-6), 0, 5*pow(10, -7), "cos_1703", 300);

    instance.boundary_disc(12,  13, 7);
    instance.boundary_disc(15,  44, 11);
    instance.boundary_disc(83,  48, 9);
    instance.boundary_disc(72,  10, 9);
    instance.boundary_disc(58,  40, 8);
    instance.boundary_disc(122, 25, 11);

    instance.open();
    instance.initialize(1);
    instance.run();
    cout << "equilibrated" << endl;
    instance.write_u("cos");

    instance.define_sources(108, 45);
    instance.initialize_C(10, 25 , 0, 0.005);

    mat C_in = instance.ADE(T);

    for (int t=0; t < T; t++){
            C_in(t, 0) = 0.0001*cos(4*2*3.14159*t/T)*cos(4*2*3.14159*t/T);
    }

    instance.clear_g();
    instance.ADE_back(T, C_in, "cos", T);
    instance.write_u("cos_back");
    return 0;
  }
