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
    int T  = 750000;
    MainClass instance(Nx, Ny, 2, 0.50 + 6*pow(10,-5), 5*pow(10,-5), 0, 5*pow(10, -7), "step", 300);

    instance.boundary_disc(12,  13, 7);
    instance.boundary_disc(15,  44, 13);
    instance.boundary_disc(60,  33, 8);
    instance.boundary_disc(83,  53, 12);
    instance.boundary_disc(72,  10, 9);
    instance.boundary_disc(91,  22, 11);
    instance.boundary_disc(122, 25, 14);
    instance.boundary_disc(70, 32, 10);

    instance.open();
    instance.initialize(1);
    instance.run();
    cout << "equilibrated" << endl;
    instance.write_u("step");
    instance.define_sources(108, 45);
    instance.initialize_C(10, 25 , 0, 1);

    mat C_in = instance.ADE(T);

    for (int t=0; t < T; t++){
            C_in(t, 0) = 1;
    }

    instance.clear_g();
    instance.ADE_back(T, C_in);
    instance.write_u("step_back");
    return 0;
  }
