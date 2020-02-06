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
    int Nx = 150;  //atoi(argv[1]);
    int Ny = 64;   //atoi(argv[2]);
    int T  = 400000;
    double tau_g = 0.50 + 6*pow(10,-5);

    MainClass instance(Nx, Ny, 2, tau_g, 5*pow(10,-6), 0, 5*pow(10, -7), "060220_to", 200);

    instance.boundary_disc(134,  14, 5);
    instance.boundary_disc(33,   44, 7);
    instance.boundary_disc(115,  27, 9);
    instance.boundary_disc(16,   19, 3);
    instance.boundary_disc(58,   16, 4);
    instance.boundary_disc(101,  13, 5);
    instance.boundary_disc(42,   37, 6);
    instance.boundary_disc(130,  26, 7);

    instance.open();
    instance.initialize(1);
    instance.run();
    cout << "equilibrated" << endl;

    instance.initialize_C(14, 26 , 0, 100);
    mat C_in = instance.ADE(T);

    double D_factor = 1.1;
    double new_tau_g = tau_g/D_factor + 0.5*(D_factor - 1);
    MainClass instance_back(Nx, Ny, 2, new_tau_g, -D_factor*5*pow(10,-6), 0, 5*pow(10, -7), "060220_from", 200);

    instance_back.boundary_disc(134,  14, 5);
    instance_back.boundary_disc(33,   44, 7);
    instance_back.boundary_disc(115,  27, 9);
    instance_back.boundary_disc(16,   19, 3);
    instance_back.boundary_disc(58,   16, 4);
    instance_back.boundary_disc(101,  13, 5);
    instance_back.boundary_disc(42,   37, 6);
    instance_back.boundary_disc(130,  26, 7);

    instance_back.open();
    instance_back.initialize(1);
    instance_back.run();
    cout << "equilibrated" << endl;

    instance_back.initialize_C(70, 50 , 0, 100);
    mat C_in_2 = instance_back.ADE(T);
    instance_back.write_u();
    return 0;
  }
