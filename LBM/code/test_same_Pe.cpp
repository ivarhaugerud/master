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
    int T  = 10000;

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
    instance.test_mass_diffusion();
    instance.run();
    cout << "equilibrated" << endl;

    instance.initialize_C(14, 32 , 0, 100);
    mat C_in = instance.ADE(T);

    instance.clear_g();

    instance.change_D(1.1);
    instance.change_F(1.1*5*pow(10,-7), 0);

    instance.initialize(1);
    instance.run();
    cout << "equilibrated" << endl;

    instance.initialize_C(40, 18 , 0, 100);
    instance.ADE_back(T, C_in);
    instance.write_u();
    return 0;
  }
