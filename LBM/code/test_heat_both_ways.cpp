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
    double T_inject = 1;
    double T_background = 0.001*T_inject;
    int Time = 6000000;
    MainClass instance(Nx, Ny, 2, 0.50 + 6*pow(10,-5), 3.5*pow(10,-7), 0, 1.81*pow(10, -5), "0303heat_F10-7", 300);

    instance.boundary_disc(12,  13, 7);
    instance.boundary_disc(15,  44, 11);
    instance.boundary_disc(83,  53, 12);
    instance.boundary_disc(72,  10, 9);
    instance.boundary_disc(58,  40, 8);
    instance.boundary_disc(122, 25, 11);

    /*
    instance.open();
    instance.initialize(1);
    instance.run();

    instance.heat_fluid(T_background);
    instance.initialize_C(32, 32 , 0, T_inject);

    mat C_in = instance.ADE_heat(Time, T_background, "to");
    instance.write_u("heat");
    */
    instance.open();
    instance.initialize(1);

    instance.change_F(-3.5*pow(10,-7), 0);
    instance.run();

    instance.clear_g();
    instance.heat_fluid(T_background);
    instance.initialize_C(100, 36 , 0, T_inject);
    instance.ADE_heat(Time, T_background, "from");

    return 0;
  }