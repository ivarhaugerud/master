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
  	int Nx = 102;  //atoi(argv[1]);
  	int Ny = 64;   //atoi(argv[2]);
    int T  = 2000000;

    MainClass instance(Nx, Ny, 2, 0.50 + 2*pow(10,-4), 5*pow(10,-6), 0, 6*pow(10, -4), "1307_reciproc_5_oscls", 900);

    instance.boundary_disc(12,  13, 7);
    instance.boundary_disc(15,  44, 11);
    instance.boundary_disc(83,  48, 9);
    instance.boundary_disc(72,  10, 9);
    instance.boundary_disc(58,  40, 8);
    instance.boundary_disc(16,  25, 4);
    //instance.boundary_disc(68,  25, 6);


    instance.open();
    instance.initialize(1);
    instance.run();
    cout << "equilibrated" << endl;

    instance.initialize_C(10, 25 , 0, 0.001);

    mat C_in = instance.ADE(T);

    instance.clear_g();
    instance.initialize_C(90, 40 , 0, 0.001);
    instance.ADE_back_no_source(T, "back");
    return 0;
  }
