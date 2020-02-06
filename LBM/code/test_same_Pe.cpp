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

    MainClass instance(Nx, Ny, 2, 0.50 + 6*pow(10,-5), 5*pow(10,-7), 0, 5*pow(10, -5), "0602reciproc", 100);

    //instance.boundary_disc(12,  13, 7);
    //instance.boundary_disc(15,  44, 7);
    instance.boundary_disc(60,  33, 8);
    //instance.boundary_disc(83,  53, 7);
    //instance.boundary_disc(72,  10, 9);
    //instance.boundary_disc(91,  22, 7);
    //instance.boundary_disc(122, 25, 7);

    instance.open();
    instance.initialize(1);
    instance.run();
    cout << "equilibrated" << endl;
    instance.write_u();

    /*
    instance.initialize_C(32, 32 , 0, 100);
    mat C_in = instance.ADE(60000);


    instance.clear_g();
    instance.initialize_C(34, 34 , 0, 100);

    double D_factor = 1.03;
    instance.change_D(D_factor);
    instance.change_F(D_factor*5*pow(10,-6), 0);

    instance.ADE_back(60000, C_in);
    */
    return 0;
  }