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
    int Nx = 512;
    int Ny = 128;
    int No = 17;
    //                 Nx  Ny tau  tau_g               Fx           Fy    tolerance  file_name           datafiles
    MainClass instance(Nx, Ny, 2, 0.50 + 2*pow(10,-4), 5*pow(10,-8), 0, pow(10, -9), "benchmark_testing", 900);
    

    for (int n=0; n < No; ++n){
        int R = rand()%30;
        instance.boundary_disc(rand()%(Nx-R),  rand()%(Ny-R), R);}
    
    instance.boundary_disc(32, 32, 10);
    instance.open();
    instance.initialize(1);
    instance.run();

    cout << "equilibrated" << endl;
    instance.write_u("benchmark_pousielle");

    //instance.test_mass_cons();
    //instance.test_mass_diffusion();
    //instance.test_pousielle();
    //instance.initialize_C(10, 25 , 0, 0.001);
    //mat C_in = instance.ADE(T);

    return 0;
  }
