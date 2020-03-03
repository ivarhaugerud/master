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
    int T  = 230000;
    int T_back = 1.2*T
    MainClass instance(Nx, Ny, 2, 0.50 + 6*pow(10,-5), 5*pow(10,-6), 0, 5*pow(10, -7), "peak", 300);

    instance.boundary_disc(12,  13, 7);
    instance.boundary_disc(15,  44, 11);
    instance.boundary_disc(83,  48, 9);
    instance.boundary_disc(72,  10, 9);
    instance.boundary_disc(58,  40, 8);
    instance.boundary_disc(122, 25, 11);

    instance.open();
    instance.initialize(1);
    instance.run();
    instance.write_u("peak");

    for (int i = 0; i < Ny; i++)
    {
    instance.define_sources(100, i);
    }

    instance.initialize_C(25, 25 , 0, 100);
    mat C_in = instance.ADE(T);

    //no cutoff 
    instance.ADE_back(T_back, C_in, "1_0");

    //cutoff 0.5
    instance.clear_g();

    double max_val = C_in.max();
    double cut_off_portion = 0.5;
    double cut_off = max_val*cut_off_portion;

    for (int y=0; y < Ny; y++){
        for (int t=0; t < T; t++){
            if (C_in(t, y) < cut_off){
                C_in(t, y) = 0;
            }
        }
    }

    instance.ADE_back(T_back, C_in, "0_5");

    //cutoff 0.8
    instance.clear_g();

    max_val = C_in.max();
    cut_off_portion = 0.8;
    cut_off = max_val*cut_off_portion;

    for (int y=0; y < Ny; y++){
        for (int t=0; t < T; t++){
            if (C_in(t, y) < cut_off){
                C_in(t, y) = 0;
            }
        }
    }

    instance.ADE_back(T_back, C_in, "0_8");

    //cutoff 0.9
    instance.clear_g();

    max_val = C_in.max();
    cut_off_portion = 0.9;
    cut_off = max_val*cut_off_portion;

    for (int y=0; y < Ny; y++){
        for (int t=0; t < T; t++){
            if (C_in(t, y) < cut_off){
                C_in(t, y) = 0;
            }
        }
    }

    instance.ADE_back(T_back, C_in, "0_9");

    // exponential
    instance.clear_g();
    double total_sum = 0;
    for (int y=0; y < Ny; y++){
        for (int t=0; t < T; t++){
            if (C_in(t, y) < cut_off){
                C_in(t, y) = exp(C_in(t, y));
                total_sum += C_in(t, y);
            }
        }
    }
    for (int y=0; y < Ny; y++){
        for (int t=0; t < T; t++){
            C_in(t, y) = 100*C_in(t, y)/total_sum;
        }}

    instance.ADE_back(T_back, C_in, "exponential");
    return 0;
  }
