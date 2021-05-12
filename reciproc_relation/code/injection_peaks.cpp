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
    int T  = 250000;
    int T_back = 1.25*T;
    MainClass instance(Nx, Ny, 2, 0.50 + 6*pow(10,-5), 5*pow(10,-6), 0, 5*pow(10, -7), "peak_1104_2", 300);

    instance.boundary_disc(12,  13, 7);
    instance.boundary_disc(15,  44, 11);
    instance.boundary_disc(83,  48, 9);
    instance.boundary_disc(72,  10, 9);
    instance.boundary_disc(58,  40, 8);
    instance.boundary_disc(122, 25, 11);

    instance.open();
    instance.initialize(1);
    instance.run();
    instance.write_u("peak_1104");
    instance.change_D(50);

    for (int i = 0; i < Ny; i++)
    {
    instance.define_sources(100, i);
    }

    instance.initialize_C(25, 25 , 0, 0.1);
    mat C_in = instance.ADE(T);
    instance.clear_g();
    C_in = C_in*1000;

    for (int y=0; y < Ny; y++){
        for (int t=0; t < 0.1*T; t++){
            C_in(t, y) = 0;}}

    //no cutoff 
    
    instance.ADE_back(T_back, C_in, "1_0", T);

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

    instance.ADE_back(T_back, C_in, "0_5", T);

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

    instance.ADE_back(T_back, C_in, "0_8", T);

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

    instance.ADE_back(T_back, C_in, "0_9", T);

    //maximum single point
    instance.clear_g();
    double global_heighest = 0;
    int    t_index = 0;
    int    y_index = 0;

    for (int i = 0; i < Ny; i++){
        for (int t = 0; t < T; t++){
            if (C_in(t, i) > global_heighest)
            {global_heighest = C_in(t, i);
             t_index = t;
             y_index = i;}}}

    Mat<double> C_out;
    C_out = Mat<double>(T, Ny);
    C_out(t_index, y_index) = global_heighest;

    instance.ADE_back(T_back, C_out, "single_maxima", T);

    //maxima of each point along line
    instance.clear_g();
    double current_heighest = 0;
    int    index_of_heighest = 0;

    for (int i = 0; i < Ny; i++){
        current_heighest = 0;
        for (int t = 0; t < T; t++){
            if (C_in(t, i) > current_heighest)
            {current_heighest = C_in(t, i);
             index_of_heighest = t;}}

     for (int t = 0; t < T; t++){
        if ( t != index_of_heighest){
            C_in(t, i) = 0;}}
    }

    instance.ADE_back(T_back, C_in, "maxima", T);
    return 0;
  }
