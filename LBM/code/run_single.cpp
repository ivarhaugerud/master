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
    MainClass instance(10, 10, 2, 0, 0, 0.000001, "filename", 5);
    //           (file name, matrix_size)

    instance.initialize(1);
    instance.run();
    return 0;
  }
