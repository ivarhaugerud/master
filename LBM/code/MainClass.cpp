#include "MainClass.hpp"

MainClass::MainClass()
{}

//all variabels we neeed for the class, in addition to the two functions for the energy and wave function
MainClass::MainClass(int NX, int NY, double TAU, double TAU_G, double FX, double FY, double tolerence, string FILENAME, int amount_of_data)
{
  //the file name for saving data, and number of lines when saving data
  filename   = FILENAME;
  data_lines = amount_of_data;

  Nx  = NX;
  Ny  = NY;
  tau = TAU;
  tau_g = TAU_G;

  //advection
  f       = Cube<double>(Nx, Ny, 9);
  f_prev  = Cube<double>(Nx, Ny, 9);
  f_star  = Cube<double>(Nx, Ny, 9);
  u       = Cube<double>(Nx, Ny, 2);
  prev_u  = Cube<double>(Nx, Ny, 2);

  //diffusion
  g       = Cube<double>(Nx, Ny, 9);
  g_star  = Cube<double>(Nx, Ny, 9);

  rho      = Mat<double>(Nx, Ny);
  C        = Mat<double>(Nx, Ny);
  F        = Col<double>(2);

  //use delta_t = 1
  alpha = 1 - 1/tau;
  beta  = 1/tau;
  gamma = 3*(1 - 1/(2*tau));
  delta = (1 - 1/(2*tau));
  eta   = 1-1/tau_g;
  zeta  = 1/tau_g;
  tol   = tolerence;

  F(0) = FX;
  F(1) = FY;
}

void MainClass::initialize(double rho)
{
  {for (int k = 0; k < rest.size(); k++){
    x = get<0>(rest[k]);
    y = get<1>(rest[k]);
    f(x, y, 0) = 4*rho/9;
    f(x, y, 1) = rho/9;
    f(x, y, 2) = rho/9;
    f(x, y, 3) = rho/9;
    f(x, y, 4) = rho/9;
    f(x, y, 5) = rho/36;
    f(x, y, 6) = rho/36;
    f(x, y, 7) = rho/36;
    f(x, y, 8) = rho/36;
    }
  }
}

void MainClass::initialize_other(int x, int y, int i, double rho)
{f(x, y, i) = rho;}

void MainClass::initialize_C(int x, int y, int i, double rho)
{
  g(x, y, i) = rho;
}

void MainClass::change_D(double D_factor)
{
  double new_tau_g = 0.5 + D_factor*(tau_g - 0.5);
  cout << new_tau_g << " " << tau_g << endl;
  eta   = 1-1/new_tau_g;
  zeta  = 1/new_tau_g;

}

void MainClass::change_F(double FX, double FY)
{
  F(0) = FX;
  F(1) = FY;
}

void MainClass::set_boundary()
{
  for (int i = 0; i < Nx; i ++)
  {
    boundary.emplace_back(i, 0);
    boundary.emplace_back(i, Ny-1);
  }
}

void MainClass::define_sources(int x, int y)
{
  source.emplace_back(x, y);
}

void MainClass::boundary_disc(int x, int y, double R)
{
  bool in_boundary;
  int X;
  int Y;

  for (int i = 0; i < Nx; i++){
  for (int j = 0; j < Ny; j++){
    in_boundary = false;

    if (sqrt( (x-i)*(x-i) + (y-j)*(y-j) ) < R){
        for (int k = 0; k < boundary.size(); k++){
          X = get<0>(boundary[k]);
          Y = get<1>(boundary[k]);

          if ((i == X) && (j == Y)){
            in_boundary = true;
            continue;}
          }
        if (not in_boundary)
        {boundary.emplace_back(i,j);}}}}
}

//gets indices of all lattice points not in a boundary
void MainClass::open()
{
  bool in_boundary;

  for (int i = 0; i < Nx; i++){
  for (int j = 0; j < Ny; j++){
    in_boundary = false;
    for (int k = 0; k < boundary.size(); k++){
       x = get<0>(boundary[k]);
       y = get<1>(boundary[k]);

       if ((i == x) && (j == y)){
        in_boundary = true;
        continue;}}

    if (not in_boundary)
    {rest.emplace_back(i,j);} 
}}}

void MainClass::clear_g()
{
  for (int x = 0; x < Nx; x++)
    {for (int y = 0; y < Ny; y++)
      {g(x,y,0) = 0;
       g(x,y,1) = 0;
       g(x,y,2) = 0;
       g(x,y,3) = 0;
       g(x,y,4) = 0;
       g(x,y,5) = 0;
       g(x,y,6) = 0;
       g(x,y,7) = 0;
       g(x,y,8) = 0;

       C(x,y) = 0;
      }
  }
}

void MainClass::heat_fluid(double wall_T)
{
  for (int s = 0; s < rest.size(); s++)
    {x = get<0>(rest[s]);
     y = get<1>(rest[s]);
     initialize_C(x, y, 0, 4*wall_T/9);
     initialize_C(x, y, 1,   wall_T/9);
     initialize_C(x, y, 2,   wall_T/9);
     initialize_C(x, y, 3,   wall_T/9);
     initialize_C(x, y, 4,   wall_T/9);
     initialize_C(x, y, 5,   wall_T/36);
     initialize_C(x, y, 6,   wall_T/36);
     initialize_C(x, y, 7,   wall_T/36);
     initialize_C(x, y, 8,   wall_T/36);
     }
}

void MainClass::run()
{
  bool equil = false;
  int counter = 0;
  double sum_difference = 0;
  double sum = 0;
  //double current_max_u;

  while (not equil)
  //for (int t = 0; t < 20000; t++)
    {for (int k = 0; k < rest.size(); k++)
      {x = get<0>(rest[k]);
       y = get<1>(rest[k]);

      rho(x, y)  = f(x, y, 0) + f(x, y, 1) + f(x, y, 2) + f(x, y, 3) + f(x, y, 4) + f(x, y, 5) + f(x, y, 6) + f(x, y, 7) + f(x, y, 8); 
      u(x, y, 0) = (f(x, y, 1) - f(x, y, 3) + f(x, y, 5) - f(x, y, 6) - f(x, y, 7) + f(x, y, 8))/rho(x,y) + F(0)/(2*rho(x,y));
      u(x, y, 1) = (f(x, y, 2) - f(x, y, 4) + f(x, y, 5) + f(x, y, 6) - f(x, y, 7) - f(x, y, 8))/rho(x,y) + F(1)/(2*rho(x,y));
      three_u_squared = 3*u(x, y, 0)*u(x, y, 0) + 3*u(x, y, 1)*u(x, y, 1);
      FU = 3*(F(0)*u(x,y,0) + F(1)*u(x,y,1));

      sum_difference += (u(x,y,0)-prev_u(x,y,0))*(u(x,y,0)-prev_u(x,y,0)) + (u(x,y,1)-prev_u(x,y,1))*(u(x,y,1)-prev_u(x,y,1));
      sum            += u(x,y,0)*u(x,y,0) + u(x,y,1)*u(x,y,1);

      f(x, y, 0) = f(x, y, 0)*alpha + beta*rho(x, y)*(2 - three_u_squared)*2/9 - delta*FU*4/9;
      f(x, y, 1) = f(x, y, 1)*alpha + beta*rho(x, y)*(2 + 6*u(x,y,0) + 9*u(x,y,0)*u(x,y,0) - three_u_squared)/18 + delta*(3*F(0) + 9*F(0)*u(x,y,0) - FU)/9;
      f(x, y, 2) = f(x, y, 2)*alpha + beta*rho(x, y)*(2 + 6*u(x,y,1) + 9*u(x,y,1)*u(x,y,1) - three_u_squared)/18 + delta*(3*F(1) + 9*F(1)*u(x,y,1) - FU)/9;
      f(x, y, 3) = f(x, y, 3)*alpha + beta*rho(x, y)*(2 - 6*u(x,y,0) + 9*u(x,y,0)*u(x,y,0) - three_u_squared)/18 + delta*(-3*F(0) + 9*F(0)*u(x,y,0) - FU)/9;
      f(x, y, 4) = f(x, y, 4)*alpha + beta*rho(x, y)*(2 - 6*u(x,y,1) + 9*u(x,y,1)*u(x,y,1) - three_u_squared)/18 + delta*(-3*F(1) + 9*F(1)*u(x,y,1) - FU)/9;
      f(x, y, 5) = f(x, y, 5)*alpha + beta*rho(x, y)*(1 + 3*(u(x,y,0)+u(x,y,1)) + 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36 + delta*( 3*( F(0) + F(1)) + 9*(F(0)  + F(1))*( u(x,y,0)+u(x,y,1)) - FU)/36;
      f(x, y, 6) = f(x, y, 6)*alpha + beta*rho(x, y)*(1 - 3*(u(x,y,0)-u(x,y,1)) - 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36 + delta*( 3*(-F(0) + F(1)) + 9*(-F(0) + F(1))*(-u(x,y,0)+u(x,y,1)) - FU)/36;
      f(x, y, 7) = f(x, y, 7)*alpha + beta*rho(x, y)*(1 - 3*(u(x,y,0)+u(x,y,1)) + 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36 + delta*(-3*( F(0) + F(1)) + 9*(F(0)  + F(1))*( u(x,y,0)+u(x,y,1)) - FU)/36;
      f(x, y, 8) = f(x, y, 8)*alpha + beta*rho(x, y)*(1 + 3*(u(x,y,0)-u(x,y,1)) - 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36 + delta*( 3*( F(0) - F(1)) + 9*(F(0)  - F(1))*( u(x,y,0)-u(x,y,1)) - FU)/36;

      x_next = x+1;
      x_prev = x-1;
      y_next = y+1;
      y_prev = y-1;

      if (x == 0)
        x_prev = Nx-1;
      if (x == Nx-1)
        x_next = 0;
      if (y == 0)
        y_prev = Ny-1;
      if (y == Ny-1)
        y_next = 0;

      f_star(x, y, 0)      = f(x, y, 0);
      f_star(x_next, y, 1) = f(x, y, 1);
      f_star(x, y_next, 2) = f(x, y, 2);
      f_star(x_prev, y, 3) = f(x, y, 3);
      f_star(x, y_prev, 4) = f(x, y, 4);
      f_star(x_next, y_next, 5) = f(x, y, 5);
      f_star(x_prev, y_next, 6) = f(x, y, 6);
      f_star(x_prev, y_prev, 7) = f(x, y, 7);
      f_star(x_next, y_prev, 8) = f(x, y, 8);
    }
    
    f_prev = f_star;
    for (int i = 0; i < boundary.size(); i++)
      {
        x = get<0>(boundary[i]);
        y = get<1>(boundary[i]);

        if (f_prev(x,y,0)+f_prev(x,y,1)+f_prev(x,y,2)+f_prev(x,y,3)+f_prev(x,y,4)+f_prev(x,y,5)+f_prev(x,y,6)+f_prev(x,y,7)+f_prev(x,y,8) > 0.000000001)
        {
        f(x,y,1) = f_prev(x,y,3);
        f(x,y,2) = f_prev(x,y,4);
        f(x,y,3) = f_prev(x,y,1);
        f(x,y,4) = f_prev(x,y,2);
        f(x,y,5) = f_prev(x,y,7);
        f(x,y,6) = f_prev(x,y,8);
        f(x,y,7) = f_prev(x,y,5);
        f(x,y,8) = f_prev(x,y,6);

        x_next = x+1;
        x_prev = x-1;
        y_next = y+1;
        y_prev = y-1;

        if (x == 0)
          x_prev = Nx-1;
        if (x == Nx-1)
          x_next = 0;
        if (y == 0)
          y_prev = Ny-1;
        if (y == Ny-1)
          y_next = 0;

        f_star(x, y, 0)      = f(x, y, 0);
        f_star(x_next, y, 1) = f(x, y, 1);
        f_star(x, y_next, 2) = f(x, y, 2);
        f_star(x_prev, y, 3) = f(x, y, 3);
        f_star(x, y_prev, 4) = f(x, y, 4);
        f_star(x_next, y_next, 5) = f(x, y, 5);
        f_star(x_prev, y_next, 6) = f(x, y, 6);
        f_star(x_prev, y_prev, 7) = f(x, y, 7);
        f_star(x_next, y_prev, 8) = f(x, y, 8);

        f_star(x,y,0) = 0;
        f_star(x,y,1) = 0;
        f_star(x,y,2) = 0;
        f_star(x,y,3) = 0;
        f_star(x,y,4) = 0;
        f_star(x,y,5) = 0;
        f_star(x,y,6) = 0;
        f_star(x,y,7) = 0;
        f_star(x,y,8) = 0;
      }}
      
      //cout << sqrt(sum_difference/sum) << endl;
  if (sqrt(sum_difference/sum) <  tol)
      {equil = true;
      cout << abs(u - prev_u).max() << " " << current_max_u << " " << counter << endl;}
  prev_u = u;
  f = f_star;
  counter += 1;
  }
  //cout << "number of steps for equilibration: " << counter << endl;
}


void MainClass::update_g()
{
    g(x, y, 0) = g(x, y, 0)*eta + zeta*C(x, y)*(2 - three_u_squared)*2/9;
    g(x, y, 1) = g(x, y, 1)*eta + zeta*C(x, y)*(2 + 6*u(x,y,0) + 9*u(x,y,0)*u(x,y,0) - three_u_squared)/18;
    g(x, y, 2) = g(x, y, 2)*eta + zeta*C(x, y)*(2 + 6*u(x,y,1) + 9*u(x,y,1)*u(x,y,1) - three_u_squared)/18;
    g(x, y, 3) = g(x, y, 3)*eta + zeta*C(x, y)*(2 - 6*u(x,y,0) + 9*u(x,y,0)*u(x,y,0) - three_u_squared)/18;
    g(x, y, 4) = g(x, y, 4)*eta + zeta*C(x, y)*(2 - 6*u(x,y,1) + 9*u(x,y,1)*u(x,y,1) - three_u_squared)/18;
    g(x, y, 5) = g(x, y, 5)*eta + zeta*C(x, y)*(1 + 3*(u(x,y,0)+u(x,y,1)) + 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
    g(x, y, 6) = g(x, y, 6)*eta + zeta*C(x, y)*(1 - 3*(u(x,y,0)-u(x,y,1)) - 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
    g(x, y, 7) = g(x, y, 7)*eta + zeta*C(x, y)*(1 - 3*(u(x,y,0)+u(x,y,1)) + 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
    g(x, y, 8) = g(x, y, 8)*eta + zeta*C(x, y)*(1 + 3*(u(x,y,0)-u(x,y,1)) - 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
}

mat MainClass::ADE(int T)
{
  int data_divide = T/data_lines;
  int counter = 0;

  Mat<double> C_in;
  C_in = Mat<double>(T, source.size());

  for (int t = 0; t < T; t++)
    {//open
      for (int k = 0; k < rest.size(); k++)
      {x = get<0>(rest[k]);
       y = get<1>(rest[k]);

      C(x, y)  = g(x, y, 0) + g(x, y, 1) + g(x, y, 2) + g(x, y, 3) + g(x, y, 4) + g(x, y, 5) + g(x, y, 6) + g(x, y, 7) + g(x, y, 8); 
      three_u_squared = 3*(u(x, y, 0)*u(x, y, 0) + u(x, y, 1)*u(x, y, 1));
      update_g();

      x_next = x+1;
      x_prev = x-1;
      y_next = y+1;
      y_prev = y-1;

      if (x == 0)
        x_prev = Nx-1;
      if (x == Nx-1)
        x_next = 0;
      if (y == 0)
        y_prev = Ny-1;
      if (y == Ny-1)
        y_next = 0;

      g_star(x, y,      0) = g(x, y, 0);
      g_star(x_next, y, 1) = g(x, y, 1);
      g_star(x, y_next, 2) = g(x, y, 2);
      g_star(x_prev, y, 3) = g(x, y, 3);
      g_star(x, y_prev, 4) = g(x, y, 4);

      g_star(x_next, y_next, 5) = g(x, y, 5);
      g_star(x_prev, y_next, 6) = g(x, y, 6);
      g_star(x_prev, y_prev, 7) = g(x, y, 7);
      g_star(x_next, y_prev, 8) = g(x, y, 8);
    }

    //boundary
    for (int i = 0; i < boundary.size(); i++)
      {x = get<0>(boundary[i]);
       y = get<1>(boundary[i]);

        x_next = x+1;
        x_prev = x-1;
        y_next = y+1;
        y_prev = y-1;

        if (x == 0)
          x_prev = Nx-1;
        if (x == Nx-1)
          x_next = 0;
        if (y == 0)
          y_prev = Ny-1;
        if (y == Ny-1)
          y_next = 0;

        g_star(x_next, y, 1) = g_star(x,y,3);
        g_star(x, y_next, 2) = g_star(x,y,4);
        g_star(x_prev, y, 3) = g_star(x,y,1);
        g_star(x, y_prev, 4) = g_star(x,y,2);

        g_star(x_next, y_next, 5) = g_star(x,y,7);
        g_star(x_prev, y_next, 6) = g_star(x,y,8);
        g_star(x_prev, y_prev, 7) = g_star(x,y,5);
        g_star(x_next, y_prev, 8) = g_star(x,y,6);
        
        g_star(x,y,0) = 0;
        g_star(x,y,1) = 0;
        g_star(x,y,2) = 0;
        g_star(x,y,3) = 0;
        g_star(x,y,4) = 0;
        g_star(x,y,5) = 0;
        g_star(x,y,6) = 0;
        g_star(x,y,7) = 0;
        g_star(x,y,8) = 0;

      }

  //drains  
  for (int j = 0; j < source.size(); j++){
      x = get<0>(source[j]);
      y = get<1>(source[j]);
      C_in(t, j) = C(x,y);}
  
  g = g_star;
  if (t%data_divide == 0)
    {write_C(counter, "front");
     counter += 1;
    }
  }
  write_source(C_in, T, "front");
  return C_in;
}


void MainClass::ADE_back(int T, mat C_in)
{
  int data_divide = T/data_lines;
  int counter = 0;
  Mat<double> C_out;
  C_out = Mat<double>(T, source.size());

  for (int t = 0; t < T; t++)
    {//open
      for (int k = 0; k < rest.size(); k++)
      {x = get<0>(rest[k]);
       y = get<1>(rest[k]);

      C(x, y)  = g(x, y, 0) + g(x, y, 1) + g(x, y, 2) + g(x, y, 3) + g(x, y, 4) + g(x, y, 5) + g(x, y, 6) + g(x, y, 7) + g(x, y, 8); 
      three_u_squared = 3*(u(x, y, 0)*u(x, y, 0) + u(x, y, 1)*u(x, y, 1));

      g(x, y, 0) = g(x, y, 0)*eta + zeta*C(x, y)*(2 - three_u_squared)*2/9;
      g(x, y, 1) = g(x, y, 1)*eta + zeta*C(x, y)*(2 - 6*u(x,y,0) + 9*u(x,y,0)*u(x,y,0) - three_u_squared)/18;
      g(x, y, 2) = g(x, y, 2)*eta + zeta*C(x, y)*(2 - 6*u(x,y,1) + 9*u(x,y,1)*u(x,y,1) - three_u_squared)/18;
      g(x, y, 3) = g(x, y, 3)*eta + zeta*C(x, y)*(2 + 6*u(x,y,0) + 9*u(x,y,0)*u(x,y,0) - three_u_squared)/18;
      g(x, y, 4) = g(x, y, 4)*eta + zeta*C(x, y)*(2 + 6*u(x,y,1) + 9*u(x,y,1)*u(x,y,1) - three_u_squared)/18;
      g(x, y, 5) = g(x, y, 5)*eta + zeta*C(x, y)*(1 - 3*(u(x,y,0)+u(x,y,1)) + 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
      g(x, y, 6) = g(x, y, 6)*eta + zeta*C(x, y)*(1 + 3*(u(x,y,0)-u(x,y,1)) - 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
      g(x, y, 7) = g(x, y, 7)*eta + zeta*C(x, y)*(1 + 3*(u(x,y,0)+u(x,y,1)) + 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
      g(x, y, 8) = g(x, y, 8)*eta + zeta*C(x, y)*(1 - 3*(u(x,y,0)-u(x,y,1)) - 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;

      x_next = x+1;
      x_prev = x-1;
      y_next = y+1;
      y_prev = y-1;

      if (x == 0)
        x_prev = Nx-1;
      if (x == Nx-1)
        x_next = 0;
      if (y == 0)
        y_prev = Ny-1;
      if (y == Ny-1)
        y_next = 0;

      g_star(x, y,      0) = g(x, y, 0);
      g_star(x_next, y, 1) = g(x, y, 1);
      g_star(x, y_next, 2) = g(x, y, 2);
      g_star(x_prev, y, 3) = g(x, y, 3);
      g_star(x, y_prev, 4) = g(x, y, 4);

      g_star(x_next, y_next, 5) = g(x, y, 5);
      g_star(x_prev, y_next, 6) = g(x, y, 6);
      g_star(x_prev, y_prev, 7) = g(x, y, 7);
      g_star(x_next, y_prev, 8) = g(x, y, 8);
    }

    //boundary
    for (int i = 0; i < boundary.size(); i++)
      {x = get<0>(boundary[i]);
       y = get<1>(boundary[i]);

        x_next = x+1;
        x_prev = x-1;
        y_next = y+1;
        y_prev = y-1;

        if (x == 0)
          x_prev = Nx-1;
        if (x == Nx-1)
          x_next = 0;
        if (y == 0)
          y_prev = Ny-1;
        if (y == Ny-1)
          y_next = 0;

        g_star(x_next, y, 1) = g_star(x,y,3);
        g_star(x, y_next, 2) = g_star(x,y,4);
        g_star(x_prev, y, 3) = g_star(x,y,1);
        g_star(x, y_prev, 4) = g_star(x,y,2);

        g_star(x_next, y_next, 5) = g_star(x,y,7);
        g_star(x_prev, y_next, 6) = g_star(x,y,8);
        g_star(x_prev, y_prev, 7) = g_star(x,y,5);
        g_star(x_next, y_prev, 8) = g_star(x,y,6);

        g_star(x,y,0) = 0;
        g_star(x,y,1) = 0;
        g_star(x,y,2) = 0;
        g_star(x,y,3) = 0;
        g_star(x,y,4) = 0;
        g_star(x,y,5) = 0;
        g_star(x,y,6) = 0;
        g_star(x,y,7) = 0;
        g_star(x,y,8) = 0;
      }

  //sources  
  
   for (int j = 0; j < source.size(); j++){
      x = get<0>(source[j]);
      y = get<1>(source[j]);
      g_star(x,y,0) += 4*C_in(T-t-1, j)/9;
      g_star(x,y,1) +=   C_in(T-t-1, j)/9;
      g_star(x,y,2) +=   C_in(T-t-1, j)/9;
      g_star(x,y,3) +=   C_in(T-t-1, j)/9;
      g_star(x,y,4) +=   C_in(T-t-1, j)/9;
      g_star(x,y,5) +=   C_in(T-t-1, j)/36;
      g_star(x,y,6) +=   C_in(T-t-1, j)/36;
      g_star(x,y,7) +=   C_in(T-t-1, j)/36;
      g_star(x,y,8) +=   C_in(T-t-1, j)/36;
      }
  
  g = g_star;
  if (t%data_divide == 0)
    {
      write_C(counter, "back");
     counter += 1;
     cout << counter << " " << data_lines << endl;}
  }
}

mat MainClass::ADE_heat(int T, double wall_T)
{
  int data_divide = T/data_lines;
  int counter = 0;

  Mat<double> C_in;
  C_in = Mat<double>(T, source.size());

  for (int t = 0; t < T; t++)
    {//open
      for (int k = 0; k < rest.size(); k++)
      {x = get<0>(rest[k]);
       y = get<1>(rest[k]);

      C(x, y)  = g(x, y, 0) + g(x, y, 1) + g(x, y, 2) + g(x, y, 3) + g(x, y, 4) + g(x, y, 5) + g(x, y, 6) + g(x, y, 7) + g(x, y, 8); 
      three_u_squared = 3*(u(x, y, 0)*u(x, y, 0) + u(x, y, 1)*u(x, y, 1));

      g(x, y, 0) = g(x, y, 0)*eta + zeta*C(x, y)*(2 - three_u_squared)*2/9;
      g(x, y, 1) = g(x, y, 1)*eta + zeta*C(x, y)*(2 + 6*u(x,y,0) + 9*u(x,y,0)*u(x,y,0) - three_u_squared)/18;
      g(x, y, 2) = g(x, y, 2)*eta + zeta*C(x, y)*(2 + 6*u(x,y,1) + 9*u(x,y,1)*u(x,y,1) - three_u_squared)/18;
      g(x, y, 3) = g(x, y, 3)*eta + zeta*C(x, y)*(2 - 6*u(x,y,0) + 9*u(x,y,0)*u(x,y,0) - three_u_squared)/18;
      g(x, y, 4) = g(x, y, 4)*eta + zeta*C(x, y)*(2 - 6*u(x,y,1) + 9*u(x,y,1)*u(x,y,1) - three_u_squared)/18;
      g(x, y, 5) = g(x, y, 5)*eta + zeta*C(x, y)*(1 + 3*(u(x,y,0)+u(x,y,1)) + 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
      g(x, y, 6) = g(x, y, 6)*eta + zeta*C(x, y)*(1 - 3*(u(x,y,0)-u(x,y,1)) - 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
      g(x, y, 7) = g(x, y, 7)*eta + zeta*C(x, y)*(1 - 3*(u(x,y,0)+u(x,y,1)) + 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
      g(x, y, 8) = g(x, y, 8)*eta + zeta*C(x, y)*(1 + 3*(u(x,y,0)-u(x,y,1)) - 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;

      x_next = x+1;
      x_prev = x-1;
      y_next = y+1;
      y_prev = y-1;

      if (x == 0)
        x_prev = Nx-1;
      if (x == Nx-1)
        x_next = 0;
      if (y == 0)
        y_prev = Ny-1;
      if (y == Ny-1)
        y_next = 0;

      g_star(x, y,      0) = g(x, y, 0);
      g_star(x_next, y, 1) = g(x, y, 1);
      g_star(x, y_next, 2) = g(x, y, 2);
      g_star(x_prev, y, 3) = g(x, y, 3);
      g_star(x, y_prev, 4) = g(x, y, 4);

      g_star(x_next, y_next, 5) = g(x, y, 5);
      g_star(x_prev, y_next, 6) = g(x, y, 6);
      g_star(x_prev, y_prev, 7) = g(x, y, 7);
      g_star(x_next, y_prev, 8) = g(x, y, 8);
    }

    //boundary
    for (int i = 0; i < boundary.size(); i++)
      {x = get<0>(boundary[i]);
       y = get<1>(boundary[i]);

        x_next = x+1;
        x_prev = x-1;
        y_next = y+1;
        y_prev = y-1;

        if (x == 0)
          x_prev = Nx-1;
        if (x == Nx-1)
          x_next = 0;
        if (y == 0)
          y_prev = Ny-1;
        if (y == Ny-1)
          y_next = 0;

        g_star(x_next, y, 1) = -g_star(x,y,3) + 2*wall_T/9;
        g_star(x, y_next, 2) = -g_star(x,y,4) + 2*wall_T/9;
        g_star(x_prev, y, 3) = -g_star(x,y,1) + 2*wall_T/9;
        g_star(x, y_prev, 4) = -g_star(x,y,2) + 2*wall_T/9;

        g_star(x_next, y_next, 5) = -g_star(x,y,7) + wall_T/18;
        g_star(x_prev, y_next, 6) = -g_star(x,y,8) + wall_T/18;
        g_star(x_prev, y_prev, 7) = -g_star(x,y,5) + wall_T/18;
        g_star(x_next, y_prev, 8) = -g_star(x,y,6) + wall_T/18;
        
        g_star(x,y,0) = 0;
        g_star(x,y,1) = 0;
        g_star(x,y,2) = 0;
        g_star(x,y,3) = 0;
        g_star(x,y,4) = 0;
        g_star(x,y,5) = 0;
        g_star(x,y,6) = 0;
        g_star(x,y,7) = 0;
        g_star(x,y,8) = 0;
      }

  //drains  
  for (int j = 0; j < source.size(); j++){
      x = get<0>(source[j]);
      y = get<1>(source[j]);
      C_in(t, j) = C(x,y);}
  
  g = g_star;
  if (t%data_divide == 0)
    {write_C(counter, "front");
     counter += 1;
    }
  }
  write_source(C_in, T, "front");
  return C_in;
}

  void MainClass::write_u(string name)
  {
    ofstream outfile("../data/" + filename + "_" + name + "_u.txt");
    if (!outfile.is_open())
     cout<<"Could not open file" << endl;
    for (int i = 0; i < Nx; i++)
      {
        for (int j = 0; j < Ny; j++)
        {
          outfile << u(i, j, 0) << " "; 
        }
      }
    outfile << "\n"; 
    for (int i = 0; i < Nx; i++)
    {
      for (int j = 0; j < Ny; j++)
      {
        outfile << u(i, j, 1) << " "; 
      }
    }
  }


void MainClass::write_C(int T, string filename2)
{
  ofstream outfile("../data/" + filename + "_C_" + to_string(T) + "_" + filename2 + ".txt");
  if (!outfile.is_open())
  cout<<"Could not open file" << endl;
  for (int i = 0; i < Nx; i++)
    {
      for (int j = 0; j < Ny; j++){
      outfile << C(i, j) << " "; 
    }
  }
}

void MainClass::write_source(mat data, int T, string filename2)
{
  ofstream outfile("../data/" + filename + "_source_" + filename2 + ".txt");
  if (!outfile.is_open())
  cout<<"Could not open file" << endl;
  for (int j = 0; j < source.size(); j++)
    {
      for (int t = 0; t < T; t++){
      outfile << data(t, j) << " "; 
    }
    outfile << "\n";
  }
}

void MainClass::test_mass_cons()
  {
    double initial_mass = 0;
    for (int x = 0; x < Nx; x++)
    {
      for (int y = 0; y < Ny; y++)
      {
        initial_mass +=  f(x, y, 0) + f(x, y, 1) + f(x, y, 2) + f(x, y, 3) + f(x, y, 4) + f(x, y, 5) + f(x, y, 6) + f(x, y, 7) + f(x, y, 8); 
      }
    }

    run();

    double final_mass = 0;
    for (int x = 0; x < Nx; x++)
    {
      for (int y = 0; y < Ny; y++)
      {
        final_mass +=  f(x, y, 0) + f(x, y, 1) + f(x, y, 2) + f(x, y, 3) + f(x, y, 4) + f(x, y, 5) + f(x, y, 6) + f(x, y, 7) + f(x, y, 8); 
      }
    }

  cout << initial_mass << " " << final_mass << " " << initial_mass-final_mass << endl;
  if (abs(initial_mass-final_mass) < 0.0001*initial_mass)
  {cout << "MASS IS CONSERVED" << endl;}
  else
    {cout << "FUCK! MASS IS NOT CONSERVED!" << endl;}
  }

void MainClass::test_mass_diffusion()
  {
    double initial_mass = 0;
    initialize_C(16, 32, 0, 100);
    for (int x = 0; x < Nx; x++)
    {for (int y = 0; y < Ny; y++)
      {initial_mass +=  g(x, y, 0) + g(x, y, 1) + g(x, y, 2) + g(x, y, 3) + g(x, y, 4) + g(x, y, 5) + g(x, y, 6) + g(x, y, 7) + g(x, y, 8);}
    }

    ADE(500);

    double final_mass = 0;
    for (int x = 0; x < Nx; x++)
    {for (int y = 0; y < Ny; y++)
      {final_mass +=  g(x, y, 0) + g(x, y, 1) + g(x, y, 2) + g(x, y, 3) + g(x, y, 4) + g(x, y, 5) + g(x, y, 6) + g(x, y, 7) + g(x, y, 8);}
    }
  cout << initial_mass << " " << final_mass << " " << initial_mass-final_mass << endl;
  if (abs(initial_mass-final_mass) < 0.0001*initial_mass)
  {cout << "CONCENTRATION IS CONSERVED" << endl;}
  else
    {cout << "FUCK! CONCENTRATION IS NOT CONSERVED!" << endl;}
  }
