#include "MainClass.hpp"

MainClass::MainClass()
{}

//all variabels we neeed for the class, in addition to the two functions for the energy and wave function
MainClass::MainClass(int NX, int NY, double TAU, double TAU_G, double FX, double FY, double tolerence, string filename, int amount_of_data)
{
  //the file name for saving data, and number of lines when saving data
  filename   = filename;
  data_lines = amount_of_data;

  Nx  = NX;
  Ny  = NY;
  tau = TAU;
  tau_g = TAU_G;

  //advection
  f       = Cube<double>(Nx, Ny, 9);
  f_prev  = Cube<double>(Nx, Ny, 9);
  f_star  = Cube<double>(Nx, Ny, 9);
  f_eq    = Cube<double>(Nx, Ny, 9);
  S       = Cube<double>(Nx, Ny, 9);
  u       = Cube<double>(Nx, Ny, 2);
  prev_u  = Cube<double>(Nx, Ny, 2);

  //diffusion
  g       = Cube<double>(Nx, Ny, 9);
  g_eq    = Cube<double>(Nx, Ny, 9);
  g_star  = Cube<double>(Nx, Ny, 9);
  g_prev  = Cube<double>(Nx, Ny, 9);

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
{
  f(x, y, i) = rho;
}

void MainClass::initialize_C(int x, int y, int i, double rho)
{
  g(x, y, i) = rho;
}

void MainClass::set_boundary()
{
  for (int i = 0; i < Nx; i ++)
  {
    boundary.emplace_back(i, 0);
    boundary.emplace_back(i, Ny-1);
  }
}

void MainClass::boundary_disc(int x, int y, double R)
{
  for (int i = 0; i < Nx; i++)
    {for (int j = 0; j < Ny; j++)
      {if (sqrt( (x-i)*(x-i) + (y-j)*(y-j) ) < R)
      {boundary.emplace_back(i, j);}}}
}

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

void MainClass::run()
{
  bool equil = false;
  int counter = 0;
  double current_max_u;

  while (not equil)
    {for (int k = 0; k < rest.size(); k ++)
      {x = get<0>(rest[k]);
       y = get<1>(rest[k]);

      rho(x, y)  = f(x, y, 0) + f(x, y, 1) + f(x, y, 2) + f(x, y, 3) + f(x, y, 4) + f(x, y, 5) + f(x, y, 6) + f(x, y, 7) + f(x, y, 8); 
      u(x, y, 0) = (f(x, y, 1) - f(x, y, 3) + f(x, y, 5) - f(x, y, 6) - f(x, y, 7) + f(x, y, 8))/rho(x,y) + F(0)/(2*rho(x,y));
      u(x, y, 1) = (f(x, y, 2) - f(x, y, 4) + f(x, y, 5) + f(x, y, 6) - f(x, y, 7) - f(x, y, 8))/rho(x,y) + F(1)/(2*rho(x,y));
      three_u_squared = u(x, y, 0)*u(x, y, 0) + u(x, y, 1)*u(x, y, 1);

      f_eq(x, y, 0) = rho(x, y)*(2 - three_u_squared)*2/9;
      f_eq(x, y, 1) = rho(x, y)*(2 + 6*u(x,y,0) + 9*u(x,y,0)*u(x,y,0) - three_u_squared)/18;
      f_eq(x, y, 2) = rho(x, y)*(2 + 6*u(x,y,1) + 9*u(x,y,1)*u(x,y,1) - three_u_squared)/18;
      f_eq(x, y, 3) = rho(x, y)*(2 - 6*u(x,y,0) + 9*u(x,y,0)*u(x,y,0) - three_u_squared)/18;
      f_eq(x, y, 4) = rho(x, y)*(2 - 6*u(x,y,1) + 9*u(x,y,1)*u(x,y,1) - three_u_squared)/18;
      f_eq(x, y, 5) = rho(x, y)*(1 + 3*(u(x,y,0)+u(x,y,1)) + 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
      f_eq(x, y, 6) = rho(x, y)*(1 - 3*(u(x,y,0)-u(x,y,1)) - 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
      f_eq(x, y, 7) = rho(x, y)*(1 - 3*(u(x,y,0)+u(x,y,1)) + 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
      f_eq(x, y, 8) = rho(x, y)*(1 + 3*(u(x,y,0)-u(x,y,1)) - 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;

      FU = 3*(F(0)*u(x,y,0) + F(1)*u(x,y,1));
      S(x, y, 0) = -FU;
      S(x, y, 1) =  3*F(0) + 9*F(0)*u(x,y,0) - FU;
      S(x, y, 2) =  3*F(1) + 9*F(1)*u(x,y,1) - FU;
      S(x, y, 3) = -3*F(0) + 9*F(0)*u(x,y,0) - FU;
      S(x, y, 4) = -3*F(1) + 9*F(1)*u(x,y,1) - FU;
      S(x, y, 5) =  3*( F(0) + F(1)) + 9*(F(0)  + F(1))*( u(x,y,0)+u(x,y,1)) - FU;
      S(x, y, 6) =  3*(-F(0) + F(1)) + 9*(-F(0) + F(1))*(-u(x,y,0)+u(x,y,1)) - FU;
      S(x, y, 7) = -3*( F(0) + F(1)) + 9*(F(0)  + F(1))*( u(x,y,0)+u(x,y,1)) - FU;
      S(x, y, 8) =  3*( F(0) - F(1)) + 9*(F(0)  - F(1))*( u(x,y,0)-u(x,y,1)) - FU;

      f(x, y, 0) = f(x, y, 0)*alpha + f_eq(x, y, 0)*beta + delta*S(x,y,0)*4/9;
      f(x, y, 1) = f(x, y, 1)*alpha + f_eq(x, y, 1)*beta + delta*S(x,y,1)*1/9;
      f(x, y, 2) = f(x, y, 2)*alpha + f_eq(x, y, 2)*beta + delta*S(x,y,2)*1/9;
      f(x, y, 3) = f(x, y, 3)*alpha + f_eq(x, y, 3)*beta + delta*S(x,y,3)*1/9;
      f(x, y, 4) = f(x, y, 4)*alpha + f_eq(x, y, 4)*beta + delta*S(x,y,4)*1/9;
      f(x, y, 5) = f(x, y, 5)*alpha + f_eq(x, y, 5)*beta + delta*S(x,y,5)*1/36;
      f(x, y, 6) = f(x, y, 6)*alpha + f_eq(x, y, 6)*beta + delta*S(x,y,6)*1/36;
      f(x, y, 7) = f(x, y, 7)*alpha + f_eq(x, y, 7)*beta + delta*S(x,y,7)*1/36;
      f(x, y, 8) = f(x, y, 8)*alpha + f_eq(x, y, 8)*beta + delta*S(x,y,8)*1/36;

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
      }}
      
  current_max_u = u.max();

  if (abs(u - prev_u).max() <  tol*current_max_u)
      {equil = true;
      cout << abs(u - prev_u).max() << " " << current_max_u << endl << " " << counter;}
  prev_u = u;
  f = f_star;
  counter += 1;
  }
  cout << counter << endl;
}


void MainClass::ADE(int T)
{
  for (int t = 0; t < T; t++)
    {for (int k = 0; k < rest.size(); k++)
      {x = get<0>(rest[k]);
       y = get<1>(rest[k]);

      C(x, y)  = g(x, y, 0) + g(x, y, 1) + g(x, y, 2) + g(x, y, 3) + g(x, y, 4) + g(x, y, 5) + g(x, y, 6) + g(x, y, 7) + g(x, y, 8); 

      three_u_squared = 3*(u(x, y, 0)*u(x, y, 0) + u(x, y, 1)*u(x, y, 1));

      g_eq(x, y, 0) = C(x, y)*(2 - three_u_squared)*2/9;
      g_eq(x, y, 1) = C(x, y)*(2 + 6*u(x,y,0) + 9*u(x,y,0)*u(x,y,0) - three_u_squared)/18;
      g_eq(x, y, 2) = C(x, y)*(2 + 6*u(x,y,1) + 9*u(x,y,1)*u(x,y,1) - three_u_squared)/18;
      g_eq(x, y, 3) = C(x, y)*(2 - 6*u(x,y,0) + 9*u(x,y,0)*u(x,y,0) - three_u_squared)/18;
      g_eq(x, y, 4) = C(x, y)*(2 - 6*u(x,y,1) + 9*u(x,y,1)*u(x,y,1) - three_u_squared)/18;
      g_eq(x, y, 5) = C(x, y)*(1 + 3*(u(x,y,0)+u(x,y,1)) + 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
      g_eq(x, y, 6) = C(x, y)*(1 - 3*(u(x,y,0)-u(x,y,1)) - 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
      g_eq(x, y, 7) = C(x, y)*(1 - 3*(u(x,y,0)+u(x,y,1)) + 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;
      g_eq(x, y, 8) = C(x, y)*(1 + 3*(u(x,y,0)-u(x,y,1)) - 9*u(x,y,0)*u(x,y,1) + three_u_squared)/36;

      g(x, y, 0) = g(x, y, 0)*eta + g_eq(x, y, 0)*zeta;
      g(x, y, 1) = g(x, y, 1)*eta + g_eq(x, y, 1)*zeta;
      g(x, y, 2) = g(x, y, 2)*eta + g_eq(x, y, 2)*zeta;
      g(x, y, 3) = g(x, y, 3)*eta + g_eq(x, y, 3)*zeta;
      g(x, y, 4) = g(x, y, 4)*eta + g_eq(x, y, 4)*zeta;
      g(x, y, 5) = g(x, y, 5)*eta + g_eq(x, y, 5)*zeta;
      g(x, y, 6) = g(x, y, 6)*eta + g_eq(x, y, 6)*zeta;
      g(x, y, 7) = g(x, y, 7)*eta + g_eq(x, y, 7)*zeta;
      g(x, y, 8) = g(x, y, 8)*eta + g_eq(x, y, 8)*zeta;

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
    g_prev = g_star;

    for (int i = 0; i < boundary.size(); i++){
      x = get<0>(boundary[i]);
      y = get<1>(boundary[i]);

      C(x,y) = g_prev(x,y,0) + g_prev(x,y,1) + g_prev(x,y,2) + g_prev(x,y,3) + g_prev(x,y,4) + g_prev(x,y,5) + g_prev(x,y,6) + g_prev(x,y,7) + g_prev(x,y,8);
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

      if (g_prev(x,y,1) != 0)
        {g_star(x_prev,y,3) = -g_prev(x,y,1) + 2*C(x,y)/9;}

      if (g_prev(x,y,2) != 0)
        {g_star(x,y_prev,4) = -g_prev(x,y,2) + 2*C(x,y)/9;}

      if (g_prev(x,y,3) != 0)
        {g_star(x_next,y,1) = -g_prev(x,y,3) + 2*C(x,y)/9;}

      if (g_prev(x,y,4) != 0)
        {g_star(x,y_next,2) = -g_prev(x,y,4) + 2*C(x,y)/9;}

       if (g_prev(x,y,5) != 0)
        {g_star(x_prev,y_prev,7) = -g_prev(x,y,5) + 2*C(x,y)/36;}

      if (g_prev(x,y,6) != 0)
        {g_star(x_next,y_prev,8) = -g_prev(x,y,6) + 2*C(x,y)/36;}

      if (g_prev(x,y,7) != 0)
        {g_star(x_next,y_next,5) = -g_prev(x,y,7) + 2*C(x,y)/36;}

      if (g_prev(x,y,8) != 0)
        {g_star(x_prev,y_next,6) = -g_prev(x,y,8) + 2*C(x,y)/36;}
  }
  g = g_star;
  }
}
  void MainClass::write_u()
  {
    ofstream outfile("../data/final_vel.txt");
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


void MainClass::write_C()
{
  ofstream outfile("../data/final_C.txt");
  if (!outfile.is_open())
  cout<<"Could not open file" << endl;
  for (int i = 0; i < Nx; i++)
    {
      for (int j = 0; j < Ny; j++){
      outfile << C(i, j) << " "; 
    }
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
        final_mass +=  rho(x,y); 
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
    for (int x = 0; x < Nx; x++)
    {
      for (int y = 0; y < Ny; y++)
      {
        initial_mass +=  g(x, y, 0) + g(x, y, 1) + g(x, y, 2) + g(x, y, 3) + g(x, y, 4) + g(x, y, 5) + g(x, y, 6) + g(x, y, 7) + g(x, y, 8); 
      }
  }

    ADE(5);

    double final_mass = 0;
    for (int x = 0; x < Nx; x++)
    {
      for (int y = 0; y < Ny; y++)
      {
        final_mass +=  g(x, y, 0) + g(x, y, 1) + g(x, y, 2) + g(x, y, 3) + g(x, y, 4) + g(x, y, 5) + g(x, y, 6) + g(x, y, 7) + g(x, y, 8); 
      }
    }
  cout << initial_mass << " " << final_mass << " " << initial_mass-final_mass << endl;
  if (abs(initial_mass-final_mass) < 0.0001*initial_mass)
  {cout << "CONCENTRATION IS CONSERVED" << endl;}
  else
    {cout << "FUCK! CONCENTRATION IS NOT CONSERVED!" << endl;}
  }
