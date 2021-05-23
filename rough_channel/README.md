
# Dispersion in a rough channel


## Current directory

The file ```flow.py``` is the main file used for this section. It solves the stationary Navier--Stokes equations, and the Brenner equations for the produced velocity field. The file is commented on to make the code more readable.
To solve the equations, it program uses the files ```bcs.py``` and ```mesh_box.py```. The first of which gives the periodic boundary conditions and define discrete mesh functions to decide if a lattice point is on the boundary or not.  
The file ```mesh_box.py``` creates the unit cell boundary for the rough geometry. Furthermore, three other geometries have been included; the two next fractal orders of the same geometry and a third unit cell that is not vertically symmetric.
The last file ```pure_diff.py``` does the same as ```flow.py```, but without flow. It, therefore, solves the pure diffusion equation using Brenner theory in the presented geometries.

## Plot
All of the code used to plot and analyze the generated data is contained. 
