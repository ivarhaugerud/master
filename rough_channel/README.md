
# Dispersion in a rough channel


## Current directory

The file ```flow.py``` is the main file used for this section. It solves the stationary Navier--Stokes equations, and the Brenner equations for the produced velocity field. The file is commented to make the code more readable.
To achieve this it uses the files ```bcs.py``` and ```mesh_box.py```. The first of which gives the periodic boundary conditions, and define descrete mesh functions to decide if a lattice point is on the boudnary or not.
The file ```mesh_box.py``` creates the unit cell boundary for the rough geometry. Furthermore, three other geometries have been included; the two next fractal orders of the same geometry, and a third unit cell which is not vertically symmetric.

## Plot
All of the code used to plot and analyze the generated data is contained. 
