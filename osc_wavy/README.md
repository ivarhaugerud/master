
# Dispersion in a channel with sinusoidal varying boundary and an oscillating external driving force


## Current directory

The file ```OscNs.py``` is the main file used for this section. It solves the Navier--Stokes equations, and the time-dependent Brenner equations for the produced velocity field using the finite element method. The file is commented on to make the code more readable.

## RW_simulations
Code used to perform random walk simulations, which is used to verify that the file ```OscNs.py``` is working correctly.

## finished plots
Code used to analyze and plot the generated data.

## finite_element
A one-dimensional finite element solver implemented from the ground up. Some files are used to bench-mark the solver for various difficulties of ordinary differential equations. Furthermore, some files are used to solve the equations for Brenner theory and calculate the second-order effective diffusion coefficient.

## sympy
A symbolic python package is used to verify the analytic solutions of the equations, by plugging the found solution into the original equation, and check if the left- and right-hand side is identical.

## trash code
Code which ended up not being used for any of the results

## ugradu
Folder for the plotting and analysis of results where the non-linear term in the Navier-Stokes equations is not ignored



