# EMDOT code

The EMDOT code has been created to perform hydrodynamical simulations of the atmospheric mass loss from irradiated exoplanets. THe code description can be found at  <<insert citation>>. In the following we describe the code organization and give a basi example of usage.

## Directories and files

The main directory of the code consists of the following elements:
* the main code file `Hydro_ioniz.f90`;
* the bash script `run_program.sh` that takes care of the compilation and the execution of the code
* the `src` directory, where all the code modules are stored.

The `src` directory contains three major sudirectories:
* the `/src/utils` folder contains the python3 files dedicated for the creatioin of the input interface:
* the `/src/modules` folder contains all the `.f90` files for all the subroutines of the code;
* the `/src/mod` folder stores the `.mod` modules interface files.

In the `/src/modules` subdirectory the code's modules are subdivided as follows:
* `/src/modules/files_IO` : subroutines for the input/output management (read input parameters, load initial conditions, write the simulation output);
* `/src/modules/flux` : library with the implemented numerical flux functions and wavespeed estimates;
* `/src/modules/functions` : various useful functions (energy-dependent cross sections, gravitational field, primitive-conservative variable conversion);
* `/src/modules/init` : subroutines for the initialization of the code (allocate global vectors and variables, define spatial and energy grid, pre-evaluate some vectors, set initial conditions);
* `/src/modules/nonlinear_system_solver` : subroutines from MINPACK (https://www.netlib.org/minpack/) and definition of the photoionization equilibrium system;
* `/src/modules/radiation` : subroutines related to the radiative part of the code (incident spectrum, photoheating and photoionization rates, temperature-dependent reaction rates, solution of the ionization equilibrium system);
* `/src/modules/states`: subroutines for the hydrodynamical reconstruction step (PLM, ESWENO3), definition of the boundary conditions and evaluation of the source term;
* `/src/modules/time_step` : evaluation of the right hand side of the Runge-Kutta integrator.


##  Installing and using the code

The code can be directly downloaded from the Github page or, in alternative, the repository can be cloned via 
    git clone https://github.com/AndreaCaldiroli/ExoCode2021

