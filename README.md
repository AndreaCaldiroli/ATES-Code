# EMDOT code

The EMDOT code has been created to perform hydrodynamical simulations of the atmospheric mass loss from irradiated exoplanets. For a detailed description of the code, we refer to << **Citation** >>. In the following we describe the code organization and give a basi example of usage.

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


##  Installing

### Requirements

The only software requirements to run the code are the `gfortran` compiler and a basic installation of `python3`. The only library requested outside the Standard Library of python is`numpy`.

### Installation

The code can be directly downloaded from the Github page or, in alternative, the repository can be cloned via 

    git clone https://github.com/AndreaCaldiroli/ExoCode2021
    
Once exctracted, it is necessary to give execution permission to the `run_program.sh` file:

    chmod +x run_program.sh
    

## Using the code

In order to run the code, the bash file must be executed: `./run_program.sh`. 

The user is prompted to a window in which the parameters of the system to be simulated and the numerical scheme to be used must be provided. See << **Citation** >> for a detailed explanation of such parameters. If a system is not available in the precompiled archive (which is stored in `/src/utils/params_table.txt`), it is possible, after filling all the fields, to add it to the default list for later simulations by using the `Add planet` button. 

The code proceed with the execution of the code after pressing the `Done` button. In the terminal, the current iteration number and the fractional variation of the momentum over the selected domain of interest, i.e.:

   <img src="https://render.githubusercontent.com/render/math?math=\Delta \dot{M} = \dfrac{\max\dot{M} - \min\dot{M}}{\min\dot{M}} \quad for \quad r>r_{esc}">

By default, the code writes the current output of the simulations on two files, which are saved in the `/output/` directory. The `Hydro_ioniz.txt` file stores the hydrodynamical variables, which are saved in column vectors in the following order:
1. radial distance (in unit of the planetary radius)
2. mass density (in unit of the proton mass)
3. velocity (in unit of the scale velocity - see << **Citation** >>)
4. pressure (in CGS units)
5. Temperature (in Kelvin)
6. Radiative heating rate (in CGS units)
7. Radiative cooling rate (in CGS units)
8. Heating efficiency (adimensional)
The ionization profiles are saved in the `Ion_species.txt` file. The columns of the file correspond to the number densities of HI, HII, HeI, HeII, HeIII in cm^{-3}.
If the `Load IC` flag is active in the input window, the code automatically chooses the last saved `Hydro_ioniz.txt` and `Ion_species.txt`files in the `/output/` directory and copies them onto two new files named, by default,`Hydro_ioniz_IC.txt` and `Ion_species_IC.txt`, which are loaded by the code. for the writing/reading formats consult the `/src/modules/file_IO/load_IC.f90` and `/src/modules/file_IO/write_output.f90` files.

Finally, a python script called `python_plots.py` is also available for live plots. 





