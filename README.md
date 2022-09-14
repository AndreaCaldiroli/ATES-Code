# The ATES code

The ATES code has been created to perform hydrodynamical simulations of the atmospheric mass loss from irradiated exoplanets. For a detailed description of the code, we refer to [[1]](#1) In the following we describe the code organization and how to run.
For any question or if you notice any bug please write an email to <andrea.caldiroli@univie.ac.at>

## Requirements

The code can be compiled with both `gfortran` (tested successfully in version 9.3.0 and newer) and `ifort` (tested on the 2021.2.0 and the 2021.5.0 versions). For the compiler choice, see below.
A basic installation of `python3` is required. The following libraries are used: `numpy,tkinter,os,shutil,matplotlib,sys,time`.

## Installation

The code doesn't require any special installation, and can be directly downloaded from the Github page or, in alternative, the repository can be cloned via 

    git clone https://github.com/AndreaCaldiroli/ATES-Code

    
## Directories and files

The main directory (`$MAIN`) of the code consists of the following elements:
* the main code file `$MAIN/ATES_main.f90`;
* the bash script `$MAIN/run_ATES.sh` that takes care of the compilation and the execution of the code;
* the `$MAIN/src` directory, where all the code modules are stored.
* the `$MAIN/ATES_plots.py` python3 file for live plots.
* the `$MAIN/eta_approx.py` python3 file with the approximate function of the effective efficiency from Appendix A in [[2]](#2).

The `$MAIN/src` directory contains three major sudirectories:
* the `$MAIN/src/utils` folder contains the python3 files dedicated for the creatioin of the input interface;
* the `$MAIN/src/modules` folder contains all the `.f90` files for all the subroutines of the code;
* the `$MAIN/src/mod` folder stores the `.mod` files.

In the `$MAIN/src/modules` subdirectory, the code's modules are subdivided as follows:
* `$MAIN/src/modules/files_IO` : subroutines for the input/output management (read input parameters, load initial conditions, write the simulation output);
* `$MAIN/src/modules/flux` : library with the implemented numerical flux functions and wavespeed estimates;
* `$MAIN/src/modules/functions` : various useful functions;
* `$MAIN/src/modules/init` : subroutines for the initialization of the code (allocate global vectors and variables, set initial conditions);
* `$MAIN/src/modules/nonlinear_system_solver` : subroutines from MINPACK (https://www.netlib.org/minpack/) and definition of the photoionization equilibrium system;
* `$MAIN/src/modules/post_process` : subroutine for the post processing;
* `$MAIN/src/modules/radiation` : subroutines related to the radiation;
* `$MAIN/src/modules/states`: subroutines for the hydrodynamical reconstruction step, boundary conditions and source terms;
* `$MAIN/src/modules/time_step` : evaluation of the right hand side of the Runge-Kutta integrator.


## Using the code

Once exctracted, it is necessary to give execution permission to the `$MAIN/run_program.sh` file:

    chmod +x $MAIN/run_ATES.sh
    
In order to run the code, the bash file must be executed. By default, ATES is compiled with gfortran. In the terminal, it is sufficient to execute:
   
    .$MAIN/run_ATES.sh

To force the use of the `ifort` compiler, run the following command:

    .$MAIN/run_ATES.sh --ifort

The user is asked to insert the physical parameters of the system to be simulated. See [[1]](#1) for a detailed explanation of such parameters. If a system is not available in the precompiled archive (which is stored in `$MAIN/src/utils/params_table.txt`), it is possible to add it to the default list for later simulations by using the `Add planet` button. 

The code is executed by pressing the `Done` button. In the terminal, the current iteration number and the fractional variation of the momentum over the selected domain of interest, i.e.:
   
   <img src="https://render.githubusercontent.com/render/math?math=\dfrac{\Delta \dot{M} }{\dot{M}} := \dfrac{\max\dot{M} - \min\dot{M}}{\min\dot{M}} \quad \text{for} \quad r>r_{esc}">

For planetary simulations, as explained in [[1]](#1), it is suggested to use the PLM reconstruction procedure when starting the simulation from general initial conditions and stop the simulation manually when <img src="https://render.githubusercontent.com/render/math?math=\Delta \dot{M}/\dot{M} \lesssim 0.5-1">. Then, restart the simulation using the previous outputs as initial condition (see below) and using the WENO3 reconstruction method instead.



## Output files

The code writes the current output of the simulations on two file saved in the `$MAIN/output` directory. The `$MAIN/output/Hydro_ioniz.txt` file stores the hydrodynamical variables, which are saved in column vectors in the following order:
1. radial distance (in unit of the planetary radius)
2. mass density (in unit of the proton mass)
3. velocity (in unit of the scale velocity - see [[1]](#1))
4. pressure (in CGS units)
5. Temperature (in Kelvin)
6. Radiative heating rate (in CGS units)
7. Radiative cooling rate (in CGS units)
8. Heating efficiency (adimensional)


The ionization profiles are saved in the `$MAIN/output/Ion_species.txt` file. The columns of the file correspond to the number densities of HI, HII, HeI, HeII, HeIII in <img src="https://render.githubusercontent.com/render/math?math=\text{cm}^{-3}">.

The post-processed profile are written on the `$MAIN/output/Hydro_ioniz_adv.txt` and `$MAIN/output/Ion_species_adv.txt` files. The data are formatted as the `$MAIN/output/Hydro_ioniz.txt` and `$MAIN/output/Ion_species.txt` files.

If the `Load IC` flag is active in the input window, the code automatically chooses the last saved `$MAIN/output/Hydro_ioniz.txt` and `$MAIN/output/Ion_species.txt`files in the `$MAIN/output` directory and copies them onto two new files named, by default,`$MAIN/output/Hydro_ioniz_IC.txt` and `$MAIN/output/Ion_species_IC.txt`, which are loaded by the code. For the writing/reading formats consult the `$MAIN/src/modules/file_IO/load_IC.f90` and `$MAIN/src/modules/file_IO/write_output.f90` files.

## Plotting results

The `$MAIN/ATES_plots.py` file can be used to plot the current status of the simulation or to follow the evolution of the profiles with a live animation. The script can be executed with the following syntax:

    python3 $MAIN/ATES_plots.py --live n
    
The `--live n` arguments are optional, and can therefore be omitted. If so, the content of the current `$MAIN/output/Hydro_ioniz.txt` and `$MAIN/output/Ion_species.txt` is plotted. If only the `--live` flag is used, the figure is updated by default every 4 seconds with the content of the current output files (which ATES, by defaults, overwrites every 1000th temporal iteration). To set the time update interval, specify the `n` argument with the desired number of seconds between the updates. Finally, a second figure with the post-processed profiles is created if the corresponding files (`$MAIN/output/Hydro_ioniz_adv.txt`and `$MAIN/output/Ion_species_adv.txt`) are found in the `$MAIN/output` directory.

## Approximate effective efficiency function

The file `$MAIN/eta_approx.py` contains the approximate expression for the effective efficiency presented in  [[2]](#2). The file can be simply run as:

    python3 $MAIN/eta_approx.py

The user must provide the planetary parameters directly through the terminal window. The approximate values of the effective efficiency and the mass loss rate are printed as outputs.

## References
<a id="1">[1]</a> 
Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I., Rauscher, E., 2021, "Irradiation-driven escape of primordial planetary atmospheres I. The ATES photoionization hydrodynamics code", A&A, 655, A30 (2021).

<a id="2">[2]</a> 
Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I., Rauscher, E., 2021, "Irradiation-driven escape of primordial planetary atmospheres II. Evaporation efficiency of sub-Neptunes through hot Jupiters", A&A, 663, A122, (2022).

