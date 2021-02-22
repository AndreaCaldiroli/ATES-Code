# EMDOT code

The EMDOT code has been created to perform hydrodynamical simulations of the atmospheric mass loss from irradiated exoplanets. THe code description can be found at << insert citation >>. In the following we describe the code organization and give a basi example of usage.

## Directories and files

The main directory of the code consists of the following elements:
* the main code file `<Hydro_ioniz.f90>`;
* the bash script `<run_program.sh>` that takes care of the compilation and the execution of the code
* the `<src>` directory, where all the code modules are stored.

The `<src>` directory contains three major sudirectories:
* the `</src/utils>` folder contains the python3 files dedicated for the creatioin of the input interface:
* the `</src/modules>` folder contains all the `<.f90>` files for all the subroutines of the code;
* the `</src/mod>` folder stores the `<.mod>` modules interface files.




