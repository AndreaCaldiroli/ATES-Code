#!/bin/bash

# Executable bash script to tun the full code        

#------- Directories -------#

# Current directory
PWD=$(pwd)

# Main program directory
DIR_MAIN="$PWD" 

# Source files directory
DIR_SRC="$PWD/src" 

# Mod files directory
DIR_MOD="$DIR_SRC/mod" 

# Utilities directory
DIR_UTILS="$DIR_SRC/utils"

# Modules directories
DIR_FILES="$DIR_SRC/modules/files_IO"
DIR_FLUX="$DIR_SRC/modules/flux"
DIR_FUNC="$DIR_SRC/modules/functions"
DIR_INIT="$DIR_SRC/modules/init"
DIR_NLSOLVE="$DIR_SRC/modules/nonlinear_system_solver"
DIR_RAD="$DIR_SRC/modules/radiation"
DIR_STAT="$DIR_SRC/modules/states"
DIR_TIME="$DIR_SRC/modules/time_step"

#------- Cleaning old files -------#

# Clean old executable 
rm -f $DIR_MAIN/*.x

# Clean .mod files
rm -f $DIR_MOD/*.mod

#------- Call python interface to create input file -------#

TABLE_FILE="$DIR_UTILS/params_table.txt"

if [ ! -f "$TABLE_FILE" ]; then
      python3 -W ignore "$DIR_UTILS/gen_file.py"
      mv "params_table.txt" "$DIR_UTILS/params_table.txt"
fi

python3 -W ignore "$DIR_UTILS/create_input.py"

#------- Check input parameters -------#

INPUT_FILE="$DIR_MAIN/input.inp"

# Check if input parameters file exists
echo "Searching for the input parameters file..."
if [ -f "$INPUT_FILE" ]; then

      echo "Input file found. Proceeding..."

else  # Abort if no input.inp exists

      echo "Missing input parameter file."
      echo "Please create an input.inp file"
      exit 1
fi


# Create output directory if it doesn't exists
if [ ! -d "$DIR_MAIN/output" ]; then
      mkdir "$DIR_MAIN/output"
fi

#------- Executing fortran file -------#

echo "
------- STARTING SIMULATION -------
"

# Define string with order of compilation
str=" gfortran -J"$DIR_MOD" -I"$DIR_MOD" -fopenmp -g -fbacktrace \
      $DIR_INIT/parameters.f90\
      $DIR_FILES/input_read.f90\
      $DIR_FILES/load_IC.f90\
      $DIR_FILES/write_output.f90\
      $DIR_FUNC/grav_field.f90\
      $DIR_FUNC/cross_sec.f90\
      $DIR_FUNC/UW_conversions.f90\
      $DIR_NLSOLVE/dogleg.f90\
      $DIR_NLSOLVE/enorm.f90\
      $DIR_NLSOLVE/hybrd1.f90\
      $DIR_NLSOLVE/qform.f90\
      $DIR_NLSOLVE/r1mpyq.f90\
      $DIR_NLSOLVE/System_HeH.f90\
      $DIR_NLSOLVE/dpmpar.f90\
      $DIR_NLSOLVE/fdjac1.f90\
      $DIR_NLSOLVE/hybrd.f90\
      $DIR_NLSOLVE/qrfac.f90\
      $DIR_NLSOLVE/r1updt.f90\
      $DIR_NLSOLVE/System_H.f90\
      $DIR_RAD/J_inc.f90\
      $DIR_RAD/Cool_coeff.f90\
      $DIR_RAD/ionization_equilibrium.f90\
      $DIR_STAT/Apply_BC.f90\
      $DIR_STAT/PLM_rec.f90\
      $DIR_STAT/WENO3_rec.f90\
      $DIR_STAT/Source.f90\
      $DIR_STAT/Reconstruction.f90\
      $DIR_FLUX/speed_estimate_HLLC.f90\
      $DIR_FLUX/speed_estimate_ROE.f90\
      $DIR_FLUX/Num_Fluxes.f90\
      $DIR_TIME/RK_rhs.f90\
      $DIR_INIT/define_grid.f90\
      $DIR_INIT/set_energy_vectors.f90\
      $DIR_INIT/set_gravity_grid.f90\
      $DIR_INIT/set_IC.f90\
      $DIR_INIT/init.f90"
      
# Add main file to execution
str="$str $DIR_MAIN/Hydro_ioniz.f90"

# Add output file to string of execution
compile="$str -o $DIR_MAIN/Hydro_ioniz.x"

# Compile program and modules
eval $compile

# Exectution string
exec_str="$DIR_MAIN/Hydro_ioniz.x"

# Exectute code
eval $exec_str

# Print when execution is over
echo "
----- FORTRAN EVALUATION DONE -----"

#------- Python plots -------#

# Run python_plots
#python3 -W ignore python_plots.py 



