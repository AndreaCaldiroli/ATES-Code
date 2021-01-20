#!/bin/bash

# Executable bash script to tun the full code        #
# v.1.0 22-07-2020

############# Directories ###############

# Current directory
PWD=$(pwd)

# Source files directory
DIR_SRC="$PWD/src/modules" 

# Mod files directory
DIR_MOD="$PWD/mod" 

# Main program directory
DIR_MAIN="$PWD"    

# Nonlinear solver directory
DIR_NLSOLVE="$DIR_SRC/nonlinear_system_solver"

# All utilities
DIR_UTILS="$PWD/src/utils/"

# Directory with Salz data
DIR_SALZ="$DIR_UTILS/Salz_data"

############# Cleaning old files ###############

# Clean old executable 
rm -f $DIR_MAIN/*.x

# Clean .mod files
rm -f $DIR_MOD/*.mod

############# Parameters Table ###############

# Planetary parameters table file
TABLE_FILE="$DIR_UTILS/planet_parameters.txt"

# Check if planetary table exists, otherwise create it 
echo "Searching for the table of parameters..."
if [ -f "$TABLE_FILE" ]; then

      echo "Table found. Proceeding..."
      
else

      echo "No table found. Proceding in creating table of parameters..."
      python3 "$DIR_UTILS/planets_data.py"
      mv planet_parameters.txt $DIR_UTILS
      echo "Done."
      
fi

############# Old parameter.f90 file ###############

# Check if a parameters.f90 exists and if you want to reuse the existing one
if [ -f "$DIR_SRC"/parameters.f90 ]; then
      echo "Do you want to restart the setup using the existing parameters file? ( y | n )"
      read params_ans
      
else
      params_ans="n"
fi


# Choose the procedure option (use old parameters.f90 or not)
if [ $params_ans == "n" ]; then
      
      # Clean old Salz data
      rm -f $DIR_MAIN/*_Salz.dat
      
      ############# Choice of system ###############

      # Choose which line to read
      cat $DIR_UTILS/N_header.txt 

      # Read input number of planet
      echo "Choose the value for N:"
      read Nline


      ############# Exctracting parameters ###############

      if [[ ! $Nline == 0 ]]; then  # Use one of the tabulated planets

            # Increase by 18 to skip the header
            Nline=$((Nline+18))
            
            # Remove old temporary file
            rm -f $DIR_UTILS/temp.txt
            
            # Save line to temporary file
            sed "${Nline}q;d" $TABLE_FILE > $DIR_UTILS/temp.txt

            # Extract parameters
            read -r p_name R0 Mp T0 LX LEUV a JXUV Mtilde atilde b0 rho < $DIR_UTILS/temp.txt
                       
            # Choose which luminosities to use
            i=0
            while [ $i == 0 ]
            do
                  echo "Do you want to use the tabulated luminosities [Ojanen (2017)] ( y | n )? "
                  read ans
                  
                  if [ $ans == "n" ]; then # Use user-given luminosities
                  
                        # Read user-given luminosities
                        echo "Insert value of the log10 of X-ray luminosity: "
                        read LX
                        echo "Insert value of the log10 of EUV luminosity: "
                        read LEUV
                        i=1
                        
                        # Write new luminosities to temporary file                        
                        L_str_temp="$LX   $LEUV"
                        L_str_temp="$L_str_temp $a"
                        echo "$LX" > $DIR_UTILS/temp_L.txt
                        echo "$LEUV"  >> $DIR_UTILS/temp_L.txt
                        echo "$a" >> $DIR_UTILS/temp_L.txt
                        
                        # Evaluate and read new XUV flux through fortran file
                        par_temp_str="gfortran $DIR_UTILS/L_temp.f90 -o $DIR_UTILS/L_temp.x"
                        eval $par_temp_str
                        
                        exec_str="$DIR_UTILS/L_temp.x"
                        eval $exec_str < $DIR_UTILS/temp_L.txt
                        
                        mv "$DIR_MAIN/F_temp.txt" "$DIR_UTILS/F_temp.txt"
                        
                        # Read new JXUV value
                        read JXUV < $DIR_UTILS/F_temp.txt
                        
                        # Remove excess files
                        rm $DIR_UTILS/temp_L.txt
                        rm $DIR_UTILS/L_temp.x
                        rm $DIR_UTILS/F_temp.txt
                        
                  elif [ $ans == "y" ]; then # Use tabulated luminosities
                        
                        i=1
                  
                  else # Invalid answer. Ask again
                        echo "Insert a valid answer"
                  fi
            
            done 
            
      else  # Manually insert parameters

            echo "Insert the physical parameters of the planet:"
            echo " "
            echo "Planet name:"
            read p_name
            
            # Compile fortran program to input the data
            par_temp_str="gfortran ${DIR_UTILS}parameters_input.f90 -o ${DIR_UTILS}temp_input.x"
            eval $par_temp_str
            
            # Execute fortran program to input the data
            exec_str="${DIR_UTILS}temp_input.x"
            eval $exec_str
            
            # Move temporary file to the correct location
            mv "$DIR_MAIN/temp_input_params.txt" "$DIR_UTILS/temp_input_params.txt"
            
            #Remove fortran executable
            rm /$DIR_UTILS/temp_input.x
            read -r R0 Mp T0 LX LEUV a JXUV Mtilde atilde < $DIR_UTILS/temp_input_params.txt

            # Remove temporary reading file
            rm /$DIR_UTILS/temp_input_params.txt
            
            # Write temporary file with planet parameters
            par_temp_str="$p_name"
            par_temp_str="$par_temp_str	$R0"
            par_temp_str="$par_temp_str	$Mp"
            par_temp_str="$par_temp_str	$T0"
            par_temp_str="$par_temp_str	$LX"
            par_temp_str="$par_temp_str	$LEUV"
            par_temp_str="$par_temp_str	$a"
            par_temp_str="$par_temp_str	$JXUV"
            par_temp_str="$par_temp_str	$Mtilde"
            par_temp_str="$par_temp_str	$atilde"
            par_temp_str="$par_temp_str	$b0"
            par_temp_str="$par_temp_str	$rho"
            
            echo "$par_temp_str" > $DIR_UTILS/temp.txt
            
      fi



      ############# Writing parameters file ###############


      # Write temporary file with second part of parameters.f90 file
      # Leave the spaces in this order to obtain an ordered ouput file
      echo "      ! Planetary parameters
      real*8 ::  n0 = 1.0d14			! Density at lower boundary (cm^-3)
      real*8 ::  R0 = $R0                 ! Planetary radius (cm)
      real*8 ::  Mp = $Mp 		      ! Planet mass (g)
      real*8 ::  T0 = $T0 			! Temperature at lower boundary (K)
      real*8 ::  Lx = $LX        	      ! Log10 X-Ray luminosity (erg/s)
      real*8 ::  LEUV = $LEUV             ! Log10 EUV luminosity (erg/s)
      real*8 ::  J_XUV = $JXUV  	      ! Bolometric star XUV flux (erg/cm^2*s)
      real*8 ::  Mrapp = $Mtilde 		! Ratio M_star/M_p
      real*8 ::  atilde = $atilde   	! Orbital radius in unit of R0 (a/R0)
      
      contains
      
      end module global_parameters" > $DIR_UTILS/parameters_fill.f90

      # Cat with template file
      cat $DIR_UTILS/parameters_temp.f90 $DIR_UTILS/parameters_fill.f90 >  $DIR_SRC/parameters.f90

      # Remove temporary filling file
      rm $DIR_UTILS/parameters_fill.f90


      ############# Retrieve Salz Data ###############

      # In the case of custom planet or of a planet missing in Salz table, 
      # use a temporary plotting planet (GJ1214b)
      Salz_temp="${p_name}_Salz.dat"
      
      if [ -f "$DIR_SALZ/$Salz_temp" ]; then
      
            echo "Retrieving Salz data"
            Salz_file=$Salz_temp
            
      else
            echo "No Salz data for this planet. Using WASP-77Ab b data for plot."
            Salz_file="WASP-77Ab_Salz.dat"
      fi
      
      # Copy file to main directory
      cp "$DIR_SALZ/$Salz_file" "$DIR_MAIN/$Salz_file"

      echo "Done."


      ############# Modify python file ###############

      # Get line with the radius definition
      R0_line=$(awk '/R0 =/{ print NR; exit }' $DIR_MAIN/python_plots.py)

      # Correct radius in python file
      sed -i'.bak' "$R0_line s/.*/R0 =$R0/" $DIR_MAIN/python_plots.py

      # Get line with Salz data loading
      Salz_line=$(awk '/_Salz.dat/{ print NR; exit }' $DIR_MAIN/python_plots.py)

      #Replace with exact filename
      sed -i'.bak' "$Salz_line s/.*/A = np.genfromtxt('$Salz_file')/" $DIR_MAIN/python_plots.py

      # Get line with figure title
      subtitle_line=$(awk '/fig.suptitle/{ print NR; exit }' $DIR_MAIN/python_plots.py)

      # Replace figure title 
      sed -i'.bak' "$subtitle_line s/.*/fig.suptitle('$p_name')/" $DIR_MAIN/python_plots.py

      # Clear backup file
      rm python_plots.py.bak

else  # If using the old parameters.f90 file

      # Extract parameters from the old temporary file
      read -r p_name R0 Mp T0 LX LEUV a JXUV Mtilde atilde b0 rho < $DIR_UTILS/temp.txt
      echo " Using previously written parameters file for $p_name"      
fi


############# Echo planet name ###############

echo "
--- Simulating the atmosphere of $p_name ---
"
 
############# Executing fortran file ###############

echo "Starting Fortran simulation."

# Initialize string of commands with ignored files. 
# This ensures that these files are compiled first 
str=" gfortran -J"$DIR_MOD" -I"$DIR_MOD" -fopenmp -g -fbacktrace \
      $DIR_SRC/parameters.f90\
      $DIR_SRC/J_inc.f90\
      $DIR_SRC/sigma.f90\
      $DIR_SRC/sigma_HeI.f90\
      $DIR_SRC/set_energy_vectors.f90
      $DIR_SRC/define_grid.f90\
      $DIR_SRC/grav_field.f90\
      $DIR_SRC/set_IC.f90\
      $DIR_SRC/load_IC.f90\
      $DIR_SRC/UW_conversions.f90\
      $DIR_SRC/Apply_BC.f90\
      $DIR_SRC/init.f90\
      $DIR_SRC/speed_estimate_HLLC.f90\
      $DIR_SRC/speed_estimate_ROE.f90\
      $DIR_SRC/Num_Fluxes.f90\
      $DIR_SRC/Source.f90\
      $DIR_SRC/RK_rhs.f90\
      $DIR_SRC/Cool_coeff.f90\
      $DIR_NLSOLVE/system_fcn.f90\
      $DIR_SRC/ionization_equilibrium.f90\
      $DIR_SRC/WENO3_rec.f90"


# List of fortran files containing modules ignoring files already compiled
str_ignore="$DIR_SRC/parameters.f90"
str_ignore="$str_ignore\|$DIR_SRC/define_grid.f90"
str_ignore="$str_ignore\|$DIR_SRC/J_inc.f90"
str_ignore="$str_ignore\|$DIR_SRC/sigma.f90"
str_ignore="$str_ignore\|$DIR_SRC/sigma_HeI.f90"
str_ignore="$str_ignore\|$DIR_SRC/set_energy_vectors.f90"
str_ignore="$str_ignore\|$DIR_SRC/grav_field.f90"
str_ignore="$str_ignore\|$DIR_SRC/set_IC.f90"
str_ignore="$str_ignore\|$DIR_SRC/load_IC.f90"
str_ignore="$str_ignore\|$DIR_SRC/UW_conversions.f90"
str_ignore="$str_ignore\|$DIR_SRC/Apply_BC.f90"
str_ignore="$str_ignore\|$DIR_SRC/init.f90"
str_ignore="$str_ignore\|$DIR_SRC/speed_estimate_HLLC.f90"
str_ignore="$str_ignore\|$DIR_SRC/speed_estimate_ROE.f90"
str_ignore="$str_ignore\|$DIR_SRC/Num_Fluxes.f90"
str_ignore="$str_ignore\|$DIR_SRC/Source.f90"
str_ignore="$str_ignore\|$DIR_SRC/RK_rhs.f90"
str_ignore="$str_ignore\|$DIR_SRC/Cool_coeff.f90"
str_ignore="$str_ignore\|$DIR_NLSOLVE/system_fcn.f900"
str_ignore="$str_ignore\|$DIR_SRC/ionization_equilibrium.f90"
str_ignore="$str_ignore\|$DIR_SRC/WENO3_rec.f90"



# Get list of all other files
files=$(ls $DIR_SRC/*.f90 | grep -v $str_ignore)


# Loop on modules files
for i in $files 
do
      str="$str $i "
done


# List of fortran files from the MINPACK library
files_nl=$(ls $DIR_NLSOLVE/*.f90 | grep -v "system_fcn.f90")

# Loop on MINPACK files
for i in $files_nl
do
      str="$str $i "
done

# Main file in program
main_file=$(ls | sed -n 's/\.f90$//p')

# Rename main file to match the planet name if the file has different name
if [[ ! "$main_file.f90" == "Hydro_ioniz_$p_name.f90" ]]; then

      mv "$main_file.f90" "Hydro_ioniz_$p_name.f90"
      
fi

# Find the new name of the main file without extension
main_file=$(ls | sed -n 's/\.f90$//p')

# Add main file to execution
str="$str $DIR_MAIN/$main_file.f90"

# Add output file to string of execution
compile="$str -o $DIR_MAIN/$main_file.x"

# Compile program and modules
eval $compile

# Exectution string
exec_str="$DIR_MAIN/$main_file.x"

# Exectute code
eval $exec_str

# Print when execution is over
echo "Fortran evaluation done."


############# Python plots ###############

# Run python_plots
#python3 -W ignore python_plots.py 









