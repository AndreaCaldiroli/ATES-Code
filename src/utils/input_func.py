# Library of useful functions for the input interface
import numpy as np
import tkinter as tk
import glob
import os
from shutil import copyfile

# Constants
MJ = 1.898e30           # Jupiter mass (g)
RJ = 6.9911e9           # Jupiter radius (cm) 
AU = 1.495978707e13     # Astronomical unit (cm)
Msun = 1.989e33         # Sun mass (g)
Gc = 6.67259e-8         # Gravitational constant
mu = 1.673e-24          # Proton mass
kb = 1.38e-16           # Boltzmann constant

# Update current flux
def upd_F(*args):
      
      a_str = glob.ent_list[4].get()
      if(len(a_str) == 0):
            a = 0.0
      else:
            a = float(a_str)
      
      # Get LX in double precision
      # if empty set to 0
      LX_str = glob.ent_list[6].get()
      if (len(LX_str) == 0):
            LX = 0.0
      else: 
            LX = float(LX_str)
      
      # If LX entry is disable, set 10**LX to 0
      if (glob.ent_list[6]["state"] == "disabled"):
            LX = -100000.0
                  
      # Get LEUV in double precision
      # if empty set to 0
      LEUV_str = glob.ent_list[7].get()
      if (len(LEUV_str) == 0):
            LEUV = 0.0
      else: 
            LEUV = float(LEUV_str)
      
      # Evaluate XUV flux value
      if(len(a_str) == 0):
            FXUV = ''
      else:
            FXUV = (10.0**LX + 10.0**LEUV)/(4.0*np.pi*a*a*AU*AU)
      
      # Delete current FXUV before updating
      glob.ent_list[8].delete(0,'end')
      
      # Insert FXUV value
      glob.ent_list[8].insert(0,FXUV)
	      

#-------------------------------------------

# Reset function
def reset_func():
      
      # Clear old inputs (both Combobox and entries
      for cb in glob.Cbox_list:
            cb.set('')
      for e in glob.ent_list:
            e.delete(0,'end')
      
      # Insert default fields
      glob.Cbox_list[0].set("Choose a planet..")
      glob.Cbox_list[1].set("Mixed")
      glob.Cbox_list[2].set("HLLC")
      glob.Cbox_list[3].set("WENO3")
      glob.ent_list[0].insert(0,14)
      glob.ent_list[9].insert(0,2.0)
      glob.ent_list[10].insert(0, 1.0/12.0)
      glob.Cbox_list[4].set("Mdot/4")
      glob.ent_list[12].insert(0,-1.0)
      
      # Update the alpha entry
      upd_alpha()
	
#-------------------------------------------

# Fill parameters if planet name is selected
def FillParams(*args):

      # Delete current inputs
      for e in glob.ent_list:
            e.delete(0,'end')
      
      # Insert default entries
      glob.ent_list[0].insert(0,14)
      glob.ent_list[9].insert(0,2.0)
      glob.ent_list[10].insert(0, 1.0/12.0)
      glob.ent_list[12].insert(0,-1.0)
	
      # Get number of selected planet in list
      p_num = glob.Cbox_list[0].current()
      
      # Insert planetary parameters
      glob.ent_list[1].insert(0,glob.Rp[p_num])
      glob.ent_list[2].insert(0,glob.Mp[p_num])
      glob.ent_list[3].insert(0,glob.T0[p_num])
      glob.ent_list[4].insert(0,glob.a[p_num])
      glob.ent_list[5].insert(0,glob.Ms[p_num])
      glob.ent_list[6].insert(0,glob.LX[p_num])
      glob.ent_list[7].insert(0,glob.LEUV[p_num])
      glob.ent_list[8].insert(0,glob.LEUV[p_num])
      
      # Update XUV flux value
      upd_F()




#-------------------------------------------

# Update planetary labels (r_RL,b0,log(phi))
def upd_lbl(*args):

      # Get number of selected planet in list
      p_num = glob.Cbox_list[0].current()
      
      # Get physical parameters            
      Rp = float(glob.ent_list[1].get())*RJ
      Mp = float(glob.ent_list[2].get())*MJ
      T0 = float(glob.ent_list[3].get())
      a  = float(glob.ent_list[4].get())*AU
      Ms = float(glob.ent_list[5].get())*Msun
      
      if (Rp != 0 and Mp != 0 and T0 != 0):
            # Evaluate Roche Lobe distance
            r_max = (3.0*Ms/Mp)**(-1.0/3.0)*(a/Rp)
            
            # Evaluate log of grav. pot. at surface       
            lgphi = np.log10(Gc*Mp/Rp)
            
            # Evaluate Jeans escape parameter
            b0 = (Gc*Mp*mu)/(kb*T0*Rp)        
                    
            # Evaluate density
            rho = Mp/(4.0/3.0*np.pi*Rp**3.0)        
      else:
            r_max = 0.0
            lgphi = 1.0
            b0    = 0.0
            rho   = 0.0
            
      
             
      glob.resc_lbl_var.set(glob.note[13] + "\n (" + "r_RL = "\
                         "{:7.4f}".format(r_max) + ")")
      
      glob.b0_lbl_var.set(" b0 = " + "{:7.2f}".format(b0))
      glob.lgphi_lbl_var.set("lg(-phi) = " + "{:5.2f}".format(lgphi))
      glob.rho_lbl_var.set("rho = " + "{:5.2f}".format(rho) + "g cm^-3")


#-------------------------------------------


# Update current flux
def upd_LEUV(*args):
      
      
      
      # Get LX state
      if (glob.ent_list[6]["state"] == "disabled"):
      
            # Evaluate orbital distance value
            a_str = glob.ent_list[4].get()
            if(len(a_str) == 0):
                  a = 0.0
            else:
                  a = float(a_str)*AU
            
            # Evaluate XUV flux value
            FXUV_str = glob.ent_list[8].get()
            if (len(FXUV_str) == 0):
                  FXUV = 0.0
                  LEUV = 1.0
            else: 
                  FXUV = float(FXUV_str)
                  
                  # Evaluate EUV flux   
                  LEUV = 4.0*np.pi*a*a*FXUV
                     
            

            # Delete current LEUV before updating
            glob.ent_list[7].delete(0,'end')
            
            # Update LEV value
            glob.ent_list[7].insert(0,np.log10(LEUV))
            
            
#-------------------------------------------           

# Activate and deactivate LX field
def cbEUV(*args):

      # Get button state
      st = glob.EUV_var.get()

      if (st == 1):
            glob.ent_list[6].delete(0,'end')
            glob.ent_list[6].insert(0,-100000.0)   
            glob.ent_list[6].config(state = tk.DISABLED)   
            
            upd_LEUV()
      if (st == 0):
            
            # Get number of selected planet in list
            p_num = glob.Cbox_list[0].current()

            glob.ent_list[6].config(state = tk.NORMAL)  
            glob.ent_list[6].delete(0,'end')
            glob.ent_list[7].delete(0,'end')
            glob.ent_list[6].insert(0,glob.LX[p_num])
            glob.ent_list[7].insert(0,glob.LEUV[p_num])
            
            upd_F()
            
            
#-------------------------------------------           

# Function for the done button
def save_func(*args):

      # Set variable for IC loading
      ans_IC = 'n'
      if (glob.IC_var.get() == 1):
            ans_IC = 'y'
            # Rename old output files for IC usage
            cdir = os.getcwd()
            cdir = cdir + "/output"
            old_hydro = os.path.join(cdir, 'Hydro_ioniz.txt')
            new_hydro = os.path.join(cdir, 'Hydro_ioniz_IC.txt')
            old_ioniz = os.path.join(cdir, 'Ion_species.txt')
            new_ioniz = os.path.join(cdir, 'Ion_species_IC.txt')      
            copyfile(old_hydro, new_hydro)
            copyfile(old_ioniz, new_ioniz)
            
      # Retrieve all fields needed to write the input.inp file
      glob.write_ans.set('y')
      for j in np.arange(1,7,1):
            lk = len(glob.ent_list[j].get())
            if (lk == 0):
                  glob.write_ans.set('n')
                  break
                  
      # Write only if all fields are complete
      if (glob.write_ans.get() == 'y'):
      
            # Write to file
            with open('input.inp', 'w') as f:
                      
                  f.write(("%s\n%s\n%s\n%s\n"
                           "%7.1E\n%5.3f\n%5.3f\n"
                           "%6.1f\n%6.4f\n%5.3f\n"
                           "%6.3f\n%6.3f\n"
                           "%4.1f\n%8.6f\n%s\n"
                           "%s\n%6.3f\n%5.3f")
                           
                          %(glob.Cbox_list[0].get(),
                            glob.Cbox_list[1].get(),
                            glob.Cbox_list[2].get(),
                            glob.Cbox_list[3].get(),
                            float(glob.ent_list[0].get()),
                            float(glob.ent_list[1].get()),
                            float(glob.ent_list[2].get()),
                            float(glob.ent_list[3].get()),
                            float(glob.ent_list[4].get()),
                            float(glob.ent_list[5].get()),
                            float(glob.ent_list[6].get()),
                            float(glob.ent_list[7].get()),
                            float(glob.ent_list[9].get()),
                            float(glob.ent_list[10].get()),
                            ans_IC,
                            glob.Cbox_list[4].get(),
                            float(glob.ent_list[11].get()),
                            float(glob.ent_list[12].get())
                            )
                            )
                           
            
      else:
            print('Some fields are not complete')       


#-------------------------------------------

# Add planet to table in file
def add_planet(*args):
      
      # Get table path
      cdir = os.getcwd()
      tfile = cdir + "/src/utils/params_table.txt"
	
      for j in glob.ent_list:
            if(len(j.get()) == 0):
                  print("Complete all fields before saving")
                  return

      # Retrieve data written in table
      pl_name = glob.Cbox_list[0].get()
      Rp = float(glob.ent_list[1].get())
      Mp = float(glob.ent_list[2].get())
      T0 = float(glob.ent_list[3].get())
      a  = float(glob.ent_list[4].get())
      Ms = float(glob.ent_list[5].get())
      LX = float(glob.ent_list[6].get())
      LE = float(glob.ent_list[7].get())
      
      with open(tfile, 'a') as f:
     
            f.write(("%s\t%5.3f\t%6.3f\t"
                     "%6.1f\t%6.4f\t%5.3f\t"
                     "%6.3f\t%6.3f\n")
                     
                    %(pl_name,Rp,Mp,
                      T0,a,Ms,LX,LE))
      
      print("Planet added to table.")


#-------------------------------------------

# Update value of alppha based on the 3D method chosenb

def upd_alpha(*args):
	
	app_m = glob.Cbox_list[4].get()
	
	# Activate alpha entry field
	glob.ent_list[11].delete(0,'end')
	glob.ent_list[11].insert(0,0)
	glob.ent_list[11].config(state = tk.DISABLED)  
	
	if app_m == "alpha = ":
	
		glob.ent_list[11].config(state = tk.NORMAL)  
		glob.ent_list[11].delete(0,'end')
		glob.ent_list[11].insert(0,4)
		
		
		
		
	
































