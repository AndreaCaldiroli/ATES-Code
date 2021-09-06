import glob
import numpy as np
import tkinter as tk
import input_func as fc
from tkinter import ttk
import os

# Set global variables
glob.init()



#----- Constants -----#
MJ = 1.898e30           # Jupiter mass (g)
RJ = 6.9911e9           # Jupiter radius (cm) 
AU = 1.495978707e13     # Astronomical unit (cm)
Msun = 1.989e33         # Sun mass (g)
Gc = 6.67259e-8         # Gravitational constant
mu = 1.673e-24          # Proton mass
kb = 1.38e-16           # Boltzmann constant

#----- Data vectors for the table -----#

# Grid_types
grid_type = ["Uniform","Mixed","Stretched"]

# Numerical fluxes
num_flux = ["HLLC", 
            "ROE", 
            "LLF"]

# Reconstruction Methods
rec_meth = ["PLM", "WENO3"]

# 3D approximation method
glob.appx_meth = ["Mdot/4","Rate/2 + Mdot/2",
			"Rate/4 + Mdot","alpha = "]



#------------------
cwd = os.getcwd()
fil = cwd + '/src/utils/params_table.txt'

# Load lists
glob.pl_names,glob.Rp,\
              glob.Mp,\
              glob.T0,\
              glob.a,\
              glob.Ms,\
              glob.LX,\
              glob.LEUV = np.loadtxt(fil,\
              dtype = 'U9,f,f,f,f,f,f,f', unpack = True)
        
#------------------            


# Header
glob.note = ["Planet: ",
             "Grid type ",
             "Numerical flux ",
             "Reconstruction method ",
             "Log number density at R0 [cm^-3] ",
             "Planetary radius [R_J] ",
             "Planetary mass [M_J] ",
             "Equilibrium temperature [K] ",
             "Orbital distance [AU] ",
             "Parent star mass [M_sun] ",
             "Log10 of x-ray Luminosity ",
             "Log10 of EUV Luminosity ",
             "XUV flux [erg cm^-2 s^-1]",
             "Escape radius ",
             "He/H number ratio ",
             "2D aprrox. method",
             "alpha",
             "Spectral index"]



#----- Input interface -----#

# Open window
window = tk.Tk()
window.title(" INPUT PARAMETERS")

# Global frame
g_frame = tk.Frame(window)
g_frame.grid(row = 0, column = 0)    

# Number of input parameters
N_rows = len(glob.note)

# Create frames
fld_frame = tk.Frame(g_frame)
but_frame = tk.Frame(g_frame)
chk_frame = tk.Frame(g_frame)

# Position frames
fld_frame.grid(row = 0, column = 0)
but_frame.grid(row = 0, column = 1)
chk_frame.grid(row = 1, column = 0)

g_frame.grid_columnconfigure(0, weight = 1)
g_frame.grid_columnconfigure(1, weight = 1)
g_frame.grid_rowconfigure(0, weight = 1)
g_frame.grid_rowconfigure(1, weight = 1)


# Variables initialization
glob.b0_lbl_var    = tk.StringVar()
glob.lgphi_lbl_var = tk.StringVar()
glob.rho_lbl_var   = tk.StringVar()
glob.resc_lbl_var  = tk.StringVar()
glob.b0_lbl_var.set('b0 = ')
glob.lgphi_lbl_var.set('log(-phi) = ') 
glob.rho_lbl_var.set('rho = ')
glob.resc_lbl_var.set(glob.note[13])

glob.EUV_var = tk.IntVar()
glob.IC_var  = tk.IntVar()
glob.EUV_var.set(0)
glob.IC_var.set(0) 

# Fill labels
b0_lbl    = tk.Label(but_frame,text = glob.b0_lbl_var.get())
lgphi_lbl = tk.Label(but_frame,text = glob.lgphi_lbl_var.get())
rho_lbl   = tk.Label(but_frame,text = glob.rho_lbl_var.get())

           
# Pack labels         
b0_lbl.pack()
lgphi_lbl.pack()
rho_lbl.pack()


#---- Insert labels in fields frame ----#

for i in range(N_rows):
      if (i != 13):
            lbl = tk.Label(fld_frame, text = glob.note[i])
            lbl.grid(row = i, column = 0, sticky = "e", pady = 5) 
      fld_frame.grid_rowconfigure(i, weight = 1)
fld_frame.grid_columnconfigure(0, weight = 1)      
fld_frame.grid_columnconfigure(1, weight = 1)      


resc_lbl = tk.Label(fld_frame,text = glob.resc_lbl_var.get())
resc_lbl.grid(row = 13, column = 0, sticky = "e", pady = 5) 

#---- Create empty Entries and Comboboxes ----#

# List of planets
Plist = ttk.Combobox(fld_frame, values = glob.pl_names.tolist())
Plist.grid(row = 0, column = 1, sticky = "w", pady = 5)

# Multiple choice for grid type
Gridlist = ttk.Combobox(fld_frame, values = grid_type)
Gridlist.grid(row = 1, column = 1, sticky = "w", pady = 5)

# Multiple choice for numerical flux
Fluxlist = ttk.Combobox(fld_frame, values = num_flux)
Fluxlist.grid(row = 2, column = 1, sticky = "w", pady = 5)

# Multiple choice for reconstruction method
Reclist = ttk.Combobox(fld_frame, values = rec_meth)
Reclist.grid(row = 3, column = 1, sticky = "w", pady = 5)

# n0 entry box
ent_n0 = tk.Entry(fld_frame)
ent_n0.grid(row = 4, column = 1, sticky = "w", pady = 5)

# Rp entry box
ent_Rp = tk.Entry(fld_frame)
ent_Rp.grid(row = 5, column = 1, sticky = "w", pady = 5)

# Mp entry box
ent_Mp = tk.Entry(fld_frame)
ent_Mp.grid(row = 6, column = 1, sticky = "w", pady = 5)
      
# T0 entry box
ent_T0 = tk.Entry(fld_frame)
ent_T0.grid(row = 7, column = 1, sticky = "w", pady = 5)

# a entry box
ent_a = tk.Entry(fld_frame)
ent_a.grid(row = 8, column = 1, sticky = "w", pady = 5)

# Ms entry box
ent_Ms = tk.Entry(fld_frame)
ent_Ms.grid(row = 9, column = 1, sticky = "w", pady = 5)

# LX entry box
ent_LX = tk.Entry(fld_frame)
ent_LX.grid(row = 10, column = 1, sticky = "w", pady = 5)

# LEUV entry box
ent_LEUV = tk.Entry(fld_frame)
ent_LEUV.grid(row = 11, column = 1, sticky = "w", pady = 5)

# FXUV entry box
ent_FXUV = tk.Entry(fld_frame)
ent_FXUV.grid(row = 12, column = 1, sticky = "w", pady = 5)

# r_escape entry box
ent_resc = tk.Entry(fld_frame)
ent_resc.grid(row = 13, column = 1, sticky = "w", pady = 5)

# He/H entry box
ent_HeH = tk.Entry(fld_frame)
ent_HeH.grid(row = 14, column = 1, sticky = "w", pady = 5)

# 3D approximation method
ApxMth = ttk.Combobox(fld_frame, values = glob.appx_meth)
ApxMth.grid(row = 15, column = 1, sticky = "w", pady = 4)

# 2D approximation method
ent_alpha = tk.Entry(fld_frame)
ent_alpha.grid(row = 16, column = 1, sticky = "w", pady = 4)

# Power-law spectrum index
ent_PLind = tk.Entry(fld_frame)
ent_PLind .grid(row = 17, column = 1, sticky = "w", pady = 4)

# Create list of Comboboxes and Entries
glob.Cbox_list = [Plist,Gridlist,Fluxlist,Reclist,ApxMth]
glob.ent_list  = [ent_n0,	# 0
			ent_Rp,	# 1
			ent_Mp,	# 2
			ent_T0,	# 3
			ent_a,	# 4
                  ent_Ms,	# 5 
                  ent_LX,	# 6
                  ent_LEUV,	# 7
                  ent_FXUV,	# 8
                  ent_resc,	# 9
                  ent_HeH,	# 10
                  ent_alpha,	# 11
                  ent_PLind]	# 12

# Update fluxes
fc.upd_F()

#---- Create checkbuttons ----#

# EUV checkboxes
def cbEUV_upd(*args):
      fc.cbEUV()
cEUV = tk.Checkbutton(chk_frame, 
                      text = 'Only EUV', 
                      variable = glob.EUV_var,
                      command = cbEUV_upd)

# IC loading checkboxes

cIC = tk.Checkbutton(chk_frame, 
                     text = 'Load IC', 
                     variable = glob.IC_var)

cEUV.pack(side = tk.LEFT)
cIC.pack(side = tk.RIGHT)


#---- Create buttons and button functions ----#

# Close window function
def close_func(*args):
      window.destroy()

# Done functions
glob.write_ans = tk.StringVar()
def done_func(*args):
      fc.save_func()
      
      # Rename file
      if os.path.isfile("input_temp.inp"):                  
            os.remove("input_temp.inp")
                  
      close_func()
      
B_done = tk.Button(but_frame,text = "Done" , 
                   command = done_func, 
                   height = 2, width = 10)
B_rest = tk.Button(but_frame,text = "Reset",
                   command = fc.reset_func,
                   height = 2, width = 10)
B_clos = tk.Button(but_frame,text = "Close",
                   command = close_func, 
                   height = 2, width = 10)
B_addP = tk.Button(but_frame,text = "Add Planet", 
                   command = fc.add_planet,
                   height = 2, width = 10)

# Pack buttons
B_done.pack()
B_rest.pack()
B_clos.pack()
B_addP.pack()


#---- Various functions ----#

# Update labels
def upd_labels(*args):
      
      lk = 0
      for j in np.arange(1,5,1):
            if(len(glob.ent_list[j].get()) == 0):
                  lk += 1 
                  
      # Update only if all fields are complete
      if (lk == 0):      
            fc.upd_lbl()
            b0_lbl["text"]    = glob.b0_lbl_var.get()
            lgphi_lbl["text"] = glob.lgphi_lbl_var.get()
            rho_lbl["text"]   = glob.rho_lbl_var.get()
            resc_lbl["text"]  = glob.resc_lbl_var.get()

      
# Update entries if planet is selected
def Cbox_select(*args):
      fc.FillParams()
      upd_labels()



#---- Check if an input.inp file exists and load ----#

if os.path.isfile("input.inp") or os.path.isfile("input_temp.inp"):
      
      # Rename file
      if os.path.isfile("input.inp"):
            os.rename("input.inp", "input_temp.inp")
          
      # Clear current Comboboxes and entries
      for cb in glob.Cbox_list:
            cb.set('')
      
      for e in glob.ent_list:
            e.delete(0,'end')
          
      # Remove blank space function
      def rem_n(b):
            a = b.rstrip('\n')
            return a
            
      with open("input_temp.inp",'r') as fi:
            glob.Cbox_list[0].insert(0,rem_n(fi.readline()))
            glob.Cbox_list[1].insert(0,rem_n(fi.readline()))
            glob.Cbox_list[2].insert(0,rem_n(fi.readline()))
            glob.Cbox_list[3].insert(0,rem_n(fi.readline()))
            glob.ent_list[0].insert(0,rem_n(fi.readline()))
            glob.ent_list[1].insert(0,rem_n(fi.readline()))
            glob.ent_list[2].insert(0,rem_n(fi.readline()))
            glob.ent_list[3].insert(0,rem_n(fi.readline()))
            glob.ent_list[4].insert(0,rem_n(fi.readline()))
            glob.ent_list[5].insert(0,rem_n(fi.readline()))
            glob.ent_list[6].insert(0,rem_n(fi.readline()))
            glob.ent_list[7].insert(0,rem_n(fi.readline()))
            glob.ent_list[9].insert(0,rem_n(fi.readline()))
            glob.ent_list[10].insert(0,rem_n(fi.readline()))
            fi.readline()
            glob.Cbox_list[4].insert(0,rem_n(fi.readline()))
            glob.ent_list[11].insert(0,rem_n(fi.readline()))
            glob.ent_list[12].insert(0,rem_n(fi.readline()))
      
      # Update current flux
      fc.upd_F()
      
      # Update labels
      upd_labels()
      
      # Update the alpha box status
      fc.upd_alpha()

else: # If no input.inp exists, fill default entries

      fc.reset_func()


#---- Bind functions to comboboxes and entries ----#

glob.Cbox_list[0].bind("<<ComboboxSelected>>", Cbox_select)
glob.ent_list[1].bind("<Any-KeyRelease>", upd_labels)
glob.ent_list[2].bind("<Any-KeyRelease>", upd_labels)
glob.ent_list[3].bind("<Any-KeyRelease>", upd_labels)
glob.ent_list[4].bind("<Any-KeyRelease>", upd_labels)
glob.ent_list[5].bind("<Any-KeyRelease>", upd_labels)
glob.ent_list[6].bind("<Any-KeyRelease>", fc.upd_F)   # LX
glob.ent_list[7].bind("<Any-KeyRelease>", fc.upd_F)   # LEUV
glob.ent_list[8].bind("<Any-KeyRelease>", fc.upd_LEUV)
glob.Cbox_list[4].bind("<<ComboboxSelected>>", fc.upd_alpha)


window.columnconfigure(0, weight=1)
window.rowconfigure(0, weight=1)
# End of window code
window.mainloop()


