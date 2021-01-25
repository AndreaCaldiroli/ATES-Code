import numpy as np
import tkinter as tk
from tkinter import ttk

# Corot 2 b non ha Rp sul catalogo!

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
grid_type = ["Mixed","Stretched"]

# Numerical fluxes
num_flux = ["HLLC", 
            "ROE", 
            "LLF",
            "LF"]

# Reconstruction Methods
rec_meth = ["PLM", "WENO3"]

            
# Planet names
pl_names = ["WASP-12b",  # WASP-12b
            "GJ3470b" ,  # GJ3470b
            "WASP-80b",  # WASP-80b
            "HD149026b", # HD149026b
            "HAT-P-11b", # HAT-P-11b
            "HD209458b", # HD209458b
            "55Cnce",    # 55Cnce
            "GJ1214b",   # GJ1214b
            "GJ436b",    # GJ436b
            "HD189733b", # HD189733b
            "HD97658b",  # HD97658b
            "WASP-77b",  # WASP-77b
            "WASP-43b",  # WASP-43b
            "Corot2b",   # Corot2b
            "WASP-8b",   # WASP-8b
            "WASP-10b",  # WASP-10b
            "HAT-P-2b",  # HAT-P-2b
            "HAT-P-20b", # HAT-P-20b
            "WASP-38b",  # WASP-38b
            "WASP-18b",  # WASP-18b
            "55Cncb"]    # 55Cncb
      
      
# Planet radius
Rp = [1.8,   # WASP-12b
      0.37,  # GJ3470b
      0.95,  # WASP-80b
      0.65,  # HD149026b
      0.42,  # HAT-P-11b
      1.4,   # HD209458b
      0.19,  # 55Cnce
      0.24,  # GJ1214b
      0.38,  # GJ436b
      1.1,   # HD189733b
      0.21,  # HD97658b
      1.2,   # WASP-77b
      0.93,  # WASP-43b
      1.5,   # Corot2b
      1.0,   # WASP-8b
      1.1,   # WASP-10b
      1.2,   # HAT-P-2b
      0.87,  # HAT-P-20b
      1.1,   # WASP-38b
      1.3,   # WASP-18b
      1.466] # 55Cncb



# Planet mass
Mp = [1.4,    # WASP-12b
      0.044,  # GJ3470b
      0.55,   # WASP-80b
      0.36,   # HD149026b
      0.083,  # HAT-P-11b
      0.69,   # HD209458b
      0.026,  # 55Cnce
      0.020,  # GJ1214b
      0.073,  # GJ436b
      1.1,    # HD189733b
      0.025,  # HD97658b
      1.8,    # WASP-77b
      1.8,    # WASP-43b
      3.3,    # Corot2b
      2.2,    # WASP-8b
      3.2,    # WASP-10b
      8.9,    # HAT-P-2b
      7.3,    # HAT-P-20b
      2.7,    # WASP-38b
      10.2,   # WASP-18b
      0.84]    # 55Cncb


# Planet surface temperature
T0 = [2900.0,  # WASP-12b
      650.0,   # GJ3470b
      800.0,   # WASP-80b
      1440.0,  # HD149026b
      850.0,   # HAT-P-11b
      1320.0,  # HD209458b
      1950.0,  # 55Cnce
      550.0,   # GJ1214b
      650.0,   # GJ436b
      1200.0,  # HD189733b
      750.0,   # HD97658b
      1650.0,  # WASP-77b
      1350.0,  # WASP-43b
      1550.0,  # Corot2b
      950.0,   # WASP-8b
      950.0,   # WASP-10b
      1700.0,  # HAT-P-2b
      950.0,   # HAT-P-20b
      1250.0,  # WASP-38b
      2400.0,  # WASP-18b
      700.0]   # 55Cncb


# Log10 stellar X-ray luminosity
LX = [27.58, # WASP-12b
      27.63, # GJ3470b
      27.85, # WASP-80b
      28.60, # HD149026b
      27.55, # HAT-P-11b
      26.40, # HD209458b
      26.65, # 55Cnce
      25.91, # GJ1214b
      25.96, # GJ436b
      28.18, # HD189733b
      27.22, # HD97658b
      28.13, # WASP-77b
      27.88, # WASP-43b
      29.32, # Corot2b
      28.45, # WASP-8b
      28.09, # WASP-10b
      28.91, # HAT-P-2b
      28.00, # HAT-P-20b
      28.04, # WASP-38b
      26.82, # WASP-18b
      26.65] # 55Cncb



# Log10 stellar X-ray luminosity
LEUV = [28.35, # WASP-12b
        28.37, # GJ3470b
        28.46, # WASP-80b
        28.80, # HD149026b
        28.33, # HAT-P-11b
        27.84, # HD209458b
        27.66, # 55Cnce
        26.61, # GJ1214b
        27.14, # GJ436b
        28.61, # HD189733b
        28.19, # HD97658b
        28.59, # WASP-77b
        28.48, # WASP-43b
        29.13, # Corot2b
        28.73, # WASP-8b
        28.57, # WASP-10b
        28.94, # HAT-P-2b
        28.53, # HAT-P-20b
        28.55, # WASP-38b
        28.02, # WASP-18b
        27.06] # 55Cncb


# Planet orbital distance
a = [ 0.023, # WASP-12b
      0.036, # GJ3470b
      0.034, # WASP-80b
      0.043, # HD149026b
      0.053, # HAT-P-11b
      0.047, # HD209458b
      0.015, # 55Cnce
      0.014, # GJ1214b
      0.029, # GJ436b
      0.031, # HD189733b
      0.080, # HD97658b
      0.024, # WASP-77b
      0.014, # WASP-43b
      0.028, # Corot2b
      0.080, # WASP-8b
      0.038, # WASP-10b
      0.068, # HAT-P-2b
      0.036, # HAT-P-20b
      0.076, # WASP-38b
      0.020, # WASP-18b
      0.113] # 55Cncb



# Stellar masses
Ms = [1.434, # WASP-12b
      0.510, # GJ3470b
      0.580, # WASP-80b
      1.300, # HD149026b
      0.809, # HAT-P-11b
      1.148, # HD209458b
      1.015, # 55Cnce
      0.150, # GJ1214b
      0.452, # GJ436b
      0.800, # HD189733b
      0.850, # HD97658b
      1.002, # WASP-77b
      0.717, # WASP-43b
      0.970, # Corot2b
      1.033, # WASP-8b
      0.710, # WASP-10b
      1.340, # HAT-P-2b
      0.756, # HAT-P-20b
      1.216, # WASP-38b
      1.240, # WASP-18b
      1.015] # 55Cncb

# Jeans escape parameter
b0 = [(Gc*Mp[j]*mu*MJ)/(kb*T0[j]*Rp[j]*RJ) for j in range(len(Mp))]

# Density
rho = [Mp[j]*MJ/(4./3.*np.pi*(Rp[j]*RJ)**3.) for j in range(len(Mp))]


# Header
note = ["Planet: ",
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
        "He/H number ratio "]



#----- Input interface -----#

# Open window
window = tk.Tk()
window.title(" INPUT_PARAMETERS")


# Number of input parameters
N_rows = len(note)
fields_frame = tk.Frame(window)
fields_frame.grid(row = 0, column = 0)

for i in range(N_rows):
      
      # Insert labels
      lbl = tk.Label(fields_frame, text = note[i])
      lbl.grid(row = i, column = 0, sticky = "e", pady = 5)

def reset_fields():
      
      
      # Delete existing Entry boxes
      for ch in fields_frame.winfo_children():
            if (ch.winfo_class() == "Entry"):
                  ch.delete(0,'end')
      
      def pl_selected (self):
      
            # Delete specific fields
            ent_Rp.delete(0,'end') 
            ent_Mp.delete(0,'end')
            ent_T0.delete(0,'end')
            ent_a.delete(0,'end')
            ent_Ms.delete(0,'end')
            ent_LX.delete(0,'end')
            ent_LEUV.delete(0,'end')
            ent_FXUV.delete(0,'end')
            
            
            # Get number in list of selected planet
            p_num = Plist.current()
            
            # Set all the fields for the selected planet
            ent_Rp.insert(0,Rp[p_num])
            ent_Mp.insert(0,Mp[p_num])
            ent_T0.insert(0,T0[p_num])
            ent_a.insert(0,a[p_num])
            ent_Ms.insert(0,Ms[p_num])
            ent_LX.insert(0,LX[p_num])
            ent_LEUV.insert(0,LEUV[p_num])
            F_XUV = (10.0**LX[p_num]+10.0**LEUV[p_num])/\
                    (4.0*np.pi*a[p_num]*a[p_num]*AU*AU)
            ent_FXUV.insert(0,F_XUV)
      
      
      # Update flux if luminosities are updated
      def upd_F(*args):
            
            # Delete current FXUV before updating
            ent_FXUV.delete(0,'end')
            
            # Get LX in double precision
            # if empty set to 0
            if (len(ent_LX.get()) == 0):
                  LX = 0.0
            else: 
                  LX = float(ent_LX.get())
            
            # Get LEUV in double precision
            # if empty set to 0
            if (len(ent_LEUV.get()) == 0):
                  LEUV = 0.0
            else: 
                  LEUV = float(ent_LEUV.get())
      
            # If LX entry is disable, set 10**LX to 0
            if (ent_LX["state"] == "disabled"):
                  LX = -100000.0
            
            # Get current orbital distance 
            a_p = float(ent_a.get())*AU
                 
            # Evaluate XUV flux 
            F_XUV = (10.0**LX+10.0**LEUV)/(4.0*np.pi*a_p*a_p)
            
            # Insert FXUV value
            ent_FXUV.insert(0,F_XUV)      
            
      # Create empty entry boxes
      ent_n0 = tk.Entry(fields_frame)
      ent_n0.grid(row = 4, column = 1, sticky = "w", pady = 5)
      ent_n0.insert(0,14)

      ent_Rp = tk.Entry(fields_frame)
      ent_Rp.grid(row = 5, column = 1, sticky = "w", pady = 5)

      ent_Mp = tk.Entry(fields_frame)
      ent_Mp.grid(row = 6, column = 1, sticky = "w", pady = 5)

      ent_T0 = tk.Entry(fields_frame)
      ent_T0.grid(row = 7, column = 1, sticky = "w", pady = 5)

      ent_a = tk.Entry(fields_frame)
      ent_a.grid(row = 8, column = 1, sticky = "w", pady = 5)

      ent_Ms = tk.Entry(fields_frame)
      ent_Ms.grid(row = 9, column = 1, sticky = "w", pady = 5)

      ent_LX = tk.Entry(fields_frame)
      ent_LX.grid(row = 10, column = 1, sticky = "w", pady = 5)
      ent_LX.bind("<Any-KeyRelease>",upd_F)

      ent_LEUV = tk.Entry(fields_frame)
      ent_LEUV.grid(row = 11, column = 1, sticky = "w", pady = 5)
      ent_LEUV.bind("<Any-KeyRelease>",upd_F)
      
      ent_FXUV = tk.Entry(fields_frame)
      ent_FXUV.grid(row = 12, column = 1, sticky = "w", pady = 5)
      
      
      ent_resc = tk.Entry(fields_frame)
      ent_resc.grid(row = 13, column = 1, sticky = "w", pady = 5)
      ent_resc.insert(0,2.0)

      ent_HeH = tk.Entry(fields_frame)
      ent_HeH.grid(row = 14, column = 1, sticky = "w", pady = 5)
      ent_HeH.insert(0, 1.0/12.0)

      # Multiple choice for planet name
      Plist = ttk.Combobox(fields_frame, values = pl_names)
      Plist.set("Choose a planet..")
      Plist.grid(row = 0, column = 1, sticky = "w", pady = 5)
      Plist.bind("<<ComboboxSelected>>", pl_selected)
      
      # Multiple choice for grid type
      Gridlist = ttk.Combobox(fields_frame, values = grid_type)
      Gridlist.set("Mixed")
      Gridlist.grid(row = 1, column = 1, sticky = "w", pady = 5)

      # Multiple choice for numerical flux
      Fluxlist = ttk.Combobox(fields_frame, values = num_flux)
      Fluxlist.set("HLLC")
      Fluxlist.grid(row = 2, column = 1, sticky = "w", pady = 5)

      # Multiple choice for reconstruction method
      Reclist = ttk.Combobox(fields_frame, values = rec_meth)
      Reclist.set("WENO3")
      Reclist.grid(row = 3, column = 1, sticky = "w", pady = 5)
      
      list_box = [Gridlist, Fluxlist, Reclist,\
                  ent_n0,ent_Rp, ent_Mp, ent_T0,ent_a,\
                  ent_Ms,ent_LX,ent_LEUV,ent_resc,ent_HeH,Plist]
      return list_box
            
# Create table with default values
list_box = reset_fields()

# Frame of checkboxes
check_frame = tk.Frame(window)
check_frame.grid(row = 1, column = 0)

# Checkboxes

# Disable LX entry function
def cbEUV(*args):
      st = EUV_var.get()
      if (st == 1):
            list_box[9].delete(0,'end')
            list_box[9].insert(0,-100000.0)   
            list_box[9].config(state = tk.DISABLED)   
      if (st == 0):
            list_box[9].config(state = tk.NORMAL)  
            p_num = list_box[13].current()
            list_box[9].delete(0,'end')
            list_box[9].insert(0,LX[p_num])
            
# Checkbox for EUV luminosity
EUV_var = tk.IntVar()
c1 = tk.Checkbutton(check_frame, 
                    text = 'Only EUV', 
                    variable = EUV_var,
                    command = cbEUV)
c1.pack(side = tk.LEFT)

# Checkbutton for IC loading
IC_var = tk.IntVar()
c2 = tk.Checkbutton(check_frame, 
                    text='Load IC', 
                    variable = IC_var)
c2.pack(side = tk.RIGHT)


# Buttons frame
but_frame = tk.Frame(window)
but_frame.grid(row = 0, column = 1)

# Done button      
      
# Initialize input variables
grid = tk.StringVar()
flux = tk.StringVar()
Recm = tk.StringVar()
n0   = tk.DoubleVar()
Rp_inp = tk.DoubleVar()
Mp_inp = tk.DoubleVar()
T0_inp = tk.DoubleVar()
a_inp  = tk.DoubleVar()
Ms_inp = tk.DoubleVar()
LX_inp = tk.DoubleVar()
LEUV_inp = tk.DoubleVar()
resc_inp = tk.DoubleVar()
HeH_inp = tk.DoubleVar()
def close_func(*args):
      window.destroy()
      
def done_func(list_box):
      # Retrieve all fields needed to write the input.inp file
      lk = 0
      for j in np.arange(4,10,1):
            lk = lk + len(list_box[j].get())
            
      if(lk != 0):
            grid.set(list_box[0].get())
            flux.set(list_box[1].get())
            Recm.set(list_box[2].get())
            n0.set(float(list_box[3].get()))
            Rp_inp.set(float(list_box[4].get()))
            Mp_inp.set(float(list_box[5].get()))
            T0_inp.set(float(list_box[6].get()))
            a_inp.set(float(list_box[7].get()))
            Ms_inp.set(float(list_box[8].get()))
            LX_inp.set(float(list_box[9].get()))
            LEUV_inp.set(float(list_box[10].get()))
            resc_inp.set(float(list_box[11].get()))
            HeH_inp.set(float(list_box[12].get()))
      
      close_func()
      
   
      
      
B_ok = tk.Button(but_frame, 
                 text ="Done",
                 command = lambda: done_func(list_box))
B_ok.pack()


# Reset button     
def reset_func(*args):
      for ent in fields_frame.winfo_children():
            if (ent.winfo_class() != "Label"):
                  ent.destroy()
      reset_fields()  
B_reset = tk.Button(but_frame, 
                    text ="Reset",
                    command = reset_func)
B_reset.pack()

# Close Button

B_close = tk.Button(but_frame, 
                 text = "Close",
                 command = close_func)
B_close.pack()

#------------ TODO ----------#
# cambiare flusso/luminosità
# Leggere campi e scrivere file input.inp
# Vedere se funziona anche cambiando i vari parametri dei pianeti
# Metti anche log(phi) e b0

# End of window code
window.mainloop()


# Set variable for IC loading
ans_IC = 'n'
if (IC_var.get() == 1):
      ans_IC = 'y'


        
# Write to file
with open('input.inp', 'w') as f:
                
      f.write(("%s\n%s\n%s\n%7.1E\n%5.3f\n%5.3f\n"
               "%6.1f\n%6.4f\n%5.3f\n%6.3f\n%6.3f\n"
               "%4.1f\n%8.6f\n%s")
               
              %(grid.get(),
                flux.get(),
                Recm.get(),
                n0.get(),
                Rp_inp.get(),
                Mp_inp.get(),
                T0_inp.get(),
                a_inp.get(),
                Ms_inp.get(),
                LX_inp.get(),
                LEUV_inp.get(),
                resc_inp.get(),
                HeH_inp.get(),
                ans_IC))

