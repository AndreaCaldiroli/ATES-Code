import glob
import numpy as np
import ATES_interface_functions as AIF
import os

try:                        
	import tkinter as tk    # python 3 
	from tkinter import ttk
	from tkinter import font
except ImportError:		# python 2
	import Tkinter as tk
	from Tkinter import ttk

# ----- Setting of global variables ----- #

# Load global variables
glob.init()

# Load lists of pre-tabulated planetary parameters
cwd = os.getcwd()
glob.f_table = cwd + '/src/utils/params_table.txt'
glob.pl_names,glob.pl_params['Rp'],\
              glob.pl_params['Mp'],\
              glob.pl_params['T0'],\
              glob.pl_params['a'], \
              glob.pl_params['Ms'],\
              glob.pl_params['LX'],\
              glob.pl_params['LEUV'] = np.loadtxt(glob.f_table,\
              dtype = 'U9,f,f,f,f,f,f,f', unpack = True)
              
# --- Various dimensions
fontsize = 13

# Borders
left_b = 40
top_b  = 30 

# Widgets widths and heights
lbl_w  = 220
lbl_h  = 30
ent_w  = 140
ent_h  = 30
chk_w  = 150
but_w  = 0.5*0.95*lbl_w
but_h  = 75

# Spacings
xspace = 10
yspace = 10

# Column start position
ys  = top_b
xcol1 = left_b
xcol2 = xcol1 + 0.9*lbl_w 
xcol3 = xcol2 + ent_w + 7.0*xspace # Some xspace added manually
xcol4 = xcol3 + lbl_w + xspace
xcol5 = xcol4 + ent_w + xspace
xcol6 = xcol5 + ent_w + xspace # Some xspace added manually
xcol7 = xcol6 + lbl_w + xspace

lblframe_w = lbl_w + ent_w + 30 

# Default properties for labels etc
def_font   = 'Times'
def_lbl    = dict(font = (def_font,fontsize), anchor = 'e')
def_lbl_w  = dict(font = (def_font,fontsize), anchor = 'w')
def_lbl_we = dict(font = (def_font,12), anchor = 'c')
def_but    = dict(font = (def_font,fontsize), anchor = 'center')
def_check  = dict(font = (def_font,fontsize), anchor = 'w')

# -------------------------------------------------------- #

# ----- Create window and frames ----- #

# Create window
window = tk.Tk()
window.title("ATES - ATmospheric EScape - Input parameters")
window.geometry("1200x480")
tk_font = font.nametofont("TkDefaultFont") 
tk_font.configure(family = def_font, size = 13)
tk_font.actual()

# ----- Create LabelFrames

planet_labelframe = ttk.LabelFrame(window,text = "Planet parameters")
planet_labelframe.place(x = left_b - 20, y = 10, width = lblframe_w, height = 270) 
stellar_labelframe = ttk.LabelFrame(window,text = "Stellar parameters")
stellar_labelframe.place(x = xcol3 - 20, y = 10, width = lblframe_w + 120, height = 270) 


# ----- Fill Planet frame ----- #

if "Planet_params" in glob.frame_list:

	# List of planets choice
	lbl = tk.Label(window, text = "Planet name: ", **def_lbl)
	lbl.place(x = xcol1, y = ys, width = 0.8*lbl_w, height = lbl_h)
	Plist = ttk.Combobox(window, values = glob.pl_names.tolist())
	Plist.place(x = xcol2, y = ys ,width = ent_w, height = ent_h)
	ys += lbl_h + yspace
	
	# ---

	# Log number density
	lbl = tk.Label(window, text = u"Log10 of n\u2080 [cm\u207B\u00B3]", **def_lbl)
	lbl.place(x = xcol1, y = ys, width = 0.8*lbl_w, height = lbl_h)
	ent_n0 = tk.Entry(window)
	ent_n0.place(x = xcol2, y = ys ,width = ent_w, height = ent_h)
	ys += lbl_h + yspace
	
	# ---

	# Planet radius
	lbl = tk.Label(window, text = "Planetary radius [R_J]", **def_lbl)
	lbl.place(x = xcol1, y = ys, width = 0.8*lbl_w, height = lbl_h)
	ent_Rp = tk.Entry(window)
	ent_Rp.place(x = xcol2, y = ys ,width = ent_w, height = ent_h)
	ys += lbl_h + yspace

	# ---

	# Planet mass
	lbl = tk.Label(window, text = "Planetary mass [M_J]", **def_lbl)
	lbl.place(x = xcol1, y = ys, width = 0.8*lbl_w, height = lbl_h)
	ent_Mp = tk.Entry(window)
	ent_Mp.place(x = xcol2, y = ys ,width = ent_w, height = ent_h)
	ys += lbl_h + yspace

	# ---

	# Equilibrium temperature 
	lbl = tk.Label(window, text = "Eq. temperature [K]", **def_lbl)
	lbl.place(x = xcol1, y = ys, width = 0.8*lbl_w, height = lbl_h)
	ent_T0 = tk.Entry(window)
	ent_T0.place(x = xcol2, y = ys ,width = ent_w, height = ent_h)
	ys += lbl_h + yspace

	# ---

	# Orbital distance 
	lbl = tk.Label(window, text = "Orbital distance [AU]", **def_lbl)
	lbl.place(x = xcol1, y = ys, width = 0.8*lbl_w, height = lbl_h)
	ent_a = tk.Entry(window)
	ent_a.place(x = xcol2, y = ys ,width = ent_w, height = ent_h)
	ys += lbl_h + yspace + 50
	
	# ------------------------------------------
	
	# Save position of the next frame
	yother = ys - 30

	# Other parameters labelframe
	other_labelframe = ttk.LabelFrame(window, text = "Other parameters")
	other_labelframe.place(x = left_b - 20, y = yother, width = lblframe_w, height = 160)
				
	# ------------------------------------------
				
	# Escape radius
	lbl = tk.Label(window, text = 'Escape Radius: ', **def_lbl)	 
	lbl.place(x = xcol1, y = ys, width = 0.8*lbl_w, height = 1.25*lbl_h)
	ent_resc = tk.Entry(window)
	ent_resc.place(x = xcol2, y = ys ,width = ent_w, height = ent_h)
	ys += lbl_h + yspace
	
	# ---
	
	# He/H number ratios
	lbl = tk.Label(window, text = "He/H number ratio", **def_lbl)
	lbl.place(x = xcol1, y = ys, width = 0.8*lbl_w, height = lbl_h)
	ent_heh = tk.Entry(window)
	ent_heh.place(x = xcol2, y = ys ,width = ent_w, height = ent_h)
	ys += lbl_h + yspace
	
	# ---
	
	# 2D approximation method
	y2D = ys	# Save for later use
	lbl = tk.Label(window, text = "2D approx. method", **def_lbl)
	lbl.place(x = xcol1, y = ys, width = 0.8*lbl_w, height = lbl_h)
	appxmth = ttk.Combobox(window, values = glob.appx_meth)
	appxmth.place(x = xcol2, y = ys ,width = ent_w, height = ent_h)
	ys += lbl_h + yspace

# ----- Fill Stellar frame ----- #

if "Stellar_params" in glob.frame_list:

	# Reset y position
	ys = top_b
	
	# ---
		
	# Star mass
	lbl = tk.Label(window, text = "Parent star mass [M_sun]", **def_lbl)
	lbl.place(x = xcol3, y = ys, width = lbl_w, height = lbl_h)
	ent_Ms = tk.Entry(window)
	ent_Ms.place(x = xcol4, y = ys ,width = ent_w, height = ent_h)
	ys += lbl_h + yspace
	
	# ---

	# Spectrum type
	lbl = tk.Label(window, text = "Spectrum type:", **def_lbl)
	lbl.place(x = xcol3, y = ys, width = lbl_w, height = lbl_h)
	spectrum = ttk.Combobox(window, values = glob.sp_type)
	spectrum.place(x = xcol4, y = ys, width = ent_w, height = ent_h)
	
	# Only EUV checkbox
	check_EUV = tk.Checkbutton(window,text = 'Only EUV', **def_check)
	check_EUV.place(x = xcol5, y = ys)
	ys += lbl_h + yspace
	
	# ---
	
	# Property of spectrum (PL index, file or energy)
	lbl_spectrum = tk.Label(window, text = '', **def_lbl)
	lbl_spectrum.place(x = xcol3, y = ys, width = lbl_w, height = lbl_h)
	ent_spec_prop = tk.Entry(window)
	ent_spec_prop.place(x = xcol4, y = ys, width = ent_w, height = ent_h)
	
	# ---
	
	# Browse button for spectrum file
	browse_but = tk.Button(window, text = "Browse..", font = (def_font,10), anchor = 'w')
	browse_but.place(x = xcol5, y = ys)
	ys += lbl_h + yspace

	# ---
	
	yenergy = ys - 30	# Save for later use
	
	# ---

	# X-ray luminosity
	#	ys    = y2D
	lbl = tk.Label(window, text = "Log10 of X-ray luminosity", **def_lbl) 
	lbl.place(x = xcol3, y = ys, width = lbl_w, height = lbl_h)
	ent_LX = tk.Entry(window)
	ent_LX.place(x = xcol4, y = ys, width = ent_w, height = ent_h)
	ys += lbl_h + yspace

	# ---

	# EUV-ray luminosity
	lbl = tk.Label(window, text = "Log10 of EUV luminosity", **def_lbl) 
	lbl.place(x = xcol3, y = ys, width = lbl_w, height = lbl_h)
	ent_LEUV = tk.Entry(window)
	ent_LEUV.place(x = xcol4, y = ys, width = ent_w, height = ent_h)
	ys += lbl_h + yspace

	# ---

	# XUV Flux
	lbl = tk.Label(window, text = u"XUV flux [erg cm\u207B\u00B2 s\u207B\u00B9]", **def_lbl) 
	lbl.place(x = xcol3, y = ys, width = lbl_w, height = lbl_h)
	ent_flux = tk.Entry(window)
	ent_flux.place(x = xcol4, y = ys, width = ent_w, height = ent_h)
	ys += lbl_h + yspace + 50

# ----- Fill Numerical parameters frame ----- #

if "Numerical_params" in glob.frame_list:

	# Derived paramters labelframe 
	numparams_labelframe = ttk.LabelFrame(window,text = "Numerical parameters")
	numparams_labelframe.place(x = xcol3 - 20, y = yother , width = lbl_w + 90, height = 160)
				
	# ---

	ygrid = ys # Save for later
		
	# Spatial grid
	lbl = tk.Label(window, text = "Grid type", **def_lbl)
	lbl.place(x = xcol3 , y = ys, width = lbl_w*0.55, height = lbl_h)
	grid = ttk.Combobox(window, values = glob.grid_type)
	grid.place(x = xcol3 + lbl_w*0.6, y = ys, width = 0.75*ent_w, height = ent_h)
	ys += lbl_h + yspace

	# ---
		
	# Numerical flux
	lbl = tk.Label(window, text = "Numerical flux", **def_lbl)
	lbl.place(x = xcol3 , y = ys, width = lbl_w*0.55, height = lbl_h)
	numflux = ttk.Combobox(window, values = glob.num_flux)
	numflux.place(x = xcol3 + lbl_w*0.6, y = ys, width = 0.75*ent_w, height = ent_h)
	ys += lbl_h + yspace

	# ---
		
	# Reconstruction method
	lbl = tk.Label(window, text = "Reconstruction", **def_lbl)
	lbl.place(x = xcol3  , y = ys, width = lbl_w*0.55, height = lbl_h)
	reconst = ttk.Combobox(window, values = glob.rec_meth)
	reconst.place(x = xcol3 + lbl_w*0.6, y = ys, width = 0.75*ent_w, height = ent_h)
	ys += lbl_h + yspace

# ----- Fill Tick options frame ----- #

if "Tick_options" in glob.frame_list:

	# Derived paramters labelframe 
	Tickopt_labelframe = ttk.LabelFrame(window,text = "Model options")
	Tickopt_labelframe.place(x = xcol4 + 70, y = yother, width = lbl_w - 30, height = 160)

	# Retrieve value of ygrid
	ys = ygrid - 10

	# Include He23S checkbutton
	check_He23S = tk.Checkbutton(window, text = u'Include HeI(2\u00B3S)', **def_check)
	check_He23S.place(x = xcol4 + 80, y = ys, width = chk_w, height = ent_h)
	ys += ent_h
	
	# Load IC checkbutton
	check_loadIC = tk.Checkbutton(window, text = 'Load IC', **def_check)
	check_loadIC.place(x = xcol4 + 80, y = ys, width = chk_w, height = ent_h)
	ys += ent_h

	# Only post processing
	check_onlyPP = tk.Checkbutton(window, text = 'Only post-proc.', **def_check)
	check_onlyPP.place(x = xcol4 + 80, y = ys, width = chk_w, height = ent_h)
	ys += ent_h

	# Force start
	check_force = tk.Checkbutton(window, text = 'Force start', **def_check)
	check_force.place(x = xcol4 + 80, y = ys, width = chk_w, height = ent_h)
	ys += ent_h

# ----- Fill derived parameters frame ----- #

if "Derived_parameters" in glob.frame_list:

	# Align with energy box
	ys = top_b

	# ---
		
	# Derived paramters labelframe 
	derparams_labelframe = ttk.LabelFrame(window,text = "Derived parameters")
	derparams_labelframe.place(x = xcol6 - 20, y = 10, width = 0.95*lbl_w, height = 150)

	# ---
		
	# Planet density
	lbl_rho = tk.Label(window, text = glob.empty_lbl['rho'], **def_lbl_w) 
	lbl_rho.place(x = xcol6, y = ys, width = 0.8*lbl_w, height = lbl_h)
	ys += lbl_h
	
	# ---
		
	# Jeans escape parameter
	lbl_b0 = tk.Label(window, text = glob.empty_lbl['b0'], **def_lbl_w) 
	lbl_b0.place(x = xcol6, y = ys, width = 0.8*lbl_w, height = lbl_h)
	ys += lbl_h
		
	# ---
			
	# log gravitational potential at surface
	lbl_phi = tk.Label(window, text = glob.empty_lbl['phi'], **def_lbl_w) 
	lbl_phi.place(x = xcol6, y = ys, width = 0.8*lbl_w, height = lbl_h)
	ys += lbl_h

	# ---
			
	# Roche Lobe distance
	lbl_resc = tk.Label(window, text = glob.empty_lbl['rochel'], **def_lbl_w) 
	lbl_resc.place(x = xcol6, y = ys, width = 0.8*lbl_w, height = lbl_h)
	ys += lbl_h + yspace
	
# ----- Info about energy bands ----- #

if "Bands_info" in glob.frame_list:

	# Go a little below
	ys += 5

	# Bands info labelframe 
	ebands_labelframe = ttk.LabelFrame(window,text = "Energy bands")
	ebands_labelframe.place(x = xcol6 - 20, y = ys, width = 0.95*lbl_w, height = 115)
	ys += 20

	# ---

	# EUV band info
	lbl_euv = tk.Label(window, text = glob.EUV_band_lbl, **def_lbl_we) 
	lbl_euv.place(x = xcol6 - 15, y = ys, width = 0.9*lbl_w, height = 1.5*lbl_h)
	ys += 1.5*lbl_h

	# Xray band info
	lbl_xray = tk.Label(window, text = glob.Xray_band_lbl, **def_lbl_we) 
	lbl_xray.place(x = xcol6 - 15, y = ys, width = 0.9*lbl_w, height = 1.5*lbl_h)
	ys += lbl_h

# ----- Fill buttons frame ----- #

if "Buttons" in glob.frame_list:

	# Align with LX entry
	ys = yother + 10

	# Add planet button - add to existing table
	add_but = tk.Button(window, text = "Add planet..", **def_but)
	add_but.place(x = xcol6 - 20, y = ys, width = but_w, height = but_h)
	
	# ---
		
	# Close button
	close_but = tk.Button(window, text = "Close",  **def_but)
	close_but.place(x = xcol6 - 20 + but_w, y = ys, width = but_w, height = but_h)
	ys += but_h #+ yspace
	
	# ---
		
	# Reset button
	reset_but = tk.Button(window, text = "Reset", **def_but)
	reset_but.place(x = xcol6 - 20, y = ys, width = but_w, height = but_h)
	
	# ---
		
	# Start button
	start_but = tk.Button(window, text = "Start", **def_but)
	start_but.place(x = xcol6 - 20 + but_w, y = ys, width = but_w, height = but_h)

# -------------------------------------------------------- #

# ----- Create dictionary of entries, comboboxes, checkbuttons and labels

glob.widgets['window']    = window
glob.widgets['planets']   = Plist 
glob.widgets['n0'] 	  	  = ent_n0 
glob.widgets['Rp'] 	  	  = ent_Rp 
glob.widgets['Mp'] 	  	  = ent_Mp 
glob.widgets['T0'] 	  	  = ent_T0 
glob.widgets['a'] 	  	  = ent_a 
glob.widgets['resc']      = ent_resc 
glob.widgets['heh']       = ent_heh 
glob.widgets['appxmth']   = appxmth 
glob.widgets['ICload']    = check_loadIC 
glob.widgets['He23S']	  =	check_He23S
glob.widgets['onlyPP']	  = check_onlyPP
glob.widgets['force']	  =	check_force
glob.widgets['Ms'] 	  	  = ent_Ms 
glob.widgets['spectrum']  = spectrum 	
glob.widgets['EUVonly']   = check_EUV  
glob.widgets['spec_prop'] = ent_spec_prop
glob.widgets['lbl_spec']  = lbl_spectrum
glob.widgets['browse']    = browse_but 
glob.widgets['LX'] 	  	  = ent_LX 
glob.widgets['LEUV']      = ent_LEUV 
glob.widgets['flux']      = ent_flux
glob.widgets['grid']      = grid 
glob.widgets['numflux']   = numflux 
glob.widgets['reconst']   = reconst 
glob.widgets['onlyPP']    = check_onlyPP
glob.widgets['force']     = check_force
glob.widgets['rho'] 	  = lbl_rho 
glob.widgets['b0']        = lbl_b0 
glob.widgets['phi']       = lbl_phi 
glob.widgets['resc_l']    = lbl_resc
glob.widgets['add_but']   = add_but 
glob.widgets['close_but'] = close_but  
glob.widgets['reset_but'] = reset_but 
glob.widgets['start_but'] = start_but 

# -------------------------------------------------------- #

# ----- Initialize and assign default values to some fields
AIF.init_func()

# --- Assign callbacks functions
glob.widgets['planets'].bind("<<ComboboxSelected>>", AIF.Plist_func)
glob.widgets['Rp'].bind("<Any-KeyRelease>", AIF.upd_labels)
glob.widgets['Mp'].bind("<Any-KeyRelease>", AIF.upd_labels)
glob.widgets['T0'].bind("<Any-KeyRelease>", AIF.upd_labels)
glob.widgets['a'].bind( "<Any-KeyRelease>", AIF.upd_flux)
glob.widgets['ICload'].configure(var = glob.LoadIC_var)
glob.widgets['He23S'].configure(var = glob.He23S_var)
glob.widgets['onlyPP'].configure(var = glob.onlyPP_var, command = AIF.onlyPP_func)
glob.widgets['force'].configure(var = glob.force_var, command = AIF.force_func)
glob.widgets['Ms'].bind("<Any-KeyRelease>", AIF.upd_labels)
glob.widgets['spectrum'].bind("<<ComboboxSelected>>", AIF.spectrum_func)
glob.widgets['EUVonly'].configure(var = glob.onlyEUV_var, command = AIF.onlyEUV_func)
glob.widgets['browse'].configure(command = AIF.browse_func)
glob.widgets['LEUV'].bind("<Any-KeyRelease>", AIF.upd_flux)
glob.widgets['LX'].bind("<Any-KeyRelease>", AIF.upd_flux)
glob.widgets['flux'].bind("<Any-KeyRelease>", AIF.flux_func) 
glob.widgets['add_but'].configure(command = AIF.add_func)
glob.widgets['reset_but'].configure(command = AIF.reset_func)
glob.widgets['close_but'].configure(command = AIF.close_func)
glob.widgets['start_but'].configure(command = AIF.start_func)

# End of window code
window.mainloop()
