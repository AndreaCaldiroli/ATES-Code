# Definition of global variables

try:                        
	import tkinter as tk    # python 3 
	from tkinter import ttk
except ImportError:			# python 2
	import Tkinter as tk
	from Tkinter import ttk
	
def init():
	
	# Globals allocation
	global f_table
	global widgets, frame_list, lbl_spectrum
	global pl_names,sp_type,grid_type,num_flux,rec_meth
	global pl_params_list,pl_params
	global empty_lbl, PLind_spfile_lbl,resc_lbl
	global appx_meth
	global onlyEUV_var,LoadIC_var,He23S_var,onlyPP_var,force_var
	global MJ,RJ,AU,Msun,Gc,mu,kb
	global hp_eV,c_light   
	global current_spec
	global elow,emid,ehigh,ephot
	global EUV_band_lbl,Xray_band_lbl
	
	# Widget dictionaries
	widgets = {}
	
	# List of physical parameters
	pl_names = []
	appxmth = []
	pl_params = { 'Rp' : [],
			  'Mp' : [],
			  'T0' : [],
			  'LX' : [],
			  'LEUV' : [],
			  'a' : [],
			  'Ms' : []}
	
	# lit of specific planetary parameters
	pl_params_list = ['Rp',
				'Mp',
				'T0',
				'LX',
				'LEUV',
				'a',
				'Ms']
				 
	
	# --- Predefined inputs --- #
	
	# Spectrum evaluation
	sp_type = ['Load from file..',
		     'Power-law',
		     'Monochromatic']
	
	# Grid_types
	grid_type = ['Uniform','Mixed','Stretched']
	
	# Numerical fluxes
	num_flux = ['HLLC', 
				'ROE', 
				'LLF']
	      	
	# Reconstruction Methods
	rec_meth = ['PLM', 'WENO3']
	      	
	# 3D approximation method
	appx_meth = ['Mdot/4','Rate/2 + Mdot/2',
			 	 'Rate/4 + Mdot','alpha = ']      
	
	# Empty labels
	empty_lbl = { 'rho' : u'\u03C1\u209A: ',
			  	  'b0'  : u'\u03B2\u2080 = ', 
			  	  'phi' : u'log(\u03A6\u209A) = ',
			  	  'rochel'  : 'Roche lobe: '}
	
	# List of frames
	frame_list = ['Planet_params',
			  	  'Stellar_params',
			  	  'Others_params',
				  'Bands_info',
				  'Derived_parameters',
				  'Numerical_params',
				  'Tick_options',
				  'Buttons']	
			  
	# Current spectrum type
	current_spec = []
	
	# Labels for EUV info bands
	EUV_band_lbl  = u"EUV band:\n" + \
					u"[100,912] \u212B \u2263 [13.6,124] eV"
	Xray_band_lbl = u"X-ray band:\n" + \
					u"[10,100] \u212B \u2263 [124,1240] eV"
					 
	# Define useful constants
	MJ      = 1.898e30           # Jupiter mass (g)
	RJ      = 6.9911e9           # Jupiter radius (cm) 
	AU      = 1.495978707e13     # Astronomical unit (cm)
	Msun    = 1.989e33           # Sun mass (g)
	Gc      = 6.67259e-8         # Gravitational constant
	mu      = 1.673e-24          # Proton mass
	kb      = 1.38e-16           # Boltzmann constant
	hp_eV   = 4.1357e-15	     # Planck constant in eV
	c_light = 2.99792458e10      # Speed of light in cm/s

# ----
