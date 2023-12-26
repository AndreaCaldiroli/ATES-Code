# Library of callback and other functions for ATES interface
import numpy as np
import glob
import os
from shutil import copyfile

try:                        
	import tkinter as tk    # python 3 
	from tkinter import ttk, filedialog
	from tkinter.filedialog import askopenfile
except ImportError:		# python 2
	import Tkinter as tk
	from Tkinter import ttk, filedialog
	from Tkinter.filedialog import askopenfile
	
#-------------------------------------------

# Define global tk variables
def allocate_tkvars():

	# Create checkbuttons variables
	glob.onlyEUV_var = tk.IntVar()
	glob.LoadIC_var  = tk.IntVar()
	glob.He23S_var   = tk.IntVar()
	glob.onlyPP_var  = tk.IntVar()
	glob.force_var   = tk.IntVar()
	
#-------------------------------------------

# Read word_number-th word from string_in
def get_word(string_in,word_number):
	
	# Initialize counters and string 
	count = 0 
	c_string = ''
	c_word_counter = 0
	string = string_in.strip()

	# Keep reading through string
	while count >= 0 and count <= len(string): 
		
		if count == len(string):
			out_word = c_string
			c_word_counter += 1
			# Return word
			if c_word_counter == word_number:
				return out_word	
		
		# If blank space or at the end of the string
		if string[count:count+1] != ' ':
			c_string += string[count:count+1] 
			count += 1
		else: 
			if string[count-1:count] != ' '  or count == len(string):
				# Update counters and get word			
				out_word = c_string
				c_word_counter += 1
				# Reset reading string
				c_string = ''
				# Return word
				# Update counter
				count += 1
				if c_word_counter == word_number:
					return out_word
			else:
				count += 1
				continue

#-------------------------------------------

# Update the labels with derived parameters
def upd_labels(*args):
	
	# Read physical parameters
	Rp_r = glob.widgets['Rp'].get()
	Mp_r = glob.widgets['Mp'].get()
	T0_r = glob.widgets['T0'].get()
	Ms_r = glob.widgets['Ms'].get()
	a_r  = glob.widgets['a'].get()

	# Trim input strings
	Rp_r .strip()
	Mp_r.strip()
	T0_r.strip()
	Ms_r.strip()
	a_r.strip()

	# Update rho string
	if Rp_r != '' and Mp_r != '':

		# Get floats
		Rp = float(Rp_r)*glob.RJ
		Mp = float(Mp_r)*glob.MJ
		
		if Rp != 0.0 or Mp_r != 0.0:
			rho     = Mp/(4.0/3.0*np.pi*Rp**3.0)
			rho_str = glob.empty_lbl['rho'] + '{:5.2f}'.format(rho) + u' g cm\u207B\u00B3'
			glob.widgets['rho'].configure(text = rho_str)


	# Update beta0 string
	if Rp_r != '' and Mp_r != '' and T0_r != '':

		# Get floats
		Rp = float(Rp_r)*glob.RJ
		Mp = float(Mp_r)*glob.MJ
		T0 = float(T0_r)
		
		if Rp != 0.0 and Mp_r != 0.0 and T0_r != 0.0:
			beta0     = glob.Gc*Mp*glob.mu/(Rp*glob.kb*T0)
			beta0_str = glob.empty_lbl['b0'] + '{:7.2f}'.format(beta0)
			glob.widgets['b0'].configure(text = beta0_str)


	# Update logphi string
	if Rp_r != '' and Mp_r != '':
		
		# Get floats
		Rp = float(Rp_r)*glob.RJ
		Mp = float(Mp_r)*glob.MJ
		
		if Rp != 0.0 and Mp_r != 0.0:
			phi     = np.log10(glob.Gc*Mp/Rp)
			phi_str = glob.empty_lbl['phi'] + '{:5.2f}'.format(phi) + u' erg g\u207B\u00B9'
			glob.widgets['phi'].configure(text = phi_str)


	# Update Roche lobe string
	if Rp_r != '' and Mp_r != '' and Ms_r != '' and a_r != '':

		# Get floats
		Rp = float(Rp_r)*glob.RJ
		Mp = float(Mp_r)*glob.MJ
		Ms = float(Ms_r)*glob.Msun
		a  = float(a_r)*glob.AU
		
		if Rp != 0.0 and Mp != 0.0 and Ms != 0.0 and a != 0.0:
			r_RL = (3.0*Ms/Mp)**(-1.0/3.0)*(a/Rp)
			RL_str = glob.empty_lbl['rochel'] + '{:5.2f}'.format(r_RL) + u' R\u209A'
			glob.widgets['resc_l'].configure(text = RL_str)

#-------------------------------------------

# Update the global XUV flux
def upd_flux(*args):
	
	# Integrate if numerical spectrum in active
	if glob.widgets['spectrum'].get() == glob.sp_type[0]:
		
		# Get file name
		sp_file = glob.widgets['spec_prop'].get()

		# Return if sp_file is not specified
		if sp_file == '': return
		
		# Integrate and update luminosities
		isthere_sp_file = os.path.isfile(sp_file)
		if isthere_sp_file: integrate_spectrum(sp_file)
	
	# Get the current status of the LX entry box
	LX_status = glob.widgets['LX'].cget('state')

	# Get values of LEUV and a
	LEUV_r = glob.widgets['LEUV'].get()
	a_r    = glob.widgets['a'].get()
	
	# Return if there is any empty input
	if LEUV_r == '' or a_r == '': return

	# Convert to float
	LEUV = 10.0**float(LEUV_r)
	a    = float(a_r)*glob.AU
		
	# Return if there is any empty input
	if a == 0.0: return
	
	# Get XUV luminosity if enabled
	if LX_status == 'disabled':
		LX = 0.0
	else: 
		LX_r = glob.widgets['LX'].get()
		if LX_r == '': return
		LX   = 10.0**float(LX_r)
	
	# Calculate flux
	FXUV = (LEUV + LX)/(4.0*np.pi*a*a)
	
	# Put into formatted string
	FXUV_str = '{:10.3f}'.format(FXUV)
	
	# Clear flux entry
	glob.widgets['flux'].delete(0,'end')
	
	# Insert flux value into the entry
	glob.widgets['flux'].insert(0,FXUV_str.strip())
	
	# Update also labels
	upd_labels()

#-------------------------------------------

# Reset function
def reset_func(*args):
	
	# Clear old inputs and reset checkboxes
	for field in glob.widgets:
		
		# Current element of dictionary
		c_elem = glob.widgets[field]
      
		if type(c_elem) is tk.Entry:		# Clear entries
			c_elem.delete(0,'end')
		if type(c_elem) is ttk.Combobox:	# Clear Comboboxes
			c_elem.set('')
		if type(c_elem) is tk.Checkbutton:	# Disable checkbuttons
			c_elem.deselect()
		
		
	# Reset default labels manually
	glob.widgets['planets'].set('Choose a planet..')
	glob.widgets['n0'].insert(0,'14.00')
	glob.widgets['resc'].insert(0,'2.00')
	glob.widgets['heh'].insert(0,'0.083333333') 
	glob.widgets['appxmth'].set('Rate/2 + Mdot/2') 
	glob.widgets['spectrum'].set('Choose..')
	glob.widgets['spec_prop'].insert(0,'-1.0')
	glob.widgets['grid'].set('Mixed')
	glob.widgets['numflux'].set('HLLC')
	glob.widgets['reconst'].set('PLM')
	glob.widgets['rho'].config(text = 'Planet density:')
	glob.widgets['b0'].config(text = 'beta_0:')
	glob.widgets['phi'].config(text = 'log(phi_P):')
	
#-------------------------------------------   

# Load old input file
def load_input(*args):
	
	# Look which file to open
	if os.path.isfile('input.inp'):	# From old run simulation - rename
		os.rename('input.inp', 'input_temp.inp')

	# Input file
	inp_file = 'input_temp.inp'
			
	# Go through file and read keyword by keyword	
	f = open(inp_file,'r')

	# Planet name
	line = f.readline()
	pname = get_word(line,3)
	
	# Log10 of n0
	line = f.readline()
	n0 = get_word(line,7)

	# Planet radius
	line = f.readline()
	Rp = get_word(line,4)
	
	# Planet mass
	line = f.readline()
	Mp = get_word(line,4)
	
	# Equilibrium temperature
	line = f.readline()
	T0 = get_word(line,4)
	
	# Orbital distance
	line = f.readline()
	a = get_word(line,4)
	
	# Escape radius
	line = f.readline()
	resc = get_word(line,4)
	
	# He/H number ratio
	line = f.readline()
	heh = get_word(line,4)
	
	# 2D approximate method
	line = f.readline()
	appxmth = get_word(line,4) 
	
	# Read alpha if selected
	if appxmth == 'alpha': alpha = get_word(line,6)
	
	# Parent star mass
	line = f.readline()
	Ms = get_word(line,5)
	
	# Spectrum type 
	line = f.readline()
	sp_type = get_word(line,3)
	
	# Next read properties of spectrum
	if sp_type == 'Load':			# Load from file
		line = f.readline()
		sp_file = get_word(line,3)
		isthere_sp_file = os.path.isfile(sp_file)
		print('WARNING: The numerical spectrum file specified in input.inp does not exists. ' 
		  		'Please select a new one.')
		if not isthere_sp_file: sp_file = ''
	
	if sp_type == 'Power-law':		# Power law
		line = f.readline()
		PLind = get_word(line,3)
	
	if sp_type == 'Monochromatic':	# Monochromatic
		line = f.readline()
		glob.ephot = get_word(line,4)
	
	
	# Only EUV status
	line = f.readline()
	onlyEUV = get_word(line,4)

	# If not monochromatic, read energy bands
	LX = '0.0'
	if sp_type != 'Monochromatic':
	
		if onlyEUV == 'True':
			line = f.readline()
			glob.elow = get_word(line,4)
			glob.emid = get_word(line,6)
			
			# Set others to default
			glob.ehigh = '1.24e3'
			
		if onlyEUV == 'False':
			line  = f.readline()
			glob.elow  = get_word(line,4)
			glob.emid  = get_word(line,6)
			glob.ehigh = get_word(line,8)

			# Read X-ray luminosity if included
			line = f.readline()
			LX = get_word(line,6)
			
	# LEUV luminosity
	line = f.readline()
	LEUV = get_word(line,6)
	
	# Grid type
	line = f.readline()
	grid = get_word(line,3)
	
	# Numerical flux
	line = f.readline()
	numflux = get_word(line,3)
	
	# Reconstruction scheme
	line = f.readline()
	reconst = get_word(line,3)

	# Include He23S
	line = f.readline()
	IncludeHe23S = get_word(line,3)

	# IC status
	line = f.readline()
	LoadIC = get_word(line,3)
	
	# Do only post-processing
	line = f.readline()
	do_only_pp = get_word(line,4)

	# Force start of sim.
	line = f.readline()
	force_start = get_word(line,3)

	# Close file
	f.close()
	
	# Clean pre-existing values and reset
	for field in glob.widgets:
		
		# Current element of dictionary
		c_elem = glob.widgets[field]
      
		if type(c_elem) is tk.Entry:		# Clear entries
			c_elem.delete(0,'end')
		if type(c_elem) is ttk.Combobox:	# Clear Comboboxes
			c_elem.set('')
		if type(c_elem) is tk.Checkbutton:	# Disable checkbuttons
			c_elem.deselect()
	
	# Fill with parameters that have been read
	glob.widgets['planets'].set(pname)
	glob.widgets['n0'].insert(0,n0)
	glob.widgets['Rp'].insert(0,Rp) 
	glob.widgets['Mp'].insert(0,Mp)
	glob.widgets['T0'].insert(0,T0)
	glob.widgets['a'].insert(0,a)
	glob.widgets['resc'].insert(0,resc)
	glob.widgets['heh'].insert(0,heh)
	glob.widgets['Ms'].insert(0,Ms)
	glob.widgets['LX'].insert(0,LX)
	glob.widgets['LEUV'].insert(0,LEUV)
	glob.widgets['grid'].set(grid) 
	glob.widgets['numflux'].set(numflux) 
	glob.widgets['reconst'].set(reconst) 
	
	# Take care of alpha box
	if appxmth == 'Mdot/4':
		glob.widgets['appxmth'].insert(0,glob.appx_meth[0])
	if appxmth == 'Rate/2':
		glob.widgets['appxmth'].insert(0,glob.appx_meth[1])
	if appxmth == 'Rate/4':
		glob.widgets['appxmth'].insert(0,glob.appx_meth[2])
	if appxmth == 'alpha':
		alpha_str = glob.appx_meth[3] + alpha
		glob.widgets['appxmth']. insert(0,alpha_str)
	
	# Spectrum specifications
	glob.widgets['spectrum'].set(sp_type)
	glob.current_spec = sp_type
	
	# Next set the properties of spectrum
	if sp_type == 'Load':			# Load from file
		glob.widgets['spectrum'].delete(0,'end')
		glob.widgets['spectrum'].insert(0,glob.sp_type[0])
		# Change state of spectrum property entry
		glob.widgets['spec_prop'].insert(0,sp_file)
	
	if sp_type == 'Power-law':		# Power law
	
		glob.widgets['spec_prop'].delete(0,'end')
		# Change state of spectrum property entry
		glob.widgets['spec_prop'].insert(0,PLind)
		
	if sp_type == glob.sp_type[2]:	# Monochromatic
	
		glob.widgets['spec_prop'].insert(0,glob.ephot)
	
	if onlyEUV == 'True':
		glob.onlyEUV_var.set(1)

	if IncludeHe23S == 'True': 
		glob.He23S_var.set(1)

	if LoadIC == 'True': 
		glob.LoadIC_var.set(1)

	if do_only_pp == 'True':
		glob.onlyPP_var.set(1)
	
	if force_start == 'True':
		glob.force_var.set(1)

	# Update spectrum properties
	spectrum_func() # Include flux update and labels update
	# Execute callback function for only EUV button
	onlyEUV_func()
	
#-------------------------------------------

# Initialize interface at startup
def init_func(*args):
	
	# Allocate global tk variables
	allocate_tkvars()

	# If old file is found, load it and fill
	if os.path.isfile("input.inp") or os.path.isfile("input_temp.inp"):	
		
		# Change state of spectrum property entry
		glob.widgets['spec_prop'].configure(state = tk.NORMAL)
		
		load_input()
		
	else:	# Use reset function to default entries
		
		# Define default energy bands
		glob.elow  = '13.60' 
		glob.emid  = '123.98'
		glob.ehigh = '1.24e3'

		# Disable PLindex and spectrum file entries
		glob.widgets['spec_prop'].config(state = tk.DISABLED)
		
		reset_func()

#-------------------------------------------

# Fill parameters at planet selection (Plist)
def Plist_func(*args):
	
	# Clean entries to be changed
	for c_ent in glob.pl_params_list:
		if type(c_ent) is ttk.Combobox:	# Clear Comboboxes
			continue
		glob.widgets[c_ent].delete(0,'end')
	
	# ---
	
	# Get current planet name
	planet_id = glob.widgets['planets'].current()
	
	# Fill with planetary parameters
	for c_ent in glob.pl_params_list:
		c_value = glob.pl_params[c_ent][planet_id]
		glob.widgets[c_ent].insert(0,c_value)
	
	# ---

	# Update flux
	upd_flux()
	
	# Update labels
	upd_labels()

#-------------------------------------------

# Callback for type of spectrum selection
def spectrum_func(*args):

	# Activate spectrum properties field
	glob.widgets['spec_prop'].configure(state = tk.NORMAL)

	# Get selected spectrum type
	spec_type = glob.widgets['spectrum'].get()

	if get_word(spec_type,1) != glob.current_spec:
		glob.widgets['spec_prop'].delete(0,'end')
	
	# Enable onlyEUV button
	glob.widgets['EUVonly'].configure(state = tk.NORMAL)
	
	# Redefine energy bands if not monochromatic
	if spec_type != glob.sp_type[2]:
		glob.elow  = '13.60' 
		glob.emid  = '123.98'
		glob.ehigh = '1.24e3'

	# Deselect onlyEUV button if selected
	if glob.onlyEUV_var.get() == 1:
		glob.onlyEUV_var.set(0)
		onlyEUV_func()
	
	# Activate/deactivate fields according to spectrum type
	if spec_type == glob.sp_type[0]:	# Load from file
	
		# Change label
		glob.widgets['lbl_spec'].config(text = "Spectrum file:")
		
		# Clear flux, LEUV and LX entries
		glob.widgets['LX'].delete(0,'end')
		glob.widgets['LEUV'].delete(0,'end')
		glob.widgets['flux'].delete(0,'end')
		
		# Select and deselect appropriatte entries
		glob.widgets['browse'].config(state = tk.NORMAL)

	if spec_type == glob.sp_type[1]:	# Power-law
		
		# Change label
		glob.widgets['lbl_spec'].config(text = "Power-law index:")
	
		glob.widgets['browse'].config(state = tk.DISABLED)
		
	if spec_type == glob.sp_type[2]:	# Monochromatic
			
		# Change label
		glob.widgets['lbl_spec'].config(text = "Radiation energy [eV]")
		
		# Disable PLindex and spectrum file entries
		glob.widgets['browse'].config(state = tk.DISABLED)
		
		# Check and disable onlyEUV button, run its callback
		glob.onlyEUV_var.set(1)
		glob.widgets['EUVonly'].configure(state = tk.DISABLED)
		onlyEUV_func()
	
	# Update flux
	upd_flux()
	
#-------------------------------------------

# Callback function for onlyPP tickbox
def onlyPP_func(*args):
	
	# Get current status of the button
	onlyPP_state = glob.onlyPP_var.get()

	# Disable force start if onlyPP is selected 
	if onlyPP_state == 1: # If selected 
		# Disable force
		glob.widgets['force'].config(state = tk.DISABLED)
		glob.force_var.set(0)
	else: 
		# Enable force
		glob.widgets['force'].config(state = tk.NORMAL)
	return

#-------------------------------------------

# Callback function for force tickbox
def force_func(*args):
	
	# Get current status of the button
	force_state = glob.force_var.get()
	
	# Disable force start if onlyPP is selected 
	if force_state == 1: # If selected 
		# Disable onlyPP
		glob.widgets['onlyPP'].config(state = tk.DISABLED)
		glob.onlyPP_var.set(0)
	else:
		# Enable onlyPP
		glob.widgets['onlyPP'].config(state = tk.NORMAL)
	return

#-------------------------------------------

# Callback for Only-EUV checkbutton
def onlyEUV_func(*args):
	
	# Get current status of button
	EUV_state = glob.onlyEUV_var.get()
	
	# Disable LX entry if chbutton in on
	if EUV_state == 1:	# If selected
		# Disable LX
		glob.widgets['LX'].config(state = tk.DISABLED)
	else: # If deselected
		# Enable LX
		glob.widgets['LX'].config(state = tk.NORMAL)
	# Update flux
	upd_flux()
			
#-------------------------------------------

# Numerical integration of the numerical spectrum
def integrate_spectrum(sp_file):
	
	# The file has to be formatted in two columns:
	#	   1) bin central wavelength [Angstrom] 
	#	   2) flux at the planet surface [erg cm^{-2} s^{-2} A^{-1}]
	#	Data have to be order with increasing wavelength

	# Load data
	A = np.loadtxt(sp_file)
	
	# Check if the matrix has two columns
	if len(A[0,:]) != 2: 
		# Error message
		print('Check spectrum file format!') 
		
		# Pass invalid string to function
		LEUV = 'invalid'
		LX   = 0.0
		return LEUV,LX
	
	# Convert to floats
	e_low  = float(glob.elow)
	e_mid  = float(glob.emid)
	e_high = float(glob.ehigh)
	
	# Get onlyEUV variable
	onlyEUV_status = glob.onlyEUV_var.get()
	
	# Do not proceed if energies are not ordered
	if e_low > e_mid : return
	
	if onlyEUV_status == 0:	# Only check top if X is included
		if e_mid > e_high : return
	

	# Convert wavelengths to energies
	A[:,0] = glob.c_light/(A[:,0]*1e-8)*glob.hp_eV
	
	# Convert fluxes in erg/(..*eV)
	A[:,1] = A[:,1]*glob.c_light*1e8/A[:,0]**2.0*glob.hp_eV
	
	# --- Cut out the energy and flux vectors in the two bands
	
	# Initialize flux lists 
	E_X   = []
	E_EUV = []
	FX   = []
	FEUV = []
	
	
	# Cut vectors
	j = 0		# Row counter
	E = A[j,0]	# Highest energy in the spectrum
	
	# Start from X band (if selected)
	if onlyEUV_status == 0:
		while E > e_mid:
			
			# Current energy
			E = A[j,0]
			
			if E > e_high: # If not already in the range
				j += 1
				continue
			else:
				E_X.append(A[j,0])
				FX.append(A[j,1])
				# Update counter
				j += 1
				continue

	# Go into EUV range
	while E > e_low:
		
		# Current energy
		E = A[j,0]

		# If energy is below HI ionization threshold, skip
		if E < 13.60: continue
		
		E_EUV.append(E)
		FEUV.append(A[j,1])
		# Update counter
		j += 1
		continue
	
	# Convert to numpy arrays
	E_X   = np.array(E_X)
	E_EUV = np.array(E_EUV)
	FX	= np.array(FX)
	FEUV  = np.array(FEUV)
	
	# Number of selected wavelengths
	NE_X   = len(E_X)
	NE_EUV = len(E_EUV)
	
	# Get orbital distance value
	a_r  = glob.widgets['a'].get()
	a    = float(a_r)*glob.AU
	
	# If a = 0, abort
	if a_r == '': return
	
	
	# Integrate X-ray first (if active)
	if onlyEUV_status == 0:
		
		# Initialize
		dE_X = np.zeros(NE_X)
		
		# Calculate energy cells width
		dE_X[0]        = -0.5*(E_X[1] - E_X[0])
		dE_X[1:NE_X-2] = -0.5*(E_X[2:NE_X-1] - E_X[0:NE_X-3])
		dE_X[NE_X-1]   = -0.5*(E_X[NE_X-1] - E_X[NE_X-2])
	
		# Calculate luminosity
		LX = np.sum(dE_X[:]*FX[:])*4.0*np.pi*a*a
	else:
		LX = 0.0
		
	# Continue with EUV
		
	# Initialize
	dE_EUV = np.zeros(NE_EUV)
	
	# Calculate energy cells width
	dE_EUV[0]        = -0.5*(E_EUV[1] - E_EUV[0])
	dE_EUV[1:NE_EUV-2] = -0.5*(E_EUV[2:NE_EUV-1] - E_EUV[0:NE_EUV-3])
	dE_EUV[NE_EUV-1]   = -0.5*(E_EUV[NE_EUV-1] - E_EUV[NE_EUV-2])

	# Calculate luminosity
	LEUV = np.sum(dE_EUV[:]*FEUV[:])*4.0*np.pi*a*a
	
	# -----
	
	# Convert to log10
	if onlyEUV_status == 0: LX = np.log10(LX)
	LEUV = np.log10(LEUV)
	
	# Clean LEUV and LX fields
	glob.widgets['LX'].delete(0,'end')
	glob.widgets['LEUV'].delete(0,'end')
	
	# Insert into the specific fields
	glob.widgets['LX'].insert(0,'{:5.2f}'.format(LX))
	glob.widgets['LEUV'].insert(0,'{:5.2f}'.format(LEUV))

#-------------------------------------------

# Callback for browse button
def browse_func(*args):

	# Spectrum file
	sp_file = filedialog.askopenfile(title = 'Select spectrum file..')
	
	# Exit if not selected
	if sp_file == None : return
	
	# Clean the Spectrum File entry
	glob.widgets['spec_prop'].delete(0,'end')
	
	# Fill the Spectrum File entry with the selected file 
	glob.widgets['spec_prop'].insert(0,sp_file.name)
	
	# Update flux by integration 
	integrate_spectrum(sp_file)

	# Update flux
	upd_flux()

#-------------------------------------------

# Update luminosity if only EUV flux is selected
def flux_func(*args):
	
	# Get state of onlyEUV button 
	onlyEUV_state = glob.onlyEUV_var.get()
	
	# Quit if both bands are selected
	if onlyEUV_state == 0:	return
	
	# Quit if numerical spectrum is selected
	if glob.widgets['spectrum'].get() == glob.sp_type[0]: return
	
	# Get data
	a_r    = glob.widgets['a'].get()
	FXUV_r = glob.widgets['flux'].get()

	# Quit if entries are empty
	if a_r == '' or FXUV_r == '': return
	
	# Convert to floats
	a    = float(a_r)*glob.AU
	FXUV = float(FXUV_r)

	# Calculate new EUV luminosity
	LEUV = FXUV*4.0*np.pi*a*a
	
	# Convert to log
	lgLEUV = np.log10(LEUV)
	
	# Clean LEUV entry
	glob.widgets['LEUV'].delete(0,'end')
	
	# Inser new value
	glob.widgets['LEUV'].insert(0,'{:5.2f}'.format(lgLEUV))

#-------------------------------------------

# Add planet to planet_table file
def add_func(*args):

	# Retrieve data for table
	pl_name = glob.widgets['planets'].get()
	Rp_r    = glob.widgets['Rp'].get()
	Mp_r    = glob.widgets['Mp'].get()
	T0_r    = glob.widgets['T0'].get()
	a_r     = glob.widgets['a'].get()
	Ms_r    = glob.widgets['Ms'].get()
	LX_r    = glob.widgets['LX'].get()
	LEUV_r  = glob.widgets['LEUV'].get()

	# Write warning message if not all fields are full
	if pl_name == '' or \
	   	Rp_r    == '' or \
	   	Mp_r    == '' or \
	   	T0_r    == '' or \
	   	a_r     == '' or \
	   	Ms_r    == '' or \
	   	LX_r    == '' or \
	   	LEUV_r  == '': 
		print('WARNING: Fill all fields before adding to table')
		return

	# Get onlyEUV variable
	onlyEUV_status = glob.onlyEUV_var.get()

	# Set X luminosity to zero if only EUV is checked
	if onlyEUV_status == 1: LX = 0.0

	# Convert to real
	Rp    = float(Rp_r)
	Mp    = float(Mp_r)
	T0    = float(T0_r)
	a     = float(a_r)
	Ms    = float(Ms_r)
	LX    = float(LX_r)
	LEUV  = float(LEUV_r)

	
	
	# Append data to table
	with open(glob.f_table, 'a') as f:
		f.write(("%s\t%5.3f\t%6.3f\t%6.1f\t%6.4f\t%5.3f\t%6.3f\t%6.3f\n")
			%(pl_name,Rp,Mp,T0,a,Ms,LX,LEUV))
      
	print("Planet added to table.")

	return

#-------------------------------------------

# Close and do not save
def close_func():
	
	# Destroy window
	glob.widgets['window'].destroy()

#-------------------------------------------

# Start function: Write datas to file
def start_func(*args):

	# Open input file and clean contents
	open('input_temp.inp', 'w').close()

	# Open file with append method
	with open('input_temp.inp', 'a') as f:

		# --- Write data to file 
		
		# Planet name
		p_name = glob.widgets['planets'].get()
		f.write('%s' %('Planet name: ' + p_name.strip()))
		
		# Log10 of n0
		n0 = glob.widgets['n0'].get()
		f.write('%s' %('\nLog10 lower boundary number density [cm^-3]: ' + n0.strip()))

		# Planet radius
		Rp = glob.widgets['Rp'].get()
		f.write('%s' %('\nPlanet radius [R_J]: ' + Rp.strip()))

		# Planet mass
		Mp = glob.widgets['Mp'].get()
		f.write('%s' %('\nPlanet mass [M_J]: ' + Mp.strip()))

		# Equilibrium temperature
		T0 = glob.widgets['T0'].get()
		f.write('%s' %('\nEquilibrium temperature [K]: ' + T0.strip()))

		# Orbital distance
		a = glob.widgets['a'].get()
		f.write('%s' %('\nOrbital distance [AU]: ' + a.strip()))

		# Escape radius
		resc = glob.widgets['resc'].get()
		f.write('%s' %('\nEscape radius [R_p]: ' + resc.strip()))
		
		# He/H number ratio
		heh = glob.widgets['heh'].get()
		f.write('%s' %('\nHe/H number ratio: ' + heh.strip()))
		
		# 2D approximate method
		appxmth = glob.widgets['appxmth'].get()
		f.write('%s' %('\n2D approximate method: ' + appxmth.strip()))
		
		# Parent star mass
		Ms = glob.widgets['Ms'].get()
		f.write('%s' %('\nParent star mass [M_sun]: ' + Ms.strip()))
		
		# Spectrum type 
		sp_type = glob.widgets['spectrum'].get()
		f.write('%s' %('\nSpectrum type: ' + sp_type.strip()))
		
		# Next print properties of spectrum
		if sp_type == glob.sp_type[0]:	# Load from file
			sp_file = glob.widgets['spec_prop'].get()
			f.write('%s' %('\nSpectrum file: ' + sp_file.strip()))	
		
		elif sp_type == glob.sp_type[1]:	# Power law
			PLind = glob.widgets['spec_prop'].get()
			f.write('%s' %('\nPower-law index: ' + PLind.strip()))	
		
		elif sp_type == glob.sp_type[2]:	# Monochromatic
			elow = glob.widgets['spec_prop'].get()
			f.write('%s' %('\nPhoton energy [eV]: ' + elow.strip()))
		else: # No valid spectrum file selected
			print('ERROR: Please select a valid spectrum type')
			
         # Close the input file and remove it before returning
			f.close()
			os.remove('input_temp.inp')
			return


		# Only EUV status
		OnlyEUV = glob.onlyEUV_var.get()
		if OnlyEUV == 0:
			f.write('%s' %('\nUse only EUV? ' + 'False'))
		if OnlyEUV == 1:
			f.write('%s' %('\nUse only EUV? ' + 'True'))
		
		
		# If not monochromatic, print energy bands
		if sp_type != glob.sp_type[2]:
			# Get energies
			e_low  = glob.elow
			e_mid  = glob.emid
			e_high = glob.ehigh
		
			# Do not print e_high if onlyEUV is selected
			if OnlyEUV == 1:
				f.write('%s' %('\n[E_low,E_mid] = [ ' + 
					e_low.strip() + ' - ' +
					e_mid.strip() + ' ]' ))
			# EUV + X
			if OnlyEUV == 0:
				f.write('%s' %('\n[E_low,E_mid,E_high] = [ ' + 
					e_low.strip()  + ' - ' +
					e_mid.strip()  + ' - ' +
					e_high.strip() + ' ]' ))
		
		
		# Print X-ray luminosity if included
		if OnlyEUV == 0:
			LX = glob.widgets['LX'].get()
			f.write('%s' %('\nLog10 of X-ray luminosity [erg/s]: ' + LX.strip()))
		
		# LEUV luminosity
		LEUV = glob.widgets['LEUV'].get()
		f.write('%s' %('\nLog10 of EUV luminosity [erg/s]: ' + LEUV.strip()))
		
		# Grid type 
		grid = glob.widgets['grid'].get()
		f.write('%s' %('\nGrid type: ' + grid.strip()))
		
		# Numerical flux
		numflux = glob.widgets['numflux'].get()
		f.write('%s' %('\nNumerical flux: ' + numflux.strip()))
		
		# Reconstruction scheme
		reconst = glob.widgets['reconst'].get()
		f.write('%s' %('\nReconstruction scheme: ' + reconst.strip()))

		include_He23S = glob.He23S_var.get()
		if include_He23S == 0:
			f.write('%s' %('\nInclude He23S? ' + 'False'))
		if include_He23S == 1:
			f.write('%s' %('\nInclude He23S? ' + 'True'))
		
		# Load IC status
		LoadIC = glob.LoadIC_var.get()
		if LoadIC == 0:
			f.write('%s' %('\nLoad IC? ' + 'False'))
		if LoadIC == 1:
			f.write('%s' %('\nLoad IC? ' + 'True'))
			
			# Create file copies if true
			cdir = os.getcwd()
			cdir = cdir + "/output"
			old_hydro = os.path.join(cdir, 'Hydro_ioniz.txt')
			new_hydro = os.path.join(cdir, 'Hydro_ioniz_IC.txt')
			old_ioniz = os.path.join(cdir, 'Ion_species.txt')
			new_ioniz = os.path.join(cdir, 'Ion_species_IC.txt')      
			copyfile(old_hydro, new_hydro)
			copyfile(old_ioniz, new_ioniz)

		do_only_PP = glob.onlyPP_var.get()
		if do_only_PP == 0:
			f.write('%s' %('\nDo only PP: ' + 'False'))
		if do_only_PP == 1:
			f.write('%s' %('\nDo only PP: ' + 'True'))

		force_start = glob.force_var.get()
		if force_start == 0:
			f.write('%s' %('\nForce start: ' + 'False'))
		if force_start == 1:
			f.write('%s' %('\nForce start: ' + 'True'))

	# Close file
	f.close()
	
	# Rename file only if started
	os.rename('input_temp.inp', 'input.inp')

	
	# Close window
	close_func()
