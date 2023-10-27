import numpy as np
import matplotlib.pyplot as plt
from scipy.special import wofz
from astropy.convolution import convolve
import time

start = time.time()

# ----- CONSTANTS ----- #

# Physicals constants
kb  = 1.380649e-23		# Boltzmann constant [J/K]
G   = 6.67e-11			# Gravitational constant [m3/kg/s2]
mp  = 1.672623e-27		# Proton mass [Kg]
me  = 9.109384e-31		# Electron mass [Kg]
mD  = 3.344497e-27		# Deuterium mass [Kg]
c_light = 2.99792458e8	        # Speed of light [m/s]
AU  = 1.495978707e11		# Astronomical unit [m]
E0  = 8.854188e-12		# Vacuum permittivity [F/m]
h   = 6.626070e-34		# Planck constant [J*sec]
ht  = 1.054572e-34		# Reduced Planck constant (h slash) [J*sec]
e   = -1.602176e-19		# Electron charge [C]
mHe = 6.64648157e-27		# Helium mass [Kg]
RJ  = 6.9911e7				# Jupiter radius
MJ  = 1.898e27				# Jupiter mass
R_sun = 6.96000000e8 	                # Sun radius [m]
M_sun = 1.989e30			# Sun mass [kg]
gray  = '#a0a0a0'			# Color gray for line plot

# ----- NEEDED USER INPUTS ----- #

# Fill below here the necessary inputs to calculate the transmission spectra
#
# path       = path to the folder where ATES simulation has been performed
# Input_file = ATES' auto-generated input file
# Hydro_file = ATES' hydro output (with path)
# Ioniz_file = ATES' ionization output (with path)
# fig_name   = Name of the output figure (leave empty for not saving the figure)
# abs_file   = Name of the output file with absorption data(leave empty for not saving the file)

path = '.'  # ATES' files destination folder
Input_file = path + '/input.inp'
Hydro_file = path + '/output/Hydro_ioniz_adv.txt'
Ioniz_file = path + '/output/Ion_species_adv.txt'
fig_name_hei = ''
fig_name_lya = ''
abs_file   = ''

# Data not in input_file
R_star    = 0.44*R_sun  # Stellar radius
Instr_res_HeTR = 8e4       # Instrument resolution CARMENES: 80,000 -- GIANO-B: 50,000
Instr_res_HI = 5e4       # Instrument resolution HST-STIS 100 - 100,000
# Planet rotation period [days]
rot_period = 4.88 # [days]

# ------------------------------ #

# Number of discretization points
Grid_Number = 200

# Select range of considered wavelength and number of wavelengths
lmin_HeTR = 10828.2
lmax_HeTR = 10831.2
number_lambda_HeTR = 201

lmin_HI = 1214.1
lmax_HI = 1217.0
number_lambda_HI = 201

# ------------------------------ #

# Wavelengths in air for HeI metastable transitions (from NIST) [m]
l_He3_1 = 10830.33977e-10 
l_He3_2 = 10830.25010e-10 
l_He3_3 = 10829.09114e-10 
# Wavelengths in vacuum for HI and Deuterium Ly-alpha (from NIST) [m]
lA  = 1215.6701e-10
lD = 1215.3379e-10

# Wave frequencies for HeI metastable transitions [s^-1]
nu_He3_1 = c_light/l_He3_1
nu_He3_2 = c_light/l_He3_2
nu_He3_3 = c_light/l_He3_3
# Wave frequencies for HI and D transitions [s^-1]
nu_HI = c_light/lA
nu_D = c_light/lD

# Oscillation strenghts for He
f10830_34 = 2.9958e-1
f10830_25 = 1.7974e-1
f10829_09 = 5.9902e-2
# Oscillation strenghts for HI and Deuterium
f_la = 4.1641e-1
f_D  = 4.1630e-1

# Einstein coeff. for Helium
A12_HeTR = 1.0216e7
# Einstein coeff. for HI and D
A12_HI = 4.6986e8
A12_D = 4.6999e8

# Common constant for Faddeeva integrals
Fadd_const = np.sqrt(np.pi)*e**2.0/(4.0*np.pi*E0*me*c_light) 

# ----- USER DEFINED FUNCTIONS ----- #

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

# ------------------------- #

# Read useful parameters from the input file of ATES
with open(Input_file,'r') as f:

	data = f.readline()
	num = 1
	while data:
		data = f.readline()
		if num == 2:  Rp = float(get_word(data,4))*RJ
		if num == 3:  Mp = float(get_word(data,4))*MJ
		if num == 4:  T0 = float(get_word(data,4))
		if num == 5:  a_orb = float(get_word(data,4))*AU
		if num == 9:  Mstar = float(get_word(data,5))*M_sun
		num += 1
		
f.close()

# Load profiles
r,rho,v,p,T,heat,cool = np.loadtxt(Hydro_file, unpack = True)
r,nhi,nhii,nhei,nheii,nheiii,nheiTR = np.loadtxt(Ioniz_file, unpack = True)

# Save inverted profiles
r_I 	   = -np.flip(r)
T_I 	   =  np.flip(T)
v_I 	   = -np.flip(v)
nheiTR_I =  np.flip(nheiTR)
nH = nhi/(1 + 2.2e-5)
nH_I = np.flip(nH)
nD = nH*(2.2e-5)
nD_I = np.flip(nD)
Rib  = min(r[-1],R_star/Rp)

# Areas
A_star    = np.pi*R_star**2.0
A_planet  = np.pi*Rp**2.0
A_atm     = np.pi*(Rib*Rp)**2.0

# Instrument parameters
FWHM_HeTR 			= 1.0830e-6/Instr_res_HeTR # FWHM for convolution at 10830 Angstrom
sigma_gauss_HeTR = FWHM_HeTR/(2.0*np.sqrt(2.0*np.log(2.0)))
FWHM_HI 			= lA/Instr_res_HI # FWHM for convolution at 1215.67 Angstrom
sigma_gauss_HI = FWHM_HI/(2.0*np.sqrt(2.0*np.log(2.0)))

# ------------------------- #

# Construct grid 
r_grid = np.array([Rib**(j/(Grid_Number-1)) for j in range(Grid_Number)]) 

# Save wavelengths for figure title
lmin_lbl_HeTR = lmin_HeTR
lmax_lbl_HeTR = lmax_HeTR
lmin_lbl_HI = lmin_HI
lmax_lbl_HI = lmax_HI

# Convert to meter
lmin_HeTR = lmin_HeTR*1e-10
lmax_HeTR = lmax_HeTR*1e-10
lmin_HI = lmin_HI*1e-10
lmax_HI = lmax_HI*1e-10

# Generate vector of wavelength
l_onde_HeTR = np.linspace(lmin_HeTR, lmax_HeTR, number_lambda_HeTR)
l_plot_HeTR = l_onde_HeTR*1.0e10
l_onde_HI = np.linspace(lmin_HI, lmax_HI, number_lambda_HI)
l_plot_HI = l_onde_HI*1.0e10

# ------------------------- #

# Create double vectors 

# Radius
data_r = np.append(r_I*Rp, r*Rp)

# Temperature
data_T = np.append(T_I, T)

# Velocity
data_v = np.append(v_I, v)*1e-2

# HeI triplet density [m^3]
data_nheiTR = np.concatenate((nheiTR_I*1.0e6, nheiTR*1.0e6))

# HI and D density [m^3]
data_nHI = np.concatenate((nH_I*1.0e6, nH*1.0e6))
data_nD = np.concatenate((nD_I*1.0e6, nD*1.0e6))

# ------------------------- #

# Initialize temporary variables
exp_tau_HeTR = np.zeros((Grid_Number, number_lambda_HeTR))
prob_temp_HeTR = np.zeros((Grid_Number, number_lambda_HeTR))
prob_tot_HeTR = np.zeros(number_lambda_HeTR)
# --
exp_tau_HD = np.zeros((Grid_Number, number_lambda_HI))
prob_temp_HD = np.zeros((Grid_Number, number_lambda_HI))
prob_tot_HD = np.zeros(number_lambda_HI)


# ----- LOOP ON GRID POINTS ----- #
for p in range(Grid_Number):

	# Center distance
	r_temp = r_grid[p]*Rp

	# My version of x_value + x_coordinate
	x_LOS = np.array([])	# Initialize as empty

	# Get index of current inner points
	data_arg = np.where( (abs(data_r) >= r_temp))[0]
	r_LOS    = data_r[data_arg]
	
	for rad in r_LOS: 
		
		# Get sign of radius
		Rad_sign = np.sign(rad)
		
		# Calculate distance
		dist     = np.sqrt(rad**2.0 - r_temp**2.0)*Rad_sign
			
		# Append to vector
		x_LOS = np.append(x_LOS,dist)
	
	# Get number of points on current LOS
	len_x = len(x_LOS)
	
	# Get x-grid spacing
	data_dx = np.abs(x_LOS[1:] - x_LOS[:-1])
	
	# ------------------------ #

	# Get values of temperature and velocity corresponding to the selected points 
	T_LOS = data_T[data_arg]
	v_th_HeTR = np.sqrt(2.0*kb*T_LOS/mHe)
	v_th_HI = np.sqrt(2.0*kb*T_LOS/mp)
	v_th_D = np.sqrt(2.0*kb*T_LOS/mD)
	n_HeTR = data_nheiTR[data_arg]
	n_HI = data_nHI[data_arg]
	n_D = data_nD[data_arg]
	v_LOS = data_v[data_arg]
	v_x = x_LOS*v_LOS/r_LOS

	# ------------------------ #

	# Doppler shifts parameters
	Dnu_He3_1 = (nu_He3_1*v_th_HeTR)/c_light 
	Dnu_He3_2 = (nu_He3_2*v_th_HeTR)/c_light 
	Dnu_He3_3 = (nu_He3_3*v_th_HeTR)/c_light

	Dnu_HI = (nu_HI*v_th_HI)/c_light
	Dnu_D = (nu_D*v_th_D)/c_light 

	# Arguments of Voigt Function 
	a_He3_1 = A12_HeTR/(4.0*np.pi*Dnu_He3_1)
	a_He3_2 = A12_HeTR/(4.0*np.pi*Dnu_He3_2)
	a_He3_3 = A12_HeTR/(4.0*np.pi*Dnu_He3_3)

	a_HI = A12_HI/(4.0*np.pi*Dnu_HI)
	a_D = A12_D/(4.0*np.pi*Dnu_D)
	
	# ------------------------ #
	
	for l_idx,l in enumerate(l_onde_HeTR):

		# Arguments of Voigt Function
		X_He3_2 = (c_light/l - nu_He3_2)/Dnu_He3_2[:]
		X_He3_3 = (c_light/l - nu_He3_3)/Dnu_He3_3[:]
		X_He3_1 = (c_light/l - nu_He3_1)/Dnu_He3_1[:]

		# ----- Calculate absorption integrals via Faddeeva method ----- #
		
		# First line HeI3
		arg_Fadd_1 = X_He3_1[:] - v_x[:]/v_th_HeTR[:] + 1j*a_He3_1[:]
		Voigt_1  = f10830_34*Fadd_const/Dnu_He3_1[:]*wofz(arg_Fadd_1).real
			
		# Second line HeI3
		arg_Fadd_2 = X_He3_2[:] - v_x[:]/v_th_HeTR[:] + 1j*a_He3_2[:]
		Voigt_2  = f10830_25*Fadd_const/Dnu_He3_2[:]*wofz(arg_Fadd_2).real

		# Third line HeI3
		arg_Fadd_3 = X_He3_3[:] - v_x[:]/v_th_HeTR[:] + 1j*a_He3_3[:]
		Voigt_3  = f10829_09*Fadd_const/Dnu_He3_3[:]*wofz(arg_Fadd_3).real

		# Calculate integrands
		I_He3_1 = n_HeTR[:]*Voigt_1[:]
		I_He3_2 = n_HeTR[:]*Voigt_2[:]
		I_He3_3 = n_HeTR[:]*Voigt_3[:]	

		# Calculate optical depth via trapezoids
		tau_not_sum_HeTR = data_dx/2.0*(
		   I_He3_1[:-1] + I_He3_1[1:]  + 
		   I_He3_2[:-1] + I_He3_2[1:]  +
		   I_He3_3[:-1] + I_He3_3[1:])
		tau_v_HeTR = sum(tau_not_sum_HeTR)
		
		# Calculate transmission probability and append to matrix
		exp_tau_HeTR[p,l_idx] = np.exp(-tau_v_HeTR)
		
# ----- End loop for HeI metastable triplet ----- #


	for l_idx,l in enumerate(l_onde_HI):

		# Arguments of Voigt Function
		X_HI = (c_light/l - nu_HI)/Dnu_HI[:]
		X_D = (c_light/l - nu_D)/Dnu_D[:]

		# ----- Calculate absorption integrals via Faddeeva method ----- #
		
		# Hydrogen
		arg_Fadd_HI = X_HI[:] - v_x[:]/v_th_HI[:] + 1j*a_HI[:]
		Voigt_HI  = f_la*Fadd_const/Dnu_HI[:]*wofz(arg_Fadd_HI).real
			
		# Deuterium
		arg_Fadd_D = X_D[:] - v_x[:]/v_th_D[:] + 1j*a_D[:]
		Voigt_D  = f_D*Fadd_const/Dnu_D[:]*wofz(arg_Fadd_D).real


		# Calculate integrands
		I_HI = n_HI[:]*Voigt_HI[:]
		I_D = n_D[:]*Voigt_D[:]

		# Calculate optical depth via trapezoids
		tau_not_sum_HD = data_dx/2.0*(
		   I_HI[:-1] + I_HI[1:]  + 
		   I_D[:-1] + I_D[1:])
		tau_v_HD = sum(tau_not_sum_HD)
		
		# Calculate transmission probability and append to matrix
		exp_tau_HD[p,l_idx] = np.exp(-tau_v_HD)
		
# ----- End loop for Hydrogen and Deuterium ----- #



# Integral over the planet's projected area metastable HeI triplet
for l in range(number_lambda_HeTR):
	prob_tot_HeTR[l] = np.trapz(x = r_grid, y = 2.0*exp_tau_HeTR[:,l]*r_grid)*A_planet/(A_atm - A_planet)

# Do geometric average with star area
avg_prob_HeTR = ((A_star - A_atm) + (A_atm - A_planet)*prob_tot_HeTR[:])/A_star

# Normalization to continuum
avg_prob_HeTR = avg_prob_HeTR[:]*A_star/(A_star - A_planet)


# Integral over the planet's projected area Hydrogen and Deuterium
for l in range(number_lambda_HI):
	prob_tot_HD[l] = np.trapz(x = r_grid, y = 2.0*exp_tau_HD[:,l]*r_grid)*A_planet/(A_atm - A_planet)
	               
# Do geometric average with star area
avg_prob_HD = ((A_star - A_atm) + (A_atm - A_planet)*prob_tot_HD[:])/A_star

# Normalization to continuum
avg_prob_HD = avg_prob_HD[:]*A_star/(A_star - A_planet)

# -------------------------------------------------------------------- # 

# Create vector of wavelengths to convolve metastable HeI triplet profile
l_range_HeTR = l_onde_HeTR[-1] - l_onde_HeTR[0]
v_gauss_HeTR = np.linspace(-l_range_HeTR*0.5, l_range_HeTR*0.5, number_lambda_HeTR)

# Do convolution with gaussian nIR instrument Resolution
gaussian_HeTR = np.exp(-0.5*(v_gauss_HeTR[:]/sigma_gauss_HeTR)**2.0)   # Normalized at 1
convolved_avg_prob_HeTR = convolve(avg_prob_HeTR, gaussian_HeTR, boundary = 'extend')

# -----

# Create vector of wavelengths to convolve Hydrogen and Deuterium profile
l_range_HD = l_onde_HI[-1] - l_onde_HI[0]
v_gauss_HD = np.linspace(-l_range_HD*0.5, l_range_HD*0.5, number_lambda_HI)

# Do convolution with gaussian UV instrument Resolution
gaussian_HD = np.exp(-0.5*(v_gauss_HD[:]/sigma_gauss_HI)**2.0)   # Normalized at 1
convolved_avg_prob_HD = convolve(avg_prob_HD, gaussian_HD, boundary = 'extend')

# ------------------------- #

# Do convolution including planet rotation (for metastable HeI triplet)

transit_depth = (Rp/R_star)**2.0
hHeTR = 1.0 - convolved_avg_prob_HeTR.min()
if hHeTR > (Rib**2.0 - 1.0)*transit_depth:
    Reff_HeTR = Rp*Rib
else:
    Reff_HeTR = Rp*np.sqrt((hHeTR + transit_depth)/transit_depth)
    
v_ang     	  = 2.0*np.pi/(60.0*60.0*24.0*rot_period)
v_rot_HeTR 	  = v_ang*Reff_HeTR  # [m/sec]
dl_gauss_HeTR 	  = 10830.0*1e-10*(v_rot_HeTR/c_light)
sigma_rot_HeTR 	  = dl_gauss_HeTR/(2.0*np.sqrt(2.0*np.log(2.0)))
gaussian_vrot_HeTR = np.exp(-0.5*(v_gauss_HeTR[:]/sigma_rot_HeTR)**2.0)
convolved_rot_prob_HeTR = convolve(convolved_avg_prob_HeTR, gaussian_vrot_HeTR, boundary = 'extend')

# -----

# Do convolution including planet rotation (for Hydrogen and Deuterium)

hHD = 1.0 - convolved_avg_prob_HD.min()
if hHD > (Rib**2-1)*transit_depth:
    Reff_HD = Rp*Rib
else:
    Reff_HD = Rp*np.sqrt((hHD + transit_depth)/transit_depth)
    
v_rot_HD 	  = v_ang*Reff_HD  # [m/sec]
dl_gauss_HD 	  = 1215.0*1e-10*(v_rot_HD/c_light)
sigma_rot_HD 	  = dl_gauss_HD/(2.0*np.sqrt(2.0*np.log(2.0)))
gaussian_vrot_HD = np.exp(-0.5*(v_gauss_HD[:]/sigma_rot_HD)**2.0)
convolved_rot_prob_HD = convolve(convolved_avg_prob_HD, gaussian_vrot_HD, boundary = 'extend')

# Transmission minima

Tl_HeI3 = (1.0 - avg_prob_HeTR.min())*100.0
Tl_HeI3_conv = (1.0 - convolved_avg_prob_HeTR.min())*100.0
Tl_HeI3_conv_rot = (1.0 - convolved_rot_prob_HeTR.min())*100.0

Tl_HI = (1.0 - avg_prob_HD.min())*100.0
Tl_HI_conv = (1.0 - convolved_avg_prob_HD.min())*100.0
Tl_HI_conv_rot = (1.0 - convolved_rot_prob_HD.min())*100.0

# ----- Setup of the figure ----- #

plt.figure(figsize = (8,7))

##### Figure metastable HeI triplet #####

# Line plots
plt.plot([l_He3_1*1.0e10,l_He3_1*1.0e10], [0.0, 1.1], '--', color = gray)
plt.plot([l_He3_2*1.0e10,l_He3_2*1.0e10], [0.0, 1.1], '--', color = gray)
plt.plot([l_He3_3*1.0e10,l_He3_3*1.0e10], [0.0, 1.1], '--', color = gray)

# Plot the curves of transmission
plt.plot(l_plot_HeTR, avg_prob_HeTR, '--', label = r'Theoretical T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_HeI3, 2)))
plt.plot(l_plot_HeTR, convolved_avg_prob_HeTR, '-.', label = 'Instrument conv. T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_HeI3_conv,2)))
plt.plot(l_plot_HeTR, convolved_rot_prob_HeTR, label = 'Planet rot. + Inst. conv T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_HeI3_conv_rot,2)))

# Axis setup
plt.xlabel(r"Wavelength [$\AA{}$]", fontsize = 15)
plt.ylabel(r"T$_{\lambda}$", fontsize = 15)
plt.xlim([lmin_HeTR*1.0e10, lmax_HeTR*1.0e10])
plt.ylim([0.99*avg_prob_HeTR.min(),1.02*avg_prob_HeTR.max()])	
plt.legend(loc = 'best', labelspacing = 1)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

# Construct the title based on inputs
title = r'Avg. Transm. Prob. with $' + str(Grid_Number) + \
         '\\times' + str(Grid_Number) + '$ grid points' + \
         ' -- ' + str(number_lambda_HeTR) + ' pt in $\lambda$ -- $ ' + \
         str(lmin_lbl_HeTR) + ' < \lambda < ' + str(lmax_lbl_HeTR) + '~ \AA{}$'

plt.title(title)

# Print out max transmission probability
print('\n ----- Transmission probability at peak metastable HeI triplet ----- \n')
print(' - Theoretical: ', Tl_HeI3, '%')
print(' - Instrument convolution: ', Tl_HeI3_conv, '%')
print(' - Planet rot. + Inst. convolution: ', Tl_HeI3_conv_rot, '%')
print('\n')

# Save figure
if len(fig_name_hei) > 0 : plt.savefig(fig_name_hei)

plt.show(block=False)

##### Figure Hydrogen and Deuterium #####

plt.figure(figsize=(8,7))
# Line plots
plt.plot([lA*1.0e10,lA*1.0e10], [0.0, 1.2], '--', color = gray)
plt.plot([lD*1.0e10,lD*1.0e10], [0.0, 1.2], '--', color = gray)

# Plot the curves of transmission
plt.plot(l_plot_HI, avg_prob_HD, '--', label = 'Theoretical T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_HI,2)))
plt.plot(l_plot_HI, convolved_avg_prob_HD, '-.', label = 'Instrument conv. T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_HI_conv,2)))
plt.plot(l_plot_HI, convolved_rot_prob_HD, label = 'Planet rot. + Inst. conv. T$_{{\lambda}}$ = {} $\%$'.format(round(Tl_HI_conv_rot,2)))

dlm_ism = (1.0 - 3e4/c_light)*1.0e10
dlp_ism = (1.0 + 3e4/c_light)*1.0e10
plt.fill_between(np.arange(lA*dlm_ism,lA*dlp_ism,0.01),0, 1.2, \
                 color = gray, hatch = '//', alpha = 0.1,
                 label = r'ISM absorption $\pm 30$ km/s')

# Axis setup
plt.xlabel(r"Wavelength [$\AA{}$]", fontsize = 15)
plt.ylabel(r"T$_{\lambda}$", fontsize = 15)
plt.xlim([lmin_HI*1.0e10, lmax_HI*1.0e10])
plt.ylim([0.8*avg_prob_HD.min(),1.1*avg_prob_HD.max()])	
plt.legend(loc = 'best', labelspacing = 1)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

# Construct the title based on inputs
title_HD = r'Avg. Transm. Prob. with $' + str(Grid_Number) + \
           '\\times' + str(Grid_Number) + '$ grid points' +  \
           ' -- ' + str(number_lambda_HI) + ' pt in $\lambda$ -- $ ' + \
           str(lmin_lbl_HI) + ' < \lambda < ' + str(lmax_lbl_HI) + '~ \AA{}$'

plt.title(title_HD)

# Print out max transmission probability
print('\n ----- Transmission probability at peak HI Lya ----- \n')
print(' - Theoretical: ', Tl_HI, '%')
print(' - Instrument convolution: ', Tl_HI_conv, '%')
print(' - Planet rot. + Inst. convolution: ', Tl_HI_conv_rot, '%')
print('\n')


print("--- Execution time: %s seconds ---" % (time.time() - start))

# Save figure
if len(fig_name_lya) > 0 : plt.savefig(fig_name_lya)

plt.show()




