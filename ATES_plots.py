import numpy as np
import os
import matplotlib.pyplot as plt
import time
import sys

# Parse options for animated plots
if len(sys.argv) == 1: animate = 'False'
if len(sys.argv) == 2: animate = 'True'; sec = 4
if len(sys.argv) == 3: animate = 'True'; sec = float(sys.argv[2])

# ------------------------------------------------------------------

# Useful functions

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
				
# --------------------------------------------------------------------
# Phisical constants
kb = 1.38e-16           # Boltzmann constant in CGS units
mu = 1.673e-24          # Hydrogen atom mass (g)       
gam = 5.0/3.0           # Polytropic index
RJ = 6.9911e9


# LOAD DATA  

# Planetary parameters
# Get planetary radius from the input.inp file
with open("input.inp",'r') as f:

	data = f.readline()
	num = 1
	while data:
		data = f.readline()
		if num == 2:
			R0 = float(get_word(data,4))*RJ
		if num == 8:
			appx_mth = get_word(data,4)
			if appx_mth == 'Rate/2': appx_mth = 'Rate/2 + Mdot/2'
			if appx_mth == 'Rate/4': appx_mth = 'Rate/4 + Mdot'
		num += 1
		
f.close()

# Hydro profiles
r,rho,v,p,T,heat,cool = \
      np.loadtxt('./output/Hydro_ioniz.txt',unpack = True)
rho = rho*mu

# Load our ionization profiles
r,nhi,nhii,nhei,nheii,nheiii,nheiTR = \
      np.loadtxt('./output/Ion_species.txt',unpack = True)


#--------------------------------------------------

# Adimensionalization parameters
N  = r.size               # Number of cells
vlim = 1.2*v.max()*1e-5
if v.max() < 0 :
	vlim = 15.0

# Internal energy
Eint = p[:]/(gam -1.0)

# Ion densities and fractions
nh  = nhi[:] + nhii[:]                    # Total hydrogen density
nhe = nhei[:] + nheii[:] + nheiii[:]      # Total helium density
ne  = nhii[:] + nheii[:] + 2.*nheiii[:]   # Free electron density
fhi    = nhi[:]/nh[:]       # HI fraction
fhii   = nhii[:]/nh[:]      # HII fraction
fhei   = nhei[:]/nhe[:]     # HeI fraction
fheii  = nheii[:]/nhe[:]    # HeII fraction
fheiii = nheiii[:]/nhe[:]   # HeIII fraction
fheiTR = nheiTR[:]/nhe[:]   # HeITR fraction

# Spherical momentum
mom   = 4.*np.pi*v[:]*rho[:]*r[:]**2.*R0**2.            
lgmom = np.zeros(N)           
for j in range(N):
	if mom[j] > 0:
		lgmom[j] = np.log10(mom[j])
	else:
		lgmom[j] = -20.0	

# Minimum for momentum plot
mom_inf = lgmom.min()
if mom_inf == -20: mom_inf = 0.8*lgmom[200]
	
# Correct mass flux for the current 3D approximation method
mom_out = lgmom[-20]
if appx_mth.strip() == "Rate/2 + Mdot/2":
	mom_out = mom_out + np.log10(0.5)
if appx_mth.strip() == "Mdot/4":
	mom_out = mom_out + np.log10(0.25)


#----------------------------------------------------#

#------- Hydro plot -------#

fig, ax = plt.subplots(2,4)
fig.set_size_inches(11, 7)
fig.subplots_adjust(left   = 0.05,
                    bottom = 0.08,
                    right  = 0.97,
                    top    = 0.93,
                    hspace = 0.30,
                    wspace = 0.31)  

# Density
rho_line, = ax[0,0].semilogy(r,rho)
ax[0,0].set_xlim([r[0],r[-1]])
ax[0,0].set_title('Density [g cm$^{-3}$]', fontdict={'weight':'bold'})
ax[0,0].set_xlabel('r/R$_P$')


# Velocity
v_line, = ax[0,1].semilogy(r,v*1e-5)
ax[0,1].set_xlim([r[0],r[-1]])
ax[0,1].set_ylim([1.e-3,vlim])
ax[0,1].set_title('Velocity [km s$^{-1}$]', fontdict={'weight':'bold'})
ax[0,1].set_xlabel('r/R$_P$')

# Pressure
p_line, = ax[0,2].semilogy(r,p)
ax[0,2].set_xlim([r[0],r[-1]])
ax[0,2].set_title('Pressure [erg cm$^{-3}$]', fontdict={'weight':'bold'})
ax[0,2].set_xlabel('r/R$_P$')

# Temperature
T_line, = ax[0,3].plot(r,T)
ax[0,3].set_xlim([r[0],r[-1]])
ax[0,3].set_title('Temperature [K]', fontdict={'weight':'bold'})
ax[0,3].set_xlabel('r/R$_P$')

# Momentum
mom_line, = ax[1,0].plot(r,lgmom)
ax[1,0].set_xlim([r[0],r[-1]])
ax[1,0].set_title('Log10 Momentum [g s$^{-1}$]', fontdict={'weight':'bold'})
ax[1,0].set_xlabel('r/R$_P$')
ax[1,0].set_ylim([mom_inf, 1.2*lgmom.max()])

# Ionization densities
nhi_line,    = ax[1,1].semilogy(r,nhi,label = '$n_{HI}$')
nhii_line,   = ax[1,1].semilogy(r,nhii,label = '$n_{HII}$')
nhei_line,   = ax[1,1].semilogy(r,nhei,label = '$n_{HeI}$')
nheii_line,  = ax[1,1].semilogy(r,nheii,label = '$n_{HeII}$')
nheiii_line, = ax[1,1].semilogy(r,nheiii,label = '$n_{HeIII}$')
nheiTR_line, = ax[1,1].semilogy(r,nheiTR,label = '$n_{HeI3}$', color = '#5baca7')
ax[1,1].set_title('Ion densities [cm$^{-3}$]', fontdict={'weight':'bold'})
ax[1,1].set_xlim([r[0],r[-1]])
ax[1,1].set_xlabel('r/R$_P$')
ax[1,1].legend(loc = 'upper right')


# H fractions
fhi_line,  = ax[1,2].plot(r,fhi,label = '$f_{HI}$')
fhii_line, = ax[1,2].plot(r,fhii,label = '$f_{HII}$')
ax[1,2].set_title('H fractions', fontdict={'weight':'bold'})
ax[1,2].set_xlim([r[0],r[-1]])
ax[1,2].set_ylim([0,1])
ax[1,2].set_xlabel('r/R$_P$')
ax[1,2].legend(loc = 'best')

# He fractions
fhei_line,   = ax[1,3].plot(r,fhei,label = '$f_{HeI}$')
fheii_line,  = ax[1,3].plot(r,fheii,label = '$f_{HeII}$')
fheiii_line, = ax[1,3].plot(r,fheiii,label = '$f_{HeIII}$')
fheiTR_line, = ax[1,3].plot(r,fheiTR,label = '$f_{HeI3}$', color = '#5baca7')
ax[1,3].set_title('He fractions', fontdict={'weight':'bold'})
ax[1,3].set_xlim([r[0],r[-1]])
ax[1,3].set_ylim([0,1])
ax[1,3].set_xlabel('r/R$_P$')
ax[1,3].legend(loc = 'best')


#------------------------------------------------------------

# Plot post processed data if the files are there

# Get file paths
cdir      = os.getcwd()
adv_hydro = cdir + "/output/Hydro_ioniz_adv.txt"
adv_ioniz = cdir + "/output/Ion_species_adv.txt"

# Make second figure if _adv files are there
if os.path.isfile(adv_hydro) and os.path.isfile(adv_ioniz):

	
	# Hydro profiles
	r,rho,v,p,T,heat,cool = \
		np.loadtxt('./output/Hydro_ioniz_adv.txt',unpack = True)
	rho = rho*mu

	# Load our ionization profiles
	r,nhi,nhii,nhei,nheii,nheiii,nheiTR = \
		np.loadtxt('./output/Ion_species_adv.txt',unpack = True)


	#--------------------------------------------------


	# Ion densities and fractions
	nh  = nhi[:] + nhii[:]                    # Total hydrogen density
	nhe = nhei[:] + nheii[:] + nheiii[:]      # Total helium density
	fhi    = nhi[:]/nh[:]       # HI fraction
	fhii   = nhii[:]/nh[:]      # HII fraction
	fhei   = nhei[:]/nhe[:]     # HeI fraction
	fheii  = nheii[:]/nhe[:]    # HeII fraction
	fheiii = nheiii[:]/nhe[:]   # HeIII fraction
	fheiTR = nheiTR[:]/nhe[:]   # HeITR fraction
	
	#----------------------------------------------------#
	
	# Add Post-processed profiles
	
	# Pressure
	ax[0,2].semilogy(r,p,'--')

	# Temperature
	ax[0,3].plot(r,T,'--')
	
	# Ionization densities
	ax[1,1].semilogy(r,nhi,'--',label = '$n_{HI}$',color = '#1f77b4')
	ax[1,1].semilogy(r,nhii,'--',label = '$n_{HII}$',color = '#ff7f0e')
	ax[1,1].semilogy(r,nhei,'--',label = '$n_{HeI}$',color = '#2ca02c')
	ax[1,1].semilogy(r,nheii,'--',label = '$n_{HeII}$',color = '#d62728')
	ax[1,1].semilogy(r,nheiii,'--',label = '$n_{HeIII}$',color = '#9467bd')
	ax[1,1].semilogy(r,nheiTR,'--',label = '$n_{HeI3}$', color = '#5baca7')
	
	# H fractions
	ax[1,2].plot(r,fhi,'--',label = '$f_{HI}$',color = '#1f77b4')
	ax[1,2].plot(r,fhii,'--',label = '$f_{HII}$',color = '#ff7f0e')

	# He fractions
	ax[1,3].plot(r,fhei,'--',label = '$f_{HeI}$',color = '#1f77b4')
	ax[1,3].plot(r,fheii,'--',label = '$f_{HeII}$',color = '#ff7f0e')
	ax[1,3].plot(r,fheiii,'--',label = '$f_{HeIII}$',color = '#2ca02c')
	ax[1,3].plot(r,fheiTR,'--',label = '$f_{HeI3}$', color = '#5baca7')

# Print the mass loss rate
print("2D approximation method: ", appx_mth.strip())
print("Log10 of mass-loss-rate = ", mom_out)

if animate == 'False':
	plt.show()
else:
	plt.show(block = False)

#----------------------------------------------------#


# Update plot every 5 seconds if --live option is set
if animate == 'True':
	k = 1
	while k > 0:

		# Hydro profiles
		r,rho,v,p,T,heat,cool = \
			np.loadtxt('./output/Hydro_ioniz.txt',unpack = True)
		rho = rho*mu

		# Load our ionization profiles
		r,nhi,nhii,nhei,nheii,nheiii,nheiTR = \
			np.loadtxt('./output/Ion_species.txt',unpack = True)

		# Spherical momentum
		mom   = 4.*np.pi*v[:]*rho[:]*r[:]**2.*R0**2.            
		lgmom = np.zeros(N)           
		for j in range(N):
			if mom[j] > 0:
				lgmom[j] = np.log10(mom[j])
			else:
				lgmom[j] = -20.0	
				
		# Ion densities and fractions
		nh  = nhi[:] + nhii[:]                    # Total hydrogen density
		nhe = nhei[:] + nheii[:] + nheiii[:]      # Total helium density
		ne  = nhii[:] + nheii[:] + 2.*nheiii[:]   # Free electron density
		fhi    = nhi[:]/nh[:]       			# HI fraction
		fhii   = nhii[:]/nh[:]      			# HII fraction
		fhei   = nhei[:]/nhe[:]     			# HeI fraction
		fheii  = nheii[:]/nhe[:]    			# HeII fraction
		fheiii = nheiii[:]/nhe[:]   			# HeIII fraction
		fheiTR = nheiTR[:]/nhe[:]   			# HeIII fraction
		
		# Update data for plot
		rho_line.set_ydata(rho) 
		v_line.set_ydata(v*1e-5)
		p_line.set_ydata(p)
		T_line.set_ydata(T)
		mom_line.set_ydata(lgmom)
		nhi_line.set_ydata(nhi)
		nhii_line.set_ydata(nhii)
		nhei_line.set_ydata(nhei)
		nheii_line.set_ydata(nheii)
		nheiii_line.set_ydata(nheiii)
		nheiTR_line.set_ydata(nheiii)
		fhi_line.set_ydata(fhi)
		fhii_line.set_ydata(fhii)
		fhei_line.set_ydata(fhei)
		fheii_line.set_ydata(fheii)
		fheiii_line.set_ydata(fheiii)
		fheiTR_line.set_ydata(fheiTR)
		
		# Update axis limits
		
		vlim = 1.2*v.max()*1e-5
		if v.max() < 0 :
			vlim = 15.0
		
		# Minimum for momentum plot
		mom_inf = lgmom.min()
		if mom_inf == -20: mom_inf = 0.8*lgmom[200]
		ax[0,1].set_ylim([1.e-3,vlim])
		ax[0,3].set_ylim(0.8*T.min(),1.2*T.max())
		ax[1,0].set_ylim([mom_inf, 1.2*lgmom.max()])
		
		fig.canvas.draw()
		fig.canvas.flush_events()
		
		# Wait 5 seconds before next update
		time.sleep(sec)
		k = k + 1
	
# ------------------------------------------------------------------- !


