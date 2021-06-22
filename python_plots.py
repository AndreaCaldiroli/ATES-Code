import numpy as np
import os
import matplotlib.pyplot as plt


# Phisical constants
kb = 1.38e-16           # Boltzmann constant in CGS units
mu = 1.673e-24          # Hydrogen atom mass (g)       
gam = 5.0/3.0           # Polytropic index
RJ = 6.9911e9


# LOAD DATA  

# Planetary parameters
# Get planetary radius from the input.inp file
with open("input.inp",'r') as fi:
      data = fi.readline()
      num = 1
      while data:
            data = fi.readline()
            if(num == 5):
                  R0 = float(data)*RJ
            if(num == 15):
                  appx_mth = data
            num += 1


# Hydro profiles
r,rho,v,p,T,heat,cool,eta = \
      np.loadtxt('./output/Hydro_ioniz.txt',unpack = True)
rho = rho*mu

# Load our ionization profiles
r,nhi,nhii,nhei,nheii,nheiii = \
      np.loadtxt('./output/Ion_species.txt',unpack = True)


#--------------------------------------------------

# Adimensionalization parameters
T0 = T[0]                 # Surface temperature
v0 = np.sqrt(kb*T0/mu)    # Scale velocity
N  = r.size               # Number of cells
vlim = 1.2*v.max()*v0*1e-5
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

# Spherical momentum
mom   = 4.*np.pi*v[:]*rho[:]*r[:]**2.*v0*R0**2.            
lgmom = np.zeros(N)           
for j in range(N):
	if mom[j] > 0:
		lgmom[j] = np.log10(mom[j])
	else:
		lgmom[j] = -20.0	

# Correct mass flux for the current 3D approximation method
mom_out = lgmom[-20]
if appx_mth.strip() == "Rate/2 + Mdot/2":
	mom_out = mom_out + np.log10(0.5)
if appx_mth.strip() == "Mdot/4":
	mom_out = mom_out + np.log10(0.25)


#----------------------------------------------------#

#------- Hydro plot -------#

fig, ax = plt.subplots(2,4)
plt.suptitle('Before post-processing',fontsize = 15,fontweight = 'bold')
fig.set_size_inches(11, 7)
fig.subplots_adjust(left   = 0.05,
                    bottom = 0.08,
                    right  = 0.97,
                    top    = 0.93,
                    hspace = 0.30,
                    wspace = 0.31)  

# Density
ax[0,0].semilogy(r,rho)
ax[0,0].set_xlim([r[0],r[-1]])
ax[0,0].set_title('Density [g cm$^{-3}$]', fontdict={'weight':'bold'})
ax[0,0].set_xlabel('r/R$_P$')


# Velocity
ax[0,1].semilogy(r,v*v0*1e-5)
ax[0,1].set_xlim([r[0],r[-1]])
ax[0,1].set_ylim([1.e-3,vlim])
ax[0,1].set_title('Velocity [km s$^{-1}$]', fontdict={'weight':'bold'})
ax[0,1].set_xlabel('r/R$_P$')

# Pressure
ax[0,2].semilogy(r,p)
ax[0,2].set_xlim([r[0],r[-1]])
ax[0,2].set_title('Pressure [erg cm$^{-3}$]', fontdict={'weight':'bold'})
ax[0,2].set_xlabel('r/R$_P$')

# Temperature
ax[0,3].plot(r,T)
ax[0,3].set_xlim([r[0],r[-1]])
ax[0,3].set_title('Temperature [K]', fontdict={'weight':'bold'})
ax[0,3].set_xlabel('r/R$_P$')

# Momentum
ax[1,0].plot(r,lgmom)
ax[1,0].set_xlim([r[0],r[-1]])
ax[1,0].set_title('Log10 Momentum [g s$^{-1}$]', fontdict={'weight':'bold'})
ax[1,0].set_xlabel('r/R$_P$')
ax[1,0].set_ylim([7,13])

# Ionization densities
ax[1,1].semilogy(r,nhi,label = '$n_{HI}$')
ax[1,1].semilogy(r,nhii,label = '$n_{HII}$')
ax[1,1].semilogy(r,nhei,label = '$n_{HeI}$')
ax[1,1].semilogy(r,nheii,label = '$n_{HeII}$')
ax[1,1].semilogy(r,nheiii,label = '$n_{HeIII}$')
ax[1,1].set_title('Ion densities [cm$^{-3}$]', fontdict={'weight':'bold'})
ax[1,1].set_xlim([r[0],r[-1]])
ax[1,1].set_xlabel('r/R$_P$')
ax[1,1].legend(loc = 'upper right')


# H fractions
ax[1,2].plot(r,fhi,label = '$f_{HI}$')
ax[1,2].plot(r,fhii,label = '$f_{HII}$')
ax[1,2].set_title('H fractions', fontdict={'weight':'bold'})
ax[1,2].set_xlim([r[0],r[-1]])
ax[1,2].set_ylim([0,1])
ax[1,2].set_xlabel('r/R$_P$')
ax[1,2].legend(loc = 'best')

# He fractions
ax[1,3].plot(r,fhei,label = '$f_{HeI}$')
ax[1,3].plot(r,fheii,label = '$f_{HeII}$')
ax[1,3].plot(r,fheiii,label = '$f_{HeIII}$')
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
if os.path.isfile(adv_hydro ) and os.path.isfile(adv_ioniz):

	
	# Hydro profiles
	r,rho,v,p,T,heat,cool,eta = \
		np.loadtxt('./output/Hydro_ioniz_adv.txt',unpack = True)
	rho = rho*mu

	# Load our ionization profiles
	r,nhi,nhii,nhei,nheii,nheiii = \
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

	# H fractions
	ax[1,2].plot(r,fhi,'--',label = '$f_{HI}$',color = '#1f77b4')
	ax[1,2].plot(r,fhii,'--',label = '$f_{HII}$',color = '#ff7f0e')

	# He fractions
	ax[1,3].plot(r,fhei,'--',label = '$f_{HeI}$',color = '#1f77b4')
	ax[1,3].plot(r,fheii,'--',label = '$f_{HeII}$',color = '#ff7f0e')
	ax[1,3].plot(r,fheiii,'--',label = '$f_{HeIII}$',color = '#2ca02c')





# Print the mass loss rate
print("2D approximation method: ", appx_mth.strip())
print("Log10 of mass-loss-rate = ", mom_out)

plt.show()

