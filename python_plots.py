import numpy as np
import math as mt
import matplotlib.pyplot as plt


# Phisical constants
kb = 1.38e-16           # Boltzmann constant in CGS units
mu = 1.673e-24          # Hydrogen atom mass (g)       
gam = 5.0/3.0           # Polytropic index
RJ = 6.9911e9


# LOAD DATA  

# Our hydrodynamical profiles
A = np.loadtxt('./output/Hydro_ioniz.txt')
r = A[:,0]
rho_tmp = A[:,1]
v   = A[:,2]
p   = A[:,3]
T   = A[:,4]
heat = A[:,5]
cool = A[:,6]
eta  = A[:,7]


# Planetary parameters
# Get planetary radius from the input.inp file
with open("input.inp",'r') as fi:
      data = fi.readline()
      num = 1
      while data:
            data = fi.readline()
            if(num == 5):
                  R0 = float(data)*RJ
                  break
            num += 1

T0=T[0]                 # Surface temperature
v0=mt.sqrt(kb*T0/mu)    # Scale velocity
N = r.size              # Number of grid points


# Load our ionization profiles
A=np.loadtxt('./output/Ion_species.txt')
r_ion  = A[:,0]
nhi    = A[:,1]
nhii   = A[:,2]
nhei   = A[:,3]
nheii  = A[:,4]
nheiii = A[:,5]


nh  = [nhi[j] + nhii[j] for j in range(N)] # Total hydrogen density
nhe = [nhei[j] + nheii[j] + nheiii[j] for j in range(N)] # Total helium density
ne  = [nhii[j] + nheii[j] + 2.*nheiii [j] for j in range(N)] # Free electron density

rho = [nh[j] + 4.0*nhe[j] for j in range(N)]



# Initialization of various vectors
mom = np.zeros(N)             # Radial mass loss
v_sound = np.zeros(N)         # Sound speed
rec_cool_hii = np.zeros(N)    # HII rec coefficient
rec_cool_heii = np.zeros(N)   # HeII rec coefficient
rec_cool_heiii = np.zeros(N)  # HeIII rec coefficient
brem = np.zeros(N)            # Bremsstrahlung
coex = np.zeros(N)            # Collisional excitation
reco = np.zeros(N)            # Recombination
coio = np.zeros(N)            # Collisional ionization
adiab = np.zeros(N)           # Adiabatic cooling
advec = np.zeros(N)           # Advection cooling
fr_hi    = np.zeros(N)        # HI fraction
fr_hii   = np.zeros(N)        # HII fraction
fr_hei   = np.zeros(N)        # HeI fraction
fr_heii  = np.zeros(N)        # HeII fraction
fr_heiii = np.zeros(N)        # HeIII fraction

#--------------------------------------------------

# Evaluate some physical quantities
fr_hi    = nhi/nh       # HI fraction
fr_hii   = nhii/nh      # HII fraction
fr_hei   = nhei/nhe     # HeI fraction
fr_heii  = nheii/nhe    # HeII fraction
fr_heiii = nheiii/nhe   # HeIII fraction


mom   = 4.*np.pi*v*rho*r**2.*mu*v0*R0**2.      # Our momentum
v_sound = np.sqrt(gam*T/T0)*v0      # Sound speed

#--------------------------------------------------

# Evaluate cooling contributions  

#---- Bremsstrahung ----#
# Gaunt factors   
GF_HI= 0.79464 + 0.1243*np.log10(T)
GF_He = 0.79464 + 0.1243*np.log10(T/4.)

brem = 1.426e-27*np.sqrt(T)*(GF_HI*nhii+4.0*GF_He*(nheii+nheiii))
brem = [brem[j]*ne[j] for j in range(N)]


#---- Collisional excitation ----#

coex = [7.5e-19/(1.+np.sqrt(T[j]/1.0e5))*np.exp(-118348/T[j])*nhi[j]
        + 1.1e-19*(T[j])**0.082*np.exp(-2.3e5/T[j])*nhei[j]
        +5.54e-17*T[j]**(-0.397)/(1.+np.sqrt(T[j]/1.0e5))
        *np.exp(-473638/T[j])*nheii[j] for j in range(N)]

coex = [coex[j]*ne[j] for j in range(N)]

#---- Recombination ----#

rec_cool_hii = [3.435e-30*x*(2*157807.0/x)**1.970
                  /(1.0+(2*157807.0/x/2.250)**0.376)**3.720
                  for x in T]

rec_cool_heii = [1.38e-16*x*1.26e-14*(2.0*285335.0/x)**0.750
                  for x in T]

rec_cool_heiii = [8.0*3.435e-30*x*(2.0*631515.0/x)**1.970/
                  (1.0+((2.0*631515.0/x)/2.250)**0.376)**3.720
                  for x in T]

reco = [rec_cool_hii[j]*nhii[j] + rec_cool_heii[j]*nheii[j]
            + rec_cool_heiii[j]*nheiii[j] for j in range(N)]
        
reco = [reco[j]*ne[j] for j in range(N)]

#---- Collisional ionization ----#

# Temperature in eV
TeV = np.log(T*8.61733e-5)

# Collisional ionization coefficients (Gloover 2007)
ion_coeff_HI = [np.exp(-3.271396786e1 + 1.35365560e1*th
              -5.73932875e0*th**2.0e0  + 1.56315498e0*th**3.0e0
              -2.87705600e-1*th**4.0e0 + 3.48255977e-2*th**5.0e0
              -2.63197617e-3*th**6.0e0 + 1.11954395e-4*th**7.0e0
              -2.03914985e-6*th**8.0e0) for th in TeV]

ion_coeff_HeI = [np.exp(-4.409864886e1 + 2.391596563e1*th
              -1.07532302e1*th**(2.0e0) + 3.05803875e0*th**(3.0e0)
              -5.6851189e-1*th**(4.0e0) + 6.79539123e-2*th**(5.0e0)
              -5.0090561e-3*th**(6.0e0) + 2.06723616e-4*th**(7.0e0)
              -3.64916141e-6*th**(8.0e0)) for th in TeV]

ion_coeff_HeII = [19.95*np.exp(-(2.0*631515.0/x)/2.0)*x**(-1.5)*
                  (2.0*631515.0/x)**(-1.089)/
                  (1.0+((2.0*631515.0/x)/0.553)**0.735)**1.275
                  for x in T]

# Collisional ionization
coio = [2.179e-11*ion_coeff_HI[j]*nhi[j] +
        3.94e-11*ion_coeff_HeI[j]*nhei[j] +
        kb*2.0*631515.0*ion_coeff_HeII[j]*nheii[j]
        for j in range(N)]

coio = [coio[j]*ne[j] for j in range(N)]

#---- Advection and adiabatic ----#

# Evaluate adiabatic and advection cooling contributions
for j in range(N-1):
      advec[j] = (p[j+1]*v[j+1]*(r[j+1])**2.    \
                 -p[j]*v[j]*(r[j])**2.)         \
	           /(r[j+1]-r[j])*1./(r[j])**2. 
	      
      adiab[j] = p[j]/(r[j])**2.                       \
                 *(v[j+1]*(r[j+1])**2.-v[j]*(r[j])**2) \
	           /(r[j+1]-r[j]) 
	   
adiab[N-1] = adiab[N-2]
advec[N-1] = advec[N-2]


#--------------------------------------------------------------


# Plots

# Thermodynamical quantities
fig, ax = plt.subplots(3,2)
fig.suptitle('')
# Change whitespace between figure
fig.subplots_adjust(left   = 0.12, 
                    bottom = 0.09, 
                    right  = 0.95, 
                    top    = 0.90, 
                    wspace = 0.38, 
                    hspace = 0.75)    


# Density
ax[0,0].semilogy(r,rho)
ax[0,0].set_title('Density')
ax[0,0].set_xlabel('r/r_p')
ax[0,0].set_ylabel('[cm^-3]')
ax[0,0].set_xlabel('r/r_p')
ax[0,0].set_ylabel('[cm^-3]')
ax[0,0].set_xlim([r[0],r[-1]])

# Velocity
ax[0,1].plot(r,v*v0*1e-5)
ax[0,1].set_title('Velocity')
ax[0,1].set_xlabel('r/r_p')
ax[0,1].set_ylabel('[km/s]')
ax[0,1].set_xlim([r[0],r[-1]])
ax[0,1].set_ylim([-1.0,1.2*v.max()*v0*1e-5])

# Pressure
ax[1,0].semilogy(r,p)
ax[1,0].set_title('Pressure')
ax[1,0].set_xlim([r[0],r[-1]])

# Temperature
ax[1,1].plot(r,T)
ax[1,1].set_title('Temperature')
ax[1,1].set_xlabel('r/r_p')
ax[1,1].set_ylabel('[K]')
ax[1,1].set_xlim([r[0],r[-1]])

# Momentum
ax[2,0].semilogy(r,mom)
ax[2,0].set_title('Momentum')
ax[2,0].set_xlim([r[0],r[-1]])
ax[2,0].set_xlabel('r/r_p')
ax[2,0].set_ylabel('Mdot [g/s]')

ax[2,1].plot(r,eta)
ax[2,1].set_title('Heating efficiency')
ax[2,1].set_xlabel('r/r_p')
ax[2,1].set_ylabel('eta')

# Ionization profiles and ionization fractions

plt.figure()
plt.semilogy(r,brem,'--',label='Bremsstrahlung')
plt.semilogy(r,coex,'--',label='Coll. excit.')
plt.semilogy(r,reco,'--',label='Recombination')
plt.semilogy(r,coio,'--',label='Coll. ioniz.')
plt.semilogy(r,np.abs(advec)*v0/R0,'--',label='Advection')
plt.semilogy(r,adiab*v0/R0,'--',label='Adiabatic')
plt.semilogy(r,heat,'r',label='Total Heating')
plt.semilogy(r,cool,'b',label='Total Cooling')
plt.legend(loc='upper right')
plt.xlim([r[0],r[-1]])
plt.ylim([1.e-12, 1.0e-4])



ion_fig,ax = plt.subplots(1,2)
ion_fig.subplots_adjust(left   = 0.12, 
                        bottom = 0.09, 
                        right  = 0.95, 
                        top    = 0.90, 
                        wspace = 0.38, 
                        hspace = 0.75)   
                         
# Ionization profiles
ax[0].semilogy(r,nhi,label='n_HI')
ax[0].semilogy(r,nhii,label='n_HII')
ax[0].semilogy(r,nhei,label='n_HeI')
ax[0].semilogy(r,nheii,label='n_HeII')
ax[0].semilogy(r,nheiii,label='n_HeIII')
ax[0].semilogy(r,nhi+nhii,'--k',label='n_H')
ax[0].semilogy(r,nhei+nheii+nheiii,'-.k',label='n_He')
ax[0].legend(loc='upper right')
ax[0].set_xlabel('r/r_p')
ax[0].set_ylabel('n_i [cm^-3]')
ax[0].set_xlim([r[0],r[-1]])
ax[0].set_ylim([1e0,1e15])

# Ionization fractions
ax[1].plot(r,fr_hi, label = 'f_HI')
ax[1].plot(r,fr_hii, label = 'f_HII')
ax[1].plot(r,fr_hei, label = 'f_HeI')
ax[1].plot(r,fr_heii, label = 'f_HeII')
ax[1].plot(r,fr_heiii, label = 'f_HeIII')
ax[1].legend(loc='upper right')
ax[1].set_xlabel('r/r_p')
ax[1].set_ylabel('f_i')
ax[1].set_xlim([r[0],r[-1]])
ax[1].set_ylim([0,1])

# Print the mass loss rate
print("Log10 of mass-loss-rate = ", np.log10(mom[-20]/4.0))

plt.show()
