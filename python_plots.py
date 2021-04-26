import numpy as np
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

# Salz data
r_S,rho_S,v_S,p_S,T_S,tmp,tmp,tmp,fhi_S,fhii_S,fhei_S,fheii_S,fheiii_S = \
      np.genfromtxt('gj3470b.dat',unpack = True)


#--------------------------------------------------

# Adimensionalization parameters
T0 = T[0]                 # Surface temperature
v0 = np.sqrt(kb*T0/mu)    # Scale velocity
N  = r.size               # Number of cells

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
mom   = 4.*np.pi*v[:]*rho[:]*r[:]**2.*v0*R0**2.            # Our momentum
mom_S = 4.*np.pi*(v_S[:]*1.e5)*rho_S[:]*r_S[:]**2.*R0**2.  # Salz's momentum

# Correct mass flux for the current 3D approximation method
if appx_mth.strip() == "Rate/2 + Mdot/2":
	mom_out = 0.5*mom[-20]
if appx_mth.strip() == "Mdot/4":
	mom_out = 0.25*mom[-20]


#--------------------------------------------------


# Evaluate cooling contributions  

#---- Bremsstrahung ----#

# Gaunt factors  
GF_HI = 0.79464 + 0.1243*np.log10(T)
GF_He = 0.79464 + 0.1243*np.log10(T/4.)

brem = 1.426e-27*np.sqrt(T[:])*(GF_HI[:]*nhii[:] +
                            4.0*GF_He[:]*(nheii[:]+nheiii[:]))*ne[:]


#---- Collisional excitation ----#

coex_rate_HI   = 7.5e-19/(1.0 + np.sqrt(T[:]/1.0e5))  \
                 *np.exp(-118348.0/T[:])
coex_rate_HeI  = 1.1e-19*T[:]**0.082*np.exp(-2.3e5/T[:])
coex_rate_HeII = 5.54e-17*T[:]**(-0.397)/(1.0 + np.sqrt(T[:]/1.0e5))  \
      	     *np.exp(-473638.0/T[:])

# Cooling rate
coex = (coex_rate_HI[:]*nhi[:]
      + coex_rate_HeI[:]*nhei[:]
      + coex_rate_HeII[:]*nheii[:])*ne[:]


#---- Recombination ----#

# HII
xl = 2.0*157807.0/T[:]
rec_cool_HII = 3.435e-30*T*xl[:]**1.970/              \
               (1.0+(xl[:]/2.250)**0.376)**3.720
               
# HeII
xl         = 2.0*285335.0/T[:]
rec_HeII_B = 1.26e-14*xl[:]**0.750
rec_cool_HeII = 1.38e-16*T*rec_HeII_B[:]

# HeIII
xl = 2.0*631515.0/T[:]
rec_cool_HeIII = 8.0*3.435e-30*T[:]*xl[:]**1.970/      \
                 (1.0+(xl[:]/2.250)**0.376)**3.720

reco  = (rec_cool_HII[:]*nhii[:]
       + rec_cool_HeII[:]*nheii[:]
       + rec_cool_HeIII[:]*nheiii[:])*ne[:]


#---- Collisional ionization ----#

# Temperature in eV
TeV = np.log(T*8.61733e-5)

# Collisional ionization coefficients (Gloover 2007)
ion_coeff_HI = np.exp(-3.271396786e1 + 1.35365560e1*TeV[:]           \
              -5.73932875e0*TeV[:]**2.0 + 1.56315498e0*TeV[:]**3.0   \
              -2.87705600e-1*TeV[:]**4.0 + 3.48255977e-2*TeV[:]**5.0 \
              -2.63197617e-3*TeV[:]**6.0 + 1.11954395e-4*TeV[:]**7.0 \
              -2.03914985e-6*TeV[:]**8.0)

ion_coeff_HeI = np.exp(-4.409864886e1 + 2.391596563e1*TeV[:]         \
              -1.07532302e1*TeV[:]**2.0 + 3.05803875e0*TeV[:]**3.0   \
              -5.6851189e-1*TeV[:]**4.0 + 6.79539123e-2*TeV[:]**5.0  \
              -5.0090561e-3*TeV[:]**6.0 + 2.06723616e-4*TeV[:]**7.0  \
              -3.64916141e-6*TeV[:]**8.0)

xl = 2.0*631515.0/T[:]
ion_coeff_HeII =  19.95*np.exp(-xl[:]/2.0)*T[:]**(-1.5)*             \
      	      xl[:]**(-1.089)/(1.0+(xl[:]/0.553)**0.735)**1.275

# Collisional ionization
coio = (2.179e-11*ion_coeff_HI[:]*nhi[:] +
        3.94e-11*ion_coeff_HeI[:]*nhei[:] +
        kb*2.0*631515.0*ion_coeff_HeII[:]*nheii[:])*ne[:]


#---- Advection and adiabatic ----#

# Evaluate adiabatic and advection cooling contributions
advec = np.zeros(N)
adiab = np.zeros(N)
advec[0:N-1] = np.diff(Eint*v*r**2.0)/np.diff(r)*1./r[0:N-1]**2.0
adiab[0:N-1] = np.diff(v*r**2.0)/np.diff(r)*p[0:N-1]/(r[0:N-1])**2.	   
advec[N-1] = advec[N-2]
adiab[N-1] = adiab[N-2]

# Evaluate heating and cooling per unit mass
brem_m = brem[:]/rho[:]
coex_m = coex/rho[:]
reco_m = reco/rho[:]
coio_m = coio/rho[:]
advec_m = advec/rho[:]*v0/R0
adiab_m = adiab/rho[:]*v0/R0
heat_m = heat/rho[:]
cool_m = cool/rho[:]




#----------------------------------------------------#

#------- Hydro plot -------#

fig, ax = plt.subplots(3,2)
fig.subplots_adjust(left   = 0.13,
                    bottom = 0.05,
                    right  = 0.89,
                    top    = 0.94,
                    hspace = 0.20,
                    wspace = 0.06)  

# Density
ax[0,0].semilogy(r,rho)
ax[0,0].semilogy(r_S,rho_S,'r')
ax[0,0].set_xlim([r[0],r[-1]])
ax[0,0].set_ylabel('rho [g cm$^{-3}$]')

# Velocity
ax[0,1].plot(r,v*v0*1e-5)
ax[0,1].plot(r_S,v_S,'r')
ax[0,1].set_xlim([r[0],r[-1]])
ax[0,1].set_ylim([-1.0,1.2*v.max()*v0*1e-5])
ax[0,1].yaxis.set_label_position("right")
ax[0,1].yaxis.tick_right()
ax[0,1].set_ylabel('v [km s$^{-1}$]')

# Pressure
ax[1,0].semilogy(r,p)
ax[1,0].semilogy(r_S,p_S,'r')
ax[1,0].set_xlim([r[0],r[-1]])
ax[1,0].set_ylabel('p [erg cm$^{-3}$]')

# Temperature
ax[1,1].plot(r,T)
ax[1,1].plot(r_S,T_S,'r')
ax[1,1].set_xlim([r[0],r[-1]])
ax[1,1].yaxis.set_label_position("right")
ax[1,1].yaxis.tick_right()
ax[1,1].set_ylabel('T [K]')

# Momentum
ax[2,0].plot(r,np.log10(mom))
ax[2,0].plot(r_S,np.log10(mom_S),'r')
ax[2,0].set_xlim([r[0],r[-1]])
ax[2,0].set_ylabel('Mdot [g s$^{-1}$]')
ax[2,0].set_xlabel('r/R$_P$')

ax[2,1].plot(r,eta)
ax[2,1].set_ylabel('Heat. eff.')
ax[2,1].yaxis.set_label_position("right")
ax[2,1].yaxis.tick_right()
ax[2,1].set_xlabel('r/R$_P$')


#------- HC plot -------#

# Colors definition
tomato = [255/255.0,99/255.0,71/255.0]
SteelBlue = [70/255.0, 130/255.0, 180/255.0]
dark_blue = [0, 0, 189/255.0]
DarkTurquoise = [0, 206/255.0, 209/255.0]

plt.figure()
plt.loglog(r,brem_m,'--',label = 'Brems.', color = dark_blue)
plt.loglog(r,coex_m,'-.',label = 'Coll. excit.', color = dark_blue)
plt.loglog(r,reco_m,linestyle = 'dotted',label = 'Reco.', color = dark_blue)
plt.loglog(r,coio_m,'--',label = 'Coll. ioniz.', color = dark_blue)
plt.loglog(r,advec_m,'--',label = 'Adv. cool.', color = SteelBlue)
plt.loglog(r,-advec_m,'--',label = 'Adv. heat.', color = tomato)
plt.loglog(r,adiab_m,'--',label = 'Adiab.',color = DarkTurquoise)
plt.loglog(r,heat_m,'r',label = 'Total Heating')
plt.loglog(r,cool_m,label = 'Total Cooling', color = dark_blue)
plt.legend(loc = 'upper right')
plt.xlim([r[0],r[-1]])
plt.ylim([1.e5, 5.0*max(heat_m)])
plt.xlabel('r/R$_P$')
plt.ylabel('Heat. and Cool. [erg cm$^{-3}$ s$^{-1}$')


#------- Ionization profiles -------#

fig, ax = plt.subplots(2,2)
fig.subplots_adjust(left   = 0.13,
                    bottom = 0.05,
                    right  = 0.89,
                    top    = 0.94,
                    hspace = 0.20,
                    wspace = 0.06)  

ax[0,0].semilogy(r,nhi,label='n_HI')
ax[0,0].semilogy(r,nhii,label='n_HII')
ax[0,0].semilogy(r,nhei,label='n_HeI')
ax[0,0].semilogy(r,nheii,label='n_HeII')
ax[0,0].semilogy(r,nheiii,label='n_HeIII')
ax[0,0].semilogy(r,nhi+nhii,'--k',label='n_H')
ax[0,0].semilogy(r,nhei+nheii+nheiii,'-.k',label='n_He')
ax[0,0].set_ylabel('n$_i$ [cm$^-3$]')
ax[0,0].set_xlim([r[0],r[-1]])
ax[0,0].set_ylim([1e0,1e14])

# HI
ax[0,1].plot(r_S,fhi_S,'r')
ax[0,1].plot(r,fhi)
ax[0,1].yaxis.set_label_position("right")
ax[0,1].yaxis.tick_right()
ax[0,1].set_ylabel('n$_{HI}$/n$_H$')
ax[0,1].set_xlim([r[0],r[-1]])
ax[0,1].set_ylim([0,1])

# HeI
ax[1,0].plot(r_S,fhei_S,'r')
ax[1,0].plot(r,fhei)
ax[1,0].set_ylabel('n$_{HeI}$/n$_{He}$')
ax[1,0].set_xlim([r[0],r[-1]])
ax[1,0].set_ylim([0,1])

# HeII
ax[1,1].plot(r_S,fheii_S,'r')
ax[1,1].plot(r,fheii)
ax[1,1].set_ylabel('n$_{HeII}$/n$_{He}$', rotation=270)
ax[1,1].yaxis.set_label_position("right")
ax[1,1].yaxis.tick_right()
ax[1,1].set_xlim([r[0],r[-1]])
ax[1,1].set_ylim([0,1])


# Print the mass loss rate
print("2D approximation method: ", appx_mth.strip())
print("Log10 of mass-loss-rate = ", np.log10(mom_out))



plt.show()


