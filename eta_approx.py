import numpy as np
import matplotlib.pyplot as plt


#--------------------------------------------------------------


def etaeff(Frho,Kphi):
	
		# Requires in input F,Kphi NOT in log10
		
		# Def of sigmoid function
		def sigmoid(phi,beta):
			
			phi0 = 10.0**13.22
			y = 1.0/(1.0 + (phi/phi0)**beta)
			
			return y	
		
		
		# Parameters
		Fmin = 10.0**2.0
		
		# Evaluate rescaled flux-to-density ratio
		Ftilde2 = Frho/Fmin
		
		# Coefficient and exponents
		A     = 1.682*abs(np.log10(Ftilde2))**0.2802 - 5.488
		
		
		alpha = 0.02489*(Ftilde2)**(-0.0860)\
			- 0.01007*(Ftilde2)**(-0.9543)
		
		eta0  = -0.03973*abs(np.log10(Ftilde2))**2.173 - 0.01359

		
		beta  = - 0.01799*(Ftilde2)**0.1723\
			  - 3.38751*(Ftilde2)**0.0140
	
		sigma = sigmoid(Kphi,beta) 
		 
		# Approximate value of log10(eta_eff)
		eta = A*Kphi**alpha*sigma + eta0*(1.0-sigma)

		return eta
	
#---------------------------------------------------------------

# Physical constants

MJ   = 1.898e30           # Jupiter mass (g)
RJ   = 6.9911e9           # Jupiter radius (cm) 
AU   = 1.495978707e13     # Astronomical unit (cm)
Msun = 1.989e33           # Sun mass (g)
Gc   = 6.67259e-8         # Gravitational constant
mu   = 1.673e-24          # Proton mass
kb   = 1.38e-16           # Boltzmann constant

# Jovian constants

rhoJ = MJ/(4.0/3.0*np.pi*RJ**3.0)
xiJ  = (MJ/(3.0*Msun))**(1.0/3.0)*AU/RJ
phiJ = Gc*MJ/RJ

#-------------------

# Read input data from terminal
Rp   = float(input('Planetary radius [R_J]: '))
Mp   = float(input('Planetary mass [M_J]: '))
a    = float(input('Orbital distance [AU]: '))
Ms   = float(input('Parent star mass [M_sun]: '))
LX   = float(input('Log10 of X-ray luminosity [erg/s]: '))
LEUV = float(input('Log10 of EUV luminosity [erg/s]: '))	

#-------------------

# Derived parameters
rhop = Mp/Rp**3.0*rhoJ
phi  = Mp/Rp*phiJ
xi   = (Mp/Ms)**(1.0/3.0)*a/Rp*xiJ
K    = 1.0 - 3.0/(2.0*xi) + 1.0/(2.0*xi**3.0)
FXUV = (10.0**LX + 10.0**LEUV)/(4.0*np.pi*a*a*AU*AU)
Kphi = K*phi
Frho = FXUV/rhop 

#-------------------

# Warning conditions if input values are outside of validity range:

# F/rho range
if (Frho < 10.0**2.0 or Frho > 10.0**6.0):

	print('\n !!! WARNING: F/rho value (%6.4e) outside of validity range [1e2,1e6].' 
		%(Frho),
		'\n eta_eff and Mdot values may not be accurate')

# K*phi range
if (Kphi < 10.0**12.17 or Kphi > 10.0**13.29):
	print('\n !!! WARNING: K*phi value (%5.2f) outside of validity range [12.17, 13.29].' 
		%(np.log10(Kphi)), 
		'\n eta_eff and Mdot values may not be accurate')

#-------------------
	
# Evaluate approximate eta_eff and Mdot

eta_eff  = 10.0**etaeff(Frho,Kphi)
Mdot_eff = np.log10(3.0*eta_eff*FXUV/(4.0*Gc*K*rhop))

# Print to terminal
print(' ')
print('OUTPUT:')
print(' ')
print('Value of eta_eff: %6.4e' %(eta_eff))
print('Value of log10(Mdot): %5.2f' %(Mdot_eff))










	
