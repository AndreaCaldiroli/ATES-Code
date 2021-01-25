import numpy as np
import math as mt

# Generate table with planet parameters

# Common parameters
MJ = 1.898e30 # Jupiter mass (g)
RJ = 6.9911e9 # Jupiter radius (cm) 
AU = 1.495978707e13 # Astronomical unit (cm)
Msun = 1.989e33   # Sun mass (g)
Gc = 6.67259e-8 # Costante gravitazionale
mu = 1.673e-24 # Proton mass
kb = 1.38e-16 # Boltzmann constant


# Planet names
pl_names = ["HAT-P-11b", 
            "HD97658b_F1e6_PLM", 
            "55Cnce",   
            "GJ436b",   
            "CoRoT-7b",  
            "GJ1214b",  
            "HAT-P-2b",  
            "WASP-38b",  
            "WASP-77Ab", 
            "WASP-10b",  
            "WASP-8b",   
            "WASP-80b",  
            "WASP-43b",  
            "HD209458b",
            "HD189733b",
            "WASP-69b",
            "GJ9827b",
            "LHS1140b",
		"Murray-Clay"]
      
      
# Planet radius
Rp = [0.422,   # HAT-P-11 b
      0.209,   # HD 97658 b
      0.186,   # 55 Cnc e
      0.377,   # GJ 436 b
      0.150,   # CoRoT-7 b
      0.239,   # GJ 1214 b
      1.157,   # HAT-P-2 b
      1.090,   # WASP-38 b
      1.210,   # WASP-77A b
      1.080,   # WASP-10 b
      1.038,   # WASP-8 b
      0.952,   # WASP-80 b
      0.930,   # WASP-43 b
      1.359,   # HD 209458 b
      1.138,   # HD 189733 b
      1.057,   # WASP 69 b
      0.148,   # GJ 9827 b
      0.131,   # LHS 1140 b
      1.400]   # Murray-Clay
     
     
     
Rp = [x*RJ for x in Rp]


# Planet mass
Mp = [0.0825,   # HAT-P-11 b
      0.0247,   # HD 97658 b
      0.0262,   # 55 Cnc e
      0.0726,   # GJ 436 b
      0.0156,   # CoRoT-7 b
      0.0204,   # GJ 1214 b
      8.860,    # HAT-P-2 b
      2.710,    # WASP-38 b
      1.759,    # WASP-77A b
      3.190,    # WASP-10 b
      2.137,    # WASP-8 b
      0.552,    # WASP-80 b
      1.780,    # WASP-43 b
      0.689,    # HD 209458 b
      1.140,    # HD 189733 b
      0.260,    # WASP 69 b
      0.0108,   # GJ 9827 b
      0.0212,   # LHS 1140 b
      0.700]    # Murray-Clay
     
     
Mp = [x*MJ for x in Mp]


# Planet surface temperature
T0 = [850.0,    # HAT-P-11 b
      750.0,    # HD 97658 b
      1950.0,   # 55 Cnc e
      650.0,    # GJ 436 b
      1800.0,   # CoRoT-7 b
      550.0,    # GJ 1214 b
      1700.0,   # HAT-P-2 b
      1250.0,   # WASP-38 b
      1650.0,   # WASP-77A b
      950.0,    # WASP-10 b
      950.0,    # WASP-8 b
      800.0,    # WASP-80 b
      1350.0,   # WASP-43 b
      1320.0,   # HD 209458 b
      1200.0,   # HD 189733 b
      963.0,    # WASP 69 b
      1075.0,   # GJ 9827 b
      230.0,    # LHS 1140 b
	1000.0]   # Murray-Clay

# Log10 stellar X-ray luminosity
LX = [27.9,   # HAT-P-11 b
      26.5,   # HD 97658 b
      26.1,   # 55 Cnc e
      26.0,   # GJ 436 b
      28.5,   # CoRoT-7 b
      25.8,   # GJ 1214 b
      28.9,   # HAT-P-2 b
      28.0,   # WASP-38 b
      0.00,   # WASP-77A b
      28.1,   # WASP-10 b
      28.4,   # WASP-8 b
      27.8,   # WASP-80 b
      27.9,   # WASP-43 b
      26.5,   # HD 209458 b
      0.00,   # HD 189733 b
      28.3,   # WASP 69 b
      27.4,   # GJ 9827 b
      26.2,   # LHS 1140 b
	0.00]	  # Murray-Clay
	
# Stellar EUV-luminosity
LEUV = [4.80 + 0.860*x for x in LX]


# Planet orbital distance
a = [0.0525,   # HAT-P-11 b
     0.0796,   # HD 97658 b
     0.0154,   # 55 Cnc e
     0.0287,   # GJ 436 b
     0.0172,   # CoRoT-7 b
     0.0143,   # GJ 1214 b
     0.0679,   # HAT-P-2 b
     0.0758,   # WASP-38 b
     0.0241,   # WASP-77A b
     0.0385,   # WASP-10 b
     0.0802,   # WASP-8 b
     0.0345,   # WASP-80 b
     0.0142,   # WASP-43 b
     0.0472,   # HD 209458 b
     0.0310,   # HD 189733 b
     0.0453,   # WASP 69 b
     0.0189,   # GJ 9827 b
     0.0875,   # LHS 1140 b 
     0.0500]   # Murray-Clay

LEUV[14] = np.log10(10.0**(6.0)*4.0*np.pi*a[1]*a[1]*AU*AU)
     
# Stellar XUV bolometric flux
JXUV = [(10.0**(LX[j])+10.0**(LEUV[j]))/(4.0*mt.pi*a[j]*a[j]*AU*AU) 
        for j in range(len(Mp))]

# Stellar masses
Ms = [0.809,   # HAT-P-11 b
      0.850,   # HD 97658 b
      1.015,   # 55 Cnc e
      0.452,   # GJ 436 b
      0.930,   # CoRoT-7 b
      0.150,   # GJ 1214 b
      1.340,   # HAT-P-2 b
      1.216,   # WASP-38 b
      1.002,   # WASP-77A b
      0.710,   # WASP-10 b
      1.033,   # WASP-8 b
      0.580,   # WASP-80 b
      0.717,   # WASP-43 b
      1.148,   # HD 209458 b
      0.800,   # HD 189733 b
      0.826,   # WASP 69 b
      0.614,   # GJ 9827 b
      0.146,   # LHS 1140 b
      1.000]   # Murray-Clay
      
Ms = [x*Msun for x in Ms]     


# Ratio Mstar/Mplanet
Mtilde = [Ms[j]/Mp[j] for j in range(len(Ms))]

# Ratio a/Rp
atilde = [a[j]*AU/Rp[j] for j in range(len(a))]

# Jeans escape parameter
b0 = [(Gc*Mp[j]*mu)/(kb*T0[j]*Rp[j]) for j in range(len(Mp))]

# Density
rho = [Mp[j]/(4./3.*mt.pi*Rp[j]**3.) for j in range(len(Mp))]

# Description
description = "# The following file contains the planetary parameters for the planets listed in Ojanen W. et al., ''The effect of stellar radiation on exoplanet atmospheric heating and mass loss '', AAS, 229, 245.25 (2017). \n # The table contains the following parameters:\n # Planet's name \n # Radius [cm] \n # Mass [g] \n # Surface temperature [K] \n # X-ray luminosity [10-100 Angstrom] \n # EUV luminosity [100-912 Angstrom] obtained from the scaling equation from Poppenhaeger, K., et al. 2012, A&A, 541, A26 \n # Orbital distance [AU] \n # Bolometrc XUV luminosity \n # Mtilde = M_star/M_planet \n # atilde = a/R_planet \n # Jeans escape parameter at the planet surface: b0 = (G*M_planet*mu)/(kb*T0*R_planet), with mu = 1.673e-24 [g], G = 6.67259e-8 [cm^3 g^-1 s^-2], \n # \t kb = 1.38e-16 erg/K \n # rho g cm^-3] \n \n \n"

# Header
head = ["# Planet_name",  
        "Radius",       
        "Mass",         
        "T0",           
        "L_X",          
        "L_EUV",        
        "a",            
        "J_XUV",        
        "Mtilde",       
        "atilde",       
        "b0",           
        "rho"]
        
# Write to file
with open('planet_parameters.txt', 'w') as f:

      f.write("%s" %description)

      f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"  \
              %(head[0],   # Planet_name
                head[1],   # Radius
                head[2],   # Mass
                head[3],   # T0
                head[4],   # L_X
                head[5],   # L_EUV
                head[6],   # a
                head[7],   # J_XUV
                head[8],   # Mtilde
                head[9],   # atilde
                head[10],  # b0
                head[11])) # rho
                
      for j in range(len(a)):
            f.write("%s\t%e\t%e\t%.0f\t%.3f\t%.3f\t%.5f\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\n" \
                    %(pl_names[j],  
                      Rp[j],        
                      Mp[j],        
                      T0[j],        
                      LX[j],        
                      LEUV[j],      
                      a[j],         
                      JXUV[j],      
                      Mtilde[j],    
                      atilde[j],    
                      b0[j],        
                      rho[j]))

      















