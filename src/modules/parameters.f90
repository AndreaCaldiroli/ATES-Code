      module global_parameters
      implicit none
      
      integer, parameter :: N = 500 ! Number of computational cells
      integer, parameter :: Ng = 2  ! Number of ghost cells
      integer, parameter :: Nl = 200  ! Number of points for energy integrations
      integer :: j_min
      integer :: count
      
      character(len=*), parameter :: grid_type = 'Stretched'!'Mixed'!
      character(len=*), parameter :: flux = 'HLLC'
      character(len=*), parameter :: rec_method = 'PLM'!'WENO3'!
      
      !------- Global constants -------!
      
      real*8,parameter ::  g = 1.66667 		      ! Polytropic index
      real*8,parameter ::  CFL = 0.4                  ! CFL number
      real*8,parameter ::  pi = 3.1415926536          ! pi   
      real*8,parameter ::  mu = 1.673e-24 		! Hydrogen mass (g)
      real*8,parameter ::  kb_erg = 1.38e-16 		! Boltzmann constant in CGS units
      real*8,parameter ::  Gc = 6.67259e-8 		! Gravitational constant (CGS)   
      real*8,parameter ::  HeH = 0.0833333 		! Initial He = H (=1/12)	
	real*8,parameter ::  erg2eV = 6.241509075e11 	! 1 erg meaured in eV
	real*8,parameter ::  ih = 1.0 			! Atomic number of Hydrogen
      real*8,parameter ::  ihe = 2.0 			! Atomic number of Helium
      real*8,parameter ::  hp_erg = 6.62620e-27 	! Planck constant in CGS units
      real*8,parameter ::  hp_eV = 4.1357e-15 	      ! Planck constant (eV*s)
      real*8,parameter ::  kb_eV = 8.6167e-05 	      ! Boltzmann constant (eV)
      real*8,parameter ::  e_top = 12.398e3           ! Maximum energy in spectrum (eV) 
      real*8,parameter ::  e_mid = 123.98             ! 100 Angstrom in eV
      real*8,parameter ::  e_low = 13.6               ! 912 Angstrom in eV
      real*8,parameter ::  a_tau = 0.0                ! Rate correction coefficient
      real*8,parameter ::  du_th = 1.0e-3             ! Escape momentum variation
      real*8 :: du                                    ! Initial momentum variation
      real*8,parameter ::  rho1 = (1.0 + 4.0*HeH)/(1.0 + HeH)     ! Adimensional number density at origin
     			 
            
      !------- Global vectors -------!
      
      real*8, dimension(1:Nl) :: e_v, s_hi, s_hei,s_heii, F_XUV
      real*8, dimension(1-Ng:N+Ng) :: r,r_edg,dr_j
      
      
      !------- Planetary parameters -------!

      real*8,parameter  ::  n0 = 1.0e14			! Density at lower boundary (cm^-3)
      real*8,parameter  ::  R0 = 2.950244e09          ! Planetary radius (cm)
      real*8,parameter  ::  Mp = 1.565850e29 		! Planet mass (g)
      real*8,parameter  ::  T0 = 850 			! Temperature at lower boundary (K)
      real*8,parameter  ::  Lx = 27.900        	      ! Log10 X-Ray luminosity (erg/s)
      real*8,parameter  ::  LEUV = 28.794             ! Log10 EUV luminosity (erg/s)
      real*8,parameter  ::  J_XUV = 9053.003  	      ! Bolometric star XUV flux (erg/cm^2*s)
      real*8,parameter  ::  Mrapp = 10276.214 		! Ratio M_star/M_p
      real*8,parameter  ::  atilde = 266.211   	      ! Orbital radius in unit of R0 (a/R0)
      real*8,parameter  :: r_max = (3.0*Mrapp)**(-1.0/3.0)*atilde  ! Maximum radius = Roche Lobe dimension 
      real*8,parameter  :: r_esc = 2.0
      
      
      
      !------- Normalizations -------!

      real*8,parameter :: rho0 = n0
	real*8,parameter :: v0   = sqrt(kb_erg*T0/mu)
	real*8,parameter :: t_s  = R0/v0
	real*8,parameter :: p0   = n0*mu*v0*v0
	real*8,parameter :: q0   = n0*mu*v0*v0*v0/R0	
      real*8,parameter :: b0   = (Gc*Mp*mu)/(kb_erg*T0*R0)   ! Jeans parameter at the surface
      
      
      
      contains
      
      end module global_parameters
