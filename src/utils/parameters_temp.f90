      module global_parameters
      implicit none
      
      integer, parameter :: N = 500 ! Number of computational cells
      integer, parameter :: Ng = 2  ! Number of ghost cells
      character(len=*), parameter :: flux = 'HLLC'
      character(len=*), parameter :: rec_method = 'PLM'!'WENO3'!
      
      ! Global constants
      real*8 ::  g = 1.66667 		      ! Polytropic index
      real*8 ::  eps = 1e-6               ! Constant in WENO reconstr.
      real*8 ::  pi = 3.1415926536        ! pi   
      real*8 ::  mu = 1.673d-24 		! Hydrogen mass (g)
      real*8 ::  kb_erg = 1.38d-16 		! Boltzmann constant in CGS units
      real*8 ::  Gc = 6.67259d-8 		! Gravitational constant (CGS)   
      real*8 ::  HeH = 0.0833333 		! Initial He = H (=1/12)	
	real*8 ::  erg2eV = 6.241509075d11 	! 1 erg meaured in eV
	real*8 ::  ih = 1.0 			! Atomic number of Hydrogen
      real*8 ::  ihe = 2.0 			! Atomic number of Helium
      real*8 ::  hp_erg = 6.62620d-27 	! Planck constant in CGS units
      real*8 ::  hp_eV = 4.1357d-15 	! Planck constant (eV*s)
      real*8 ::  kb_eV = 8.6167d-05 	! Boltzmann constant (eV)
      real*8 ::  e_top = 12.398d3         ! Maximum energy in spectrum (eV) 
      real*8 ::  e_mid = 123.98           ! 100 Angstrom in eV
      real*8 ::  e_low = 13.6             ! 912 Angstrom in eV
      real*8 ::  a_tau = 0.0              ! Rate correction coefficient
      
      

