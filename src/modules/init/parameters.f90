      module global_parameters
      ! Definition of global parameters and vectors
      
      implicit none
      
      integer, parameter :: N = 500   ! Number of computational cells
      integer, parameter :: Ng = 2    ! Number of ghost cells
      integer, parameter :: Nl = 200  ! Number of points for energy integrations
      integer :: j_min
      integer :: count
      
      character(len=:), allocatable :: grid_type
      character(len=:), allocatable :: flux 
      character(len=:), allocatable :: rec_method
      character(len=1) :: ans_IC
      
      !------- Global constants -------!
      
      real*8,parameter ::  g = 1.66667 		      ! Polytropic index
      real*8,parameter ::  CFL = 0.6                  ! CFL number
      real*8,parameter ::  pi = 3.1415926536          ! pi   
      real*8,parameter ::  mu = 1.673e-24 	      ! Hydrogen mass (g)
      real*8,parameter ::  kb_erg = 1.38e-16 	      ! Boltzmann constant in CGS units
      real*8,parameter ::  Gc = 6.67259e-8 	      ! Gravitational constant (CGS)   
      real*8,parameter ::  erg2eV = 6.241509075e11    ! 1 erg meaured in eV
      real*8,parameter ::  ih = 1.0 		      ! Atomic number of Hydrogen
      real*8,parameter ::  ihe = 2.0 		      ! Atomic number of Helium
      real*8,parameter ::  hp_erg = 6.62620e-27       ! Planck constant in CGS units
      real*8,parameter ::  hp_eV = 4.1357e-15 	      ! Planck constant (eV*s)
      real*8,parameter ::  kb_eV = 8.6167e-05 	      ! Boltzmann constant (eV)
      real*8,parameter ::  e_top = 12.398e3           ! Maximum energy in spectrum (eV) 
      real*8,parameter ::  e_mid = 123.98             ! 100 Angstrom in eV
      real*8,parameter ::  e_low = 13.6               ! 912 Angstrom in eV
      real*8,parameter ::  a_tau = 0.0                ! Rate correction coefficient
      real*8,parameter ::  du_th = 1.0e-3             ! Escape momentum variation
      real*8 :: du                                    ! Initial momentum variation
      real*8,parameter :: RJ = 6.9911e9
      real*8,parameter :: MJ = 1.898e30
      real*8,parameter :: Msun = 1.989e33
      real*8,parameter :: AU = 1.495978707e13
            
      !------- Global vectors -------!
      
      real*8, dimension(1:Nl) :: e_v, s_hi, s_hei,s_heii, F_XUV
      real*8, dimension(1-Ng:N+Ng) :: r,r_edg,dr_j
      real*8, dimension(1-Ng:N+Ng) :: Gphi_c,Gphi_i
      
      
      
      !------- Planetary parameters -------!

      real*8  ::  n0 		! Density at lower boundary (cm^-3)
      real*8  ::  R0            ! Planetary radius (cm)
      real*8  ::  Mp 		! Planet mass (g)
      real*8  ::  T0 		! Temperature at lower boundary (K)
      real*8  ::  Lx            ! Log10 X-Ray luminosity (erg/s)
      real*8  ::  LEUV          ! Log10 EUV luminosity (erg/s)
      real*8  ::  J_XUV  	! Bolometric star XUV flux (erg/cm^2*s)
      real*8  ::  Mstar         ! Mass of companion star
      real*8  ::  Mrapp 	! Ratio M_star/M_p
      real*8  ::  atilde    	! Orbital radius in unit of R0 (a/R0)
      real*8  ::  r_esc         ! Escape radius for constant momentum
      real*8  ::  a_orb         ! Orbital distance
      real*8  ::  r_max         ! Maximum radius = Roche Lobe dimension 
      real*8  ::  HeH		! He/H ratio	
      real*8  ::  rho1          ! Adimensional number density at origin
      
      
      !------- Normalizations -------!

      real*8 :: v0            ! Velocity normalization
      real*8 :: t_s           ! Time normalization
      real*8 :: p0            ! Pressure normalization
      real*8 :: q0            ! Scale normalization 
      real*8 :: b0            ! Jeans parameter at planet surface
     
      contains

      ! End of module      
      end module global_parameters
