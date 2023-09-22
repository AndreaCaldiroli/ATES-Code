      module global_parameters
      ! Definition of global parameters and vectors
      
      implicit none
      
      integer, parameter :: outfile = 99  ! Unit number of report file
      integer, parameter :: N = 500       ! Number of computational cells
      integer, parameter :: Ng = 2        ! Number of ghost cells
      integer, parameter :: Nl_fix = 200  ! Number of default energy bins
      integer :: Nl, NlTR                 ! Number of points for energy integrations
      integer :: N_eq                     ! Numbers of equations in NL solver
      integer :: lwa                      ! Working array length for NL solver
      integer :: info                     ! Output info variable of NL solver
      integer :: j_min
      integer :: count
      
      character(len = 9), parameter   :: inp_file = 'input.inp'
      character(len = :), allocatable :: p_name
      character(len = :), allocatable :: grid_type
      character(len = :), allocatable :: flux 
      character(len = :), allocatable :: rec_method 
      character(len = :), allocatable :: appx_mth
      character(len = :), allocatable :: sp_type
      character(len = :), allocatable :: sed_file
      
      logical :: is_mom_const  = .false.  ! Is momentum constant within tolerance
      logical :: is_zero_dt    = .false.  ! Is time derivative really zero 
      logical :: force_start   = .false.  ! Force to do first 1000 iterations
      logical :: do_only_pp    = .false.  ! Do only the post processing
      logical :: do_load_IC    = .false.  ! Load existing IC
      logical :: thereis_He    = .false.  ! Is He included in computations
      logical :: do_read_sed   = .false.  ! Is numerical SED read from file
      logical :: is_monochr    = .false.  ! Is monochromatic radiation selected
      logical :: is_PL_sed     = .false.  ! Is the SED a power law
      logical :: thereis_Xray  = .false.  ! Include only EUV band	
      logical :: thereis_HeITR = .false.  ! Include calculations for He triplet
      logical :: use_weno3     = .false.  ! Use WENO3 reconstruction
      logical :: use_plm       = .false.  ! Use PLM reconstruction

      !------- Global constants -------!
      
      ! Physical constants
      real*8,parameter ::  pi      = 3.1415926536     ! pi   
      real*8,parameter ::  kb_erg  = 1.38e-16         ! Boltzmann constant in CGS units
      real*8,parameter ::  kb_eV   = 8.6167e-05       ! Boltzmann constant (eV/K)
      real*8,parameter ::  mu      = 1.673e-24        ! Hydrogen mass (g)
      real*8,parameter ::  g       = 1.666666666667   ! Polytropic index
      real*8,parameter ::  Gc      = 6.67259e-8       ! Gravitational constant (CGS)   
      real*8,parameter ::  erg2eV  = 6.241509075e11   ! 1 erg measured in eV
      real*8,parameter ::  hp_erg  = 6.62620e-27      ! Planck constant in CGS units
      real*8,parameter ::  hp_eV   = 4.1357e-15       ! Planck constant (eV*s)
      real*8,parameter ::  c_light = 2.99792458e10    ! Speed of light in cm/s
      real*8,parameter ::  parsec  = 3.08567758147e18 ! 1 pc in cm
      real*8,parameter ::  AU      = 1.495978707e13   ! Astronomical unit
      real*8,parameter ::  RJ      = 6.9911e9         ! Jupiter radius (cm)
      real*8,parameter ::  MJ      = 1.898e30         ! Jupiter mass (g)
      real*8,parameter ::  Msun    = 1.989e33         ! Sun mass (g)
      real*8,parameter ::  R_earth = 6.3725e8         ! Earth radius (cm)
      real*8,parameter ::  M_earth = 5.9726e27        ! Earth mass (g)
      real*8,parameter ::  ih      = 1.0              ! Atomic number of Hydrogen
      real*8,parameter ::  ihe     = 2.0              ! Atomic number of Helium
	
	   ! -- - Energy constants

      ! Energy intervals
      real*8 ::  e_top
      real*8 ::  e_mid
      real*8 ::  e_low
      
      ! Threshold energies
      real*8,parameter ::  e_th_HI   = 13.6      ! Threshold for HI ionization
      real*8,parameter ::  e_th_HeI  = 24.6      ! Threshold for HeI ionization
      real*8,parameter ::  e_th_HeII = 54.4      ! Threshold for HeII ionization
      real*8,parameter ::  e_th_HeTR = 4.80      ! Threshold for HeI triplet ionization

	   ! Numerical constants
      real*8,parameter ::  CFL    = 0.6         ! CFL number
      real*8,parameter ::  du_th  = 1.0e-3      ! Escape momentum variation
      real*8,parameter ::  dtu_th = 1.0d-8      ! Threshold variation of time deriv.
      real*8           ::  du                   ! Initial momentum variation
      real*8           ::  dtu                  ! Norm of time derivative
      
      !------- Planetary parameters -------!

      real*8  ::  n0          ! Density at lower boundary (cm^-3)
      real*8  ::  R0          ! Planetary radius (cm)
      real*8  ::  Mp          ! Planet mass (g)
      real*8  ::  T0          ! Temperature at lower boundary (K)
      real*8  ::  LX          ! Log10 X-Ray luminosity (erg/s)
      real*8  ::  LEUV        ! Log10 EUV luminosity (erg/s)
      real*8  ::  J_XUV       ! Bolometric star XUV flux (erg/cm^2*s)
      real*8  ::  Mstar       ! Mass of companion star
      real*8  ::  Mrapp       ! Ratio M_star/M_p
      real*8  ::  atilde      ! Orbital radius in unit of R0 (a/R0)
      real*8  ::  r_esc       ! Escape radius for constant momentum
      real*8  ::  a_orb       ! Orbital distance
      real*8  ::  r_max       ! Maximum radius = Roche Lobe dimension 
      real*8  ::  HeH         ! He/H ratio	
      real*8  ::  rho_bc      ! Adimensional number density at origin
      real*8  ::  a_tau       ! Rate correction coefficient
      real*8  ::  PLind       ! Index of spectral power law
      real*8  ::  dp_bc       ! Boundary condition for pressure
	
      !------- Normalizations -------!

      real*8 :: v0      ! Velocity normalization
      real*8 :: t_s     ! Time normalization
      real*8 :: p0      ! Pressure normalization
      real*8 :: q0      ! Scale normalization 
      real*8 :: b0      ! Jeans parameter at planet surface

      !------- Global vectors -------!
      
      real*8, dimension(:), allocatable :: e_v, de_v
      real*8, dimension(:), allocatable :: s_hi,s_hei,s_heii,s_heiTR
      real*8, dimension(:), allocatable :: F_XUV
      real*8, dimension(1-Ng:N+Ng) :: r,r_edg,dr_j
      real*8, dimension(1-Ng:N+Ng) :: Gphi_c,Gphi_i
      
      ! NL solver vectors
      real*8, dimension(:), allocatable :: sys_sol, sys_x
      real*8, dimension(:), allocatable :: wa
      
      contains

      ! End of module      
      end module global_parameters
