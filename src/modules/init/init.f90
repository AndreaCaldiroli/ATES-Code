      module Initialization
      ! Initialize the simulation:
      ! - Set global options (thread number, time loop counter, initial
      !                       initial momentum variation)
      ! - Construct spatial grid
      ! - Construct energy grids for ionization equilibrium
      ! - Construct initial conditions (load or set)
      
      use global_parameters
      use grid_construction
      use energy_vectors_construct
      use gravity_grid_construction
      use initial_conditions
      use Conversion
      use BC_Apply
      use omp_lib
      use IC_load
      use initial_conditions
      
      implicit none 
      
      contains
      
      subroutine init(W,u,f_sp)
      ! Initialize the simulation setup

      integer :: n_omp_threads
      real*8, dimension(1-Ng:N+Ng)   :: rho,v,p,T
      real*8, dimension(1-Ng:N+Ng,6),intent(out) :: f_sp
      real*8, dimension(1-Ng:N+Ng,3),intent(out) :: W,u    
      
      write(*,*) '(init.f90) Initializing the simulation..'

      !---- Global options ----!
      
      ! Set the number of threads used 
      n_omp_threads = omp_get_max_threads()
      if (n_omp_threads .gt. 4) n_omp_threads = 4  ! Limits to a maximum of 4 threads
      call omp_set_num_threads(n_omp_threads)
      write(*,'(A12,I2,A12)') '    - Using',n_omp_threads,' OMP threads'
      
      ! Loop parameters
      count = 0  
      du  = 1.0
	   dtu = 1.0
      
      !------------------------------------------------!
      
      ! Construction of radial grid
      write(*,*) '    - Constructing the spatial grid..'
      call define_grid         
      
      !------------------------------------------------!
      
      ! Construction of energy grid
      write(*,*) '    - Precalculating energy grid, cross section and flux vector..'
      call set_energy_vectors      
      
      !------------------------------------------------!
      
      ! Pre-evaluate gravity at grid center and edges
      write(*,*) '    - Precalculating the gravitational potential..'
      call set_gravity_grid
      
      !------------------------------------------------!
      
      !---- Initial conditions ----!
      
      if (.not. do_load_IC) then

         ! Set IC to isothermal atmosphere
         write(*,*) '    - Setting the default, isothermal IC..'
      	call set_IC(W,T,f_sp)
      
      else  ! Load existing initial conditions
            
         write(*,*) '    - Loading IC from file..'
	      ! Load thermodynamic profiles
	      call load_IC(rho,v,p,T,f_sp,W)      
	                  
	      ! Change starting loop counting index
	      count = 1
    	endif
      
      !------------------------------------------------!

      ! Apply BC to initial condition
      call W_to_U(W,u)
      call Apply_BC(u,u)    

      write(*,*) '(init.f90) Done.'
      
      ! End of subroutine
      end subroutine init
      
      ! End of module
      end module Initialization
