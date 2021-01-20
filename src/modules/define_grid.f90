	module grid_construction
	! Constructs radial grid, edges and sized of grids

	use global_parameters
	use J_incident


	implicit none
	
	contains

	subroutine define_grid!(r,r_edg,dr_j)
	integer :: j,i
	integer,parameter :: N_low = 50
	integer :: N_up = N - N_low
	real*8 :: drc = 2.0e-4
	real*8 :: x0,x1
      real*8 :: f,df
      real*8 :: tol = 1.0
      real*8 :: q
  
      select case (grid_type)
      
      case ('Mixed') 
      
            !------ Salz's grid ------!
            ! Constructed with N_low uniform spaced points
            ! and N_up points in a stretched grid
            
            ! Lower ghost cells
            r(1-Ng) = 1.0 - drc
            r(2-Ng) = 1.0
            
            ! Grid centers in the uniform region
            do j = 1,N_low
                  r(j) = r(j-1) + drc
            enddo
            
            ! Solve for the region of stretched grid
            ! Solves the equation f(x) = 0 using 
            ! Newton-Raphson method with
            ! 
            !     f(x) = (1-x^N)/(1-x) - (r_up-r_low)/dr
            ! 
            ! x = stretch parameter
            ! r_up, r_low = upper and lower boundary 
            ! of the domain
            ! dr = starting grid dimension
            
            ! Initial guess
            x0 = 1.01
            
            do while(tol.ge.(1.0e-8))
                  
                  ! Evaluate function and its derivative
                  
                  ! Auxiliary constant
                  q  = (r_max-r(N_low))/drc
                  
                  ! Function
                  f  = (1.0-x0**(1.0*N_up))/(1.0-x0) - q
                  
                  ! Analytic function derivative
                  df = (f  + q - N_up*x0**(N_up-1.0))/(1.0-x0)
                   
                  ! Guess of solution
                  x1 = x0 - f/df

                  ! Evaluate tolerance
                  tol = abs(x1-x0)
                  
                  ! Update point for the next step
                  x0 = x1
                  
            ! End of while loop
            enddo
            
            ! Construct stretched grid
            do j = N_low + 1,N
                  
                  r(j) = r(j-1) + x0**(1.0*j - N_low -1)*drc
                  
            enddo
            
            ! Add ghost points at the top of the domain
            do j = 1,Ng
                  r(N+j) = 2.0*r(N+j-1) - r(N+j-2)
            enddo
            
       case ('Stretched')
        
            !------ Regular stretched grid ------!
            
            r = (/ ((r_max)**((j*1.0)/(N+Ng+0.0)), j = 1-Ng,N+Ng) /)
            
      end select
      
      
      !--- Cell edges r_{j+1/2} ---!
      
      ! Cell edges (N+2*Ng-1 points) - r_edg(j) = r_{j+1/2}
      r_edg(1-Ng:N+Ng-1) = 0.5*(r(1-Ng:N+Ng-1) + r(2-Ng:N+Ng))
      r_edg(N+Ng) = 2.0*r_edg(N+Ng-1) - r_edg(N+Ng-2)
      
      !--- Cell dimensions r_{j+1/2} - r_{j-1/2} --- !
      ! Cell size (N+2*Ng points) - dr(j) = dimension of cell j
      dr_j(2-Ng:N+Ng-1) = r_edg(3-Ng:N+Ng) - r_edg(2-Ng:N+Ng-1)
      dr_j(1-Ng) = dr_j(2-Ng)
      dr_j(N+Ng) = dr_j(N+Ng-1)
      
      
      !-----------------------------
      
      ! Select relevant domain for constant momentum
      do j = 1-Ng,N+Ng
	      if(r(j).ge.r_esc) goto 1
	enddo
		
1	j_min = j
      


	end subroutine define_grid
	
	! End of module
	end module grid_construction
