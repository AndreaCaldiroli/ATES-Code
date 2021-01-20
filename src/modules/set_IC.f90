	module initial_conditions
	! Set general initial conditions (isothermal atmosphere)
	
	use global_parameters
	use grav_func	
		
	implicit none
	
	contains
	
	subroutine set_IC(W,T,f_sp)
	!real*8, dimension(1-Ng:N+Ng), intent(in) :: r
	integer :: j
	real*8, dimension(1-Ng:N+Ng,3), intent(out) :: W
	real*8, dimension(1-Ng:N+Ng), intent(out) :: T
	real*8, dimension(1-Ng:N+Ng,5), intent(out) :: f_sp

		
	!--- Set initial conditions for thermodynamic variables ---!
	
	! Density
      W(:,1) = (/ (rho1*exp(-phi(r(j))+phi(r(0))), j = 1-Ng,N+Ng) /)
      where (W(:,1).lt.(1.0e-8)) W(:,1) = 1.0e-8
      
      ! Velocity
      W(:,2) = 0.5*(r-r(0))
      
      ! Pressure
      W(:,3) = (1.0 + 1.0d-10)*W(:,1)/rho1
      
      ! Temperature
      T = 1.0

      ! Ionized fractions      
      f_sp(:,1) = (1.0 - 1.0d-10)/(1.0 + 4.0*HeH)		!HI
      f_sp(:,2) = 1.0d-10/(1.0 + 4.0*HeH)	                  !HII
      f_sp(:,3) = HeH*(1.0 - 1.0d-10)/(1.0 + 4.0*HeH)	      !HeI
      f_sp(:,4) = 1.0d-10*HeH/(1.0 + 4.0*HeH)	            !HeII
      f_sp(:,5) = 0.0					            !HeIII
      
	
	end subroutine set_IC
	
	
	! End of module
	end module initial_conditions
