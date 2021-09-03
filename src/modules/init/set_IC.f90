	module initial_conditions
	! Set general initial conditions (isothermal atmosphere)
	
	use global_parameters
	use grav_func	
		
	implicit none
	
	contains
	
	subroutine set_IC(W,T,f_sp)
	integer :: j
	real*8 :: b0_eff
	real*8, dimension(1-Ng:N+Ng,3), intent(out) :: W
	real*8, dimension(1-Ng:N+Ng), intent(out) :: T
	real*8, dimension(1-Ng:N+Ng,5), intent(out) :: f_sp

		
	!--- Set initial conditions for thermodynamic variables ---!
	
	! Density
	b0_eff = 1.0	! Change if the planet b0 is too low - only for IC
      W(:,1) = (/ (rho1*exp(b0_eff*(-Gphi_c(j) + Gphi_c(0))), 	&
      		   j = 1-Ng,N+Ng) /)
      where (W(:,1).lt.(1.0e-8)) W(:,1) = 1.0e-8
      
      ! Velocity
      W(:,2) = 0.5*(r-r(0))
      
      ! Pressure
      W(:,3) = (1.0 + dp_bc)*W(:,1)/rho1
      
      ! Temperature
      T = 1.0

      ! Ionized fractions      
      f_sp(:,1) = (1.0 - dp_bc)/(1.0 + 4.0*HeH)		!HI
      f_sp(:,2) = dp_bc/(1.0 + 4.0*HeH)	                  !HII
      f_sp(:,3) = HeH*(1.0 - dp_bc)/(1.0 + 4.0*HeH)	      !HeI
      f_sp(:,4) = 1.0d-10*HeH/(1.0 + 4.0*HeH)	            !HeII
      f_sp(:,5) = 0.0					            !HeIII
      
	
	end subroutine set_IC
	
	
	! End of module
	end module initial_conditions
