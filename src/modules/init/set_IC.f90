	module initial_conditions
	! Set general initial conditions (isothermal atmosphere)
	
	use global_parameters
	use grav_func	
		
	implicit none
	
	contains
	
	subroutine set_IC(W,T,f_sp)
	! Subroutine to set IC for the model

	integer :: j,i_rhalf
	real*8 :: r_half, minrho
	real*8 :: b0_eff
	real*8, dimension(1-Ng:N+Ng,3), intent(out) :: W
	real*8, dimension(1-Ng:N+Ng),   intent(out) :: T
	real*8, dimension(1-Ng:N+Ng,6), intent(out) :: f_sp

		
	!--- Set initial conditions for thermodynamic variables ---!
	
	! Density
	b0_eff = 1.0d0	! Change if the planet b0 is too low - only for IC
	do	
		W(:,1) = (/ (rho_bc*exp(b0_eff*(-Gphi_c(j) + Gphi_c(0))), 	&
				j = 1-Ng,N+Ng) /)

		! Calculate minimum of density profile
		r_half = 0.5e0*(r_max + 1.0e0) 
		i_rhalf = minloc(abs(r-r_half), dim = 1)
		minrho = W(i_rhalf,1)

		if (minrho .gt. 1.0e-8) then
			b0_eff = b0_eff + 0.2
		else
			exit
		endif

	enddo

	! Write b0_eff to output
	write(*,'(A31,F5.1,A7)') '    (set_IC.f90) Using b0_eff =',b0_eff,' for IC'
	write(*,*) 

	! Fix density in outer layers
   where (W(:,1).lt.(1.0e-8)) W(:,1) = 1.0e-8
      
   ! Velocity
   W(:,2) = 0.5*(r-r(0))
      
   ! Pressure
   W(:,3) = (1.0 + dp_bc)*W(:,1)/rho_bc
      
   ! Temperature
   T = 1.0

   ! Ionized fractions      
   f_sp(:,1) = (1.0 - dp_bc)/(1.0 + 4.0*HeH)		   ! HI
   f_sp(:,2) = dp_bc/(1.0 + 4.0*HeH)	            ! HII
   f_sp(:,3) = HeH*(1.0 - dp_bc)/(1.0 + 4.0*HeH)	! HeI
   f_sp(:,4) = 1.0d-10*HeH/(1.0 + 4.0*HeH)	      ! HeII
   f_sp(:,5) = 0.0					                  ! HeIII
   f_sp(:,6) = 0.0									      ! HeITR
	
	! End of subroutine
	end subroutine set_IC
	
	! End of module
	end module initial_conditions
