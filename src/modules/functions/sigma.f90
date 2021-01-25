	module H_cross_section
	! Photoionization cross section for hydrogenic atoms 
	!     (1E-18 cm^{-2} units)
	
	use global_parameters
	
	implicit none
	
	contains
	
	! Energy-dependent cross section function
      double precision function sigma(E,Z)
      real*8, intent(in) :: Z,E
      real*8 :: E_0,eps,cut

      ! Ionization treshold
      E_0 = 13.6*Z*Z

      if (E.gt.E_0) then      ! For E > E_th
            
            ! Substitution
            eps = sqrt(E/E_0-1.0)               
            
            ! Cross section value
            sigma = 6.3/(Z*Z)*(E_0/E)**4.0        &
                   *exp(4.0-4.0*atan(eps)/eps)    &
                   /(1.0-exp(-2.0*pi/eps))
            
      else
            ! Cross section value
            sigma = 6.3/(Z*Z) 
      
      endif

      ! Correct if E = E_th
      cut = 0.99999*E_0

      if(E.lt.cut) sigma = 0.0

      end
      
      ! End of module
      end module H_cross_section
