	module Cross_sections
	! Energy-dependent photoionization cross sections (1e-18 cm^2)
	
	use global_parameters
	
	implicit none
	
	contains
	
	!----- Hydrogenic atoms -----! 
	
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
      
      !----------------------------------------
      
      !----- HeI -----! 
      
      double precision function sigma_HeI(E)
	real*8,intent(in) :: E
	real*8 :: eth

      ! Ionization treshold
	eth = 24.6*0.999

      if (E.ge.eth) then      ! if E > E_th
            sigma_HeI = 0.6935/((E*1.0d-2)**1.82+(E*1.0d-2)**3.23)
      else
	      sigma_Hei = 0.0
	endif
	
      end function sigma_HeI
      
      ! End of module
      end module Cross_sections
