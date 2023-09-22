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

      ! End of function
      end function
      
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
	
      ! End of function
      end function sigma_HeI
      
      !----------------------------------------
      
      ! HeI(23S) triplet photoionization cross section
      !	Fit of data from Norcross (1971)
      ! Note: Oklopcic calculations include HeITR only in the
      !	range [4.8-13.6] eV
	
      double precision function sigma_HeI23S(E)
      real*8, intent(in) :: E
      real*8 :: logE,logsigma
      real*8 :: a1,c1,a2,c2,a3,c3
      real*8 :: x1,x2,x3,x4,x5
      real*8 :: y3,y4,m
      
      ! Intervals of power laws
      x1 = log10(hp_eV*c_light/(2593.01*1e-8))*0.9999
      x2 = log10(hp_eV*c_light/(1655.63*1e-8))
      x3 = log10(hp_eV*c_light/(357.340*1e-8))
      x4 = log10(hp_eV*c_light/(271.940*1e-8))
      x5 = log10(hp_eV*c_light/(209.490*1e-8))*1.0001
      
      ! Coefficients of fit
      a1 = -0.8134
      a2 = -1.772
      a3 = -3.039
      c1 =  1.240
      c2 =  c1 + x2*(a1-a2)	! For continuity of the broken PL
      c3 =  5.470
      y3 =  a2*x3 + c2
      y4 =  a3*x4 + c3
      m  =  (y4-y3)/(x4-x3)
      
      ! Log of current energy
      logE = log10(E)

      if (logE .lt. x1 .or. logE .gt. x5) then
         sigma_HeI23S = 0.0
         return
      endif
     
      ! Value of log10 of cross section from the fit
      if (logE .le. x2) logsigma = a1*logE + c1
      if (logE .gt. x2 .and. logE .le. x3) logsigma = a2*logE + c2
      if (logE .gt. x3 .and. logE .lt. x4) logsigma = m*(logE - x3) + y3
      if (logE .ge. x4) logsigma = a3*logE + c3
      
      ! Return value
      sigma_HeI23S = 10.0**logsigma
      
      ! End of function
      end function sigma_HeI23S
      
      !----------------------------------------
      
      ! End of module
      end module Cross_sections