      module J_incident
      ! Energy-dependent incident spectrum
      
      use global_parameters
      
      implicit none
      
      contains
      
      double precision function J_inc(E)
      real*8, intent(in) ::  E
      real*8 :: Lrapp,J_X,J_EUV
      real*8 :: JEUVnorm,JXnorm
      real*8 :: P1
      
      ! Substitution
      P1 = PLind + 1.0
      
      ! Ratio of luminosities
      Lrapp = 10.0**(Lx-LEUV)
      
      ! X-ray Flux
      J_X = Lrapp*J_XUV/(1.0+Lrapp)
      
      ! EUV flux
      J_EUV = J_XUV/(1.0+Lrapp)
      
      ! Normalizations for different power-law index
      if (PLind.eq.(-1.0)) then
      	
      	JEUVnorm = J_EUV/log(e_mid/e_low)

      	JXnorm = J_X/log(e_top/e_mid)
      	
      	! If energy range is only in EUV
      	if (e_top.le.e_mid_XUV) then 
      		
	      	JEUVnorm = J_EUV/log(e_top/e_low)

      		JXnorm = 0.0
		endif
		
		! If energy range is only in X
      	if (e_low.ge.e_mid_XUV) then 
      		
	      	JEUVnorm = 0.0

      		JXnorm = J_X/log(e_top/e_low)
		endif
      		
      	      	
	else
      
		JEUVnorm = J_EUV*P1/(e_mid**P1 - e_low**P1)

      	JXnorm = J_X*P1/(e_top**P1 - e_mid**P1)
      	
      	
      	! If energy range is only in EUV
      	if (e_top.le.e_mid_XUV) then 
      		
	      	JEUVnorm = J_EUV*P1/(e_top**P1 - e_low**P1)

      		JXnorm = 0.0
		endif
		
		! If energy range is only in X
      	if (e_low.ge.e_mid_XUV) then 
      		
	      	JEUVnorm = 0.0

      		JXnorm = J_X*P1/(e_top**P1 - e_low**P1)
		endif
		
		print*, log10(JEUVnorm*2e-4)
		
	endif      
      
      
      ! Parametrization of incident spectrum
      if(E.lt.e_mid_XUV) then
            J_inc = JEUVnorm*E**PLind
      else
            J_inc = JXnorm*E**PLind
      endif
      
      end function
      
      ! End of module
      end module J_incident
      
