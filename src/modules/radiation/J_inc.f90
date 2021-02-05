      module J_incident
      ! Energy-dependent incident spectrum
      
      use global_parameters
      
      implicit none
      
      contains
      
      double precision function J_inc(E)
      real*8, intent(in) ::  E
      real*8 :: Lrapp,J_X,J_EUV
      
      ! Ratio of luminosities
      Lrapp = 10.0**(Lx-LEUV)
      
      ! X-ray Flux
      J_X = Lrapp*J_XUV/(1.0+Lrapp)
      
      ! EUV flux
      J_EUV = J_XUV/(1.0+Lrapp)
      
      ! Parametrization of incident spectrum
      if(E.lt.e_mid) then
            J_inc = J_EUV/log(e_mid/e_low)*1.0/E
      else
            J_inc = J_X/log(e_top/e_mid)*1.0/E
      endif
      
      end function
      
      ! End of module
      end module J_incident
      
