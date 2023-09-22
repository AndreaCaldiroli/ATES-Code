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
   Lrapp = 10.0**(LX-LEUV)
   
   ! X-ray Flux
   J_X = Lrapp*J_XUV/(1.0+Lrapp)
   
   ! EUV flux
   J_EUV = J_XUV/(1.0+Lrapp)
   
   ! Normalizations for different power-law index
   if (PLind.eq.(-1.0)) then
      
      JEUVnorm = J_EUV/log(e_mid/e_low)
      if (thereis_Xray) then
         JXnorm   = J_X/log(e_top/e_mid)
      else      		
         JXnorm   = 0.0
      endif
               
   else
   
      JEUVnorm = J_EUV*P1/(e_mid**P1 - e_low**P1)
      if (thereis_Xray) then
         JXnorm   = J_X*P1/(e_top**P1 - e_mid**P1)
      else      	
         JXnorm   = 0.0
      endif
   
   endif      
   
   ! Parametrization of incident spectrum
   if(E.lt.e_mid) then
      J_inc = JEUVnorm*E**PLind
   else
      J_inc = JXnorm*E**PLind
   endif
   
   end function
   
   ! End of module
   end module J_incident
   
