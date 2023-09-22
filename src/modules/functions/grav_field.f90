      module grav_func
      ! Expression for the gravitational field and for its derivative

      use global_parameters
      
      implicit none
      
      contains
      
      !-------------------------------------------------------!
      
      ! Gravitational field strenght (function)
      double precision function phi(r_in)
      real*8, intent(in)  :: r_in

      phi = -b0/r_in                                &
      	-b0*Mrapp/(atilde-r_in)                    &
      	-b0*(1.0+Mrapp)/(2.0*atilde**3.0)*         &
      	 (atilde*Mrapp/(1.0+Mrapp)-r_in)**2.0	
      
      ! End of function
      end function phi
      
      !-------------------------------------------------------!
      
      ! Gravitational field gradient
      double precision function Dphi(r_in)
      real*8, intent(in)  :: r_in
      
      DPhi = b0/r_in**2.0                       &
      	 -b0*Mrapp/(atilde-r_in)**2.0          &
      	 +b0*(1.0 + Mrapp)/atilde**3.0*        &
      	 (atilde*Mrapp/(1.0+Mrapp)-r_in)
      
      ! End of function
      end function Dphi
      
      ! End of module
      end module grav_func
