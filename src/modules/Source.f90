      module source_func
      ! Subroutine to evaluate the source function
      !     (gravitational + geometrical)
      
      use global_parameters
      use Conversion
      use grav_func
      
      implicit none
      
      contains
      
      subroutine source(rc,rp,rm,dr,dAp,dAm,dV,u,WR,WL,S)
      real*8, intent(in) :: rc,rp,rm,dr,dAp,dAm,dV
      real*8, intent(in) :: u(3),WR(3),WL(3) 
      real*8 :: rhoL,rhoR
      real*8 :: pC
      real*8, intent(out) :: S(3)
      
      ! Extract physical quantities

      ! Density
      rhoL = WR(1)
      rhoR = WL(1)       
      
      ! Central ressure
      pC = (g-1.0)*(u(3)-0.5*u(2)*u(2)/u(1))
      
      ! Evaluare source term
      S(1) = 0.0
      S(2) = - 0.5*(rhoL + rhoR)*(phi(rp) - phi(rm))/dr
      S(3) = 0.0
      
      ! Add source term explicitly if PLM is used
      if (rec_method.eq.'PLM') then 
            S(2) = S(2) + (dAp-dAm)/dV*pC
      endif

      
      end subroutine source
      
      
      end module source_func
