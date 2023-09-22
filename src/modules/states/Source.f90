   module source_func
   ! Subroutine to evaluate the source function
   !     (gravitational + geometrical)
   
   use global_parameters
   use Conversion
   use grav_func
   
   implicit none
   
   contains
   
   subroutine source(j,dr,dAp,dAm,dV,u,WR,WL,S)
   integer, intent(in) :: j
   real*8, intent(in) :: dr,dAp,dAm,dV
   real*8, intent(in) :: u(3),WR(3),WL(3) 
   real*8 :: rhoL,rhoR
   real*8 :: pC
   real*8, intent(out) :: S(3)
   
   !--- Extract physical quantities ---!

   ! Density
   rhoL = WR(1)
   rhoR = WL(1)       
   
   ! Central ressure
   pC = (g-1.0)*(u(3)-0.5*u(2)*u(2)/u(1))
   
   !--- Evaluare source term ---!
   
   S(1) = 0.0
   S(2) = - 0.5*(rhoL + rhoR)*(Gphi_i(j) - Gphi_i(j-1))/dr
   S(3) = 0.0
   
   ! Add source term explicitly if PLM is used
   if (use_plm) S(2) = S(2) + (dAp-dAm)/dV*pC

   ! End of subroutine
   end subroutine source
   
   ! End of module
   end module source_func
