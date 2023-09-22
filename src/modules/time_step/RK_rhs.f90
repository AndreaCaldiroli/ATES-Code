      module RK_integration
      ! Evaluate RK right hand side (convection + source)
      
      use global_parameters
      use Numerical_Fluxes
      use source_func
      
      implicit none
      
      contains
            
      subroutine RK_rhs(u_in,WL,WR,alpha,dF,S)
      real*8, dimension(1-Ng:N+Ng,3), intent(in) :: u_in
      real*8, dimension(1-Ng:N+Ng,3), intent(in) :: WL,WR
      real*8, intent(in) :: alpha 
      integer :: j 
      real*8 :: dr
      real*8 :: rp,rm
      real*8 :: dAp,dAm
      real*8 :: dV  
      real*8, dimension(3) ::  Fp,Fm
      real*8 :: dF3p
      real*8 :: pL,pR
      real*8, dimension(1-Ng:N+Ng,3), intent(out) :: dF
      real*8, dimension(1-Ng:N+Ng,3), intent(out) :: S
      
      do j = 2-Ng,N+Ng
      
         ! Substitutions
         dr = dr_j(j)
         rp = r_edg(j)
         rm = r_edg(j-1)
         dAp = rp*rp
         dAm = rm*rm
         dV = (dAp*rp - dAm*rm)/3.0
            
         ! Evaluate numerical fluxes
         if (j.eq.(2-Ng)) then
         
               ! Use flux from previous step
               call Num_flux(WL(j-1,:),WR(j-1,:),Fm,alpha,pL)
         else
         
               Fm = Fp
               pL = pR
         endif
         
         ! Evaluate flux at the right interface
         call Num_flux(WL(j,:),WR(j,:),Fp,alpha,pR)

         ! Evaluate source
         call source(j,dr,dAp,dAm,dV,    &
                     u_in(j,:),WR(j-1,:),WL(j,:),S(j,:))
      
         ! Evaluate flux differences
         dF(j,1) = (dAp*Fp(1) - dAm*Fm(1))/dV
         dF(j,2) = (dAp*Fp(2) - dAm*Fm(2))/dV 
         
         ! Correct for WENO3 discretization
         if (use_weno3)  dF(j,2) = dF(j,2) + (pR - pL)/dr
         
         dF3p  = dAp*Fp(1)*(Gphi_i(j) - Gphi_c(j))         &
               - dAm*Fm(1)*(Gphi_i(j-1) - Gphi_c(j))
         dF(j,3) = (dAp*Fp(3) - dAm*Fm(3) + dF3p)/dV 
      
      enddo
      
      ! End of subroutine
      end subroutine RK_rhs
      
      ! End of module
      end module RK_integration
