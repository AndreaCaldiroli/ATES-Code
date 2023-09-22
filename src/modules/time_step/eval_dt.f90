   module eval_time_step
    
   use global_parameters

   implicit none

   contains
        
   subroutine eval_dt(W,dt)    
   ! Subroutine to evaluate the time step according to the CFL condition
   real*8, dimension(1-Ng:N+Ng,3), intent(in) :: W
   real*8, dimension(1-Ng:N+Ng) :: rho,v,p,cs
   real*8 :: alpha
   real*8, intent(out) :: dt

   ! Extract physical variables
   rho = W(:,1)
   v   = W(:,2)
   p   = W(:,3)
    
   ! Evaluate sound speed
   cs = sqrt(g*p/rho)

   ! Maximum eigenvalue
   alpha = maxval(abs(v) + cs)

   ! Evaluate time step according to CFL condition
   dt = CFL*minval(dr_j/(abs(v) + cs))  

   ! End of subroutine 
   end subroutine eval_dt

   ! End of module
   end module eval_time_step