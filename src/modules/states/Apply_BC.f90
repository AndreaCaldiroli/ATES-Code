   module BC_Apply
   ! Implementation of boundary conditions
   
   use global_parameters
   use Conversion
   
   implicit none
   
   contains
   
   subroutine Apply_BC(u_in,u_out)
   ! Boundary conditions for conservative variables

   real*8, intent(in) :: u_in(1-Ng:N+Ng,3)
   real*8 :: W(1-Ng:N+Ng,3)
   real*8, intent(out) :: u_out(1-Ng:N+Ng,3)
   
   ! Convert to primitive variables
   call U_to_W(u_in,W)
   
   ! Apply bc to W's
   call Apply_BC_W(W,W)     
   
   ! Return to conservaive variables
   call W_to_U(W,u_out)     
   
   ! End of subroutine
   end subroutine Apply_BC
   
   !------------------------------------------!
   
   subroutine BC_component_constrho(W_in,index)
   ! Component-wise BC in case of zero-velocity gradient
   !     at the lower boundary

   real*8, intent(inout)  :: W_in(1-Ng:N+Ng,3)
   integer, intent(in)    :: index

   ! Save into output vector
   W_in(index,1) = rho_bc
   W_in(index,2) = max(W_in(1,2),0.0)
   W_in(index,3) = 1.0 + dp_bc 

   ! End of subroutine
   end subroutine BC_component_constrho

   !------------------------------------------!

   subroutine Apply_BC_W(W_in,W_out)
   ! Boundary conditions for primitive variables
   
   real*8, intent(in)  :: W_in(1-Ng:N+Ng,3)
   integer :: k
   real*8, intent(out) :: W_out(1-Ng:N+Ng,3)
   
   ! Copy input vector
   W_out = W_in
   
   ! BC with constant rho at lower boundary
   do k = 1,Ng
      call BC_component_constrho(W_out,1-k)
   enddo
         
   do k = 1,Ng
      ! Upper boundary
      W_out(N+k,:) = W_out(N,:)
      if (use_weno3) W_out(N+k,:) = 2.0*W_out(N+k-1,:) - W_out(N+k-2,:)
   enddo
   
   ! End of subroutine
   end subroutine Apply_BC_W
   
   !------------------------------------------!
   
   subroutine Rec_BC(WL_in,WR_in,WL_out,WR_out)
   ! Boundary conditions for reconstructed variables

   real*8, dimension(1-Ng:N+Ng,3), intent(in) :: WL_in, WR_in
   integer :: k
   real*8, dimension(1-Ng:N+Ng,3), intent(out) :: WL_out, WR_out
   
   WL_out = WL_in
   WR_out = WR_in
   
   ! Lower boundary
   call BC_component_constrho(WR_out,1-Ng)
   
   ! BC with constant density at lower boundary
   do k = 1,Ng
         call BC_component_constrho(WL_out,1-k)
   enddo
         
   ! Upper boundary
   do k = 1,Ng
      WR_out(N+k,:) = WR_out(N,:)
      if (use_weno3) WR_out(N+k,:) = 2.0*WR_out(N+k-1,:) - WR_out(N+k-2,:)
   enddo
   
   WL_out(N+2,:) = WL_out(N+1,:)
   if (use_weno3) WL_out(N+2,:) = 2.0*WL_out(N+1,:) - WL_out(N,:)

   ! End of subroutine
   end subroutine Rec_BC
   
   ! End of module
   end module BC_Apply
