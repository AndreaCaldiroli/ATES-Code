      module BC_Apply
      ! Implementation of boundary conditions
      
      use global_parameters
      use Conversion
      
      implicit none
      
      contains
      
      ! Boundary conditions for conservative variables
      subroutine Apply_BC(u_in,u_out)
      real*8, intent(in) :: u_in(1-Ng:N+Ng,3)
      real*8 :: W(1-Ng:N+Ng,3)
      real*8, intent(out) :: u_out(1-Ng:N+Ng,3)
      
      ! Convert to primitive variables
      call U_to_W(u_in,W)
      
      ! Apply bc to W's
      call Apply_BC_W(W,W)     
      
      ! Return to conservaive variables
      call W_to_U(W,u_out)     
      
      end subroutine Apply_BC
      
      !------------------------------------------!
      
      ! Boundary conditions for primitive variables
      subroutine Apply_BC_W(W_in,W_out)
      
      real*8, intent(in) :: W_in(1-Ng:N+Ng,3)
      real*8, intent(out) :: W_out(1-Ng:N+Ng,3)
      integer :: k
      
      
      ! Copy input vector
      W_out = W_in


      ! Apply BCs to primitive variables    
      do k = 1,Ng
            
            ! Lower boundary
            W_out(1-k,1) = rho1
            W_out(1-k,2) = W_out(1,2)
            if(W_out(1,2).le.(0.0)) W_out(1-k,2) = 0.0
            W_out(1-k,3) = 1.0 + 1.0e-10

            ! Upper boundary
            W_out(N+k,:) = W_out(N,:)
            
            if (rec_method.eq.'WENO3') then 
                  W_out(N+k,:) = 2.0*W_out(N+k-1,:) - W_out(N+k-2,:)
            endif
                  
      enddo
      
      end subroutine Apply_BC_W
      
      !------------------------------------------!
      
      ! Boundary conditions for reconstructed variables
      subroutine Rec_BC(WL_in,WR_in,WL_out,WR_out)
      real*8, dimension(1-Ng:N+Ng,3), intent(in) :: WL_in, WR_in
      integer :: k
      real*8, dimension(1-Ng:N+Ng,3), intent(out) :: WL_out, WR_out
      
      WL_out = WL_in
      WR_out = WR_in
      
      ! Lower boundary
      WR_out(-1,1) = (1.0 + 4.0*HeH)/(1.0 + HeH)
      WR_out(-1,2) = WR_out(0,2)
      if(WR_out(1,2).le.(0.0)) WR_out(-1,2) = 0.0
      WR_out(-1,3) = 1.0 + 1.0e-10
      
      WL_out(-1:0,1) = (1.0 + 4.0*HeH)/(1.0 + HeH)
      WL_out(-1:0,2) = WL_out(1,2)
      if(WL_out(1,2).le.(0.0)) WL_out(-1:0,2) = 0.0
      WL_out(-1:0,3) = 1.0 + 1.0e-10
      
      
      ! Upper boundary
     	do k = 1,Ng
      	WR_out(N+k,:) = WR_out(N,:)
      enddo
      
      WL_out(N+2,:) = WL_out(N+1,:)
      
      
      end subroutine Rec_BC
      
      
      
      end module BC_Apply
