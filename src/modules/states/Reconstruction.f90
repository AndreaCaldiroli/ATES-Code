      module Reconstruction_step
      ! Collection of reconstruction procedure subroutine
      
      use global_parameters
      use BC_Apply
      use PLM_reconstruction
      use WENO3_reconstruction

      implicit none
      
      contains
      
      subroutine Reconstruct(dt,u_in,WL_out,WR_out) 
      
      integer :: j,k
      real*8, dimension(1-Ng:N+Ng,3), intent(in) :: u_in
      real*8, intent(in) :: dt
      real*8, dimension(1-Ng:N+Ng,3) :: uL,uR,WL,WR
      real*8 :: dx
      real*8, dimension(1-Ng:N+Ng,3), intent(out) :: WL_out, WR_out

      select case (rec_method)
      
            !------------------------------------------------
            
            case ('PLM') ! Piecewise linear reconstruction

                  call PLM_rec(u_in,WL,WR) 
		       
		      ! Apply BC to reconstructed variables
                  call Rec_BC(WL,WR,WL_out,WR_out)
		      
		      !call W_to_U(WL,uL_out)
		      !call W_to_U(WR,uR_out)


		!------------------------------------------------
            
            
            case ('WENO3')    ! ESWENO3 Reconstruction
            
                  call WENO3_rec(u_in,WR,WL)
                  
                  ! Apply BC to reconstructed variables
                  call Rec_BC(WL,WR,WL_out,WR_out)
		      
                  ! Convert to conservative		
		      !call W_to_U(WL,uL_out)
		      !call W_to_U(WR,uR_out)

      
      end select
  

      end subroutine Reconstruct
      
      ! End of module
      end module Reconstruction_step
