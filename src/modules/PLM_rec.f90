      module PLM_reconstruction
      ! Piecewise linear reconstruction
      
      use global_parameters
      use Conversion
      
      implicit none            
      
      contains
      
      subroutine PLM_rec(u_in,WL_rec,WR_rec) 
      integer :: i,j,k
      real*8, dimension(1-Ng:N+Ng,3),intent(in) :: u_in
      !real*8, dimension(1-Ng:N+Ng), intent(in) :: r
      real*8 :: x(-1:1)
      real*8 :: W(-1:1,3)
      real*8, dimension(3) :: sp,sm,sd,sc
      real*8, dimension(1-Ng:N+Ng,3), intent(out)  :: WL_rec,WR_rec 
      
      do j = 0, N+1
      
            ! Extract stencil grid
            x = r(j-1:j+1)
            
            ! Convert to local primitive variables (2nd order conversion)
            do k = j-1,j+1           
                  call U_to_W_comp(u_in(k,:),W(k-j,:))
            enddo

            ! Compute derivative approximations
            sp = (W(1,:) - W(0,:))/(x(1)-x(0))
            sm = (W(0,:) - W(-1,:))/(x(0)-x(-1))
            sd = 0.5*(W(1,:) - W(-1,:))/(x(1)-x(-1))
            
            ! Compute limited slope with MC limiter
            call Minmod_MC(sp,sm,sd,sc)
            
            ! Compute recontructed boundary values
            WL_rec(j,:)   = W(0,:) + 0.5*sc*(x(1)-x(0))
            WR_rec(j-1,:) = W(0,:) - 0.5*sc*(x(0)-x(-1))
      
      enddo
      
      
      end subroutine PLM_rec
      
      !--------------------------------------------
      
      ! Minmod generalized slope limiter
      subroutine Minmod_MC(a,b,c,d)
      real*8, dimension(1,3), intent(in) :: a,b,c
      real*8 :: v_arg(3)
      real*8 :: theta = 2.0
      integer :: i
      
      real*8, dimension(1,3), intent(out) :: d

      ! Generalized Minmod limiter (Kappeli 2016)
               
      do i = 1,3
      
            ! Limiter arguments
            v_arg(1) = theta*a(1,i)
            v_arg(2) = theta*b(1,i)
            v_arg(3) = c(1,i)
      
            if (maxval(v_arg).lt.(0.0)) then
            
                  d(1,i) =maxval(v_arg)
                  
            elseif (minval(v_arg).gt.(0.0)) then
                  
                  d(1,i) = minval(v_arg)
                  
            else
                  d(1,i) = 0.0
            endif     
      enddo

      
      end subroutine Minmod_MC
      end module PLM_reconstruction
