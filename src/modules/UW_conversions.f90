      module Conversion
      ! Conversion conservative-primitive variables subroutines
      
      use global_parameters
      
      implicit none
      
      contains
      
      !-----------------------------------------------------------!
      
      ! Conservative to primitive (one component)
      subroutine U_to_W_comp(U_in,W_out)
      real*8, intent(in) :: U_in(3)
      real*8, intent(out) :: W_out(3)
      
      W_out(1) = U_in(1)
      W_out(2) = U_in(2)/U_in(1)
      W_out(3) = (g-1.0)*(U_in(3)-0.5*U_in(2)*U_in(2)/U_in(1))
 
      end subroutine U_to_W_comp
      
      !-----------------------------------------------------------!
      
      ! Conservative to primitive (one component)
      subroutine W_to_U_comp(W_in,U_out)
      real*8, intent(in) :: W_in(3)
      real*8, intent(out) :: U_out(3)
      
      U_out(1) = W_in(1)
      U_out(2) = W_in(1)*W_in(2)
      U_out(3) = 0.5*W_in(1)*W_in(2)**2.0 + W_in(3)/(g-1.0)
 
      end subroutine W_to_U_comp
      
      !-----------------------------------------------------------!
      
      ! Primitive to conservative
      subroutine W_to_U(W_in,U_out)
      real*8, intent(in) :: W_in(1-Ng:N+Ng,3)
      real*8, intent(out) :: U_out(1-Ng:N+Ng,3)
      
      U_out(:,1) = W_in(:,1)
      U_out(:,2) = W_in(:,1)*W_in(:,2)
      U_out(:,3) = 0.5*W_in(:,1)*W_in(:,2)**2.0    &
                   + W_in(:,3)/(g-1.0)
 
      end subroutine W_to_U
      
      
      !-----------------------------------------------------------!
      
      
      ! Conservative to primitive
      subroutine U_to_W(U_in,W_out)
      real*8, intent(in) :: U_in(1-Ng:N+Ng,3)
      real*8, intent(out) :: W_out(1-Ng:N+Ng,3)
      
      W_out(:,1) = U_in(:,1)
      W_out(:,2) = U_in(:,2)/U_in(:,1)
      W_out(:,3) = (g-1.0)*                   &
                   (U_in(:,3)-0.5*U_in(:,2)*U_in(:,2)/U_in(:,1))
 
      end subroutine U_to_W
        
      !-----------------------------------------------------------!
      
      
      

      end module Conversion
