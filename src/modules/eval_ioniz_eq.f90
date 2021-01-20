      module Evaluate_ioniz_eq
      ! Evaluate the ionization equilibrium equations
      
      use global_parameters
      use Conversion
      use ionization_equilibrium
      
      implicit none
      
      contains
      
      subroutine eval_ioniz_eq(u_in,f_sp_in,u_out,f_sp_out,eta,  &
 		                   heat,cool)
      
      real*8, dimension(1-Ng:N+Ng,3), intent(in) :: u_in 
      real*8, dimension(1-Ng:N+Ng,5), intent(in) :: f_sp_in 
      real*8, dimension(1-Ng:N+Ng,3) :: W
      real*8, dimension(1-Ng:N+Ng) :: rho,v,p,T
      real*8, dimension(1-Ng:N+Ng) :: nhi,nhii
      real*8, dimension(1-Ng:N+Ng) :: nhei,nheii,nheiii
      real*8, dimension(1-Ng:N+Ng) :: ne
      real*8, dimension(1-Ng:N+Ng) :: n_tot
      real*8, dimension(1-Ng:N+Ng,3), intent(out) :: u_out
      real*8, dimension(1-Ng:N+Ng,5), intent(out) :: f_sp_out
      real*8, dimension(1-Ng:N+Ng),   intent(out) :: heat,cool,eta
      
      ! Exctract primitive variables
      call U_to_W(u_in,W)
      rho = W(:,1)
      v   = W(:,2)
      p   = W(:,3)
      
	! Evaluate species densities
  	nhi    = rho*f_sp_in(:,1)
	nhii   = rho*f_sp_in(:,2)
	nhei   = rho*f_sp_in(:,3)
	nheii  = rho*f_sp_in(:,4)
	nheiii = rho*f_sp_in(:,5)
  	ne     = nhii + nheii + 2.0*nheiii
  	
  	! Total number density
  	n_tot = nhi + nhii + nhei + nheii + nheiii 
  	  	 
  	! Temperature profile
  	T = p/(n_tot + ne)
  	
	! Evaluate ionization equilibrium
	call ioniz_eq(T,rho,f_sp_in,rho,f_sp_out,heat,cool,eta)
	
	!Evaluate partial densities
  	nhi    = rho*f_sp_out(:,1)
	nhii   = rho*f_sp_out(:,2)
	nhei   = rho*f_sp_out(:,3)
	nheii  = rho*f_sp_out(:,4)
	nheiii = rho*f_sp_out(:,5)
  	ne     = nhii + nheii +2.0*nheiii
	
	! Total number density
  	n_tot = nhi + nhii + nhei + nheii + nheiii 
  	
	! Evaluate updated pressure
  	p = (n_tot + ne)*T 
  	
  	! Convert to primitive profiles
      W(:,1) = rho
      W(:,2) = v
      W(:,3) = p
      
      ! Revert to conservative 
  	call W_to_U(W,u_out)
  	
  	end subroutine eval_ioniz_eq
           
      ! End of module
      end module Evaluate_ioniz_eq
