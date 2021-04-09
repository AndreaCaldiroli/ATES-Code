	module System_implicit_adv_H
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	
	implicit none
	
	contains
	
	subroutine adv_implicit_H(N_eq,x,fvec,iflag,params)
	
	integer :: N_eq,iflag
	real*8  :: x(N_eq),fvec(N_eq)
	real*8  :: xhi_old
	real*8  :: ghi
	real*8  :: xhi,xhii,xe
	real*8  :: c1
	real*8  :: n_h
	real*8  :: ahii
      real*8  :: params(5)

	
	! Coefficients of the system
 	c1      = params(1)    ! = dr/v
 	xhi_old = params(2)    ! = nhi/nh
 	n_h     = params(3)    ! = nh
 	ghi     = params(4)    ! = P_HI  
      ahii    = params(5)    ! = rchiiB  
	
	! Substitutions
	xhi = x(1)
	xhii = 1.0 - x(1)

 	! Electron density
      xe = xhii
      
      ! System of equations      
  	fvec(1) = xhi_old + c1*(-ghi*xhi + ahii*xhii*xe*n_h) - x(1)
	
	
	end subroutine adv_implicit_H
	
	! End of module
	end module System_implicit_adv_H
