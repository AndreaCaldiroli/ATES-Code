	module System_H
	! Ionization equilibrium system with only H
	
	use global_parameters
	
	implicit none
	
	contains
	
	subroutine ion_system_H(Neq,x,fvec,iflag,params)
	
	integer :: Neq,iflag
	real*8  :: x(Neq),fvec(Neq)
	real*8  :: g_hi                   ! Photoionization rates
	real*8  :: b_hi                   ! Collisional ionization rates
	real*8  :: a_hii                  ! Recombination rates
   real*8  :: params(25)
	real*8  :: n_h,n_e
	real*8  :: n_hi,n_hii
	
	! Coefficients of the system

 	g_hi    = params(1)    ! = P_HI 
 	a_hii   = params(2)    ! = rchiiB 
 	n_h     = params(3)    ! = nh 
   b_hi    = params(4)    ! = a_ion_HI 
 	
 	! Species densities
 	n_hi  = (1.0-x(1))*n_h
 	n_hii = x(1)*n_h
 	
 	! Electron density
   n_e = n_hii
      
    ! System of equations      
  	fvec(1) = n_hi*g_hi + (n_hi*b_hi - a_hii*n_hii)*n_e

	! End of subroutine
	end subroutine ion_system_H
	
	! End of module
	end module System_H
