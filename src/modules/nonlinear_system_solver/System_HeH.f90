	module System_HeH
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	
	implicit none
	
	contains
	
	subroutine ion_system_HeH(N_eq,x,fvec,iflag,params)
	
	integer :: N_eq,iflag
	real*8  :: x(N_eq),fvec(N_eq)
	real*8  :: g_hi,g_hei,g_heii		! Photoionization rates
	real*8  :: b_hi,b_hei,b_heii		! Collisional ionization rates
	real*8  :: a_hii,a_heii,a_heiii	    ! Recombination rates
   real*8  :: params(25)
	real*8  :: n_h,n_he,n_e
	real*8  :: n_hi,n_hii
	real*8  :: n_hei,n_heii,n_heiii
	
	! Coefficients of the system

 	g_hi    = params(1)    ! = P_HI 
 	g_hei   = params(2)    ! = P_HeI
 	g_heii  = params(3)    ! = P_HeII 
 	a_hii   = params(4)    ! = rchiiB 
 	a_heii  = params(5)    ! = rcheiiB 
 	a_heiii = params(6)    ! = rcheiiiB 
 	n_h     = params(7)    ! = nh 
 	n_he    = params(8)    ! = nhe 
   b_hi    = params(9)    ! = a_ion_HI 
 	b_hei   = params(10)   ! = a_ion_HeI 
 	b_heii  = params(11)   ! = a_ion_HeII 
 	
 	! Species densities
 	n_hi    = (1.0-x(1))*n_h
 	n_hii   = x(1)*n_h
 	n_hei   = (1.0 - x(2) - x(3))*n_he 
 	n_heii  = x(2)*n_he
 	n_heiii = x(3)*n_he 
 	
 	! Electron density
   n_e = n_hii + n_heii + 2.0*n_heiii
      
      
      ! System of equations      
  	fvec(1) = n_hi*g_hi + (n_hi*b_hi - a_hii*n_hii)*n_e      
  	fvec(2) = n_hei*g_hei + (n_hei*b_hei - a_heii*n_heii)*n_e
  	fvec(3) = n_heii*g_heii + (n_heii*b_heii - a_heiii*n_heiii)*n_e
	
	return 
	
	! End of subroutine
	end subroutine ion_system_HeH
	
	! End of module
	end module System_HeH
