	module System_HeH_TR
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	
	implicit none
	
	contains
	
	subroutine ion_system_HeH_TR(Neq,x,fvec,iflag,params)
	
	integer :: Neq,iflag
	real*8  :: x(Neq),fvec(Neq)
	real*8  :: g_hi,g_hei,g_heii,g_heiTR		! Photoionization rates
	real*8  :: b_hi,b_hei,b_heii			! Collisional ionization rates
	real*8  :: A31,q13,q31a,q31b,Q31
	real*8  :: a_hii,a_heii,a_heiii,a_heiTR	! Recombination rates
   real*8  :: params(25)
	real*8  :: n_h,n_he,n_e
	real*8  :: n_hi,n_hii
	real*8  :: n_hei,n_heii,n_heiii,n_heiTR,n_heiSI
	
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
 	
 	! Triplet parameters
 	a_heiTR = params(12)   ! = rcheiTR
 	A31     = params(13)   ! = A31
 	g_heiTR = params(14)   ! = P_HeITR
 	q13     = params(15)   ! = q13
 	q31a    = params(16)   ! = q31a
 	q31b    = params(17)   ! = q31b
 	Q31     = params(18)   ! = Q31
 	
 	
 	! Species densities
 	n_hi    = (1.0-x(1))*n_h
 	n_hii   = x(1)*n_h
 	n_hei   = (1.0 - x(2) - x(3))*n_he 
 	n_heii  = x(2)*n_he
 	n_heiii = x(3)*n_he
 	n_heiSI = (1.0 - x(2) - x(3) - x(4))*n_he
 	n_heiTR = x(4)*n_he
 	
 	
 	! Electron density
   n_e = n_hii + n_heii + 2.0*n_heiii

    ! System of equations      
  	fvec(1) = n_hi*g_hi - a_hii*n_hii*n_e  
  	
  	! New equation for hei - sum of the two equations of Oklopcic
  	fvec(2) =  n_heii*(a_heiTR + a_heii)*n_e	&
  		     - n_heiSI*g_hei					      &
  		     - n_heiTR*g_heiTR
  	
  	fvec(3) = n_heii*g_heii - a_heiii*n_heiii*n_e
  	
  	fvec(4) = - n_heiTR*g_heiTR 		  	      &
  		      + n_e*( n_heii*a_heiTR   	      &
			        + n_heiSI*q13   		      &
		    	    - n_heiTR*(q31a + q31b))	   &
		       - n_heiTR*(A31 + n_hi*Q31)
	
	return 
	
	! End of subroutine
	end subroutine ion_system_HeH_TR
	
	! End of module
	end module System_HeH_TR
