	module System_implicit_adv_HeH_TR
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	
	implicit none
	
	contains
	
	subroutine adv_implicit_HeH_TR(Neq,x,fvec,iflag,params)
	
	integer :: Neq,iflag
	real*8  :: x(Neq),fvec(Neq)
	real*8  :: xhi_old,xhei_old,xheiii_old,xheiTR_old
	real*8  :: ghi,ghei,gheii,gheiTR
	real*8  :: xhi,xhii
	real*8  :: xhei,xheii,xheiii
	real*8  :: xheiTR, xheiS
	real*8  :: xe
	real*8  :: c1
	real*8  :: n_h
	real*8  :: ahii,aheii,aheiii,aheiTR
	real*8  :: ionhi,ionhei,ionheii
   real*8  :: params(25)
   real*8  :: A31,q13,q31a,q31b,Q31
	
	! Coefficients of the system
 	c1         = params(1)    ! = dr/v
 	xhi_old    = params(2)    ! = nhi/nhe 
 	xhei_old   = params(3)    ! = nheii/nhe 
 	xheiii_old = params(4)    ! = nheiii/nhe
 	n_h        = params(5)    ! = nh
 	ghi        = params(6)    ! = P_HI 
 	ghei       = params(7)    ! = P_HeI 
 	gheii      = params(8)    ! = P_HeII 
   ahii       = params(9)    ! = rchiiB  
 	aheii      = params(10)   ! = rcheiiB 
 	aheiii     = params(11)   ! = rcheiiiB  
	ionhi	     = params(12)   ! = a_ion_HI
	ionhei     = params(13)   ! = a_ion_HeI
	ionheii    = params(14)   ! = a_ion_HEII
	aheiTR     = params(15)	  ! = rcheiTR
	A31	     = params(16)   ! = A31
	gheiTR     = params(17)   ! = P_HeITR
	q13        = params(18)   ! = q13
	q31a	     = params(19)   ! = q31a
	q31b	     = params(20)   ! = q31b
	Q31 	     = params(21)   ! = Q31
	xheiTR_old = params(22)   ! = nheiTR/nh

	! Substitutions
	xhi    = x(1)
	xhii   = 1.0 - x(1)
	xhei   = x(2)
   xheii  = 1.0 - x(2) - x(3)
	xheiii = x(3)
	xheiS  = x(2) - x(4)
	xheiTR = x(4)
		
 	! Electron density
   xe = xhii + HeH*(xheii + 2.0*xheiii)
      
    ! System of equations      
  	fvec(1) =  xhi_old - xhi + c1*(  		&
  	           - (ghi + ionhi*xe*n_h)*xhi 	&
  	           +  ahii*xhii*xe*n_h) 	   
  	
  	fvec(2) =  xhei_old - xhei + c1*(  	       &
  		       xheii*(aheiTR + aheii)*xe*n_h	 &
  		     - xheiS*ghei - xheiTR*gheiTR) 
  	 	      	 	    
  	fvec(3) =  xheiii_old - xheiii + c1*(	    &
  	           (gheii + ionheii*xe*n_h)*xheii  & 
  	          - aheiii*xheiii*xe*n_h) 
  	           
	fvec(4) =    xheiTR_old - xheiTR + c1*(   &
		     - gheiTR*xheiTR				         &
		     + (xheii*aheiTR + xheiS*q13    	&
		     -  xheiTR*(q31a + q31b))*xe*n_h	&
		     - xheiTR*(A31 + xhi*Q31*n_h))
	
	! End of subroutine
	end subroutine adv_implicit_HeH_TR
	
	! End of module
	end module System_implicit_adv_HeH_TR
