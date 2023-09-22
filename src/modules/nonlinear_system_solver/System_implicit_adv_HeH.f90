	module System_implicit_adv_HeH
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	
	implicit none
	
	contains
	
	subroutine adv_implicit_HeH(N_eq,x,fvec,iflag,params)
	
	integer :: N_eq,iflag
	real*8  :: x(N_eq),fvec(N_eq)
	real*8  :: xhi_old,xhei_old,xheiii_old
	real*8  :: ghi,ghei,gheii
	real*8  :: xhi,xhii
	real*8  :: xhei,xheii,xheiii
	real*8  :: xe
	real*8  :: c1
	real*8  :: n_h
	real*8  :: ahii,aheii,aheiii
	real*8  :: ionhi,ionhei,ionheii
   real*8  :: params(25)
	
	! Coefficients of the system

 	c1         = params(1)    ! = dr/v
 	xhi_old    = params(2)    ! = nhi/nhe 
 	xhei_old   = params(3)    ! = nheii/nhe 
 	xheiii_old = params(4)    ! = nheiii/nhe
 	n_h        = params(5)    ! = nh
 	ghi        = params(6)    ! = P_HI 
 	ghei       = params(7)    ! = P_HeI = nh 
 	gheii      = params(8)    ! = P_HeII = nhe 
   ahii       = params(9)    ! = rchiiB  
 	aheii      = params(10)   ! = rcheiiB 
 	aheiii     = params(11)   ! = rcheiiiB  
	ionhi	   = params(12)   ! = a_ion_HI
	ionhei     = params(13)   ! = a_ion_HeI
	ionheii    = params(14)   ! = a_ion_HEII

	! Substitutions
	xhi    = x(1)
	xhii   = 1.0 - x(1)
	xhei   = x(2)
   xheii  = 1.0 - x(2) - x(3)
	xheiii = x(3)
		
 	! Electron density
   xe = xhii + HeH*(xheii + 2.0*xheiii)
      
   ! System of equations      
  	fvec(1) =  xhi_old - xhi   					        		&
  	        + c1*(-(ghi+ionhi)*xhi    + ahii*xhii*xe*n_h) 	   
  		     		    
  	fvec(2) =  xhei_old - xhei 						  		    &
  		    + c1*(-(ghei+ionhei)*xhei  + aheii*xheii*xe*n_h) 
  	 	      	 	    
  	fvec(3) =  xheiii_old - xheiii					  			&
  	        + c1*((gheii+ionheii)*xheii - aheiii*xheiii*xe*n_h) 
  		    
	! End of subroutine
	end subroutine adv_implicit_HeH
	
	! End of module
	end module System_implicit_adv_HeH
