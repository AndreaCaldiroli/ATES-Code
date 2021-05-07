	module equation_T
	! Ionization equilibrium system with both H and He
	
	use global_parameters
	use Cooling_Coefficients
	
	implicit none
	
	contains
	
	subroutine T_equation(N_eq,x,fvec,iflag,params)
	
	integer :: N_eq,iflag
	real*8  :: x(N_eq),fvec(N_eq)
      real*8  :: params(14)
	real*8  :: nhi,nhii
	real*8  :: nhei,nheii,nheiii
	real*8  :: ne
	real*8  :: mum,mup
	real*8  :: rhov
	real*8  :: coeff
	real*8  :: dr
	real*8  :: Told,heaold
	real*8  :: reco,coio,brem,coex,cool
	real*8  :: TT

	nhi    = params(1)
	nhii   = params(2)
	nhei   = params(3)
	nheii  = params(4)
	nheiii = params(5)
  	mup    = params(6)
      mum    = params(7)
      rhov 	 = params(8)
      coeff  = params(9)
      dr	 = params(10)
      Told   = params(11)
      heaold = params(12)
      
	
	TT = x(1)*T0
	ne = nhii + nheii + 2.0*nheiii
	
			
      ! Cooling rate
      reco  =  rec_cool_HII(TT)*nhii     & ! HII
             + rec_cool_HeII(TT)*nheii   & ! HeII
             + rec_cool_HeIII(TT)*nheiii   ! HeIII


      !-- Collisional ionization --!
      
      
      ! Cooling rate
      coio =  2.179e-11*ion_coeff_HI(TT)*nhi  	       & ! HI
            + 3.940e-11*ion_coeff_HeI(TT)*nhei 		 & ! HeI
		+ kb_erg*631515.0*ion_coeff_HeII(TT)*nheii   ! HeII

      !-- Bremsstrahlung --!
 
      ! Cooling rate
      brem = 1.426e-27*sqrt(TT)*                       &
      	 (ih**2.0*GF(TT,ih)*nhii +                 &  ! HII
      	 ihe**2.0*GF(TT,ihe)*(nheii + nheiii))        ! He
      
      !-- Collisional excitation --!
     
      ! Collisional excitation 
      coex = coex_rate_HI(TT)*nhi       &    ! HI
           + coex_rate_HeI(TT)*nhei     &    ! HeI
           + coex_rate_HeII(TT)*nheii        ! HeII

 
      ! Total cooling rate in erg/(s cm^3)
 	cool = ne*(brem + coex + reco + coio)/q0
      
      fvec(1) = mum*rhov*x(1) - mup*rhov*Told 		&
      	  - (g-1.0)*( coeff*x(1) + mup*mum*dr*(heaold - cool))
      
      	
	end subroutine T_equation
	
	! End of module
	end module equation_T
