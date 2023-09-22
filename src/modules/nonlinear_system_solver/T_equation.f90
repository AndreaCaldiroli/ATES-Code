	module equation_T
	! Equation for temperature at the steady state
	
	use global_parameters
	use utils, only : calc_ne
	use Cooling_Coefficients
	
	implicit none
	
	contains
	
	subroutine T_equation(N_T_eq,x,fvec,iflag,params)
	
	integer :: N_T_eq,iflag
	real*8  :: x(N_T_eq),fvec(N_T_eq)
   real*8  :: params(25)
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
	
	! Parameters
	nhi    = params(1)
	nhii   = params(2)
	nhei   = params(3)
	nheii  = params(4)
	nheiii = params(5)
  	mup    = params(6)
   mum    = params(7)
   rhov   = params(8)
   coeff  = params(9)
   dr	   = params(10)
   Told   = params(11)
   heaold = params(12)
      
   ! Free electron density
   if (thereis_He) then
		ne = nhii + nheii + 2.0*nheiii
	else
		ne = nhii
	endif	
      
	! Substitutions
	TT = x(1)*T0
	
	!--- Evaluate cooling rates ---!
			
    ! Cooling rate
   reco  =  rec_cool_HII_func(TT)*nhii     & ! HII
         +  rec_cool_HeII_func(TT)*nheii   & ! HeII
         +  rec_cool_HeIII_func(TT)*nheiii   ! HeIII

   !-- Collisional ionization --!
      
   ! Cooling rate
   coio =  2.179e-11*ion_coeff_HI_func(TT)*nhi  	         & ! HI
           + 3.940e-11*ion_coeff_HeI_func(TT)*nhei 		   & ! HeI
	  		  + kb_erg*631515.0*ion_coeff_HeII_func(TT)*nheii    ! HeII

   !-- Bremsstrahlung --!
 
   ! Cooling rate
   brem = 1.426e-27*sqrt(TT)*                       &
     	   (ih**2.0*GF_func(TT,ih)*nhii +             &  ! HII
      	ihe**2.0*GF_func(TT,ihe)*(nheii + nheiii))    ! He
      
   !-- Collisional excitation --!
     
   ! Collisional excitation 
   coex = coex_rate_HI_func(TT)*nhi       &    ! HI
        + coex_rate_HeI_func(TT)*nhei     &    ! HeI
        + coex_rate_HeII_func(TT)*nheii        ! HeII
 
   ! Total cooling rate in erg/(s cm^3)
 	cool = ne*(brem + coex + reco + coio)/q0

	! Equation
   fvec(1) = mum*rhov*x(1) - mup*rhov*Told 		&
           - (g-1.0)*(coeff*x(1) + mup*mum*dr*(heaold - cool))
      
   ! End of subroutine
	end subroutine T_equation
	
	! End of module
	end module equation_T
