	module post_processing
	! Subroutine to correct the output ionization profiles
	!	taking into account the ionization term
	
	use global_parameters
	use utils
	use System_implicit_adv_H
	use System_implicit_adv_HeH
	use System_implicit_adv_HeH_TR
	use Cooling_Coefficients
	use utils_ion_eq
	use output_write	
	use equation_T
	
	implicit none
	
	
	contains
	
	subroutine post_process_adv(rho,v,p,T_in,heat,cool,eta,   &
                                  nhi_in,nhii_in,		    &
                                  nhei_in,nheii_in,nheiii_in,   & 
                                  nheiTR_in)
                                  
                                  
	real*8, dimension(1-Ng:N+Ng), intent(in) :: rho,v,p,T_in
   real*8, dimension(1-Ng:N+Ng), intent(in) :: heat,cool
   real*8, dimension(1-Ng:N+Ng), intent(in) :: eta
   real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi_in,nhii_in
   real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei_in,nheii_in,   &
      							  				        nheiii_in,nheiTR_in
	
	integer i,j,k
	 
	real*8, dimension(1-Ng:N+Ng) ::  T_K,p_out,T_out     ! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng) ::  nh,nhe,ne,n_tot 
	
	! Dummy variable
	real*8, dimension(1-Ng:N+Ng) ::  dum_v

   ! Photo ionization rates
   real*8, dimension(1-Ng:N+Ng) ::  P_HI,P_HeI,P_HeII,P_HeITR

   ! Recombination coefficients
   real*8, dimension(1-Ng:N+Ng) ::  rchiiB,rcheiiB,rcheiiiB,rcheiTR
      
   real*8, dimension(1-Ng:N+Ng) ::  q13,q31a,q31b
	real*8 :: A31,Q31  
 	
 	! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng) ::  a_ion_HI,a_ion_HeI,a_ion_HeII 
      
   ! Heating, cooling
   real*8, dimension(1-Ng:N+Ng) ::  theat,tcool
      
 	! Updated species densities
   real*8, dimension(1-Ng:N+Ng) :: nhi,nhii
   real*8, dimension(1-Ng:N+Ng) :: nhei,nheii,nheiii,nheiTR,nheiS
   real*8, dimension(1-Ng:N+Ng) :: mmw
   real*8, dimension(1-Ng:N+Ng) :: nhi_w,nhii_w
   real*8, dimension(1-Ng:N+Ng) :: nhei_w,nheii_w,nheiii_w,nheiTR_w
      							  
	      
	      
   real*8 :: TT                              ! Temperature component
   real*8 :: PIR_1,PIR_15,PIR_2,PIR_TR       ! Photoionization rates
   real*8 :: deltal                          ! Optical depth
   real*8 :: dr	                           ! Grid spacing
   real*8 :: iup_1,ilo_1,        &           ! Photoheating integral variables
             iup_15,ilo_15,      &
             iup_2,ilo_2,        &
             iup_TR,ilo_TR,	   &
             iup_f,ilo_f 
   real*8 :: elo,eup                         ! Energy parameters
   real*8 :: tol,dpmpar                      ! Equilibrium system setup
   real*8 :: Hea_1 		          	         ! Heating rates
   real*8 :: brem,coex,coio,reco             ! Cooling rates
   real*8 :: iup_H,ilo_H                     ! Heating rate integral variables         
      
      
	! Substitution in the ODE solution
	real*8 :: As


   real*8, dimension(25) :: params
   real*8, dimension(12) :: paramsT
      
	real*8 :: rhop,rhom,vp,mum,mup,vm
	real*8 :: sys_sol_T(1), sys_x_T(1)
   real*8 :: wa_T(8)
      
      
   !----------------------------------------------------------!      
      
   ! Global parameters
      
   ! Numerical tolerance for system solution
   tol = sqrt(dpmpar(1))
      
   !----------------------------------!
	
	! Preliminary profiles extraction
	
	! Dimensional total number density profile and temperature
	T_K = T_in*T0
	         
   ! Initialize vectors
   nhi    = nhi_in*n0
	nhii   = nhii_in*n0
   if (thereis_He) then
		nhei   = nhei_in*n0
		nheii  = nheii_in*n0
		nheiii = nheiii_in*n0
		if (thereis_HeITR) nheiTR = nheiTR_in*n0	
	endif
	!----------------------------------!
	
	! Iterate the post processing 
	do k = 1,10	! Usually 10 gives a good convergence
	
	! Use singlet if included
	nheiS = nhei
	if (thereis_HeITR) nheiS = nhei - nheiTR
	nh  = nhi  + nhii 
	nhe = nheiS + nheii + nheiii 
	if (thereis_HeITR) nhe = nhe + nheiTR

	! Free electron density (assuming overall neutrality)
	call calc_ne(nhii,nheii,nheiii,ne) 
      
   ! Calculate the photoionization rates

	if (thereis_He) then
		call PH_heat_HHe(nhi,nhei,nheii,nheiTR,	&
					 P_HI,P_HeI,P_HeII,P_HeITR,dum_v,dum_v)
  	else
	  	call PH_heat_H(nhi,P_HI,dum_v,dum_v)
  	endif
      
   !---- Recombination rates ----!
 	
	call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii,  	&
	  			   rchiiB,rcheiiB,rcheiiiB, 			&
				   a_ion_HI,a_ion_HeI,a_ion_HeII,dum_v)

 	
 	if (thereis_HeITR) then
		call HeITR_coeffs(T_K,rcheiTR,rcheiiB,A31,q13,q31a,q31b,Q31)
		! NOTE: rcheiiB is alpha1 from Oklopcic - being overwritten

	endif
	
   !----------------------------------!
      
   ! Evolve species including the advection term in the 
   ! 	ODE form
   ! Note: we are using point values here instead of 
   ! 	volume averages; they agree up to O(dr^2)
      
   ! The ionization fraction at the inner boundary are taken 
	!	from the input vectors (completely neutral atmosphere)
     

   ! Loop to solve the differential equation
   ! It is implicitly assumed that the velocity fields does 
   !	not change by including the advection term
      
      
	      
   if (.not.thereis_He) then
	      
		do j = 2-Ng,N+Ng 
	
			! Substitutions
			dr  = (r(j) - r(j-1))*R0
			As  = dr/(v(j-1)*v0)

			! Advection coeff.
			params(1) = As 
			params(2) = nhi(j-1)/nh(j-1) 
			params(3) = nh(j) 
			params(4) = P_HI(j)
			params(5) = rchiiB(j)    
			params(6) = a_ion_HI(j)
			
			! Initial guess of solution
			sys_x(1) = nhi(j)/nh(j)     
			
			! Call hybrd1 routine (from minpack)
			call hybrd1(adv_implicit_H,N_eq,sys_x,sys_sol,   &
						tol,info,wa,lwa,params) 
			
			! Extract solution profiles	
			nhi(j)    = sys_x(1)*nh(j)
			nhii(j)   = (1.0 - sys_x(1))*nh(j)

      	enddo
      	
		   ! Force condition of zero helium
		   nhei   = 0.0
		   nheii  = 0.0
		   nheiii = 0.0			
		   nheiTR = 0.0

	else
		
		do j = 2-Ng,N+Ng 
	
			! Substitutions
			dr  = (r(j) - r(j-1))*R0
			As  = dr/(v(j-1)*v0)

			! Advection coeff.
			params(1)  = As 
			params(2)  = nhi(j-1)/nh(j-1) 
			params(3)  = nhei(j-1)/nhe(j-1)
			params(4)  = nheiii(j-1)/nhe(j-1)
			params(5)  = nh(j) 
			params(6)  = P_HI(j)
			params(7)  = P_HeI(j)
			params(8)  = P_HeII(j)
			params(9)  = rchiiB(j)  
			params(10) = rcheiiB(j)
			params(11) = rcheiiiB(j) 
			params(12) = a_ion_HI(j)
			params(13) = a_ion_HeI(j)
			params(14) = a_ion_HeII(j)

			! Add more if HeITR is present
			if (thereis_HeITR) then 
				params(15) = rcheiTR(j)
				params(16) = A31
				params(17) = P_HeITR(j)
				params(18) = q13(j)
				params(19) = q31a(j)
				params(20) = q31b(j)
				params(21) = Q31		
				params(22) = nheiTR(j-1)/nhe(j-1)
			endif 
			
			! Initial guess of solution
			sys_x(1) = nhi(j)/nh(j) 
			sys_x(2) = nhei(j)/nhe(j)
			sys_x(3) = nheiii(j)/nhe(j)
			if (thereis_HeITR) sys_x(4) = nheiTR(j)/nhe(j) 
			
			! Call hybrd1 routine (from minpack)
			if (thereis_HeITR) then 
				call hybrd1(adv_implicit_HeH_TR,N_eq,sys_x,sys_sol,   &
							tol,info,wa,lwa,params) 
			else
				call hybrd1(adv_implicit_HeH,N_eq,sys_x,sys_sol,   &
							tol,info,wa,lwa,params) 
			endif
				
			! Extract solution profiles	
			nhi(j)    = sys_x(1)*nh(j)
			nhii(j)   = (1.0 - sys_x(1))*nh(j)
			nhei(j)   = sys_x(2)*nhe(j)
			nheii(j)  = (1.0 - sys_x(2) - sys_x(3))*nhe(j)
			nheiii(j) = sys_x(3)*nhe(j)
			if (thereis_HeITR) then
				nheiTR(j) = sys_x(4)*nhe(j) 
			else
				nheiTR(j) = 0.0
			endif
			
		enddo
		
	endif ! End if thereis_He
	      
      
   !---------------------------------------------------------
      
   !------- Fix stationarity of new pressure profile -------!
      
   !---- Update densities and temperature ----!
	
	! Number densities      
   nheiS = nhei
	if (thereis_HeITR) nheiS = nhei - nheiTR
	nh  = nhi  + nhii 
	nhe = nheiS + nheii + nheiii
	if (thereis_HeITR) nhe = nhe + nheiTR

   ! Total number density
   call calc_ntot(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_tot)
	 
   ! Free electron density (assuming overall neutrality)
   call calc_ne(nhii,nheii,nheiii,ne)
      
	!----------------------------------!
	
	!---- Update photoheating rate ----!
	
	if (thereis_He) then
		call PH_heat_HHe(nhi,nhei,nheii,nheiTR,					&
						 dum_v,dum_v,dum_v,dum_v,theat,dum_v)
  	else
	  	call PH_heat_H(nhi,dum_v,theat,dum_v)
  	endif

	! Adimensionalize
	theat = theat/q0		

	!----------------------------------!
	
	!---- Solve stationary energy equation ----!
	! This procedure uses the same velocity profile and 
	!	the ionization profile after the advection correction
	
	! Initialize temperature at ghost cells
	T_out = T_K/T0
	
	! Calculate mean molecular weight
	call calc_mmw(nh,nhe,ne,mmw)

	do j = 3-Ng,N+Ng ! Start from first computational cell
		
		! Substitutions
		rhop = rho(j)
		rhom = rho(j-1)
		vm = v(j-1) 
		vp = v(j)
		dr = r(j) - r(j-1)
		mum = mmw(j-1)
		mup = mmw(j)

	 	!--- Solve equation for temperature implicitly ---!
		
		! Parameters
		paramsT(1)  = nhi(j)
	 	paramsT(2)  = nhii(j)
	 	paramsT(3)  = nhei(j)
	 	paramsT(4)  = nheii(j)
	 	paramsT(5)  = nheiii(j)
	 	paramsT(6)  = mmw(j)
	 	paramsT(7)  = mmw(j-1)
	 	paramsT(8)  = rhop*vp
	 	paramsT(9)  = mum*vp*(rhop-rhom)
	 	paramsT(10) = dr
	 	paramsT(11) = T_out(j-1)
	 	paramsT(12) = theat(j) 
	 	
	 	! Initial guess of solution
		sys_x_T(1) = T_out(j) 
		
		! Call hybrd1 routine (from minpack)
		call hybrd1(T_equation,1,sys_x_T,sys_sol_T,   &
		            tol,info,wa_T,8,paramsT) 
		
		! Extract solution profiles	
		T_out(j) = sys_x_T(1) 	
	 	
	enddo
	
	! Update pressure and temperature
	p_out = (n_tot + ne)/n0*T_out
	T_K = T_out*T0
      
	enddo ! End loop on post processing 
      
   ! ---------------------------- ! 
      
	!---- Update cooling rates ----!	
	
	call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii,  	&
	  			   dum_v,dum_v,dum_v, 			&
				   dum_v,dum_v,dum_v,tcool)

	! Adimensionalize
	tcool = tcool/q0

	!----------------------------------!
 
   ! Adimensionalize ion densities before writing
   nhi_w    = nhi/n0
   nhii_w   = nhii/n0
   nhei_w   = nhei/n0
   nheii_w  = nheii/n0
   nheiii_w = nheiii/n0
   nheiTR_w = nheiTR/n0

   !----------------------------------!
      
   ! Write updated thermodynamic and ionization profiles
   call write_output(rho,v,p_out,T_out,theat,tcool,eta,    	  &
                     nhi_w,nhii_w,nhei_w,nheii_w,nheiii_w,	  &
                     nheiTR_w,'ad')

	! End of subroutine
	end subroutine post_process_adv
	
	! End of module
	end module post_processing
