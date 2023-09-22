   module ionization_equilibrium
	! Evaluate the ionization structure and the heating and cooling functions for a given temperature 
	
	use global_parameters
   use utils
   use utils_ion_eq
   use System_HeH                ! Equilibrium equations
	use System_HeH_TR 
   use System_H
	
   implicit none

	contains 
	
	subroutine ioniz_eq(T_in,n_in,f_sp_in,n_out, &
      	 		  f_sp_out,heat_out,cool_out,q)
      	 		  
	integer :: j

	real*8, dimension(1-Ng:N+Ng),   intent(in) :: T_in,n_in
	real*8, dimension(1-Ng:N+Ng,6), intent(in) :: f_sp_in
	 
	real*8, dimension(1-Ng:N+Ng) ::  T_K      ! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng) ::  nh,nhi,nhii,                   & ! Species densities
	                                 nhe,nhei,nheii,nheiii,nheiTR,  &	
	                                 ne,n_in_dim
	
   ! Photo ionization rates
   real*8, dimension(1-Ng:N+Ng) ::  P_HI,P_HeI,P_HeII,P_HeITR
                       	
   ! Heating, cooling
   real*8, dimension(1-Ng:N+Ng) ::  heat,cool    
                 
   ! Recombination coefficients
   real*8, dimension(1-Ng:N+Ng) ::  rchiiB,rcheiiB,rcheiiiB,rcheiTR 
              
	real*8, dimension(1-Ng:N+Ng) :: q13,q31a,q31b
	real*8 :: A31,Q31
								   
   ! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng) ::  a_ion_HI,a_ion_HeI,a_ion_HeII    
	
	! Equilibrium system setup      
   real*8 :: tol,dpmpar                      
   real*8, dimension(25) :: params
      
   ! Output density
   real*8, dimension(1-Ng:N+Ng),intent(out) :: n_out 
      
   ! Output heating,cooling and absorbed energy 
   real*8, dimension(1-Ng:N+Ng),intent(out) :: heat_out,cool_out,q 
      
   ! Output species fractions          
   real*8, dimension(1-Ng:N+Ng,6),intent(out) :: f_sp_out
	
   !----------------------------------------------------------!      
   ! Global parameters
      
   ! Numerical tolerance for system solution
   tol = sqrt(dpmpar(1))

	!----------------------------------!
	
	! Preliminary profiles exctraction
	
	! Dimensional total number density profile and temperature
	n_in_dim = n_in*n0
	T_K      = T_in*T0
		
	! Extract species profiles
	nhi    = f_sp_in(:,1)*n_in_dim    ! HI
	nhii   = f_sp_in(:,2)*n_in_dim    ! HII
	if (thereis_He) then

		nhei   = f_sp_in(:,3)*n_in_dim    ! HeI
		nheii  = f_sp_in(:,4)*n_in_dim    ! HeII
		nheiii = f_sp_in(:,5)*n_in_dim    ! HeIII
		nheiTR = f_sp_in(:,6)*n_in_dim    ! HeITR
	
	else 
		
		! Enforce condition of zero helium
	    nhei   = 0.0
        nheii  = 0.0
        nheiii = 0.0
        nheiTR = 0.0

    endif
	
	! Total number densities      
	nh  = nhi  + nhii
	nhe = nhei + nheii + nheiii
	
	! Free electron density (assuming overall neutrality)
	call calc_ne(nhii,nheii,nheiii,ne)

	!----------------------------------!
      
    !---- Photoionization and photoheating ----!
      
	if (thereis_He) then
      	call PH_heat_HHe(nhi,nhei,nheii,nheiTR,	&
      			     P_HI,P_HeI,P_HeII,P_HeITR,heat,q)
	else
		call PH_heat_H(nhi,P_HI,heat,q)
	endif
	
	
   !----------------------------------!
      
   ! Evaluate cooling rates and recombination/collisional 
   ! 	ionization rates
      
    call eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii,  	&
			   	   	rchiiB,rcheiiB,rcheiiiB, 			&
			    	a_ion_HI,a_ion_HeI,a_ion_HeII,cool)
      
	if (thereis_HeITR) then
		call HeITR_coeffs(T_K,rcheiTR,rcheiiB,A31,q13,q31a,q31b,Q31)
		! NOTE: rcheiiB is alpha1 from Oklopcic - being overwritten
	endif
	
   !----------------------------------!
      
   ! Ionization equilibrium system solution	

	if (.not.thereis_He) then ! If no helium

		do j = N+Ng,1-Ng,-1
		
			! Ionization equilibrium system setup
			params(1) = P_HI(j)
			params(2) = rchiiB(j)
			params(3) = nh(j)
			params(4) = a_ion_HI(j)

			 ! Initial guess
			if (count.le.0) then
				if(r(j).le.(1.5))then
					sys_x(1) = r(j)-0.5
				else
					sys_x(1) = 1.0
				endif     
			else

				sys_x(1) = nhii(j)/nh(j)

			endif

		 	! Call hybrd1 routine (from minpack)
			call hybrd1(ion_system_H,N_eq,sys_x,sys_sol,   &
				      tol,info,wa,lwa,params) 

			! Extract solution profiles	
			nhi(j)    = nh(j)*(1.0 - sys_x(1))
			nhii(j)   = nh(j)*sys_x(1)

		enddo
		
		nhei   = 0.0
		nheii  = 0.0
		nheiii = 0.0
		nheiTR = 0.0

	else  ! If there is helium
	
		do j = N+Ng,1-Ng,-1
		
			! System coefficients
			params(1)  = P_HI(j)
			params(2)  = P_HeI(j)
			params(3)  = P_HeII(j)
			params(4)  = rchiiB(j)
			params(5)  = rcheiiB(j)
			params(6)  = rcheiiiB(j)
			params(7)  = nh(j)
			params(8)  = nhe(j)
			params(9)  = a_ion_HI(j) 
			params(10) = a_ion_HeI(j) 
			params(11) = a_ion_HeII(j)  
			
			! Add more if HeITR is present
			if (thereis_HeITR) then 
				params(12) = rcheiTR(j)
				params(13) = A31
				params(14) = P_HeITR(j)
				params(15) = q13(j)
				params(16) = q31a(j)
				params(17) = q31b(j)
				params(18) = Q31
			endif

			! Initial guess
			if (count .eq. 0) then   
				if (j .eq. N+Ng) then
					sys_x(1) = 1.0
					sys_x(2) = 1.0 
					sys_x(3) = 1.0 
					if (thereis_HeITR) sys_x(4) = 0.01
				else
					sys_x(1) = nhii(j+1)/nh(j+1)
					sys_x(2) = nheii(j+1)/nhe(j+1)
					sys_x(3) = nheiii(j+1)/nhe(j+1)
					if (thereis_HeITR) &
						sys_x(4) = nheiTR(j+1)/nhe(j+1) 
				endif 
			else
				sys_x(1) = nhii(j)/nh(j)
				sys_x(2) = nheii(j)/nhe(j)
				sys_x(3) = nheiii(j)/nhe(j)
				if (thereis_HeITR) sys_x(4) = nheiTR(j)/nhe(j) 
			endif
			
		 	! Call hybrd1 routine (from minpack)
			if (thereis_HeITR) then
				call hybrd1(ion_system_HeH_TR,N_eq,sys_x,sys_sol,   &
						    tol,info,wa,lwa,params) 
			else
				call hybrd1(ion_system_HeH,N_eq,sys_x,sys_sol,      &
                            tol,info,wa,lwa,params)
			endif
			
			! Extract solution profiles	
			nhi(j)    = nh(j)*(1.0 - sys_x(1))
			nhii(j)   = nh(j)*sys_x(1)
			nhei(j)   = nhe(j)*(1.0 - sys_x(2) - sys_x(3))
			nheii(j)  = nhe(j)*sys_x(2)
			nheiii(j) = nhe(j)*sys_x(3)
			if (thereis_HeITR) nheiTR(j) = nhe(j)*sys_x(4)
			
		enddo
	
	endif

	
	! Density with atomic numbers
   call calc_rho(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_out)    

   ! Abundancies profiles
   f_sp_out(:,1) = nhi/n_out
   f_sp_out(:,2) = nhii/n_out
   f_sp_out(:,3) = nhei/n_out
   f_sp_out(:,4) = nheii/n_out
   f_sp_out(:,5) = nheiii/n_out
   f_sp_out(:,6) = nheiTR/n_out
  
   ! Adimensional number density profile
	n_out = n_out/n0

   ! Adimensional heating and cooling rates
   heat_out = heat/q0
   cool_out = cool/q0
      
   ! Adjust value of pressure boundary condition
   dp_bc = (nhii(1-Ng) + nheii(1-Ng) + 2.0*nheiii(1-Ng))/n0

	! End of subroutine 
	end subroutine ioniz_eq
	
	! End of module
	end module ionization_equilibrium         
