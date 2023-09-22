   module utils_ion_eq

   use global_parameters
   use utils
   use Cooling_Coefficients      ! Various functions for cooling coefficients	
   use omp_lib                   ! OMP libraries
	
	! Move here photoionization and photoheating
	! Two subroutines for H and He+H
	
   implicit none
	
	contains
	
	! ------------------------------------------------------------- !

	subroutine PH_heat_H(nhi,P_HI,heat,q)
	! Computes photoionization rates and heating rates for
	!	an atmosphere composed of H and He
	
	integer :: i,j
	
	real*8, dimension(1-Ng:N+Ng),intent(in) :: nhi    

	! Dummy zero 
	real*8, dimension(1-Ng:N+Ng), parameter :: nhei   = 0.0, nheii  = 0.0
	real*8, dimension(1-Ng:N+Ng), parameter :: nheiii = 0.0, nheiTR = 0.0
	real*8, dimension(1-Ng:N+Ng) :: N15, N2, NTR

   real*8 :: dr	                ! Grid spacing
   real*8 :: PIR_1		            ! Photoionization rates
   real*8 :: Hea_1 		        ! Heating rates  
	real*8 :: q_abs           	! Absorbed energy
	
	! Integral variables
	real*8, dimension(Nl) :: tauE			
	real*8, dimension(Nl) :: int_f,int_1,int_q,int_H

	! Column densities                                 
   real*8, dimension(1-Ng:N+Ng) ::  N1 
	
   ! Photo ionization rates
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  P_HI 
      
   ! Heating efficiency 
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  q
      
   ! Heating rate
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  heat
      
	!----------------------------------!
	
	! Evaluate the column density
	call calc_column_dens(nhi,nhei,nheii,nheiTR,N1,N15,N2,NTR)

   ! Evaluate photoionization rates and photoheating rates
	do j = 1-Ng,N+Ng   

		! Initialization of integrands
		Hea_1  = 0.0
		PIR_1  = 0.0
		q_abs  = 0.0

		tauE = s_hi*N1(j)*1.0e-18	

		! Initial integrands
		int_f = F_XUV*exp(-tauE)/(1.0 + a_tau*tauE)
		int_H = int_f*(1.0-e_th_HI/e_v)*s_hi*nhi(j)
		int_1 = int_f*s_hi/e_v
		int_q = int_f*s_hi*nhi(j)

		! Value of integrals
		Hea_1 = sum(int_H*de_v)	
		PIR_1 = sum(int_1*de_v)	
		q_abs = sum(int_q*de_v)
		
		! Multiply for the dimensional coefficient 
		heat(j)   = Hea_1*1.0e-18	
		P_HI(j)   = PIR_1*1.0e-18*erg2eV	
		q(j)      = Hea_1/q_abs
	
	enddo
	
	end subroutine PH_heat_H

	! ------------------------------------------------------------- !
	
	subroutine PH_heat_HHe(nhi,nhei,nheii,nheiTR,	&
				     P_HI,P_HeI,P_HeII,P_HeITR,heat,q)
	! Computes photoionization rates and heating rates for
	!	an atmosphere composed of H and He
	
	integer :: i,j
	
	real*8, dimension(1-Ng:N+Ng),intent(in) :: nhi,nhei,nheii
	real*8, dimension(1-Ng:N+Ng),intent(in) :: nheiTR
	real*8, dimension(1-Ng:N+Ng) :: nheiS
   real*8, dimension(1-Ng:N+Ng) ::  N1,N15,N2,NTR
	
   real*8 :: PIR_1,PIR_15,PIR_2,PIR_TR     ! Photoionization rates
   real*8 :: Hea_1 		          	    ! Heating rates  
   real*8 :: q_abs                		    ! Absorbed energy
	
	! Integral variables
	real*8, dimension(Nl) :: tauE			
	real*8, dimension(Nl) :: int_f,int_1,int_15,int_2,int_TR
	real*8, dimension(Nl) :: int_q,int_H
	
   ! Photo ionization rates
	real*8, dimension(1-Ng:N+Ng), intent(out) ::  P_HI
   real*8, dimension(1-Ng:N+Ng), intent(out) ::  P_HeI,P_HeII,P_HeITR 
      
   ! Heating efficiency 
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  q
      
   ! Heating rate
	real*8, dimension(1-Ng:N+Ng),intent(out) ::  heat
      
	!----------------------------------!
	
	! Use nheiS as variable
	nheiS = nhei
	if (thereis_HeITR) nheiS = nhei - nheiTR

	! Evaluate the column density
	call calc_column_dens(nhi,nheiS,nheii,nheiTR,N1,N15,N2,NTR)
      
    !----------------------------------!
	!$OMP PARALLEL DO & 
	!$OMP SHARED ( P_HI,P_HeI,P_HeII,P_HeITR,heat,q )       	&
	!$OMP PRIVATE ( Hea_1,PIR_1,PIR_15,PIR_2,PIR_TR,        	&
	!$OMP           int_1,int_15,int_2,int_TR,int_f,int_H, 	&
	!$OMP           int_q,q_abs,tauE,j) 

	do j = 1-Ng,N+Ng 

		Hea_1  = 0.0
      PIR_1  = 0.0
      PIR_15 = 0.0
      PIR_2  = 0.0
      PIR_TR = 0.0
      q_abs  = 0.0

		! Calculate optical depth
		tauE = (s_hi*N1(j) + s_hei*N15(j) + s_heii*N2(j))*1.0e-18
		if (thereis_HeITR) tauE = tauE + s_heiTR*NTR(j)*1.0e-18

		! Calculate photoionization integrals
		int_f  = F_XUV*exp(-tauE)/(1.0 + a_tau*tauE)
		int_1  = int_f*s_hi/e_v
		int_15 = int_f*s_hei/e_v
		int_2  = int_f*s_heii/e_v
		if (thereis_HeITR) int_TR =  int_f*s_heiTR/e_v
		
		! Photoheating integral
		int_H  =  int_f*(                                 &
		         (1.0-e_th_HI/e_v)*s_hi*nhi(j) +          &
		         (1.0-e_th_HeI/e_v)*s_hei*nheiS(j) +      &
		         (1.0-e_th_HeII/e_v)*s_heii*nheii(j))
		
		! Absorbed energy integral
		int_q  =  int_f*(s_hi*nhi(j)  +     &
						 s_hei*nheiS(j)   +     &
						 s_heii*nheii(j))

		! Midpoint rule for integrals
		PIR_1  = sum(int_1*de_v)	
		PIR_15 = sum(int_15*de_v)	
		PIR_2  = sum(int_2*de_v)
		if(thereis_HeITR) PIR_TR  = sum(int_TR*de_v)	
		Hea_1  = sum(int_H*de_v)	
		q_abs  = sum(int_q*de_v)

		!$OMP CRITICAL
		! Save into vector
    	P_HI(j)    = PIR_1*1.0e-18*erg2eV	
    	P_HeI(j)   = PIR_15*1.0e-18*erg2eV
		P_HeII(j)  = PIR_2*1.0e-18*erg2eV
		P_HeITR(j) = PIR_TR*1.0e-18*erg2eV
		heat(j)    = Hea_1*1.0e-18	
		q(j)       = Hea_1/q_abs
		!$OMP END CRITICAL

	enddo
	!$OMP END PARALLEL DO

	! End of subroutine 
	end subroutine PH_heat_HHe
	
	!----------------------------------!
	
	subroutine eval_cool(T_K,nhi,nhii,nhei,nheii,nheiii,  	&
				   rchiiB,rcheiiB,rcheiiiB, 			&
				   a_ion_HI,a_ion_HeI,a_ion_HeII,cool)

	! Evaluate the cooling rate contributions to energy and 
	!	rate equations
	
	integer :: j
	
	real*8, dimension(1-Ng:N+Ng),intent(in)  :: nhi,nhii,          & 
	                                 		    nhei,nheii,nheiii
	                                 
	! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng),intent(in) ::  T_K      

   real*8, dimension(1-Ng:N+Ng) :: brem,coex,coio,reco  ! Cooling rates
	real*8, dimension(1-Ng:N+Ng) :: ne		  			 ! Electron number density
	real*8, dimension(1-Ng:N+Ng) :: GF_H,GF_He			 ! Gaunt factors
	
   ! Recombination rate coefficients
   real*8, dimension(1-Ng:N+Ng),intent(out) :: rchiiB,	 &
												rcheiiB, &
												rcheiiiB 
	
	! Recombination cooling coefficients
	real*8, dimension(1-Ng:N+Ng) :: coeff_rec_cool_HII,  &
									coeff_rec_cool_HeII, &
									coeff_rec_cool_HeIII
	
   ! Ionization coefficients
   real*8, dimension(1-Ng:N+Ng),intent(out) ::  a_ion_HI,	&
      							   				 a_ion_HeI, &
												 a_ion_HeII  

	real*8, dimension(1-Ng:N+Ng) :: coeff_coex_rate_HI,    &
	 								coeff_coex_rate_HeI,   &
									coeff_coex_rate_HeII
									 
	! Heating, cooling
	real*8, dimension(1-Ng:N+Ng),intent(out) ::  cool 
									 
   ! Free electron density
   call calc_ne(nhii,nheii,nheiii,ne)	
	

	!-- Recombination --!

	! Rate coefficients
	call rec_HII_B(T_K,rchiiB)      ! HII
	call rec_HeII_B(T_K,rcheiiB)    ! HeII
	call rec_HeIII_B(T_K,rcheiiiB)  ! HeIII
	
	! Cooling rate coefficients
	call rec_cool_HII(T_K,coeff_rec_cool_HII)
	call rec_cool_HeII(T_K,coeff_rec_cool_HeII)
	call rec_cool_HeIII(T_K,coeff_rec_cool_HeIII)

	! Cooling rate
	reco  = coeff_rec_cool_HII*nhii     & ! HII
		   + coeff_rec_cool_HeII*nheii   & ! HeII
		   + coeff_rec_cool_HeIII*nheiii   ! HeIII

	!-- Collisional ionization --!
	
	! Rate coefficients
	call ion_coeff_HI(T_K,a_ion_HI)      ! HI
	call ion_coeff_HeI(T_K,a_ion_HeI)    ! HeI
	call ion_coeff_HeII(T_K,a_ion_HeII)  ! HeII
	
	! Cooling rate
	coio =  2.179e-11*a_ion_HI*nhi  	 & ! HI
		  + 3.940e-11*a_ion_HeI*nhei 	 & ! HeI
		  + 8.715e-11*a_ion_HeII*nheii   ! HeII
	
	!-- Bremsstrahlung --!
	
	! Gaunt factors
	call GF(T_K,ih,GF_H)
	call GF(T_K,ihe,GF_He)

	! Cooling rate
	brem = 1.426e-27*sqrt(T_K)*             &
		   (ih**2.0*GF_H*nhii +              & ! HII
			ihe**2.0*GF_He*(nheii + nheiii))    ! He

	!-- Collisional excitation --!

	! Rate coefficients
	call coex_rate_HI(T_K,coeff_coex_rate_HI) 		! HI
	call coex_rate_HeI(T_K,coeff_coex_rate_HeI)   	! HeI
	call coex_rate_HeII(T_K,coeff_coex_rate_HeII)  	! HeII

	! Cooling rate 
	coex = coeff_coex_rate_HI*nhi       &    ! HI
		  + coeff_coex_rate_HeI*nhei     &    ! HeI
		  + coeff_coex_rate_HeII*nheii        ! HeII
	
	! Total cooling rate
	cool = ne*(brem + coex + reco + coio)

	! End of subroutine
	end subroutine eval_cool
	
	!----------------------------------!
	
	! Various coefficients for HeI triplet chemistry
	subroutine HeITR_coeffs(T_K,rcheiTR,rcheii,A31,q13,q31a,q31b,Q31)
	
	! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng),intent(in) ::  T_K 
	
	real*8, intent(out) :: A31,Q31	
	real*8, dimension(1-Ng:N+Ng), intent(out) :: rcheiTR,rcheii,   &
								   q13,q31a,q31b
								   
	call rec_HeII_23S(T_K,rcheiTR)
	call rec_HeII_11S(T_K,rcheii)
	call coex_HeI_1S_23S(T_K,q13)
	call coex_HeI_23S_21S(T_K,q31a)
	call coex_HeI_23S_21P(T_K,q31b)
	A31 = 1.272e-4
	Q31 = 5.00e-10

   ! End of subroutine
	end subroutine HeITR_coeffs
	
	! End of module
	end module utils_ion_eq
