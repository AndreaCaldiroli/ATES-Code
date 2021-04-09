	module post_processing
	! Subroutine to correct the output ionization profiles
	!	taking into account the ionization term
	
	use global_parameters
	use System_implicit_adv_HeH	
	use System_implicit_adv_H
	use Cooling_Coefficients
	use output_write
	
	implicit none
	
	
	contains
	
	subroutine post_process_adv(rho,v,p,T_in,heat,cool,eta,    &
                                  nhi_in,nhii_in,		     &
                                  nhei_in,nheii_in,nheiii_in)
                                  
                                  
	real*8, dimension(1-Ng:N+Ng), intent(in) :: rho,v,p,T_in
      real*8, dimension(1-Ng:N+Ng), intent(in) :: heat,cool
      real*8, dimension(1-Ng:N+Ng), intent(in) :: eta
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi_in,nhii_in
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei_in,nheii_in, &
      							  nheiii_in
	
	
	integer :: N_eq    ! Numbers of equations
	integer :: lwa     ! Working array length
	integer i,j,k,info
	 
	real*8, dimension(1-Ng:N+Ng) ::  T,T_K      ! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng) ::  nh,nhe,ne,n_tot 
	
	! Column densities                                 
      real*8, dimension(1-Ng:N+Ng) ::  N1,N15,N2               
      
      ! Photo ionization rates
      real*8, dimension(1-Ng:N+Ng) ::  P_HI,P_HeI,P_HeII 

      ! Recombination coefficients
      real*8, dimension(1-Ng:N+Ng) ::  rchiiB,rcheiiB,rcheiiiB  
 
 	! Updated species densities
      real*8, dimension(1-Ng:N+Ng) :: nhi,nhii
      real*8, dimension(1-Ng:N+Ng) :: nhei,nheii,nheiii
      
      real*8, dimension(1-Ng:N+Ng) :: nhi_w,nhii_w
      real*8, dimension(1-Ng:N+Ng) :: nhei_w,nheii_w,nheiii_w
      							  
	      
	      
      real*8 :: TT                              ! Temperature component
      real*8 :: PIR_1,PIR_15,PIR_2              ! Photoionization rates
      real*8 :: deltal                          ! Optical depth
      real*8 :: dr	                        ! Grid spacing
      real*8 :: iup_1,ilo_1,        &           ! Photoheating integral variables
                iup_15,ilo_15,      &
                iup_2,ilo_2,        &
                iup_f,ilo_f 
      real*8 :: elo,eup                         ! Energy parameters
      real*8 :: e_th_HI,e_th_HeI,e_th_HeII      ! Treshold energies
      real*8 :: tol,dpmpar                      ! Equilibrium system setup


	! Substitution in the ODE solution
	real*8 :: r_p,r_m
	real*8 :: As


      real*8, dimension(:), allocatable :: sys_sol, sys_x
      real*8, dimension(:), allocatable :: wa
      real*8, dimension(:), allocatable :: params

	
	!----------------------------------------------------------!
      
      ! Allocate variables according to composition
      if (HeH.gt.(0.0)) then   
              
            N_eq = 3
            allocate (params(11))
            
      else
            
            N_eq = 1
            allocate (params(5))
            
      endif
      
      lwa  = (N_eq*(3*N_eq+13))/2
      allocate (sys_sol(N_eq))
      allocate (sys_x(N_eq))
      allocate (wa(lwa))
      
      
      !----------------------------------------------------------!      
      
      ! Global parameters
      
      ! Numerical tolerance for system solution
      tol = sqrt(dpmpar(1))
		
	! Limits of integration
	e_th_HI   = 13.6	 ! Treshold for HI ionization
	e_th_HeI  = 24.6   ! Treshold for HeI ionization
	e_th_HeII = 54.4   ! Treshold for HeII ionization
      
      !----------------------------------!
	
	! Preliminary profiles exctraction
	
	! Dimensional total number density profile and temperature
	T_K = T_in*T0
	         
      ! Initialize vectors
      nhi    = nhi_in*n0
      nhii   = nhii_in*n0
      nhei   = nhei_in*n0
      nheii  = nheii_in*n0
      nheiii = nheiii_in*n0	
      
	! Free electron density (assuming overall neutrality)
	nh  = nhi  + nhii 
	nhe = nhei + nheii + nheiii 
	ne  = nhii + nheii + 2.0*nheiii 
      
      
      do j = N+Ng-1,2-Ng,-1	
	
	      ! Spacing
	      dr = 0.5*(r(j+1)-r(j-1))*R0
	      
	      ! Evaluate new column densities by integration
	      N1(j)  = N1(j+1)  + nhi(j)*dr         ! HI
	      N15(j) = N15(j+1) + nhei(j)*dr	  ! HeI 
	      N2(j)  = N2(j+1)  + nheii(j)*dr       ! HeII 

	enddo
	
	! Boundary values
	N1(1-Ng)  = N1(2-Ng)   
	N15(1-Ng) = N15(2-Ng)  
	N2(1-Ng)  = N2(2-Ng)
	
	N1(N+Ng)  = N1(N+Ng-1)   
	N15(N+Ng) = N15(N+Ng-1)  
	N2(N+Ng)  = N2(N+Ng-1)   
	

      !----------------------------------!
      
      !---- Photoionization and photoheating ----!
	
	!$OMP PARALLEL DO & 
	!$OMP SHARED ( P_HI,P_HeI,P_HeII)           &
	!$OMP PRIVATE ( PIR_1,PIR_15,PIR_2,         &
	!$OMP           ilo_1,ilo_15,ilo_2,ilo_f,   &
	!$OMP           iup_1,iup_15,iup_2,iup_f,   &
	!$OMP           elo,eup,deltal,i,j) 
      
    
	do j = 1-Ng,N+Ng   

	! Initialization of integrands
      PIR_1  = 0.0
      PIR_15 = 0.0
      PIR_2  = 0.0
           
	elo    = e_v(1)		!Initial lower boundary of integration
	
	! Initial optical depth 
      deltal = (s_hi(1)*N1(j)                               &
              + s_hei(1)*N15(j)                             &       
              + s_heii(1)*N2(j))*1.0d-18	

      ilo_f =  F_XUV(1)*exp(-deltal)/(1.0+a_tau*deltal)

      	               
	ilo_1  =  ilo_f*s_hi(1)/elo
	ilo_15 =  ilo_f*s_hei(1)/elo
	ilo_2  =  ilo_f*s_heii(1)/elo

	                            	            
               do i = 2,Nl
               
               ! Upper boundary of i-th interval of integration 
               eup    = e_v(i)	

               ! Total optical depth at heigth r
		   deltal = (s_hi(i)*N1(j)                    &
		            +s_hei(i)*N15(j)                  &
               	      +s_heii(i)*N2(j))*1.0e-18	
     		   
     		   
     		   ! Upper integrand for photoionization rate
     		   iup_f = F_XUV(i)*exp(-deltal)/(1.0+a_tau*deltal)
     		   
     		   iup_1  = iup_f*s_hi(i)/eup 
	         iup_15 = iup_f*s_hei(i)/eup      
     		   iup_2  = iup_f*s_heii(i)/eup

               
               ! Evaluate photoheating integral (eV/s)
		   PIR_1  = PIR_1                                     &
      	            +0.5*(iup_1+ilo_1)*(eup-elo) 
      	   PIR_15 = PIR_15                                    &
      	            +0.5*(iup_15+ilo_15)*(eup-elo)  
      	   PIR_2  = PIR_2                                     &
      	            +0.5*(iup_2+ilo_2)*(eup-elo)          
      	                   
                            
     		   ! Renew lower boundary and integrand
               elo = eup             
               ilo_1 = iup_1
               ilo_15 = iup_15
               ilo_2 = iup_2
               
               enddo


      ! Multiply for the dimensional coefficient 
      !$OMP CRITICAL
      P_HI(j)   = PIR_1*1.0e-18*erg2eV	
      P_HeI(j)  = PIR_15*1.0e-18*erg2eV
	P_HeII(j) = PIR_2*1.0e-18*erg2eV
	!$OMP END CRITICAL
      
      enddo
	!$OMP END PARALLEL DO
      
      
      !---- Recombination rates ----!
 	
      do j = 1-Ng,N+Ng
	
	      ! Temperature substitution
	      TT = T_K(j)
	      
	      !-- Recombination --!
	      
	      ! Coefficients 	
	      rchiiB(j)   = rec_HII_B(TT)     ! HII
            rcheiiB(j)  = rec_HeII_B(TT)    ! HeII
	      rcheiiiB(j) = rec_HeIII_B(TT)   ! HeIII
                       
      enddo
 
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
      
      do j = 2-Ng,N+Ng ! 
      
      	! Substitutions
      	dr  = (r(j) - r(j-1))*R0
		As  = dr/(v(j-1)*v0)
	      
	      if (HeH.eq.(0.0)) then
	      
	      	! Advection coeff.
			params(1)  = As 
		 	params(2)  = nhi(j-1)/nh(j-1) 
		 	params(3)  = nh(j) 
		 	params(4)  = P_HI(j)
			params(5)  = rchiiB(j)    
			
			! Initial guess of solution
			sys_x(1) = nhi(j)/nh(j)     
			
			! Call hybrd1 routine (from minpack)
		      call hybrd1(adv_implicit_H,N_eq,sys_x,sys_sol,   &
		                  tol,info,wa,lwa,params) 
		      
		      ! Exctract solution profiles	
		      nhi(j)    = sys_x(1)*nh(j)
		      nhii(j)   = (1.0 - sys_x(1))*nh(j)
		      nhei(j)   = 0.0
		      nheii(j)  = 0.0
		      nheiii(j) = 0.0			
	      	
	      else
	      
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
			
			! Initial guess of solution
			sys_x(1) = nhi(j)/nh(j) 
			sys_x(2) = nhei(j)/nhe(j)
			sys_x(3) = nheiii(j)/nhe(j)      
			
			! Call hybrd1 routine (from minpack)
		      call hybrd1(adv_implicit_HeH,N_eq,sys_x,sys_sol,   &
		                  tol,info,wa,lwa,params) 
		      
		      ! Exctract solution profiles	
		      nhi(j)    = sys_x(1)*nh(j)
		      nhii(j)   = (1.0 - sys_x(1))*nh(j)
		      nhei(j)   = sys_x(2)*nhe(j)
		      nheii(j)  = (1.0 - sys_x(2) - sys_x(3))*nhe(j)
		      nheiii(j) = sys_x(3)*nhe(j)
			
	      endif ! End if on HeH
	      
	enddo
      
      
      !----------------------------------!
      
      ! Update ionization-dependent vectors before writing
      
      ! Total number density
      n_tot = nhi + nhii + nhei + nheii + nheiii
      
      ! Electron density 
      ne  = nhii + nheii + 2.0*nheiii
      
      ! Temperature profile
      T = p*n0/(n_tot + ne)   
      
      ! Adimensionalize before writing
      nhi_w    = nhi/n0
      nhii_w   = nhii/n0
      nhei_w   = nhei/n0
      nheii_w  = nheii/n0
      nheiii_w = nheiii/n0
      
      !----------------------------------!
      
      ! Write updated thermodynamic and ionization profiles
      call write_output(rho,v,p,T,heat,cool,eta,    		        &
                        nhi_w,nhii_w,nhei_w,nheii_w,nheiii_w,'ad')
                                   
                                   	
	end subroutine post_process_adv
	
	
	! End of module
	end module post_processing
