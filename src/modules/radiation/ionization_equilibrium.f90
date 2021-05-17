	module ionization_equilibrium
	! Evaluate the ionization structure and the heating and cooling functions for a given temperature 
	
	use global_parameters
      use Cooling_Coefficients      ! Various functions for cooling coefficients	
	use System_HeH                ! Equilibrium equations
	use System_H
	use omp_lib                   ! OMP libraries
	
	implicit none

	   
	contains 
	
	subroutine ioniz_eq(T_in,n_in,f_sp_in,n_out, &
      	 		  f_sp_out,heat_out,cool_out,q)
      	 		  
                             
	integer :: N_eq    ! Numbers of equations
	integer :: lwa     ! Working array length
	integer i,j,k,info

	real*8, dimension(1-Ng:N+Ng), intent(in)   :: T_in,n_in
	real*8, dimension(1-Ng:N+Ng,5), intent(in) :: f_sp_in
	 
	real*8, dimension(1-Ng:N+Ng) ::  T_K      ! Dimensional temperature
	real*8, dimension(1-Ng:N+Ng) ::  nh,nhi,nhii,                   & ! Species densities
	                                 nhe,nhei,nheii,nheiii,         &	
	                                 nhii_sol,nhei_sol,nheii_sol,   &
	                                 ne,n_in_dim
	
	! Column densities                                 
      real*8, dimension(1-Ng:N+Ng) ::  N1,N15,N2               
      
      ! Photo ionization rates
      real*8, dimension(1-Ng:N+Ng) ::  P_HI,P_HeI,P_HeII 
                       	
      ! Heating, cooling
      real*8, dimension(1-Ng:N+Ng) ::  heat,cool    
                 
      ! Recombination coefficients
      real*8, dimension(1-Ng:N+Ng) ::  rchiiB,rcheiiB,rcheiiiB  
              
      ! Ionization coefficients
      real*8, dimension(1-Ng:N+Ng) ::  a_ion_HI,a_ion_HeI,a_ion_HeII    
	      
	      
      real*8 :: TT                              ! Temperature component
      real*8 :: PIR_1,PIR_15,PIR_2              ! Photoionization rates
      real*8 :: Hea_1 		          	      ! Heating rates
      real*8 :: brem,coex,coio,reco             ! Cooling rates
      real*8 :: deltal                          ! Optical depth
      real*8 :: dr	                        ! Grid spacing
      real*8 :: iup_H,ilo_H                     ! Heating rate integral variables         
      real*8 :: iup_1,ilo_1,        &           ! Photoheating integral variables
                iup_15,ilo_15,      &
                iup_2,ilo_2,        &
                iup_f,ilo_f 
      real*8 :: ilo_q,iup_q,q_abs               ! Absorbed energy
      real*8 :: elo,eup                         ! Energy parameters
      real*8 :: e_th_HI,e_th_HeI,e_th_HeII      ! Treshold energies
      real*8 :: tol,dpmpar                      ! Equilibrium system setup

      real*8, dimension(:), allocatable :: sys_sol, sys_x
      real*8, dimension(:), allocatable :: wa
      real*8, dimension(14) :: params
      
      
      ! Output density
      real*8, dimension(1-Ng:N+Ng),intent(out) :: n_out 
      
      ! Output heating,cooling and absorbed energy 
      real*8, dimension(1-Ng:N+Ng),intent(out) :: heat_out,cool_out,q 
      
      ! Output species fractions          
      real*8, dimension(1-Ng:N+Ng,5),intent(out) :: f_sp_out

      
      !----------------------------------------------------------!
      
      ! Allocate variables according to composition
      if (HeH.gt.(0.0)) then   
              
            N_eq = 3
            
      else
            
            N_eq = 1
            
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
	n_in_dim = n_in*n0
	T_K      = T_in*T0
		
	! Exctract species profiles
	nhi    = f_sp_in(:,1)*n_in_dim    ! HI
	nhii   = f_sp_in(:,2)*n_in_dim    ! HII
	nhei   = f_sp_in(:,3)*n_in_dim    ! HeI
	nheii  = f_sp_in(:,4)*n_in_dim    ! HeII
	nheiii = f_sp_in(:,5)*n_in_dim    ! HeIII
	
	if (HeH.eq.(0.0)) then
	      
	      nhei   = 0.0
            nheii  = 0.0
            nheiii = 0.0

      endif
      
	! Free electron density (assuming overall neutrality)
	nh  = nhi  + nhii
	nhe = nhei + nheii + nheiii
	ne  = nhii + nheii + 2.0*nheiii
	
	!----------------------------------!
		
	N1(N+Ng)  = dr_j(N+Ng)*R0*nhi(N+Ng)
	N15(N+Ng) = dr_j(N+Ng)*R0*nhei(N+Ng)
	N2(N+Ng)  = dr_j(N+Ng)*R0*nheii(N+Ng)
	
	
      do j = N+Ng-1,1-Ng,-1	
	
	      ! Spacing
	      dr = dr_j(j)*R0
	      
	      ! Evaluate new column densities by integration
	      N1(j)  = N1(j+1)  + nhi(j)*dr         ! HI
	      N15(j) = N15(j+1) + nhei(j)*dr	  ! HeI 
	      N2(j)  = N2(j+1)  + nheii(j)*dr       ! HeII 

	enddo 
	

      !----------------------------------!
      
      !---- Photoionization and photoheating ----!
	
	!$OMP PARALLEL DO & 
	!$OMP SHARED ( P_HI,P_HeI,P_HeII,heat,q )       &
	!$OMP PRIVATE ( Hea_1,PIR_1,PIR_15,PIR_2,       &
	!$OMP           ilo_1,ilo_15,ilo_2,ilo_f,ilo_H, &
	!$OMP           iup_1,iup_15,iup_2,iup_f,iup_H, &
	!$OMP           ilo_q,iup_q,q_abs,              &
	!$OMP           elo,eup,deltal,i,j) 
      
    
	do j = 1-Ng,N+Ng   

	! Initialization of integrands
	Hea_1  = 0.0
      PIR_1  = 0.0
      PIR_15 = 0.0
      PIR_2  = 0.0
      q_abs  = 0.0
           
	elo    = e_v(1)		!Initial lower boundary of integration
	
	! Initial optical depth 
      deltal = (s_hi(1)*N1(j)                               &
              + s_hei(1)*N15(j)                             &       
              + s_heii(1)*N2(j))*1.0e-18	

      ilo_f =  F_XUV(1)*exp(-deltal)/(1.0+a_tau*deltal)


	ilo_H  = ilo_f*(                                          &
	         nhi(j)*(1.0-e_th_HI/elo)*s_hi(1)    +      &
	         nhei(j)*(1.0-e_th_HeI/elo)*s_hei(1) +      &
	         nheii(j)*(1.0-e_th_HeII/elo)*s_heii(1))
      	               
	ilo_1  =  ilo_f*s_hi(1)/elo
	ilo_15 =  ilo_f*s_hei(1)/elo
	ilo_2  =  ilo_f*s_heii(1)/elo
	         	            
	ilo_q  =  ilo_f*(nhi(j)*s_hi(1)   +    &
	                 nhei(j)*s_hei(1) +    &
	                 nheii(j)*s_heii(1))
	                            	            
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
     		        	
     		   iup_q  =  iup_f*(nhi(j)*s_hi(i)   +    &
	                          nhei(j)*s_hei(i) +    &
	                          nheii(j)*s_heii(i))
	          
	                               		     
     		   ! Upper integrand for photoheating rate 
     	         iup_H = iup_f*(                                    &
	                 nhi(j)*(1.0-e_th_HI/eup)*s_hi(i) +         &
	                 nhei(j)*(1.0-e_th_HeI/eup)*s_hei(i) +      &
	                 nheii(j)*(1.0-e_th_HeII/eup)*s_heii(i))
               
               
               ! Evaluate photoheating integral (eV/s)
		   PIR_1  = PIR_1                                     &
      	            +0.5*(iup_1+ilo_1)*(eup-elo) 
      	   PIR_15 = PIR_15                                    &
      	            +0.5*(iup_15+ilo_15)*(eup-elo)  
      	   PIR_2  = PIR_2                                     &
      	            +0.5*(iup_2+ilo_2)*(eup-elo)          
      	                   
               ! Evaluate photoheating integral (eV/s)
		   Hea_1  = Hea_1                                     &
      	            +0.5*(iup_H+ilo_H)*(eup-elo) 
               
               ! Evaluate total absorbed energy 
               q_abs = q_abs + 0.5*(ilo_q+iup_q)*(eup-elo)
                            
     		   ! Renew lower boundary and integrand
               elo = eup             
               ilo_1 = iup_1
               ilo_15 = iup_15
               ilo_2 = iup_2
               ilo_H = iup_H         
               ilo_q = iup_q
               
               enddo


      ! Multiply for the dimensional coefficient 
      !$OMP CRITICAL
      heat(j)   = Hea_1*1.0e-18	
      P_HI(j)   = PIR_1*1.0e-18*erg2eV	
      P_HeI(j)  = PIR_15*1.0e-18*erg2eV
	P_HeII(j) = PIR_2*1.0e-18*erg2eV
	q(j)      = Hea_1/q_abs
	!$OMP END CRITICAL
	
	enddo
	!$OMP END PARALLEL DO

      !----------------------------------!
      
      !---- Cooling rates ----!
	
	!$OMP PARALLEL DO                               &
	!$OMP SHARED ( rchiiB,rcheiiB,rcheiiiB,         &
	!$OMP          a_ion_HI,a_ion_HeI,a_ion_HeII,   &
	!$OMP          cool )                           &
	!$OMP PRIVATE ( TT,j,brem,coex,reco,coio )
 	
      do j = 1-Ng,N+Ng
	
	      ! Temperature substitution
	      TT = T_K(j)
	      
	      !-- Recombination --!
	      
	      ! Coefficients 	
	      rchiiB(j)   = rec_HII_B(TT)     ! HII
            rcheiiB(j)  = rec_HeII_B(TT)    ! HeII
	      rcheiiiB(j) = rec_HeIII_B(TT)   ! HeIII
            
            
            ! Cooling rate
            reco  =  rec_cool_HII(TT)*nhii(j)     & ! HII
                   + rec_cool_HeII(TT)*nheii(j)   & ! HeII
                   + rec_cool_HeIII(TT)*nheiii(j)   ! HeIII


            !-- Collisional ionization --!
            
            ! Rate coefficients
            a_ion_HI(j)   = ion_coeff_HI(TT)      ! HI
            a_ion_HeI(j)  = ion_coeff_HeI(TT)     ! HeI
            a_ion_HeII(j) = ion_coeff_HeII(TT)    ! HeII
            
            
            ! Cooling rate
            coio =  2.179e-11*a_ion_HI(j)*nhi(j)  		   & ! HI
                  + 3.940e-11*a_ion_HeI(j)*nhei(j) 		   & ! HeI
		      + kb_erg*631515.0*a_ion_HeII(j)*nheii(j)   ! HeII

            !-- Bremsstrahlung --!
       
            ! Cooling rate
            brem = 1.426e-27*sqrt(TT)*                         &
            	 (ih**2.0*GF(TT,ih)*nhii(j) +                & ! HII
            	 ihe**2.0*GF(TT,ihe)*(nheii(j)+nheiii(j)) )    ! He
            
            !-- Collisional excitation --!
           
            ! Collisional excitation 
            coex = coex_rate_HI(TT)*nhi(j)       &    ! HI
                 + coex_rate_HeI(TT)*nhei(j)     &    ! HeI
                 + coex_rate_HeII(TT)*nheii(j)        ! HeII

       
            ! Total cooling rate in erg/(s cm^3)
            !$OMP CRITICAL
       	cool(j) = ne(j)*(brem + coex + reco + coio)
       	!$OMP END CRITICAL
 	
      enddo
 
      !$OMP END PARALLEL DO
      

      !----------------------------------!
      
      ! Ionization equilibrium system solution	
      do j = N+Ng,1-Ng,-1
      
            if (HeH.eq.(0.0)) then ! If no helium
            
	            ! Ionization equilibrium system setup
                  params(1) = P_HI(j)
                  params(2) = rchiiB(j)
                  params(3) = nh(j)
                  params(4) = a_ion_HI(j)

                   ! Initial guess
                  if(count.le.0) then            
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

                  ! Exctract solution profiles	
	            nhi(j)    = nh(j)*(1.0 - sys_x(1))
	            nhii(j)   = nh(j)*sys_x(1)
	            nhei(j)   = 0.0
	            nheii(j)  = 0.0
	            nheiii(j) = 0.0

            else  ! If there is helium

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
                  
                  ! Initial guess
                  if(count.le.0) then            
                        if(r(j).le.(1.5))then
                  	      sys_x(1) = r(j)-0.5
                  	      sys_x(2) = r(j)-0.5
                  	      sys_x(3) = r(j)-0.5
                        else
                  	      sys_x(1) = 1.0
                  	      sys_x(2) = 1.0
                  	      sys_x(3) = 1.0    
                        endif     
                  else

		            sys_x(1) = nhii(j)/nh(j)
		            sys_x(2) = nheii(j)/nhe(j)
		            sys_x(3) = nheiii(j)/nhe(j)

                  endif

             	! Call hybrd1 routine (from minpack)
                  call hybrd1(ion_system_HeH,N_eq,sys_x,sys_sol,   &
                              tol,info,wa,lwa,params) 

                  
	            ! Exctract solution profiles	
	            nhi(j)    = nh(j)*(1.0 - sys_x(1))
	            nhii(j)   = nh(j)*sys_x(1)
	            nhei(j)   = nhe(j)*(1.0 - sys_x(2) - sys_x(3))
	            nheii(j)  = nhe(j)*sys_x(2)
	            nheiii(j) = nhe(j)*sys_x(3)
	      
	      endif
	
      enddo

	
	! Density with atomic numbers
      n_out = nhi + nhii  + 4.0*(nhei + nheii + nheiii)      

       	     
      ! Abundacies profiles
      f_sp_out(:,1) = nhi/n_out
      f_sp_out(:,2) = nhii/n_out
      f_sp_out(:,3) = nhei/n_out
      f_sp_out(:,4) = nheii/n_out
      f_sp_out(:,5) = nheiii/n_out
      
  
      ! Adimensional number density profile
	n_out = n_out/n0

      ! Adimensional heating and cooling rates
      heat_out = heat/q0
      cool_out = cool/q0
      
      ! Adjust value of pressure bounday condition
      dp_bc = (nhii(1-Ng) + nheii(1-Ng) + 2.0*nheiii(1-Ng))/n0
	 
	end subroutine ioniz_eq
	
	! End of module
	end module ionization_equilibrium         
