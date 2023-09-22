   module energy_vectors_construct
   ! Construct energy grid for the computation of radiative
	!      contributions and radiative equilibrium
      
   use global_parameters
   use sed_reader
   use J_incident
   use Cross_sections 
	
   implicit none
      
   contains
      
   subroutine set_energy_vectors
	! Subroutine to construct the energy grid

   integer :: i,j
	integer, parameter :: num_HI = 50
	integer, parameter :: num_HeI = 50
	integer, parameter :: num_HeII = 50
	integer, parameter :: num_X = 50
	
	real*8 :: e_min, e_max
	real*8, dimension(:), allocatable :: dum_E,dum_F,dum_dE
	real*8 :: temp1,temp2
	
	! Set the energy band of the spectrum
	if (is_PL_sed) then 
		
		! Add points for He triplet if included
		NlTR = 0
		if (thereis_HeITR) NlTR = 20
		Nl = Nl_fix + NlTR
		
		! Allocate vectors
		allocate(e_v(Nl),de_v(Nl))
		allocate(s_hi(Nl), s_hei(Nl),s_heii(Nl))
		allocate(F_XUV(Nl))
				
		! --- Construct energy grid --- ! 
				
		! Points between 4.8 eV and 13.6 eV if HeI triplet is included
		if (thereis_HeITR) then
			do j = 1,NlTR
				e_v(j) = e_th_HeTR*	&
					  (e_th_HI/e_th_HeTR)**((j-1.0)/(NlTR*1.0))
		  	enddo		
		endif

		! Pure-HI grid [13.6 eV, 24.6 eV]
		e_min = e_th_HI
		e_max = e_th_HeI
		do j = 1,num_HI
			e_v(NlTR + j) = 	&
				e_min*(e_max/e_min)**((j-1.0)/(num_HI*1.0)) 
		enddo 
		 
		! HI+HeI grid [24.6 eV, 54.4 eV]
		e_min = e_th_HeI
		e_max = e_th_HeII
		do j = 1,num_HeI
			e_v(NlTR + num_HI + j) = 	&
				e_min*(e_max/e_min)**((j-1.0)/(num_HeI*1.0)) 
		enddo 
		 
		! HI+HeI+HeII grid	[54.4 eV,124 eV] 
		e_min = e_th_HeII
		e_max = e_mid
		do j = 1,num_HeII
			e_v(NlTR + num_HI + num_HeI + j) = 	&
				e_min*(e_max/e_min)**((j-1.0)/(num_HeII*1.0)) 
		enddo 
		 
		! X-ray grid [124 eV, e_XUV_top]	 
		e_min = e_mid
		e_max = e_top
		do j = 1,num_X
			e_v(NlTR + num_HI + num_HeI + num_HeII + j) = 	&
				e_min*(e_max/e_min)**((j-1.0)/(num_X*1.0))
		enddo 

		! Construct bin width
		de_v(1)      = 0.5*(e_v(2) - e_v(1))
		de_v(2:Nl-1) = 0.5*(e_v(3:Nl) - e_v(1:Nl-2))
		de_v(Nl)     = 0.5*(e_v(Nl) - e_v(Nl-1))
	
	else if (is_monochr) then	! If monochromatic radiation
		
		! Set only one wavelength
		Nl = 1
		
		! Allocate vectors
		allocate(e_v(Nl),de_v(Nl))
		allocate(s_hi(Nl), s_hei(Nl),s_heii(Nl))
		allocate(F_XUV(Nl))
		
		! Define the only energy value according to the input		
		e_v(1)  = e_low
		de_v(1) = 1.0
		
		! Define the total flux at this energy
		F_XUV(1) = 10.0**LEUV/(4.0*pi*a_orb**2.0)
	
	else if (do_read_sed) then    ! If read sed
	 
		! Read the SED file
	 	call read_sed
 	
		! Flip energy vectors
		allocate(dum_E(Nl),dum_F(Nl),dum_dE(Nl))
		
		dum_E  = e_v
		dum_F  = F_XUV
		dum_dE = de_v
		
		do i = 1,Nl
			e_v(i)   = dum_E(Nl-i+1)
			F_XUV(i) = dum_F(Nl-i+1)
			de_v(i)  = dum_dE(Nl-i+1)
		enddo
	
	 	deallocate(dum_E,dum_F,dum_dE)
	endif
	

	! Calculate the integrated value of the flux
	J_XUV  = (10.0**LX + 10.0**LEUV)/(4.0*pi*a_orb**2.0)
      	
	! HI photoioiniz. cross section 
	s_hi = (/ (sigma(e_v(i),ih), i = 1,Nl) /)
	
	! HeI photoioiniz. cross section                         
	s_hei = (/ (sigma_hei(e_v(i)), i = 1,Nl) /)
	
	! HeII photoioiniz. cross section                       
	s_heii = (/ (sigma(e_v(i),ihe), i = 1,Nl) /)
	
	! HeI triplet photoioiniz. cross section                       
	s_heiTR = (/ (sigma_HeI23S(e_v(i)), i = 1,Nl) /)
	
	! --------------------------------------------------

	! Incident flux    
	if (is_PL_sed) then 

		if (appx_mth.eq.'Rate/2 + Mdot/2') then

			F_XUV  = (/ (0.5e0*J_inc(e_v(i)), i = 1,Nl) /)

		elseif (appx_mth.eq.'Rate/4 + Mdot') then

			F_XUV  = (/ (0.25e0*J_inc(e_v(i)), i = 1,Nl) /)

		else

			F_XUV  = (/ (J_inc(e_v(i)), i = 1,Nl) /)

		endif      
		
	else 
	
		if (appx_mth.eq.'Rate/2 + Mdot/2') then

			F_XUV  = 0.5e0*F_XUV

		elseif (appx_mth.eq.'Rate/4 + Mdot') then

			F_XUV  = 0.25e0*F_XUV

		endif 
	
	endif      
	
	! End of subroutine 
    end subroutine set_energy_vectors
      
    ! End of module
    end module energy_vectors_construct

