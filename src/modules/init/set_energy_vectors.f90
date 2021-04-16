      module energy_vectors_construct
      ! Construct energy grid for the computation of radiative
      !      contributions and radiative equilibrium
      
      use global_parameters
      use J_incident
      use Cross_sections 
	
      implicit none
      
      contains
      
      subroutine set_energy_vectors
      integer :: i

	! Energy grid
	e_v = (/ ( e_low*(e_top/e_low)**((i-1.0)/(Nl-1.0)), i =1,Nl) /) 
	 
	! HI photoioiniz. cross section 
	s_hi = (/ (sigma(e_v(i),ih), i = 1,Nl) /)
	
	! HeI photoioiniz. cross section                         
	s_hei = (/ (sigma_hei(e_v(i)), i = 1,Nl) /)
	
	! HeII photoioiniz. cross section                       
	s_heii = (/ (sigma(e_v(i),ihe), i = 1,Nl) /)
	
	! Incident flux         
	if (appx_mth.eq.'Rate/2 + Mdot/2') then

      	F_XUV  = (/ (0.5*J_inc(e_v(i)), i = 1,Nl) /)
      
      elseif (appx_mth.eq.'Rate/4 + Mdot') then
      
      	F_XUV  = (/ (0.25*J_inc(e_v(i)), i = 1,Nl) /)
	
	else
	
		F_XUV  = (/ (J_inc(e_v(i)), i = 1,Nl) /)
	
	endif      
      
      
      end subroutine set_energy_vectors
      
      
      ! End of module
      end module energy_vectors_construct

