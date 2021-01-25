	module HeI_cross_section
	! Energy-dependent photoionization cross section for HeI
	!     (1E-18 cm^{-2} units)
	implicit none
	
	contains

      ! Energy-dependent photoionization cross section for HeI
      double precision function sigma_HeI(E)
	real*8,intent(in) :: E
	real*8 :: eth

      ! Ionization treshold
	eth = 24.6*0.999

      if (E.ge.eth) then      ! if E > E_th
            sigma_HeI = 0.6935/((E*1.0d-2)**1.82+(E*1.0d-2)**3.23)
      else
	      sigma_Hei = 0.0
	endif
	
      end function sigma_HeI
      
      ! End of module
      end module HeI_cross_section
