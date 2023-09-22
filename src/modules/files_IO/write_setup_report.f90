	module setup_report
	! Contains subroutine to write the setup of the current 
	!	simulation to file
	
	use global_parameters
	
	implicit none
	contains
	
	! ------------------------------ !
	
	subroutine write_setup_report
	! Write summary of the current simulation setup to file

	write(*,*) '(write_setup_report.f90) Writing the setup report on ATES.out..'
	
	write(outfile,*) '######## Simulation for ', p_name, ' ########'
	write(outfile,*) ' '
	write(outfile,*) ' ----- Planetary parameters ----- '
	write(outfile,*) ' '
	write(outfile,1)  &
      ' - Planet mass: ', Mp/MJ, ' [M_J], ', Mp/M_earth, ' [M_earth]'
	write(outfile,2)  &
      ' - Planet radius: ', R0/RJ, ' [R_J], ', R0/R_earth, ' [R_earth]'
	write(outfile,3) ' - Orbital distance: ', a_orb/AU, ' [AU]'
	write(outfile,4) ' - Equilibrium temperature: ', T0, ' [K]'
	write(outfile,5) ' - Jeans parameter (beta_0): ', b0
	write(outfile,15) & 
      ' - Log surface gravitational potential: ', log10(Gc*Mp/R0), ' erg/g'	 
	write(outfile,*) ' '
	write(outfile,*) ' ----- Stellar parameters ----- '
	write(outfile,*) ' '
	write(outfile,6) ' - Star mass: ', Mstar/Msun, ' [M_sun]'
	write(outfile,7) ' - Log EUV luminosity: ', LEUV, ' [erg/s]'
	write(outfile,8) ' - Log X-ray luminosity: ', LX, ' [erg/s]'
	write(outfile,9)  &
      ' - Log XUV flux at the planet distance: ', log10(J_XUV), ' [erg/(cm^2 s)]'
	write(outfile,*)	
	write(outfile,*) ' ----- Simulation setup parameters -----'
	write(outfile,*)
	write(outfile,10) &
      ' - Upper boundary of the domain: ', r_max, ' [R_p]'
	if (thereis_He) then
		write(outfile,11) & 
         ' - Simulating a H/He atmosphere with',' He/H ratio of ', HeH
		if (thereis_HeITR) 	&
			write(outfile,*) '- Including helium triplet chemistry'
	else 
		write(outfile,*) '- Simulating a pure H atmosphere'
	endif
	if (do_read_sed) &
		write(outfile,*) & 
         ' - Spectrum read from external file: ', sed_file
	if (is_PL_sed) &
      	write(outfile,12)  & 
            ' - Using power-law spectrum ', 'with index ', PLind
	if (is_monochr) & 
		write(outfile,13)  &
         ' - Using monochromatic radiation with energy ', e_low
	if (appx_mth.eq.'alpha') then 
		write(outfile,14) ' - 2D approximation used: alpha =, with alpha = ',a_tau
	else
		write(outfile,*) '- 2D approximation used: ', appx_mth
	endif
	write(outfile,*) 
	write(outfile,*) '----- Numerical parameters -----'	
	write(outfile,*)
	write(outfile,*) '- Grid type: ', grid_type
	write(outfile,*) '- Numerical flux: ', flux
	write(outfile,*) '- Reconstruction method: ', rec_method
	write(outfile,*) 
	if (.not.do_load_IC) &
		write(outfile,*) '----- Starting a new simulation ----- ' 
	if (do_load_IC) 		&
		write(outfile,*) '----- Continuing existing simulation ----- '
	if (do_only_pp)		&
		write(outfile,*) '- Evaluating post processing only'
	if (force_start) &
		write(outfile,*) '- Forcing the simulation for the first 1000 steps'
	   
1	format(A16,F5.3,A8,F8.5,A10)
2	format(A18,F5.3,A8,F8.5,A10)	
3	format(A21,F6.4,A5)	
4	format(A28,F6.1,A4)	
5	format(A29,F7.2)
6	format(A14,F5.3,A8)	
7	format(A23,F6.3,A8)	
8	format(A25,F6.3,A8)	
9	format(A40,F5.3,A15)	
10	format(A33,F5.2,A6)	
11	format(A36,A15,F8.6)		
12	format(A28,A11,F5.2)	
13 format(A45,F7.2)	
14	format(A48,F6.3)		
15 format(A40,F5.2,A6)

	write(*,*) '(write_setup_report.f90) Done.'

	end subroutine write_setup_report
	
	
	
	! End of module
	end module setup_report
