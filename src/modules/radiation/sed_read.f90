	module sed_reader
	! Module containing the subroutine for reading the numerical SED
		
	use global_parameters

	implicit none

	contains

	!---------------------------------------------------!
	
	subroutine read_sed
	! Subroutine to read the numerical SED from file.
	!	The file has to be formatted in two columns:
	!	   1) bin central wavelength [Angstrom] 
	!	   2) flux at the planet surface [erg cm^{-2} s^{-2} A^{-1}]
	!	Data have to be order with increasing wavelength
	
	integer :: io
	integer :: j,skip
	real*8 :: e_top_read
	real*8 :: dum_w,dum_f,dum_e
	real*8 :: LEUV_int,LX_int,Df_int
	real*8, dimension(:), allocatable :: wave_c 

   write(*,*) '(sed_read.f90) Reading the numerical spetrum..'
   
	! Open file to read number of lines to be skipped
	!	according to the selected energy interval
	open(unit = 1, file = sed_file, status = 'old')
	
		! Extend to 4.8 eV if He triplet is included
		if (thereis_HeITR) e_low = e_th_HeTR 

		! Set highest energy in SED
		e_top_read = e_top
		if (.not.thereis_Xray) e_top_read = e_mid 

		! Initialize counters
		Nl   = 0
		skip = 0
		
		! Read through file
		do 
			! Read wavelength and flux
			read(1,*,iostat = io) dum_w,dum_f
			
			! Convert wavelength to eV
			dum_e = hp_eV*c_light/(dum_w*1e-8)
			
			! Add one line to skip if current energy > e_top_read
			!	and keep reading
			if (dum_e .gt. e_top_read) then
				skip = skip + 1
				cycle
			endif
			
			! if current energy < e_low quit loop
			if (dum_e .lt. e_low .or. io .lt. 0) exit
			
			! Update number of selected points
			Nl = Nl + 1 	
			
		enddo
		
		! Return to the beginning of the file
		rewind (unit = 1)

		! Allocate energy vectors
		allocate(wave_c(Nl))
		allocate(F_XUV(Nl))
		allocate(e_v(Nl),de_v(Nl))

		! Read lines to be skipped
		do j = 1,skip; read(1,*) dum_w,dum_f; enddo

		! Read lines of desired energy interval
		do j = 1,Nl
			read(1,*) wave_c(j),F_XUV(j)
		enddo

	! Close SED file
	close(1)
	write(*,*) ' (sed_read.f90) Done.'

	! Convert wavelength to energy
	e_v   = hp_eV*c_light/(wave_c*1e-8)
	
	! Calculate energy bins width 
	de_v(1)      = -0.5*(e_v(2) - e_v(1))
	de_v(2:Nl-1) = -0.5*(e_v(3:Nl) - e_v(1:Nl-2))
	de_v(Nl)     = -0.5*(e_v(Nl) - e_v(Nl-1))

	! Rescale from F_lambda to F_E
	F_XUV = F_XUV*hp_eV*c_light*1.0e8/e_v**2.0e0
	
	! Calculate LEUV and LX luminosities
	LEUV_int  = 0.0
	LX_int    =	0.0	
	do j = 1,Nl
		
		! Skip if energy is below ionization threshold
		if (e_v(j) .lt. e_th_HI) exit
		
		! Integrate F(E) w.r.t. energy
		Df_int = de_v(j)*F_XUV(j)*4.0e0*pi*a_orb*a_orb
					
		! Update separately EUV and X-ray luminosities
		if (e_v(j) .lt. e_mid) then
			LEUV_int = LEUV_int + Df_int
			else
			LX_int = LX_int + Df_int
		endif
		
	enddo
	
	! Convert to log10 
	LX   = log10(LX_int)
	LEUV = log10(LEUV_int)
	
	deallocate(wave_c)

	! End of subroutine
	end subroutine read_sed
	
	! End of module
	end module sed_reader	
