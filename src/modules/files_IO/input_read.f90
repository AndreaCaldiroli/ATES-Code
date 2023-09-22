   module  Read_input
   ! Read input planetary parameters adn define 
      
   use global_parameters
      
   implicit none
      
   contains
      
   subroutine input_read
   ! Subroutine to read the input file and assign names and values 
   !	to global constants
	
	character(len = :), allocatable :: str
	character(len = 250) 		    :: line
      
   ! ----- Read planetary parameters from input file ----- !
      
   ! Open file for reading 
	write(*,*) '(input_read.f90) Reading the input.inp file..'
   open(unit = 11, file = inp_file)

	! --- Go line by line and read
	
		! Planet name
		read(11,'(A)') line
		p_name = get_word(line, 3)

     		! Log10 of n0
      	read(11,'(A)') line
      	str = get_word(line, 7)
     		read(str,*) n0
    
	     	! Planet radius
	     	read(11,'(A)') line
	     	str = get_word(line, 4)
	     	read(str,*) R0 
	     	
	     	! Planet mass
	     	read(11,'(A)') line
	     	str = get_word(line, 4)
	     	read(str,*) Mp 
	     	
	     	! Equilibrium temperature
	     	read(11,'(A)') line
	     	str = get_word(line, 4)
	     	read(str,*) T0 
	     	
	     	! Orbital distance
	     	read(11,'(A)') line
	     	str = get_word(line, 4)
	     	read(str,*) a_orb 
	     	
	     	! Escape radius
	     	read(11,'(A)') line
	     	str = get_word(line, 4)
	     	read(str,*) r_esc 
	     	
	     	! He/H number ratio
	     	read(11,'(A)') line
	     	str = get_word(line, 4)
	     	read(str,*) HeH 
		if (HeH .gt. 0.0e0) thereis_He = .true. 
	     	
		! 2D approximate method
		read(11,'(A)') line
		appx_mth = get_word(line, 4)
		
		! Read alpha if selected
		if (appx_mth .eq. 'alpha') then
			str = get_word(line, 6)
			read(str,*) a_tau 
		else
			a_tau = 0.0
		endif
		
		! Correct appx_meth keywords
		if (appx_mth .eq. 'Rate/4') appx_mth = 'Rate/4 + Mdot'
		if (appx_mth .eq. 'Rate/2') appx_mth = 'Rate/2 + Mdot/2'
		
		! Parent star mass
		read(11,'(A)') line
		str = get_word(line, 5)
	     	read(str,*) Mstar 
	     	
		! Spectrum type 
		read(11,'(A)') line
		sp_type = get_word(line, 3)

		! Next read properties of spectrum
		select case (sp_type)
		
			case ('Load')			! Load from file
				read(11,'(A)') line
				sed_file = get_word(line, 3)
				do_read_sed = .true.
			
			case ('Power-law')
				read(11,'(A)') line
				str = get_word(line, 3)
				read(str,*) PLind 
				is_PL_sed = .true.
	
			case ('Monochromatic')
			
				! Set corresponding logical to true
				is_monochr = .true.
				
				! Read photon enerrgy
				read(11,'(A)') line
				str = get_word(line, 4)
				read(str,*) e_low  
				
				! Remove helium if monochromatic and
				!	photon energy lower than helium ionization threshold
				if (e_low .lt. e_th_HeI) thereis_He = .false.
			     	
		end select
		
		! Only EUV status
		read(11,'(A)') line
		str = get_word(line, 4)
		if (str .eq. 'False') thereis_Xray = .true.
		
		! If not monochromatic, read energy bands
		if (.not. is_monochr ) then 
			
			if (.not.thereis_Xray) then 
			
				! Read e_low
				read(11,'(A)') line
				str = get_word(line, 4)
			     	read(str,*) e_low  
			     	
				! Read e_mid
				str = get_word(line, 6)
				read(str,*) e_mid 
				
				! Set e_top to default
				e_top = 1.24e3 
			     	
			else
				! Read e_low
				read(11,'(A)') line
				str = get_word(line, 4)
			     	read(str,*) e_low  
			     	
				! Read e_mid
				str = get_word(line, 6)
			     	read(str,*) e_mid 
				
				! Read e_top
				str = get_word(line, 8)
			     	read(str,*) e_top 
			     	
			endif
			
		endif
			
		! Read X-ray luminosity if included
		if (thereis_Xray) then
			read(11,'(A)') line
			str = get_word(line, 6)
		     	read(str,*) LX  
		else
			LX = 0.0
		endif
		
		! Read LEUV luminosity
		read(11,'(A)') line
		str = get_word(line, 6)
		read(str,*) LEUV  
		
		! Read grid type
		read(11,'(A)') line
		grid_type = get_word(line, 3)
		
		! Read numerical flux
		read(11,'(A)') line
		flux = get_word(line, 3)
		
		! Read reconstruction scheme
		read(11,'(A)') line
		rec_method = get_word(line, 3)
		if (rec_method.eq.'WENO3') use_weno3 = .true.
		if (rec_method.eq.'PLM')   use_plm = .true.
		
		! Include He23S
		read(11,'(A)') line
		str = get_word(line, 3)
		if (str .eq. 'True')  thereis_HeITR = .true.

		! Remove HeITR chemistry if He is not included
		if (.not. thereis_He) thereis_HeITR = .false.

		! IC status
		read(11,'(A)') line
		str = get_word(line, 3)
		if (str .eq. 'True')  do_load_IC = .true.
		
		! Do only post-processing
		read(11,'(A)') line
		str = get_word(line, 4)
		if (str .eq. 'True')  then
			do_only_pp  = .true.
			force_start = .false. ! Set to false to avoid overlap
		endif

		! Force start of sim.
		read(11,'(A)') line
		str = get_word(line, 3)
		if (str .eq. 'True')  then 
			force_start = .true.
			do_only_pp  = .false. ! Set to false to avoid overlap
		endif

   close(unit = 1)
	write(*,*) '(input_read.f90) Done'

   !------ Definition of physical parameters ------!
      
   n0     = 10.0**(n0)
   R0     = R0*RJ
   Mp     = Mp*MJ
   a_orb  = a_orb*AU
   Mstar  = Mstar*Msun
   Mrapp  = Mstar/Mp
   atilde = a_orb/R0
   r_max  = (3.0*Mrapp)**(-1.0/3.0)*atilde
            
	!------ Normalization constants ------!
      
   rho_bc = (1.0 + 4.0*HeH)/(1.0 + HeH)
   v0     = sqrt(kb_erg*T0/mu)       
   t_s    = R0/v0                    
   p0     = n0*mu*v0*v0              
   q0     = n0*mu*v0*v0*v0/R0	      
   b0     = (Gc*Mp*mu)/(kb_erg*T0*R0)       
   dp_bc  = 1.0e-10
	
   !------ Allocations ------!
      
   ! Allocate variables according to composition
   if (.not.thereis_He) then 
      	N_eq = 1
   else 
		if (thereis_HeITR) then 
			N_eq = 4
		else
			N_eq = 3
		endif
	endif
	
   lwa  = (N_eq*(3*N_eq+13))/2
   allocate (sys_sol(N_eq))
   allocate (sys_x(N_eq))
   allocate (wa(lwa)) 
       
   ! End of subroutine
   end subroutine input_read
      
   ! ------------------------------------------------------- !
      
   function get_word(string_in,n_word)
	! Function to read the nth_word in the current string
	! 	"Words" are separated by spaces
	
	character(len = *), intent(in) :: string_in
	integer, intent(in) :: n_word
	
	character(len = :), allocatable  :: string
	character(len = 300) :: c_string	
	character :: p_char,c_char
	integer :: counter
	integer :: c_word_counter
	integer :: str_len
	
	character(len = :), allocatable :: get_word
	
	! Initialize counters and strings
	counter        = 1
	c_word_counter = 0
	string   = trim(string_in)
	str_len  = len(string)
	p_char = ''
	c_char = ''
	c_string = ''
	
	! Loop inside the string	
	do while (counter .ge. 0 .and. counter .le. str_len)

		! Characters
		if (counter .ge. 2) then	! Skip if its the first iteration
			p_char = string(counter - 1:counter - 1)
		endif
		c_char = string(counter:counter)
		
		! If a character is found
		if (c_char .ne. '') then	
			
			! Attach character to current string
			c_string = trim(c_string) // c_char
			
			! Update counter and continue
			counter = counter + 1 
		
			continue			
		
		else
		
			! If it's a first space after a character 		
			if (p_char .ne. '') then		
				
				! Update word counter
				c_word_counter = c_word_counter + 1
				
				! Exit from loop if word counter 
				! is equal to n_word in input				
				if (c_word_counter .eq. n_word) then 
					get_word = trim(c_string)
					return
				endif
			
				! Reset current string	
				c_string = ''
			
				! Update counter
				counter = counter + 1
				
			else	! If multiple spaces
				
				counter = counter + 1
				continue
			endif
				
		endif
		
		! If last character, return
		if (counter .eq. str_len) then 
			
			! Update string 
			c_char = string(counter:counter)
			c_string = trim(c_string) // c_char
			
			! Update word counter
			c_word_counter = c_word_counter + 1
				
			! Exit from loop if word counter 
			! is equal to n_word in input				
			if (c_word_counter .eq. n_word) then 
				get_word = trim(c_string)
				return
			endif
		endif
		
	enddo	! End of while loop
	
	! End of get_word function
	end function
      
    ! End of module
	end module Read_input
