      module  Read_input
      ! Read input planetary parameters adn define 
      
      use global_parameters
      
      implicit none
      
      contains
      
      
      subroutine input_read
      
      character(len = 15) :: grid_type_inp
      character(len = 15) :: flux_inp
      character(len = 15) :: rec_method_inp
      character(len = 15) :: p_name
      character(len = 15) :: appx_mth_inp
      
      !------ Read planetary parameters from input file ------!
      
      open(unit = 1, file = "input.inp")
            read(1,200) p_name
            read(1,200) grid_type_inp
            read(1,200) flux_inp
            read(1,200) rec_method_inp
            read(1,300) n0
            read(1,400) R0
            read(1,400) Mp
            read(1,500) T0
            read(1,900) a_orb
            read(1,400) Mstar
            read(1,600) LX
            read(1,600) LEUV
            read(1,800) r_esc
            read(1,901) HeH
            read(1,200) ans_IC
            read(1,200) appx_mth_inp
            read(1,600) a_tau
            read(1,900) PLind
      close(1)


      !------ Trim and allocate character variables ------!
      grid_type  = trim(grid_type_inp)
      flux       = trim(flux_inp)
      rec_method = trim(rec_method_inp)
      appx_mth   = trim(appx_mth_inp)
      
      !------ Definition of physical parameters ------!
      
      n0     = 10.0**(n0)
      R0     = R0*RJ
      Mp     = Mp*MJ
      a_orb  = a_orb*AU
      J_XUV  = (10.0**LX + 10.0**LEUV)/(4.0*pi*a_orb**2.0)
      Mstar  = Mstar*Msun
      Mrapp  = Mstar/Mp
      atilde = a_orb/R0
      r_max  = (3.0*Mrapp)**(-1.0/3.0)*atilde      ! Roche Lobe distance
            
      !------ Normalization constants ------!
      
      rho1 = (1.0 + 4.0*HeH)/(1.0 + HeH)
      v0   = sqrt(kb_erg*T0/mu)       
      t_s  = R0/v0                    
      p0   = n0*mu*v0*v0              
      q0   = n0*mu*v0*v0*v0/R0	      
      b0   = (Gc*Mp*mu)/(kb_erg*T0*R0)       
       
       
      !------ Reading format specification ------!
200   format (A15)
300   format (E7.1)
400   format (F5.3)
500   format (F6.1)    
600   format (F6.3)
800   format (F4.1)      
900   format (F6.4)
901   format (F8.6)  
            
      end subroutine input_read
      
      ! End of module
      end module Read_input
