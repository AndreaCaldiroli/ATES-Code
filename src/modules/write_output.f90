      module output_write
      ! Write the output to the standard output files
      
      use global_parameters
      
      contains
      
      
      subroutine write_output(rho,v,p,T,heat,cool,eta,    &
                              nhi,nhii,nhei,nheii,nheiii)
                              
      character(len=6) :: tab = '      '
      integer :: j
      real*8, dimension(1-Ng:N+Ng), intent(in) :: rho,v,p,T
      real*8, dimension(1-Ng:N+Ng), intent(in) :: heat,cool
      real*8, dimension(1-Ng:N+Ng), intent(in) :: eta
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi,nhii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei,nheii,nheiii
      
      
      !---- Write thermodynamic profiles ----!      
      open(unit = 2, file = 'Hydro_ioniz.txt')
      
            do j = 1-Ng,N+Ng
                  write(2,2) r(j),tab,          &     ! Rad. dist.
                             rho(j)*rho0,tab,   &     ! Density
                             v(j),tab,          &     ! Velocity
                             p(j)*p0,tab,       &     ! Pressure
                		     T(j)*T0,tab,       &     ! Temperature
                		     heat(j)*q0,tab,    &     ! Rad. heat.
                		     cool(j)*q0,tab,    &     ! Rad. cool.
                		     eta(j)                   ! Heat. eff.
            enddo
      close(2)
      
      
      !---- Write ionization profiles ----!  
      open(unit = 3, file = 'Ion_species.txt')
            do j = 1-Ng,N+Ng
            
                  write(3,3) r(j),tab,          &
                             nhi(j)*rho0,tab,   & ! HI
            	 	     nhii(j)*rho0,tab,  & ! HII
            		     nhei(j)*rho0,tab,  & ! HeI
            		     nheii(j)*rho0,tab, & ! HeII
             		     nheiii(j)*rho0   ! HeIII
            enddo
      close(3)
                  
                        
      !---- Format specification ----!
      
      ! Thermodynamic variables
2     format (F8.5,A,     &
              E12.5,A,    &
              E12.5,A,    &
              E12.5,A,    &
              F8.2,A,     &
              E12.5,A,    &
              E12.5,A,    &
              E12.5)             
 
      ! Ionization Profiles 
3     format (F8.5,A,4(E12.5,A),E12.5)
 
 
            
      end subroutine write_output
      
      
      ! End of module
      end module output_write
