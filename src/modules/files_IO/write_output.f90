      module output_write
      ! Write the output to the standard output files
      
      use global_parameters
      
      contains
      
      subroutine write_output(rho,v,p,T,heat,cool,eta,    &
                              nhi,nhii,nhei,nheii,nheiii,nheiTR,flag)
                              
      character(len = 2) :: flag
      integer :: j
      real*8, dimension(1-Ng:N+Ng), intent(in) :: rho,v,p,T
      real*8, dimension(1-Ng:N+Ng), intent(in) :: heat,cool
      real*8, dimension(1-Ng:N+Ng), intent(in) :: eta
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhi,nhii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nhei,nheii,nheiii
      real*8, dimension(1-Ng:N+Ng), intent(in) :: nheiTR
      
      
      !---- Write thermodynamic profiles ----! 
          
      if (flag.eq.'eq') then    
      	open(unit = 2, file = './output/Hydro_ioniz.txt')
	   else	! Change output file after postprocessing
		   open(unit = 2, file = './output/Hydro_ioniz_adv.txt')
	   endif
      
         do j = 1-Ng,N+Ng
            write(2,*) r(j),        &     ! Rad. dist.
                     rho(j)*n0,     &     ! Density
                     v(j)*v0,       &     ! Velocity
                     p(j)*p0,       &     ! Pressure
                     T(j)*T0,       &     ! Temperature
                     heat(j)*q0,    &     ! Rad. heat.
                     cool(j)*q0           ! Rad. cool.
         enddo
      close(2)
      
      !---- Write ionization profiles ----!  
      if (flag.eq.'eq') then   
      	open(unit = 3, file = './output/Ion_species.txt')
	   else	! Change output file after postprocessing
		   open(unit = 3, file = './output/Ion_species_adv.txt')
	   endif
	
      do j = 1-Ng,N+Ng
         
         write(3,*) r(j),       & ! Rad. dist.
                  nhi(j)*n0,    & ! HI
                  nhii(j)*n0,   & ! HII
                  nhei(j)*n0,   & ! HeI
                  nheii(j)*n0,  & ! HeII
                  nheiii(j)*n0, & ! HeIII
                  nheiTR(j)*n0    ! HeITR
      enddo
      close(3)
                  
      ! End of subroutine
      end subroutine write_output
      
      ! End of module
      end module output_write
