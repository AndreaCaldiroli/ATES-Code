      module IC_load
      ! Module to load previous ICs, stored in the following files:
      ! - Hydro_ioniz_IC.txt
      ! - Ion_species_IC.txt
      
      use global_parameters
      
      implicit none
      
      contains
      
      subroutine load_IC(rho,v,p,T,f_sp,W)
      
      ! Integer variables
      integer :: j

      ! Thermodynamic loading variables
      real*8, dimension(1-Ng:N+Ng) :: nhi_l,nhii_l
      real*8, dimension(1-Ng:N+Ng) :: nhei_l,nheii_l,nheiii_l,nheiTR_l
      
      ! Auxiliary temporary variable
      real*8 :: tmp 
      
      ! Output variables
      real*8, dimension(1-Ng:N+Ng),   intent(out) :: rho,v,p,T
      real*8, dimension(1-Ng:N+Ng,6), intent(out) :: f_sp
      real*8, dimension(1-Ng:N+Ng,3), intent(out) :: W
      
      
	   !-------------------------------------!
	
      ! Load thermodynamic variables
      open(unit = 1, file = 'output/Hydro_ioniz_IC.txt')
         do j = 1-Ng,N+Ng
            read(1,*) tmp, tmp, v(j), p(j), T(j), tmp, tmp
         enddo
      close(1)
      
      ! Adimensionalize
      v = v/v0
      p = p/p0
      T = T/T0
      
      ! Load ionization profiles
      open(unit = 2, file = 'output/Ion_species_IC.txt')
         do j = 1-Ng,N+Ng
         read(2,*)  r(j),      & ! Rad.
                  nhi_l(j),    & ! HI
                  nhii_l(j),   & ! HII
                  nhei_l(j),   & ! HeI
                  nheii_l(j),  & ! HeII
                  nheiii_l(j), & ! HeIII
                  nheiTR_l(j)
         enddo
      close(2)
      
      
      ! Construct mass density profile (adimensional)
      rho = (nhi_l + nhii_l + 4.0*(nhei_l + nheii_l + nheiii_l))/n0	
         
      ! Construct the ration n/rho
      f_sp(:,1) = nhi_l/(rho*n0)	    ! HI
      f_sp(:,2) = nhii_l/(rho*n0)	    ! HII
      f_sp(:,3) = nhei_l/(rho*n0)	    ! HeI
      f_sp(:,4) = nheii_l/(rho*n0)      ! HeII		
      f_sp(:,5) = nheiii_l/(rho*n0)     ! HeIII
      f_sp(:,6) = nheiTR_l/(rho*n0)
         
      ! Construct matrix of primitive profiles
      W(:,1) = rho
      W(:,2) = v
      W(:,3) = p

      ! End of subroutine
      end subroutine load_IC
      
      ! End of module
      end module IC_load
