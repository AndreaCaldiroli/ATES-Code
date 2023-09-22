   module utils
	! Collection of auxiliary subroutines
	
   use global_parameters

   implicit none

	contains

	! ------------------------------------------------------!

	subroutine calc_ne(nhii,nheii,nheiii,ne)
	! Calculate the free electron density
	
	real*8, dimension(1-Ng:N+Ng), intent(in) :: nhii
	real*8, dimension(1-Ng:N+Ng), intent(in) :: nheii,nheiii
	real*8, dimension(1-Ng:N+Ng), intent(out) :: ne
	
	if (thereis_He) then
		ne = nhii + nheii + 2.0*nheiii
	else
		ne = nhii
	endif	
	
	end subroutine calc_ne
	
	! ------------------------------------------------------!

	subroutine calc_ntot(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_tot)
	! Calculate the total atomic number density
	
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi,nhii
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhei,nheii,nheiii,nheiTR
	real*8, dimension(1-Ng:N+Ng), intent(out) :: n_tot
	
	if (thereis_He) then
		n_tot = nhi + nhii + nhei + nheii + nheiii
		if (thereis_HeITR) n_tot = n_tot + nheiTR
	else
		n_tot = nhi + nhii 
	endif

	end subroutine calc_ntot

	! ------------------------------------------------------!
	
	subroutine calc_rho(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_out)
	! Calculate the total mass density (adimensional)
	
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi,nhii
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhei,nheii,nheiii,nheiTR
	real*8, dimension(1-Ng:N+Ng), intent(out) :: n_out
	
	if (thereis_He) then
		n_out = nhi + nhii + 4.0*(nhei + nheii + nheiii)
		if (thereis_HeITR) n_out = n_out + 4.0*nheiTR
	else
		n_out = nhi + nhii 
	endif

	end subroutine calc_rho

	! ------------------------------------------------------!

	subroutine calc_column_dens(nhi,nhei,nheii,nheiTR,N1,N15,N2,NTR)
	! Calculates the column densities for given ionization profiles
	! 	by method of rectangles
    
	integer :: j
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhi
	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nhei,nheii,nheiTR
   real*8 :: dr 
	real*8, dimension(1-Ng:N+Ng), intent(out) :: N1
	real*8, dimension(1-Ng:N+Ng), intent(out) :: N15,N2,NTR
	
	! Initialize outputs
	N1  = 0.0
	N15 = 0.0
	N2  = 0.0
	NTR = 0.0
	
	! Outer point
	N1(N+Ng)  = dr_j(N+Ng)*R0*nhi(N+Ng)
	if (thereis_He) then 

		N15(N+Ng) = dr_j(N+Ng)*R0*nhei(N+Ng)
		N2(N+Ng)  = dr_j(N+Ng)*R0*nheii(N+Ng)
		if(thereis_HeITR) NTR(N+Ng) = dr_j(N+Ng)*R0*nheiTR(N+Ng)

	endif
    
	do j = N+Ng-1,1-Ng,-1	
	
	    ! Spacing
	    dr = dr_j(j)*R0
	      
	      ! Evaluate new column densities by integration
	    N1(j)  = N1(j+1)  + nhi(j)*dr         ! HI
      	
		if (thereis_He) then 
			N15(j) = N15(j+1) + nhei(j)*dr	  					  ! HeI 
			N2(j)  = N2(j+1)  + nheii(j)*dr       				  ! HeII 
			if(thereis_HeITR) NTR(j) = NTR(j+1) + nheiTR(j)*dr	  ! HeI triplet
		endif 

	enddo 
	
	! End of subroutine
	end subroutine calc_column_dens

	! ------------------------------------------------------!

	subroutine calc_mmw(nh,nhe,ne,mmw)
	! Calculate the mean molecular weight for a certain ionization profile

	real*8, dimension(1-Ng:N+Ng), intent(in)  :: nh,nhe,ne
	real*8, dimension(1-Ng:N+Ng), intent(out) :: mmw
	
	if (thereis_He) then
		mmw = (nh + 4.0*nhe)/(nh + nhe + ne)
	else
		mmw = nh/(nh + ne)
	endif

	! End of subroutine
	end subroutine calc_mmw

	! End of module
	end module utils 
