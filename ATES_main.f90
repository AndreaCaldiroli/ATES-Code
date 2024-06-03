      program Hydro_ioniz

      use global_parameters
      use Read_input
      use Initialization
      use setup_report
      use eval_time_step
      use utils
      use Conversion
      use ionization_equilibrium
      use Reconstruction_step
      use RK_integration
      use BC_Apply
      use output_write
      use post_processing
      use ionization_equilibrium
      
      implicit none
      
      ! Logical variables
      logical :: l_isnan = .false.
      
      ! Integers variables
      integer :: j,k
      
      ! Timing variables
      integer :: n_hrs 
      integer :: n_min 
      real*8  :: start, finish 
      real*8  :: exec_time
      real*8  :: n_sec
      real*8  :: dum
      
      ! Momentum variables
      real*8, dimension(1-Ng:N+Ng) :: mom
      real*8 :: mom_max,mom_min

      ! Maximum eigenvalue      
      real*8 :: alpha

      ! Temporal step
      real*8 :: dt
      
      ! Mdot value
      real*8 :: Mdot
      
      ! Vectors of thermodynamical variables
      real*8, dimension(1-Ng:N+Ng) :: rho,v,E,p,T,cs
      real*8, dimension(1-Ng:N+Ng) :: heat,cool
      real*8, dimension(1-Ng:N+Ng) :: eta 
      real*8, dimension(1-Ng:N+Ng) :: nhi,nhii 
      real*8, dimension(1-Ng:N+Ng) :: nhei,nheii,nheiii 
      real*8, dimension(1-Ng:N+Ng) :: nheiTR
      real*8, dimension(1-Ng:N+Ng) :: ne,n_tot  
      real*8, dimension(1-Ng:N+Ng,6) :: f_sp
       
      ! Conservative and primitive vectors
      real*8, dimension(1-Ng:N+Ng,3) :: u,u1,u2,u_old
      real*8, dimension(1-Ng:N+Ng,3) :: W,WL,WR
          
      ! Flux and source vectors
      real*8, dimension(1-Ng:N+Ng,3) :: dF,S
      
      !------------------------------------------------! 
      
      ! Open output report file 
      open(unit = outfile, file = 'ATES.out')
      
      !------------------------------------------------!
      
      ! Read planetary parameters from file
      call input_read
      
      !------------------------------------------------!
      
      ! Initialize simulations
      call init(W,u,f_sp)     

      !------------------------------------------------!
      
      ! Generate report of the current setup
      call write_setup_report

	!---------------------------------------------------!

	! Close outfile
	close(unit = outfile)

      !------------------------------------------------!

      !---- Start computation ----!
      
      ! Get starting time
      write(*,*) '(ATES_main.f90) Starting time integration..'
      start = omp_get_wtime()
      
      !------ Main temporal loop ------!
      do while( .not.is_mom_const .or. force_start)

            !---- Time step evaluation ----! 
            call eval_dt(W,dt)
            
            !-------------------------------------------------!
            
            !--- Thermodynamic evolution ---!
            
            ! Save previous step solution
            u_old = u
            
            ! FIRST RK STEP
            
            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            call Reconstruct(u,WL,WR) 
            
            ! Evaluate flux difference and source terms
            call RK_rhs(u,WL,WR,alpha,dF,S)
            
            u1 = u - dt*(dF - S)
             
            ! Apply boundary conditions
            call Apply_BC(u1,u1)
                      
            !----------------------------
            
            ! SECOND RK STEP

            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            call Reconstruct(u1,WL,WR) 

            ! Evaluate flux difference and source terms
            call RK_rhs(u1,WL,WR,alpha,dF,S)
                    
            ! Advance in time
            u2 = (3.0*u + u1 - dt*(dF - S))/4.0

            ! Apply boundary conditions            
            call Apply_BC(u2,u2)
            
            !----------------------------
            
            ! THIRD RK STEP

            ! Reconstruct u+_{j+1/2}, u-_{j+1/2}
            call Reconstruct(u2,WL,WR) 

            ! Evaluate flux difference and source terms
            call RK_rhs(u2,WL,WR,alpha,dF,S)
                  
            ! Advance in time
            u = (u + 2.0*(u2 - dt*(dF - S)))/3.0

            ! Apply boundary conditions
            call Apply_BC(u,u)
 		
		!------------------------------------------------!
 		
 		!---- Ionization Equilibrium ----!
 		
            ! Extract primitive variables
            call U_to_W(u,W)
            rho = W(:,1)
            v   = W(:,2)
            p   = W(:,3)
            
            ! Evaluate species densities
            nhi    = rho*f_sp(:,1)
            nhii   = rho*f_sp(:,2)
            if (thereis_He) then 
                  nhei   = rho*f_sp(:,3)
                  nheii  = rho*f_sp(:,4)
                  nheiii = rho*f_sp(:,5)
                  if (thereis_HeITR) nheiTR = rho*f_sp(:,6)
            endif
            call calc_ne(nhii,nheii,nheiii,ne)
            
            ! Total number density
            call calc_ntot(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_tot)
              	 
            ! Temperature profile
            T = p/(n_tot + ne)
            
            ! Evaluate ionization equilibrium
            call ioniz_eq(T,rho,f_sp,rho,f_sp,heat,cool,eta)
            
            !Evaluate partial densities
            nhi    = rho*f_sp(:,1)
            nhii   = rho*f_sp(:,2)
            if (thereis_He) then 
                  nhei   = rho*f_sp(:,3)
                  nheii  = rho*f_sp(:,4)
                  nheiii = rho*f_sp(:,5)
                  if (thereis_HeITR) nheiTR = rho*f_sp(:,6)
            endif
            call calc_ne(nhii,nheii,nheiii,ne)
            
            ! Total number density
            call calc_ntot(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_tot)
            
            ! Evaluate updated pressure
            p = (n_tot + ne)*T 
            
            ! Convert to primitive profiles
            W(:,1) = rho
            W(:,2) = v
            W(:,3) = p
            
            ! Revert to conservative 
            call W_to_U(W,u)
  	
  	      !------------------------------------------------!
  	      
  	      ! Evolve solution in time due to source terms
            !     using a forward Euler step
            
            u(:,3) = u(:,3) + dt*(heat - cool)

		call Apply_BC(u,u)

            !------------------------------------------------!
            
            ! Convert to physical variables and extract profiles
            call U_to_W(u,W)
            rho = W(:,1)
            v   = W(:,2)
            p   = W(:,3)
            E   = u(:,3)
            
            ! Evaluate ionized densities and electron density
            nhi    = rho*f_sp(:,1)
            nhii   = rho*f_sp(:,2)
            if (thereis_He) then 
                  nhei   = rho*f_sp(:,3)
                  nheii  = rho*f_sp(:,4)
                  nheiii = rho*f_sp(:,5)
                  if (thereis_HeITR) nheiTR = rho*f_sp(:,6)
            endif
  	      call calc_ne(nhii,nheii,nheiii,ne)
  	      
  	      ! Total number density
  	      call calc_ntot(nhi,nhii,nhei,nheii,nheiii,nheiTR,n_tot)
  	      
  	      ! Temperature profile
  	      T = p/(n_tot + ne)
 	      		
		! Evaluate momentum
		mom = rho*v*r*r
		
            !---------------------------------------------------!
        
            !--- Loop counters and escape condition ---!
            
            ! Update counter
            count = count + 1
            
            ! Maximum and minimum value for momentum
            mom_max = maxval(abs(mom(j_min:N)))
            mom_min = minval(abs(mom(j_min:N)))

            ! Evaluate relative momentum variation
            du = abs((mom_max-mom_min)/mom_min)
            
            ! Evaluate variation of time derivative          
      	u(:,2) = u(:,2) + 1.0e-16 ! To avoid division by zero          
	      
	      ! --- Infty-norm
	      dtu = max(maxval(abs(1.0-u(j_min:N,1)/u_old(j_min:N,1))), &
	      	    maxval(abs(1.0-u(j_min:N,2)/u_old(j_min:N,2))))
		dtu = max(maxval(abs(1.0-u(j_min:N,3)/u_old(j_min:N,3))), &
	      	    dtu)
            
            ! Adjust logical for loop if necessary
            is_mom_const = du .lt. du_th
	      is_zero_dt   = dtu .lt. dtu_th
            
            ! Write to standard output
            write(*,*) count,du !,dtu 
				
            !---------------------------------------------------!
            
            ! Detect NaNs
            do j = 1-Ng,N+Ng
            	do k = 1,3
            		dum = u(j,k)
            		if (dum.ne.dum) then
                              write(*,*)
            			write(*,'(A20,F8.6)') 'NaN detected at r = ',r(j)
                              write(*,'(A6,E13.6)') 'rho = ',W(j,1)*n0
                              write(*,'(A4,E13.6)') 'v = ',W(j,2)*v0/1.0e5
                              write(*,'(A4,E13.6)') 'p = ',W(j,3)*p0
                              write(*,'(A4,F8.1)')  'T = ',T(j)*T0
                              write(*,'(A6,E13.6)') 'nhi = ',nhi(j)*n0
                              write(*,'(A7,E13.6)') 'nhii = ',nhii(j)*n0
                              write(*,'(A7,E13.6)') 'nhei = ',nhei(j)*n0
                              write(*,'(A8,E13.6)') 'nheii = ',nheii(j)*n0
                              write(*,'(A9,E13.6)') 'nheiii = ',nheiii(j)*n0
                              write(*,'(A9,E13.6)') 'nheiTR = ',nheiTR(j)*n0
                              l_isnan = .true.
                        endif
            	enddo
         	enddo
            if(l_isnan) exit
            
            !---------------------------------------------------!
            
            !--- Write to files every 1000th iteration---! 
            if (mod(count,1000).eq.1) then
                  
                  ! Write thermodynamic and ionization profiles
                  call write_output(rho,v,p,T,heat,cool,eta,    &
                                    nhi,nhii,nhei,nheii,nheiii, &
                                    nheiTR,'eq')
                                   
            endif     
            
            ! Exit from temporal loop if only post processing has to be done
	      if (do_only_pp) exit

            ! Force continue for the first 1000 loops if force_start is enabled
            if (force_start) force_start = count .le. 1000
            
      !---------------------------------------------------!
            
      ! End of temporal while loop      
      enddo
      write(*,*) '(ATES_main.f90) Time integration done.'

      !---------------------------------------------------!
      
      ! Write final thermodynamic and ionization profiles
      write(*,*) '(ATES_main.f90) Writing final results to file..'
      call write_output(rho,v,p,T,heat,cool,eta,    &
                        nhi,nhii,nhei,nheii,nheiii,nheiTR,'eq')

      !---------------------------------------------------!
      
      ! Post processing to include advection
      write(*,*) '(ATES_main.f90) Starting the post processing routine..'

      call post_process_adv(rho,v,p,T,heat,cool,eta,    &
                            nhi,nhii,nhei,nheii,nheiii,nheiTR)
      
      write(*,*) '(ATES_main.f90) Post processing routine done.'                           
      
      !---------------------------------------------------!                            
                                    
      ! Get CPU time
      finish = omp_get_wtime()
      
      ! Write final execution time
      exec_time = finish - start
      n_hrs = floor(exec_time/3600.0)
      n_min = floor((exec_time - 3600.0*n_hrs)/60.0)
      n_sec = exec_time - 3600.0*n_hrs - 60.0*n_min

      write(*,*) ' '
      write(*,100) 'Execution Time = ', n_hrs,' h ', &
                                        n_min,' m ', &
                                        n_sec,' s'
100   format  (A17,I2,A3,I2,A3,F8.5,A2) 
      
      !---------------------------------------------------!
      
      ! Choose index to evaluate Mdot far enough from the top boundary  
      j = N - 20
      
      ! Evaluate steady state log of Mdot
      Mdot = log10(4.0*pi*rho(j)*v(j)*r(j)*r(j)*n0*mu*R0*R0)
      
      ! Correct for the 2D approximation used
      if (appx_mth.eq.'Rate/2 + Mdot/2') Mdot = Mdot - log10(2.0)
      if (appx_mth.eq.'Mdot/4') Mdot = Mdot - log10(4.0)
      
      
      ! Write Mdot in output
      write(*,*) ' '
      write(*,*) '----- Results -----'
      write(*,*) ' '
      write(*,*) '---> 2D approximate method: ', appx_mth
      write(*,101) ' ---> Log10 of steady-state Mdot = ', Mdot, ' g/s'
      
      ! Write Mdot to report file 
      open(unit = outfile, file = 'ATES.out', access = 'append' )
      	write(outfile,101) ' '
      	write(outfile,102) ' - Log10 of steady-state Mdot = ', Mdot, ' g/s'
      close(unit = outfile)
      
101   format (A35,F5.2,A4)
102   format (A32,F5.2,A4)

      !---------------------------------------------------!
      
      ! End of program
      end program Hydro_ioniz
