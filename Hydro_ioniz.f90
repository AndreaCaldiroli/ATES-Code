      program Hydro_ioniz

      use global_parameters
      use Read_input
      use Initialization
      use Conversion
      use ionization_equilibrium
      use Reconstruction_step
      use RK_integration
      use BC_Apply
      use output_write
      
      use ionization_equilibrium
      
      implicit none
      

      ! Integers variables
      integer :: k,j
      
      ! Timing variables
      integer :: n_hrs 
      integer :: n_min 
      real*8  :: start, finish 
      real*8  :: exec_time
      real*8  :: n_sec
      
      ! Momentum variables
      real*8, dimension(1-Ng:N+Ng) :: mom
      real*8 :: mom_max,mom_min

      ! Maximum eigenvalue      
      real*8 :: alpha

      ! Temporal step
      real*8 :: dt
 
      ! Vectors of thermodynamical variables
      real*8, dimension(1-Ng:N+Ng) :: rho,v,E,p,T,cs
      real*8, dimension(1-Ng:N+Ng) :: heat,cool
      real*8, dimension(1-Ng:N+Ng) :: eta 
      real*8, dimension(1-Ng:N+Ng) :: nhi,nhii 
      real*8, dimension(1-Ng:N+Ng) :: nhei,nheii,nheiii 
      real*8, dimension(1-Ng:N+Ng) :: ne,n_tot  
      real*8, dimension(1-Ng:N+Ng,5) :: f_sp
       
      ! Conservative and primitive vectors
      real*8, dimension(1-Ng:N+Ng,3) :: u,u1,u2 
      real*8, dimension(1-Ng:N+Ng,3) :: W,WL,WR
          
      ! Flux and source vectors
      real*8, dimension(1-Ng:N+Ng,3) :: dF,S
      
      
      !------------------------------------------------!
      
      ! Read planetary parameters from file
      call input_read
      
      !------------------------------------------------!
      
      ! Initialize simulations
      call init(W,u,f_sp)     
            
      !------------------------------------------------!
      
      !---- Start computation ----!
      
      ! Get starting time
      start = omp_get_wtime()
      
      !------ Main temporal loop ------!
      
      do while(du.ge.du_th.or.count.le.10)

            !---- Time step evaluation ----! 
            
            ! Extract phisical variables
            rho = W(:,1)
            v   = W(:,2)
            p   = W(:,3)
            
            ! Evaluate sound speed
            cs = sqrt(g*p/rho)
            
            ! Maximum eigenvalue
            alpha = maxval(abs(v) + cs)

            ! Evaluate time step according to CFL condition
            dt = CFL*minval(dr_j/(abs(v) + cs))            
            
            !-------------------------------------------------!
            
            !--- Thermodynamic evolution ---!
            
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
 		
            ! Exctract primitive variables
            call U_to_W(u,W)
            rho = W(:,1)
            v   = W(:,2)
            p   = W(:,3)
            
            ! Evaluate species densities
            nhi    = rho*f_sp(:,1)
            nhii   = rho*f_sp(:,2)
            nhei   = rho*f_sp(:,3)
            nheii  = rho*f_sp(:,4)
            nheiii = rho*f_sp(:,5)
            ne     = nhii + nheii + 2.0*nheiii
            
            ! Total number density
            n_tot = nhi + nhii + nhei + nheii + nheiii 
              	 
            ! Temperature profile
            T = p/(n_tot + ne)
            
            ! Evaluate ionization equilibrium
            call ioniz_eq(T,rho,f_sp,rho,f_sp,heat,cool,eta)
            
            !Evaluate partial densities
            nhi    = rho*f_sp(:,1)
            nhii   = rho*f_sp(:,2)
            nhei   = rho*f_sp(:,3)
            nheii  = rho*f_sp(:,4)
            nheiii = rho*f_sp(:,5)
            ne     = nhii + nheii + 2.0*nheiii
            
            ! Total number density
            n_tot = nhi + nhii + nhei + nheii + nheiii 
            
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
            
            ! Convert to physical variables and exctract profiles
            call U_to_W(u,W)
            rho = W(:,1)
            v   = W(:,2)
            p   = W(:,3)
            E   = u(:,3)
            
            ! Evaluate ionized densities and electron density
  	      nhi    = rho*f_sp(:,1)
		nhii   = rho*f_sp(:,2)
		nhei   = rho*f_sp(:,3)
		nheii  = rho*f_sp(:,4)
		nheiii = rho*f_sp(:,5)
  	      ne     = nhii + nheii + 2.0*nheiii
  	      
  	      ! Total number density
  	      n_tot = nhi + nhii + nhei + nheii + nheiii  
  	      
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
                       
            ! Write to standard output
            write(*,*) count,du
            
            !---------------------------------------------------!
            
            !--- Write to files every 100th iteration---! 
            
            if (mod(count,100).eq.1) then
                  
                  ! Write thermodynamic and ionization profiles
                  call write_output(rho,v,p,T,heat,cool,eta,    &
                                    nhi,nhii,nhei,nheii,nheiii)
                                   
            endif     
            
      !---------------------------------------------------!
            
      ! End of temporal while loop      
      enddo

      !---------------------------------------------------!
      
      ! Get CPU time
      finish = omp_get_wtime()
      
      ! Write final execution time
      exec_time = finish - start
      n_hrs = nint(exec_time/3600.0)
      n_min = nint((exec_time - 3600.0*n_hrs)/60.0)
      n_sec = exec_time - 3600.0*n_hrs - 60.0*n_min

      
      write(*,100) 'Execution Time = ', n_hrs,' h ', &
                                        n_min,' m ', &
                                        n_sec,' s'
100   format  (A17,I2,A3,I2,A3,F8.5,A2) 
      
      !---------------------------------------------------!
            
      ! End of program
      end program Hydro_ioniz
