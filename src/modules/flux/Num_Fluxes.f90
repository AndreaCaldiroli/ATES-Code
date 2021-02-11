      module Numerical_Fluxes
      ! Collection of numerical flux functions
      
      use global_parameters
      use Conversion
      use S_estimate_HLLC
      use S_estimate_ROE
  
      implicit none
      
      contains
      
      ! Subroutine for the numerical flux
      subroutine Num_flux(WL,WR,NF,alpha,p_out)
      
      real*8,intent(in) :: WL(3),WR(3),alpha
      real*8 :: uL(3),uR(3)
      real*8 :: FL(3),FR(3)
      real*8 :: rhoL,vL,pL,aL,EL,HL,SL,phiL
      real*8 :: rhoR,vR,pR,aR,ER,HR,SR,phiR
      real*8 :: S_star
      real*8 :: usL(3),usR(3)
      real*8 :: s,rho_avg,v_avg,H_avg,a_avg          
	real*8 :: l1,l2,l3
	real*8 :: l1L,l1R,l3L,l3R 
	real*8 :: v_star,aL_star,aR_star
	real*8 :: drho,dvel,dp
	real*8 :: a1,a2,a3   
	real*8, dimension(3) :: K1,K2,K3
      real*8, intent(out) :: NF(3),p_out   

      
      ! Exctract left state
      rhoL = WL(1)
      vL   = WL(2)
      pL   = WL(3)
      EL   = 0.5*rhoL*vL*vL + pL/(g-1.0)
      HL   = (EL+pL)/rhoL
      aL   = sqrt(g*pL/rhoL)   
      
      ! Exctract right state
      rhoR = WR(1)
      vR   = WR(2)
      pR   = WR(3)
      ER   = 0.5*rhoR*vR*vR + pR/(g-1.0)
      HR   = (ER+pR)/rhoR
      aR   = sqrt(g*pR/rhoR)
            
      ! Evaluate numerical flux      
      select case(flux)
      
      case ('LLF') ! Local Lax Friedrichs
            
            ! Evaluate physical left and right flux
            call Phys_flux(WL,FL)
            call Phys_flux(WR,FR)
            
            ! Maximum eigenvalue between adjacent cells
		a1 = max(abs(vL+aL),abs(vR+aR))		
		
		! Get vector of conservative variables
		call W_to_U_comp(WL,uL)
            call W_to_U_comp(WR,uR)
            
            ! Evaluate numerical flux
            NF = 0.5*(FL + FR - a1*(uR-uL))
            
            ! Output pressure
            p_out = 0.5*(pR + pL)
      
      !----------------------------------------------!      
      
      case ('LF') ! Global Lax Friedrichs
            
            ! Evaluate physical left and right flux
            call Phys_flux(uL,FL)
            call Phys_flux(uR,FR)
            
            ! Get vector of conservative variables
            call W_to_U_comp(WL,uL)
            call W_to_U_comp(WR,uR)
            
            ! Evaluate numerical flux
            NF = 0.5*(FL + FR - alpha*(uR-uL))
            
            ! Output pressure
            p_out = 0.5*(pR + pL)
            
      !----------------------------------------------!
      
      case('HLLC') ! HLLC solver
            
            call W_to_U_comp(WL,uL)
            call W_to_U_comp(WR,uR)
      
      
            ! Speed estimates
            SL = min(0.0,min(vL-aL, vR-aR))
            SR = max(0.0,max(vL+aL, vR+aR))
            phiL = rhoL*(SL-vL)
            phiR = rhoR*(SR-vR)
            S_star = (pR - pL + phiL*vL - phiR*vR)/(phiL - phiR)
            

            ! HLLC state approximation
            
            ! Left state
            usL(1) = 1.0
            usL(2) = S_star
            usL(3) = EL/rhoL + (S_star-vL)*(S_star + pL/phiL)            
            usL = rhoL*(SL-vL)/(SL-S_star)*usL
            

            ! Right state
            usR(1) = 1.0
            usR(2) = S_star
            usR(3) = ER/rhoR + (S_star-vR)*(S_star + pR/phiR)           
            usR = rhoR*(SR-vR)/(SR-S_star)*usR

            
            ! Evaluate HLLC flux 
            
            if(SL.ge.(0.0)) then
                  
                  call Phys_flux(WL,NF)
                  
                  p_out = pL
                  
            elseif(SL.lt.(0.0).and.S_star.ge.(0.0)) then
            
                  call Phys_flux(WL,NF)

                  NF = NF + SL * (usL-uL)
                  
                  p_out = pL
                  
            elseif(S_star.lt.(0.0).and.SR.ge.(0.0)) then
            
                  call Phys_flux(WR,NF)
                  NF = NF + SR * (usR-uR)
                  
                  p_out = pR
            else
            
                  call Phys_flux(WR,NF)
                  
                  p_out = pR
                  
            endif
            
      !----------------------------------------------!
      
      case ('ROE')
      
		! Construct averaged states
		
		! Auxiliary parameter
		s = sqrt(rhoL)/(sqrt(rhoL)+sqrt(rhoR))
		
		! Density
		rho_avg = sqrt(rhoL*rhoR)
		
		! Velocity
		v_avg = s*vL + (1.0-s)*vR
		
		! Entalpy
		H_avg = s*HL + (1.0-s)*HR
		
		! Sound speed
		a_avg = sqrt((g-1.0)*(H_avg-0.5*v_avg*v_avg))
      
      
      	!---------------------------------!
      	
      	! Average eigenvalues
      
      	l1 = v_avg - a_avg
      	l2 = v_avg
      	l3 = v_avg + a_avg
      	
      	!---------------------------------!
      
		! Entropy correction
		
		! Evaluate the approximate velocities
		call speed_estimate_ROE(WL,WR,v_star,aL_star,aR_star)
		

		! Intermediate eigenvalues
		l1L = vL - aL
		l1R = v_star - aL_star
		l3L = v_star + aR_star
		l3R = vR + aR
		
		! Modify the eigenvalues if a rarefaction is present
		
		! Left rarefaction
		if (l1L.lt.(0.0).and.l1R.gt.(0.0)) then
		      l1 = l1L*(l1R-l1)/(l1R-l1L)
		endif
		
		! Right rarefaction
		if (l3L.lt.(0.0).and.l3R.gt.(0.0)) then
		      l3 = l3R*(l3-l3L)/(l3R-l3L)
		endif
		
		!---------------------------------!
		
		! Averaged eigenevectors
		
		! K1
		K1(1) = 1.0 
		K1(2) = v_avg - a_avg
		K1(3) = H_avg - v_avg*a_avg
		
		! K2
		K2(1) = 1.0 
		K2(2) = v_avg
		K2(3) = 0.5*v_avg*v_avg
		
		! K3
		K3(1) = 1.0 
		K3(2) = v_avg + a_avg
		K3(3) = H_avg + v_avg*a_avg
      
      	!---------------------------------!
      
		! Evaluate the conserved-variables differences		
		drho = rhoR - rhoL
		dvel = vR - vL
		dp   = pR - PL
		    
		! Evaluate the expansion coefficients		
		a1 = 0.5/a_avg**2.0*(dp - rho_avg*a_avg*dvel)
		a2 = drho - dp/a_avg**2.0
		a3 = 0.5/a_avg**2.0*(dp + rho_avg*a_avg*dvel)
             
             
            !---------------------------------!
      
      	! Evaluate the left and right fluxes      	
      	call Phys_flux(WL,FL)
      	call Phys_flux(WR,FR)


		! Evaluate the flux at interface ( eq.[11.29] Toro )
		NF = 0.5*(FR+FL)   &
		   - 0.5*(a1*abs(l1)*K1 + a2*abs(l2)*K2 + a3*abs(l3)*K3)
	      
	      ! Output pressure
            p_out = 0.5*(pR + pL)
		   
      end select
      
      end subroutine Num_flux
      
      
      !-----------------------------------------------------------!
      
      ! Subroutine to compute the physical flux function
      subroutine Phys_flux(W,PF)
      
      real*8, intent(in)  :: W(3)
      real*8 :: rho,v,p,E
      real*8, intent(out) :: PF(3)
      
      ! Get physical variables
      rho = W(1)
      v   = W(2)
      p   = W(3)
      E   = 0.5*rho*v*v + p/(g-1.0)
            
      ! Output exact flux vector
      PF(1) = rho*v
      PF(2) = rho*v*v
      PF(3) = v*(E+p)
      
      ! Add pressure for PLM discretization
      if (rec_method.eq.'PLM') then 
            PF(2) = PF(2) + p
      endif
      
      
      end subroutine Phys_flux     
      
      ! End of module
      end module Numerical_Fluxes
