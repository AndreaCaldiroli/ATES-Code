      module S_estimate_HLLC
      ! Estimate for left, middle and right velocities
      
      use global_parameters
      
      
      implicit none
      
      contains
      
      subroutine speed_estimate_HLLC(WL,WR,SL,SR,S_star)
      real*8, intent(in) :: WL(3), WR(3)
      
      real*8 :: rhoL,uL,pL,cL
      real*8 :: rhoR,uR,pR,cR
      real*8 :: rho_bar,c_bar
      real*8 :: Q,Q_user
      real*8 :: p_min,p_max,p_star
      real*8 :: z,pLR
      real*8 :: AL,BL,AR,BR,gL,gR
      real*8 :: qL,qR,u_star
                
      real*8, intent(out) :: SL,SR,S_star
                       
      ! Exctract left state
      rhoL = WL(1)
      uL = WL(2)
      pL = WL(3)
      
      ! Exctract right state
      rhoR = WR(1)
      uR = WR(2)
      pR = WR(3)
      
      ! Sound speed
      cL = sqrt(g*pL/rhoL)
      cR = sqrt(g*pR/rhoR)
      
      ! Average density and sound speed
      rho_bar = 0.5*(rhoL+rhoR)
      c_bar   = 0.5*(cL+cR)     
      
      
      !-------------------------------------------------------------!
      
      ! Substitutions
      p_min = min(pL,pR)
      p_max = max(pL,pR)

      Q = p_max/p_min
      Q_user = 2.0

      ! Initial PVRS guess
      p_star = 0.5*(pL+pR)-0.5*(uR-uL)*rho_bar*c_bar
      u_star = 0.5*(uL+uR)-0.5*(pR-pL)/(c_bar*rho_bar)
      
      ! Construct middle state approximation
      if(Q.gt.Q_user) then      ! PVRS
      
            if(p_star.le.p_min) then      ! TRRS
                  
                  ! Subs
                  z = 0.5*(g-1.0)/g
                  pLR = (pL/pR)**z
                  
                  p_star = sqrt((cL+cR-0.5*(g-1.0)*(uR-uL))/ &
                                (cL/pL**z + cR/pR**z)) 
                  u_star = (pLR*uL/cL+uR/cR                      &
                              +2.0*(pLR-1.0)/(g-1.0))/     &
                           (pLR/cL+1.0/cR)

            elseif(p_star.ge.p_max) then  ! TSRS
                  
                  ! Subs
                  AL = 2.0/((g+1.0)*rhoL)
                  BL = (g-1.0)/(g+1.0)*pL
                  AR = 2.0/((g+1.0)*rhoR)
                  BR = (g-1.0)/(g+1.0)*pR
                  
                  p_star = max(0.0, p_star)
                  gL = sqrt(AL/(p_star + BL))
                  gR = sqrt(AR/(p_star + BR))
                  
                  p_star = (gL*pL+gR*pR-uR+uL)/(gL+gR)
                  u_star = 0.5*(uL+uR)+ &
                           0.5*((p_star-pR)*gR-(p_star-pL)*gL)

            endif
      endif
      
      
      !-------------------------------------------------------------!
      
      ! Construct speed estimates
      
      ! qL
      if(p_star.lt.pL) then
            qL = 1.0
      else
            qL = sqrt(1.0 + 0.5*(g-1.0)/g*(p_star/pL-1.0))
      endif
      
      ! qR
      if(p_star.lt.pR) then
            qR = 1.0
      else
            qR = sqrt(1.0 + 0.5*(g-1.0)/g*(p_star/pR-1.0))
      endif
      
      ! Speed estimates
      SL = uL-cL*qL
      SR = uR+cR*qR
      
      ! Use S_star approximation from Batten et al.
      S_star = (pR - PL + rhoL*uL*(SL-uL) - rhoR*uR*(SR-uR))/     &
               (rhoL*(SL-uL) - rhoR*(SR-uR))
      
      end subroutine speed_estimate_HLLC
      
      
      end module S_estimate_HLLC
