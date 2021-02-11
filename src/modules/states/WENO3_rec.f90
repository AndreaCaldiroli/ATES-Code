      module WENO3_reconstruction
      ! ESWENO3 reconstruction method from Yamaleev & Carpenter (2009) 
      !     and Mignone et al (2011)
      
      use global_parameters
      use Conversion
      
      implicit none
      
      contains
            
      subroutine WENO3_rec(u_in,WR,WL)
      real*8, dimension(1-Ng:N+Ng,3), intent(in) :: u_in
      integer :: j,k
      real*8, dimension(1-Ng:N+Ng,3) :: W,dW
	real*8, dimension(3) :: dWp,dWm
	real*8, dimension(3) :: b0,b1
	real*8, dimension(3) :: tau
	real*8, dimension(3) :: S0,S1
	real*8 :: dxj,dx2
	real*8 :: rm,rp
	real*8, dimension(1-Ng:N+Ng) :: dV,C1,C2,D1,D2
	real*8, dimension(1-Ng:N+Ng,3), intent(out) :: WL,WR
	
      ! Convert to primitive variables
      call U_to_W(u_in,W)
      
      ! Evaluate jumps at interfaces
      dW(1-Ng:N+Ng-1,:) = W(2-Ng:N+Ng,:) - W(1-Ng:N+Ng-1,:)
      dW(N+Ng,:) = 0.0
      
      ! Evaluate cell edges and cell volumes
      do j = 0, N+1
      
	      rp = 0.5*(r(j+1) + r(j))
	      rm = 0.5*(r(j) + r(j-1))
	      dV(j) = (rp*rp*rp - rm*rm*rm)
      enddo
      dV(1-Ng) = dV(2-Ng)
      dV(N+Ng) = dV(N+Ng-1)
      
      
      ! Evaluate reconstruction geometry-dependent coefficients
      do j = 1-Ng,N+Ng-1
	  	C1(j) = dV(j+1)/(dV(j) + dV(j+1))
	  	C2(j) = 1.0 - C1(j)
	enddo
	
	do j = 2-Ng,N+Ng-1  	
	  	D1(j) = dV(j+1)/(dV(j) + dV(j-1))
	  	D2(j) = dV(j-1)/(dV(j) + dV(j+1))
      enddo
      
      
      
      ! Construct smoothness indicators and the reconstructed values
      do j = 0,N+1
      
      	dWp = dW(j,:)
      	dWm = dW(j-1,:)
      	
      	dxj = 0.5*(r(j+1) - r(j-1))
      	dx2 = dxj*dxj
      	
      	do k = 1,3
      		b0(k) = dWp(k)*dWp(k) + dx2
      		b1(k) = dWm(k)*dWm(k) + dx2
      	enddo
      	
      	tau = dWp - dWm
		
      	
      	S0 = 1.0 + tau*tau/b0
      	S1 = 1.0 + tau*tau/b1
      	
      	
      	WL(j,:) = W(j,:)	&
      		  + (S0*C1(j)*dWp + D1(j)*S1*C1(j-1)*dWm) &
      		  /(S0 + D1(j)*S1)
      	
      	WR(j-1,:) = W(j,:)  &
      		  - (D2(j)*S0*C2(j)*dWp + S1*C2(j-1)*dWm) &
      		   /(D2(j)*S0 + S1)
      enddo
      
      
      end subroutine WENO3_rec
      
      
      
      
      end module WENO3_reconstruction
