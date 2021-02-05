      module Cooling_Coefficients
      ! Module containing rate coefficients and cooling rates 
      ! Included processes: bremsstrahlung, collisional ionization,
      !     collisional excitation, recombination
      
      implicit none
      
      contains
      
      !---------------------------------------------------!
      
      !--- Recombination ---!
      
      ! Case B recombination coefficient of HII
      double precision function rec_HII_B(T)
      real*8, intent(in) :: T
      real*8 :: xl
      
      xl        = 2.0*157807.0/T
      rec_HII_B = 2.753e-14*xl**1.5/(1.0+(xl/2.740)**0.407)**2.242
      
      end function rec_HII_B
      
      !--------------!
      
      ! Case B recombination coefficient of HeII
      double precision function rec_HeII_B(T)
      real*8, intent(in) :: T
      real*8 :: xl
      
      xl         = 2.0*285335.0/T
      rec_HeII_B = 1.26e-14*xl**0.750
      
      end function rec_HeII_B
      
      !--------------!
      
      ! Case B recombination coefficient of HeIII
      double precision function rec_HeIII_B(T)
      real*8, intent(in) :: T
      real*8 :: xl
      
      xl          = 2.0*631515.0/T
	rec_HeIII_B = 2.0*2.753e-14*xl**1.5/(1.0+(xl/2.740)**0.407)**2.242
      
      end function rec_HeIII_B
           
      !--------------!
      
      ! Recombination cooling rate for HII
      double precision function rec_cool_HII(T)
      real*8, intent(in) :: T
      real*8 :: xl
      
      xl = 2.0*157807.0/T
      rec_cool_HII = 3.435e-30*T*xl**1.970/         &
                     (1.0+(xl/2.250)**0.376)**3.720
      
      end function rec_cool_HII
      
      !--------------!
      
      ! Recombination cooling rate for HeII
      double precision function rec_cool_HeII(T)
      real*8, intent(in) :: T
            
      rec_cool_HeII = 1.38e-16*T*rec_HeII_B(T)
      
      end function rec_cool_HeII
      
      !--------------!
      
      ! Recombination cooling rate for HeIII
      double precision function rec_cool_HeIII(T)
      real*8, intent(in) :: T
      real*8 :: xl
      
      xl = 2.0*631515.0/T
      rec_cool_HeIII = 8.0*3.435e-30*T*xl**1.970/      &
                       (1.0+(xl/2.250)**0.376)**3.720
      
      end function rec_cool_HeIII
      
      
      !---------------------------------------------------!
      
      !--- Collisional ionization ---!
      
      
      ! Collisional ionization rate for HI
      double precision function ion_coeff_HI(T)
      real*8, intent(in) :: T
      real*8 :: th
      
      th = log(T*8.61733e-5)
      
      ion_coeff_HI = exp(-3.271396786e1 + 1.35365560e1*th       &
             -5.73932875*th**2.0    + 1.56315498*th**3.0        &
             -2.87705600e-1*th**4.0 + 3.48255977e-2*th**5.0     &
             -2.63197617e-3*th**6.0 + 1.11954395e-4*th**7.0     &
             -2.03914985e-6*th**8.0)
      
      end function ion_coeff_HI
      
      !--------------!
      
      ! Collisional ionization rate for HeI
      double precision function ion_coeff_HeI(T)
      real*8, intent(in) :: T
      real*8 :: th
      
      th = log(T*8.61733d-5)
      
      ion_coeff_HeI = exp(-4.409864886e1 + 2.391596563e1*th       &
      	 -1.07532302e1*th**2.0 + 3.05803875*th**3.0           &
      	 -5.6851189e-1*th**4.0 + 6.79539123e-2*th**5.0        &
      	 -5.0090561e-3*th**6.0 + 2.06723616e-4*th**7.0        &
      	 -3.64916141e-6*th**8.0)
      
      end function ion_coeff_HeI
            
      !--------------!
      
      ! Collisional ionization rate for HeII
      double precision function ion_coeff_HeII(T)
      real*8, intent(in) :: T
      real*8 :: xl
      
      xl = 2.0*631515.0/T
      
      ion_coeff_HeII =  19.95*exp(-xl/2.0)*T**(-1.5)*                &
      	            xl**(-1.089)/(1.0+(xl/0.553)**0.735)**1.275
      
      end function ion_coeff_HeII
      
    
      !---------------------------------------------------!
      
      !--- Bremmstrahllung ---!
      
      ! Gaunt factor function
      double precision function GF(T,Z)
      real*8,intent(in) :: T,Z
      
      if (T.lt.(Z*Z*3.2e5) ) then
            
            GF = 0.79464 + 0.1243*log10(T)     
      else
      
            GF = 2.13164 + 0.1243*log10(T)      
      endif      
      
      
      end function GF
      
      
      !---------------------------------------------------!
      
      !--- Collisional excitation ---!
      
      ! Collisional excitation rate for HI
      double precision function coex_rate_HI(T)
      real*8, intent(in) :: T
      
      coex_rate_HI = 7.5d-19/(1.0+sqrt(T/1.0e5))*exp(-118348.0/T)
      
      end function coex_rate_HI
      
      !--------------!
      
      ! Collisional excitation rate for HeI
      double precision function coex_rate_HeI(T)
      real*8, intent(in) :: T
      
      coex_rate_HeI = 1.1e-19*T**0.082*exp(-2.3e5/T)
      
      end function coex_rate_HeI
      
      !--------------!
      
      ! Collisional excitation rate for HeII
      double precision function coex_rate_HeII(T)
      real*8, intent(in) :: T
      
      coex_rate_HeII = 5.54e-17*T**(-0.397)/(1.0+sqrt(T/1.0e5))     &
      	  	                 *exp(-473638.0/T)
      
      end function coex_rate_HeII
      

      ! End of module
      end module Cooling_Coefficients
