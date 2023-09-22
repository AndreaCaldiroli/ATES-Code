   module Cooling_Coefficients
   ! Module containing rate coefficients and cooling rates 
   ! Included processes: bremsstrahlung, collisional ionization,
   !     collisional excitation, recombination
   ! 
   ! --------------------
   ! 
   ! - AC, 14/04/22: Added He metastable collisional strength 
   !	coefficients
   
   use global_parameters
   
   implicit none
   
   contains
   
   !---------------------------------------------------!
   
   !--- Recombination ---!
   
   ! Case B recombination coefficient of HII
   subroutine rec_HII_B(T,coeff_rec_HII_B)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng) :: xl
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_rec_HII_B

   xl        = 2.0*157807.0/T
   coeff_rec_HII_B = 2.753e-14*xl**1.5/(1.0+(xl/2.740)**0.407)**2.242
   
   end subroutine rec_HII_B
   
   !--------------!
   
   ! Case B recombination coefficient of HeII
   subroutine rec_HeII_B(T,coeff_rec_HeII_B)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng) :: xl
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_rec_HeII_B
   
   xl         = 2.0*285335.0/T
   coeff_rec_HeII_B = 1.26e-14*xl**0.750
   
   end subroutine rec_HeII_B
   
   !--------------!
   
   ! Case B recombination coefficient of HeIII
   subroutine rec_HeIII_B(T,coeff_rec_HeIII_B)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng) :: xl
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_rec_HeIII_B

   xl          = 2.0*631515.0/T
   coeff_rec_HeIII_B = 2.0*2.753e-14*xl**1.5/(1.0+(xl/2.740)**0.407)**2.242
   
   end subroutine rec_HeIII_B
         
   !--------------!
   
   ! Recombination cooling rate for HII
   subroutine rec_cool_HII(T,coeff_rec_cool_HII)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng) :: xl
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_rec_cool_HII

   xl = 2.0*157807.0/T
   coeff_rec_cool_HII = 3.435e-30*T*xl**1.970/           &
                  (1.0+(xl/2.250)**0.376)**3.720
   
   end subroutine rec_cool_HII
   
   !--------------!
   
   ! Recombination cooling rate for HeII
   subroutine rec_cool_HeII(T,coeff_rec_cool_HeII)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng) :: coeff_rec_HeII_B
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_rec_cool_HeII

   call rec_HeII_B(T,coeff_rec_HeII_B)
   coeff_rec_cool_HeII = 1.38e-16*T*coeff_rec_HeII_B
   
   end subroutine rec_cool_HeII
   
   !--------------!
   
   ! Recombination cooling rate for HeIII
   subroutine rec_cool_HeIII(T,coeff_rec_cool_HeIII)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng) :: xl
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_rec_cool_HeIII

   xl = 2.0*631515.0/T
   coeff_rec_cool_HeIII = 8.0*3.435e-30*T*xl**1.970/      &
                           (1.0+(xl/2.250)**0.376)**3.720
   
   end subroutine rec_cool_HeIII
   
   !---------------------------------------------------!
   
   !--- Collisional ionization ---!
   
   ! Collisional ionization rate for HI
   subroutine ion_coeff_HI(T,a_ion_coeff_HI)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng) :: th
   real*8, dimension(1-Ng:N+Ng), intent(out) :: a_ion_coeff_HI

   th = log(T*8.61733e-5)
   
   a_ion_coeff_HI = exp(-3.271396786e1 + 1.35365560e1*th       &
            -5.73932875*th**2.0    + 1.56315498*th**3.0        &
            -2.87705600e-1*th**4.0 + 3.48255977e-2*th**5.0     &
            -2.63197617e-3*th**6.0 + 1.11954395e-4*th**7.0     &
            -2.03914985e-6*th**8.0)
   
   end subroutine ion_coeff_HI
   
   !--------------!
   
   ! Collisional ionization rate for HeI
   subroutine ion_coeff_HeI(T,a_ion_coeff_HeI)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng) :: th
   real*8, dimension(1-Ng:N+Ng), intent(out) :: a_ion_coeff_HeI

   th = log(T*8.61733e-5)
   
   a_ion_coeff_HeI = exp(-4.409864886e1 + 2.391596563e1*th       &
         -1.07532302e1*th**2.0 + 3.05803875*th**3.0           &
         -5.6851189e-1*th**4.0 + 6.79539123e-2*th**5.0        &
         -5.0090561e-3*th**6.0 + 2.06723616e-4*th**7.0        &
         -3.64916141e-6*th**8.0)
   
   end subroutine ion_coeff_HeI
         
   !--------------!
   
   ! Collisional ionization rate for HeII
   subroutine ion_coeff_HeII(T,a_ion_coeff_HeII)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng) :: xl
   real*8, dimension(1-Ng:N+Ng), intent(out) :: a_ion_coeff_HeII

   xl = 2.0*631515.0/T
   
   a_ion_coeff_HeII =  19.95*exp(-xl/2.0)*T**(-1.5)*                &
                     xl**(-1.089)/(1.0+(xl/0.553)**0.735)**1.275
   
   end subroutine ion_coeff_HeII
   
   !---------------------------------------------------!
   
   !--- Bremmstrahllung ---!
   
   ! Gaunt factor function
   subroutine GF(T,Z,GF_out)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, intent(in) :: Z
   real*8 :: T_crit
   real*8, dimension(1-Ng:N+Ng), intent(out) :: GF_out
   
   T_crit = Z*Z*3.2e5
   GF_out = 0.79464 + 0.1243*log10(T)     
   where (T .gt. T_crit) GF_out = 2.13164 + 0.1243*log10(T)      

   end subroutine GF
   
   !---------------------------------------------------!
   
   !--- Collisional excitation ---!
   
   ! Collisional excitation rate for HI
   subroutine coex_rate_HI(T,coeff_coex_rate_HI)
   real*8, dimension(1-Ng:N+Ng), intent(in)  :: T
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_coex_rate_HI

   coeff_coex_rate_HI = 7.5e-19/(1.0+sqrt(T/1.0e5))*exp(-118348.0/T)
   
   end subroutine coex_rate_HI
   
   !--------------!
   
   ! Collisional excitation rate for HeI
   subroutine coex_rate_HeI(T,coeff_coex_rate_HeI)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_coex_rate_HeI

   coeff_coex_rate_HeI = 1.1e-19*T**0.082*exp(-2.3e5/T)
   
   end subroutine coex_rate_HeI
   
   !--------------!
   
   ! Collisional excitation rate for HeII
   subroutine coex_rate_HeII(T,coeff_coex_rate_HeII)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_coex_rate_HeII
   
   coeff_coex_rate_HeII = 5.54e-17*T**(-0.397)/(1.0+sqrt(T/1.0e5))     &
                           *exp(-473638.0/T)
   
   end subroutine coex_rate_HeII
   
   !---------------------------------------------------!
   
   ! Collisional excitation from HeI(1S) to HeI(23S) [q_13]
   subroutine coex_HeI_1S_23S(T,coeff_coex_HeI_1S_23S)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8 :: a,b,c,d
   real*8, dimension(1-Ng:N+Ng) :: ups
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_coex_HeI_1S_23S

   ! Collision strengths from fit of Bray (2000)
   a =  0.06703
   b = -2.673e-6
   c = -0.03584
   d = -3.880e-4
   ups = a*exp(b*T) + c*exp(d*T)

   ! Value of rate coefficient (from Oklopcic (2018) & Lampon (2020))
   coeff_coex_HeI_1S_23S = 2.10e-8*sqrt(13.60/(kb_eV*T))		&
                        *exp(-19.81/(kb_eV*T))		&
                        *ups
   end subroutine coex_HeI_1S_23S
   
   !--------------!
   
   ! Collisional excitation from HeI(23S) to HeI(21S) [q_31a]
   subroutine coex_HeI_23S_21S(T,coeff_coex_HeI_23S_21S)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8 :: a,b,c,d
   real*8, dimension(1-Ng:N+Ng) :: ups
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_coex_HeI_23S_21S

   ! Collision strengths from fit of Bray (2000) 
   a =  2.847
   b = -1.252e-5 
   c = -1.953
   d = -3.558e-4
   ups = a*exp(b*T) + c*exp(d*T)
   
   ! Value of rate coefficient (from Oklopcic (2018) & Lampon (2020))
   coeff_coex_HeI_23S_21S = 2.10e-8*sqrt(13.60/(kb_eV*T))	&
                        *exp(-0.80/(kb_eV*T))		&
                        *ups/3.0
   
   end subroutine coex_HeI_23S_21S
   
   !--------------!
   
   ! Collisional excitation from HeI(23S) to HeI(21P) [q_31b]
   subroutine coex_HeI_23S_21P(T,coeff_coex_HeI_23S_21P)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8 :: a,b,c,d
   real*8, dimension(1-Ng:N+Ng) :: ups
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_coex_HeI_23S_21P

   ! Collision strengths from fit of Bray (2000) 
   a =  1.185
   b = -4.749e-6 
   c = -0.9131
   d = -1.669e-4
   ups = a*exp(b*T) + c*exp(d*T)
   
   ! Value of rate coefficient (from Oklopcic (2018) & Lampon (2020))
   coeff_coex_HeI_23S_21P = 2.10e-8*sqrt(13.60/(kb_eV*T))	&
                        *exp(-1.40/(kb_eV*T))		&
                        *ups/3.0
   
   end subroutine coex_HeI_23S_21P
   
   !--------------!
   
   ! Recombination coefficient of HeII on 23S state
   subroutine rec_HeII_23S(T,coeff_rec_HeII_23S)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_rec_HeII_23S

   coeff_rec_HeII_23S = 2.10e-13*(T/1.0e4)**(-0.778)
   
   end subroutine rec_HeII_23S

   !--------------!
   
   ! Recombination coefficient of HeII on 23S state
   subroutine rec_HeII_11S(T,coeff_rec_HeII_11S)
   real*8, dimension(1-Ng:N+Ng), intent(in) :: T
   real*8, dimension(1-Ng:N+Ng), intent(out) :: coeff_rec_HeII_11S

   coeff_rec_HeII_11S = 1.54e-13*(T/1.0e4)**(-0.486)
   
   end subroutine rec_HeII_11S
   
   ! ---------------------------------------------------------------- !

   ! In the following there are some copies of the previous routines
   !   coded as functions - to use in calls in T_equation

   ! Case B recombination coefficient of HeII
   double precision function rec_HeII_B_func(T)
   real*8, intent(in) :: T
   real*8 :: xl
   
   xl         = 2.0*285335.0/T
   rec_HeII_B_func = 1.26e-14*xl**0.750
   
   end function rec_HeII_B_func

   !--------------!

   ! Recombination cooling rate for HII
   double precision function rec_cool_HII_func(T)
   real*8, intent(in) :: T
   real*8 :: xl
   
   xl = 2.0*157807.0/T
   rec_cool_HII_func = 3.435e-30*T*xl**1.970/           &
                  (1.0+(xl/2.250)**0.376)**3.720
   
   end function rec_cool_HII_func
   
   !--------------!
   
   ! Recombination cooling rate for HeII
   double precision function rec_cool_HeII_func(T)
   real*8, intent(in) :: T
         
   rec_cool_HeII_func = 1.38e-16*T*rec_HeII_B_func(T)
   
   end function rec_cool_HeII_func
   
   !--------------!
   
   ! Recombination cooling rate for HeIII
   double precision function rec_cool_HeIII_func(T)
   real*8, intent(in) :: T
   real*8 :: xl
   
   xl = 2.0*631515.0/T
   rec_cool_HeIII_func = 8.0*3.435e-30*T*xl**1.970/      &
                     (1.0+(xl/2.250)**0.376)**3.720
   
   end function rec_cool_HeIII_func
   
   !--------------!

   ! Collisional ionization rate for HI
   double precision function ion_coeff_HI_func(T)
   real*8, intent(in) :: T
   real*8 :: th
   
   th = log(T*8.61733e-5)
   
   ion_coeff_HI_func = exp(-3.271396786e1 + 1.35365560e1*th  &
            -5.73932875*th**2.0    + 1.56315498*th**3.0        &
            -2.87705600e-1*th**4.0 + 3.48255977e-2*th**5.0     &
            -2.63197617e-3*th**6.0 + 1.11954395e-4*th**7.0     &
            -2.03914985e-6*th**8.0)
   
   end function ion_coeff_HI_func
   
   !--------------!
   
   ! Collisional ionization rate for HeI
   double precision function ion_coeff_HeI_func(T)
   real*8, intent(in) :: T
   real*8 :: th
   
   th = log(T*8.61733e-5)
   
   ion_coeff_HeI_func = exp(-4.409864886e1 + 2.391596563e1*th  &
         -1.07532302e1*th**2.0 + 3.05803875*th**3.0           &
         -5.6851189e-1*th**4.0 + 6.79539123e-2*th**5.0        &
         -5.0090561e-3*th**6.0 + 2.06723616e-4*th**7.0        &
         -3.64916141e-6*th**8.0)
   
   end function ion_coeff_HeI_func
         
   !--------------!
   
   ! Collisional ionization rate for HeII
   double precision function ion_coeff_HeII_func(T)
   real*8, intent(in) :: T
   real*8 :: xl
   
   xl = 2.0*631515.0/T
   
   ion_coeff_HeII_func = 19.95*exp(-xl/2.0)*T**(-1.5)*                &
                        xl**(-1.089)/(1.0+(xl/0.553)**0.735)**1.275
   
   end function ion_coeff_HeII_func

   !--------------!
   
   ! Gaunt factor function
   double precision function GF_func(T,Z)
   real*8,intent(in) :: T,Z
   
   if (T.lt.(Z*Z*3.2e5) ) then
         
         GF_func = 0.79464 + 0.1243*log10(T)     
   else
   
         GF_func = 2.13164 + 0.1243*log10(T)      
   endif      
   
   end function GF_func

   !--------------!

   ! Collisional excitation rate for HI
   double precision function coex_rate_HI_func(T)
   real*8, intent(in) :: T
   
   coex_rate_HI_func = 7.5e-19/(1.0+sqrt(T/1.0e5))*exp(-118348.0/T)
   
   end function coex_rate_HI_func
   
   !--------------!
   
   ! Collisional excitation rate for HeI
   double precision function coex_rate_HeI_func(T)
   real*8, intent(in) :: T
   
   coex_rate_HeI_func = 1.1e-19*T**0.082*exp(-2.3e5/T)
   
   end function coex_rate_HeI_func
   
   !--------------!
   
   ! Collisional excitation rate for HeII
   double precision function coex_rate_HeII_func(T)
   real*8, intent(in) :: T
   
   coex_rate_HeII_func = 5.54e-17*T**(-0.397)/(1.0+sqrt(T/1.0e5))     &
                           *exp(-473638.0/T)
   
   end function coex_rate_HeII_func
   
   ! End of module
   end module Cooling_Coefficients
