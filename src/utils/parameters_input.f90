      program parameters_input
      ! Read input data and write to temporary file
      
      
      implicit none
      
      real*8 :: R0,Mp,T0,LX,LEUV,a_orb,pi,JXUV,Mstar, &
                Mtilde,atilde,RJ,MJ,AU,Ms
      character(len=6) :: tab = '      '
      character(len=1) :: ans
       
      ! Constant data
      data RJ/6.9911e9/       ! Jupiter radius (cm)
      data MJ/1.8982e30/      ! Jupiter mass (g)
      data AU/1.495978707e13/ ! Astronomical unit (cm)
      data Ms/1.98847e33/     ! Solar mass (g)
      data pi/3.1415926536d0/ ! Pi
      
      
      ! Read data
      write(*,*) "Planet radius (R_J):"
      read(*,*) R0
      R0 = R0 * RJ 
      
      write(*,*) "Planet mass (M_J):"
      read(*,*) Mp
      Mp = Mp * MJ 
      
      write(*,*) "Surface temperature (K):"
      read(*,*) T0
      
      write(*,*) "Log10 of X-ray luminosity:"
      read(*,*) LX

      write(*,*) "Use Sanz-Forcada (2011) scaling law? (y | n)"
      read(*,*) ans
      
      select case (ans)
            
            case ('y')
            LEUV = 4.80d0+0.860d0*LX
            write(*,*) "The corresponding EUV luminosity is", LEUV
            
            case ('n')
            write(*,*) "Insert log10 of X-ray luminosity:"
            read(*,*) LEUV
      
      end select
      
      
      write(*,*) "Orbital distance (AU):"
      read(*,*) a_orb
      a_orb = a_orb * AU
      
      JXUV = (10.0**LX+10.0**LEUV)/(4.0*pi*a_orb*a_orb)
      write(*,*) "JXUV flux is ", JXUV," erg cm^-2 s^-1"
      write(*,*) log10(JXUV)
      
      write(*,*) "Mass of parent star (M_sun):"
      read(*,*) Mstar
      Mstar = Mstar * Ms
      
      Mtilde = Mstar/Mp
      write(*,*) "The ratio Mstar/M_planet is", Mtilde
      
      atilde = a_orb/R0
      write(*,*) "The ratio a/Rp is", atilde
      
      
      ! Write parameters to temporary file
      open(1,file = "temp_input_params.txt")
      
100   format (ES13.6, A,    &
              ES13.6, A,    &
              F6.1, A,      &
              F5.2, A,      &
              F5.2, A,      &
              F6.4, A,      &
              F10.2,A,      &
              F10.2,A,      &
              F10.4)
                
      write(1,100)            &
            R0, tab,          &
            Mp, tab,          &
            T0, tab,          &
            LX, tab,          &
            LEUV, tab,        &
            a_orb/AU, tab,    &
            JXUV, tab,        & 
            Mtilde, tab,      &
            atilde
      close(1)

      
      end program parameters_input
      
