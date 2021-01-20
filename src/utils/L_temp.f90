      program L_temp
      ! Evaluate flux from luminosities

      implicit none
      
      real*8 :: LX,LEUV,FXUV,a,AU,pi
      
      ! Constants
      data AU/1.495978707e13/ ! Astronomical unit (cm)
      data pi/3.1415926536d0/ ! Pi
      
      read(*,*) LX
      read(*,*) LEUV
      read(*,*) a
      
      ! Evaluate XUV flux
      FXUV = (10.0**LX + 10.0**LEUV)/(4.0*pi*a*a*AU*AU)
      
      ! Format specification
100   format(F10.2)
      
      open(file = 'F_temp.txt',unit = 1)
      write(1,100) FXUV
      close(1)
      
      
      end program L_temp
