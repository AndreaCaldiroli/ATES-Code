      subroutine r1updt(m,n,s,ls,u,v,w,sing)
      integer m,n,ls
      logical sing
      double precision s(ls),u(m),v(n),w(m)

      integer i,j,jj,l,nmj,nm1
      double precision cos,cotan,giant,one,p5,p25,sin,tan,tau,temp, &
                       zero
      double precision dpmpar
      data one,p5,p25,zero /1.0d0,5.0d-1,2.5d-1,0.0d0/

      giant = dpmpar(3)

      jj = (n*(2*m - n + 1))/2 - (m - n)

      l = jj
      do 10 i = n, m
         w(i) = s(l)
         l = l + 1
   10    continue

      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 nmj = 1, nm1
         j = n - nmj
         jj = jj - (m - j + 1)
         w(j) = zero
         if (v(j) .eq. zero) go to 50

         if (dabs(v(n)) .ge. dabs(v(j))) go to 20
            cotan = v(n)/v(j)
            sin = p5/dsqrt(p25+p25*cotan**2)
            cos = sin*cotan
            tau = one
            if (dabs(cos)*giant .gt. one) tau = one/cos
            go to 30
   20    continue
            tan = v(j)/v(n)
            cos = p5/dsqrt(p25+p25*tan**2)
            sin = cos*tan
            tau = sin
   30    continue

         v(n) = sin*v(j) + cos*v(n)
         v(j) = tau

         l = jj
         do 40 i = j, m
            temp = cos*s(l) - sin*w(i)
            w(i) = sin*s(l) + cos*w(i)
            s(l) = temp
            l = l + 1
   40       continue
   50    continue
   60    continue
   70 continue

      do 80 i = 1, m
         w(i) = w(i) + v(n)*u(i)
   80    continue

      sing = .false.
      if (nm1 .lt. 1) go to 140
      do 130 j = 1, nm1
         if (w(j) .eq. zero) go to 120

         if (dabs(s(jj)) .ge. dabs(w(j))) go to 90
            cotan = s(jj)/w(j)
            sin = p5/dsqrt(p25+p25*cotan**2)
            cos = sin*cotan
            tau = one
            if (dabs(cos)*giant .gt. one) tau = one/cos
            go to 100
   90    continue
            tan = w(j)/s(jj)
            cos = p5/dsqrt(p25+p25*tan**2)
            sin = cos*tan
            tau = sin
  100    continue

         l = jj
         do 110 i = j, m
            temp = cos*s(l) + sin*w(i)
            w(i) = -sin*s(l) + cos*w(i)
            s(l) = temp
            l = l + 1
  110       continue

         w(j) = tau
  120    continue

         if (s(jj) .eq. zero) sing = .true.
         jj = jj + (m - j + 1)
  130    continue
  140 continue

      l = jj
      do 150 i = n, m
         s(l) = w(i)
         l = l + 1
  150    continue
      if (s(jj) .eq. zero) sing = .true.
      return

      end
