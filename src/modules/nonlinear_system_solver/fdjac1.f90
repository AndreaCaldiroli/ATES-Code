      subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,&
                       wa1,wa2,params)
      integer n,ldfjac,iflag,ml,mu
      double precision epsfcn
      double precision x(n),fvec(n),fjac(ldfjac,n),wa1(n),wa2(n), &
                       params(11)
      !     **********
      
      !     subroutine fdjac1
      
      !     this subroutine computes a forward-difference approximation
      !     to the n by n jacobian matrix associated with a specified
      !     problem of n functions in n variables. if the jacobian has
      !     a banded form, then function evaluations are saved by only
      !     approximating the nonzero terms.
      
      !     the subroutine statement is
      
      !       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,
      !                         wa1,wa2)
      
      !     where
      
      !       fcn is the name of the user-supplied subroutine which
      !         calculates the functions. fcn must be declared
      !         in an external statement in the user calling
      !         program, and should be written as follows.
      
      !         subroutine fcn(n,x,fvec,iflag)
      !         integer n,iflag
      !         double precision x(n),fvec(n)
      !         ----------
      !         calculate the functions at x and
      !         return this vector in fvec.
      !         ----------
      !         return
      !         end
      
      !         the value of iflag should not be changed by fcn unless
      !         the user wants to terminate execution of fdjac1.
      !         in this case set iflag to a negative integer.
      
      !       n is a positive integer input variable set to the number
      !         of functions and variables.
      
      !       x is an input array of length n.
      
      !       fvec is an input array of length n which must contain the
      !         functions evaluated at x.
      
      !       fjac is an output n by n array which contains the
      !         approximation to the jacobian matrix evaluated at x.
      
      !       ldfjac is a positive integer input variable not less than n
      !         which specifies the leading dimension of the array fjac.
      
      !       iflag is an integer variable which can be used to terminate
      !         the execution of fdjac1. see description of fcn.
      
      !       ml is a nonnegative integer input variable which specifies
      !         the number of subdiagonals within the band of the
      !         jacobian matrix. if the jacobian is not banded, set
      !         ml to at least n - 1.
      
      !       epsfcn is an input variable used in determining a suitable
      !         step length for the forward-difference approximation. this
      !         approximation assumes that the relative errors in the
      !         functions are of the order of epsfcn. if epsfcn is less
      !         than the machine precision, it is assumed that the relative
      !         errors in the functions are of the order of the machine
      !         precision.
      
      !       mu is a nonnegative integer input variable which specifies
      !         the number of superdiagonals within the band of the
      !         jacobian matrix. if the jacobian is not banded, set
      !         mu to at least n - 1.
      
      !       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at
      !         least n, then the jacobian is considered dense, and wa2 is
      !         not referenced.
      
      !     subprograms called
      
      !       minpack-supplied ... dpmpar
      
      !       fortran-supplied ... dabs,dmax1,dsqrt
      
      !     argonne national laboratory. minpack project. march 1980.
      !     burton s. garbow, kenneth e. hillstrom, jorge j. more
      
      !     **********
      integer i,j,k,msum
      double precision eps,epsmch,h,temp,zero
      double precision dpmpar
      data zero /0.0d0/
 
      !epsmch is the machine precision.
 
      epsmch = dpmpar(1)

      eps = dsqrt(dmax1(epsfcn,epsmch))
      msum = ml + mu + 1
      if (msum .lt. n) go to 40
 
         !computation of dense approximate jacobian.
 
         do 20 j = 1, n
            temp = x(j)
            h = eps*dabs(temp)
            if (h .eq. zero) h = eps
            x(j) = temp + h
            call fcn(n,x,wa1,iflag,params)
            if (iflag .lt. 0) go to 30
            x(j) = temp
            do 10 i = 1, n
               fjac(i,j) = (wa1(i) - fvec(i))/h
   10          continue
   20       continue
   30    continue
         go to 110
   40 continue
 
         !computation of banded approximate jacobian.
 
         do 90 k = 1, msum
            do 60 j = k, n, msum
               wa2(j) = x(j)
               h = eps*dabs(wa2(j))
               if (h .eq. zero) h = eps
               x(j) = wa2(j) + h
   60          continue
            call fcn(n,x,wa1,iflag,params)
            if (iflag .lt. 0) go to 100
            do 80 j = k, n, msum
               x(j) = wa2(j)
               h = eps*dabs(wa2(j))
               if (h .eq. zero) h = eps
               do 70 i = 1, n
                  fjac(i,j) = zero
                  if (i .ge. j - mu .and. i .le. j + ml) &
                     fjac(i,j) = (wa1(i) - fvec(i))/h
   70             continue
   80          continue
   90       continue
  100    continue
  110 continue
      return
 
      !last card of subroutine fdjac1.
 
      end

