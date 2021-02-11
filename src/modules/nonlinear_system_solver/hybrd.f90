      subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag, &
                       mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,&
                       qtf,wa1,wa2,wa3,wa4,params)
      integer n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr
      double precision xtol,epsfcn,factor
      double precision x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr), &
                       qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
      external fcn
      !     **********
      !           
      !     subroutine hybrd
      
      !     the purpose of hybrd is to find a zero of a system of
      !     n nonlinear functions in n variables by a modification
      !     of the powell hybrid method. the user must provide a
      !     subroutine which calculates the functions. the jacobian is
      !     then calculated by a forward-difference approximation.
      
      !     the subroutine statement is
      
      !       subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,
      !                        diag,mode,factor,nprint,info,nfev,fjac,
      !                        ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4)
      
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
      !         ---------
      !         return
      !         end
      
      !         the value of iflag should not be changed by fcn unless
      !         the user wants to terminate execution of hybrd.
      !         in this case set iflag to a negative integer.
      
      !       n is a positive integer input variable set to the number
      !         of functions and variables.
      
      !       x is an array of length n. on input x must contain
      !         an initial estimate of the solution vector. on output x
      !         contains the final estimate of the solution vector.
      
      !       fvec is an output array of length n which contains
      !         the functions evaluated at the output x.
      
      !       xtol is a nonnegative input variable. termination
      !         occurs when the relative error between two consecutive
      !         iterates is at most xtol.
      
      !       maxfev is a positive integer input variable. termination
      !         occurs when the number of calls to fcn is at least maxfev
      !         by the end of an iteration.
      
      !       ml is a nonnegative integer input variable which specifies
      !         the number of subdiagonals within the band of the
      !         jacobian matrix. if the jacobian is not banded, set
      !         ml to at least n - 1.
      
      !       mu is a nonnegative integer input variable which specifies
      !         the number of superdiagonals within the band of the
      !         jacobian matrix. if the jacobian is not banded, set
      !         mu to at least n - 1.
      
      !       epsfcn is an input variable used in determining a suitable
      !         step length for the forward-difference approximation. this
      !         approximation assumes that the relative errors in the
      !         functions are of the order of epsfcn. if epsfcn is less
      !         than the machine precision, it is assumed that the relative
      !         errors in the functions are of the order of the machine
      !         precision.
      
      !       diag is an array of length n. if mode = 1 (see
      !         below), diag is internally set. if mode = 2, diag
      !         must contain positive entries that serve as
      !         multiplicative scale factors for the variables.
      
      !       mode is an integer input variable. if mode = 1, the
      !         variables will be scaled internally. if mode = 2,
      !         the scaling is specified by the input diag. other
      !         values of mode are equivalent to mode = 1.
      
      !       factor is a positive input variable used in determining the
      !         initial step bound. this bound is set to the product of
      !         factor and the euclidean norm of diag*x if nonzero, or else
      !         to factor itself. in most cases factor should lie in the
      !         interval (.1,100.). 100. is a generally recommended value.
      
      !       nprint is an integer input variable that enables controlled
      !         printing of iterates if it is positive. in this case,
      !         fcn is called with iflag = 0 at the beginning of the first
      !         iteration and every nprint iterations thereafter and
      !         immediately prior to return, with x and fvec available
      !         for printing. if nprint is not positive, no special calls
      !         of fcn with iflag = 0 are made.
      
      !       info is an integer output variable. if the user has
      !         terminated execution, info is set to the (negative)
      !         value of iflag. see description of fcn. otherwise,
      !         info is set as follows.
      
      !         info = 0   improper input parameters.
      
      !         info = 1   relative error between two consecutive iterates
      !                    is at most xtol.
      
      !         info = 2   number of calls to fcn has reached or exceeded
      !                    maxfev.
      
      !         info = 3   xtol is too small. no further improvement in
      !                    the approximate solution x is possible.
      
      !         info = 4   iteration is not making good progress, as
      !                    measured by the improvement from the last
      !                    five jacobian evaluations.
      
      !         info = 5   iteration is not making good progress, as
      !                    measured by the improvement from the last
      !                    ten iterations.
      
      !       nfev is an integer output variable set to the number of
      !         calls to fcn.
      
      !       fjac is an output n by n array which contains the
      !         orthogonal matrix q produced by the qr factorization
      !         of the final approximate jacobian.
      
      !       ldfjac is a positive integer input variable not less than n
      !         which specifies the leading dimension of the array fjac.
      
      !       r is an output array of length lr which contains the
      !         upper triangular matrix produced by the qr factorization
      !         of the final approximate jacobian, stored rowwise.
      
      !       lr is a positive integer input variable not less than
      !         (n*(n+1))/2.
      
      !       qtf is an output array of length n which contains
      !         the vector (q transpose)*fvec.
      
      !       wa1, wa2, wa3, and wa4 are work arrays of length n.
      
      !     subprograms called
      
      !       user-supplied ...... fcn
      
      !       minpack-supplied ... dogleg,dpmpar,enorm,fdjac1,
      !                            qform,qrfac,r1mpyq,r1updt
      
      !       fortran-supplied ... dabs,dmax1,dmin1,min0,mod
      
      !     argonne national laboratory. minpack project. march 1980.
      !     burton s. garbow, kenneth e. hillstrom, jorge j. more
      
      !     **********
      integer i,iflag,iter,j,jm1,l,msum,ncfail,ncsuc,nslow1,nslow2
      integer iwa(1)
      logical jeval,sing
      double precision actred,delta,epsmch,fnorm,fnorm1,one,pnorm, &
                       prered,p1,p5,p001,p0001,ratio,sum,temp,xnorm,&
                       zero
      double precision dpmpar,enorm
      data one,p1,p5,p001,p0001,zero &
           /1.0d0,1.0d-1,5.0d-1,1.0d-3,1.0d-4,0.0d0/
 
      !epsmch is the machine precision.
 
      epsmch = dpmpar(1)

      info = 0
      iflag = 0
      nfev = 0
 
      !check the input parameters for errors.
 
      if (n .le. 0 .or. xtol .lt. zero .or. maxfev .le. 0 &
          .or. ml .lt. 0 .or. mu .lt. 0 .or. factor .le. zero &
          .or. ldfjac .lt. n .or. lr .lt. (n*(n + 1))/2) go to 300
      if (mode .ne. 2) go to 20
      do 10 j = 1, n
         if (diag(j) .le. zero) go to 300
   10    continue
   20 continue
 
      !evaluate the function at the starting point
      !and calculate its norm.
 
      iflag = 1
      call fcn(n,x,fvec,iflag,params)
      nfev = 1
      if (iflag .lt. 0) go to 300
      fnorm = enorm(n,fvec)
 
      !determine the number of calls to fcn needed to compute
      !the jacobian matrix.
 
      msum = min0(ml+mu+1,n)
 
      !initialize iteration counter and monitors.
 
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
 
      !beginning of the outer loop.
 
   30 continue
         jeval = .true.
 
         !calculate the jacobian matrix.
 
         iflag = 2
         call fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1, &
                     wa2,params)
         nfev = nfev + msum
         if (iflag .lt. 0) go to 300
 
         !compute the qr factorization of the jacobian.
 
         call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
 
         !on the first iteration and if mode is 1, scale according
         !to the norms of the columns of the initial jacobian.
 
         if (iter .ne. 1) go to 70
         if (mode .eq. 2) go to 50
         do 40 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) .eq. zero) diag(j) = one
   40       continue
   50    continue
 
         !on the first iteration, calculate the norm of the scaled x
         !and initialize the step bound delta.
 
         do 60 j = 1, n
            wa3(j) = diag(j)*x(j)
   60       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta .eq. zero) delta = factor
   70    continue
 
         !form (q transpose)*fvec and store in qtf.
 
         do 80 i = 1, n
            qtf(i) = fvec(i)
   80       continue
         do 120 j = 1, n
            if (fjac(j,j) .eq. zero) go to 110
            sum = zero
            do 90 i = j, n
               sum = sum + fjac(i,j)*qtf(i)
   90          continue
            temp = -sum/fjac(j,j)
            do 100 i = j, n
               qtf(i) = qtf(i) + fjac(i,j)*temp
  100          continue
  110       continue
  120       continue
 
         !copy the triangular factor of the qr factorization into r.
 
         sing = .false.
         do 150 j = 1, n
            l = j
            jm1 = j - 1
            if (jm1 .lt. 1) go to 140
            do 130 i = 1, jm1
               r(l) = fjac(i,j)
               l = l + n - i
  130          continue
  140       continue
            r(l) = wa1(j)
            if (wa1(j) .eq. zero) sing = .true.
  150       continue
 
         !accumulate the orthogonal factor in fjac.
 
         call qform(n,n,fjac,ldfjac,wa1)
 
         !rescale if necessary.
 
         if (mode .eq. 2) go to 170
         do 160 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  160       continue
  170    continue
 
         !beginning of the inner loop.
 
  180    continue
 
          !  if requested, call fcn to enable printing of iterates.
 
         if (nprint .le. 0) go to 190
         iflag = 0
         if (mod(iter-1,nprint) .eq. 0) call fcn(n,x,fvec,iflag,params)
         if (iflag .lt. 0) go to 300
  190    continue
 
           ! determine the direction p.
 
            call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
 
            !store the direction p and x + p. calculate the norm of p.
 
            do 200 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  200          continue
            pnorm = enorm(n,wa3)
 
            !on the first iteration, adjust the initial step bound.
 
            if (iter .eq. 1) delta = dmin1(delta,pnorm)
 
            !evaluate the function at x + p and calculate its norm.
 
            iflag = 1
            call fcn(n,wa2,wa4,iflag,params)
            nfev = nfev + 1
            if (iflag .lt. 0) go to 300
            fnorm1 = enorm(n,wa4)
 
            !compute the scaled actual reduction.
 
            actred = -one
            if (fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
 
            !compute the scaled predicted reduction.
 
            l = 1
            do 220 i = 1, n
               sum = zero
               do 210 j = i, n
                  sum = sum + r(l)*wa1(j)
                  l = l + 1
  210             continue
               wa3(i) = qtf(i) + sum
  220          continue
            temp = enorm(n,wa3)
            prered = zero
            if (temp .lt. fnorm) prered = one - (temp/fnorm)**2
 
            !compute the ratio of the actual to the predicted
            !reduction.
 
            ratio = zero
            if (prered .gt. zero) ratio = actred/prered
 
            !update the step bound.
 
            if (ratio .ge. p1) go to 230
               ncsuc = 0
               ncfail = ncfail + 1
               delta = p5*delta
               go to 240
  230       continue
               ncfail = 0
               ncsuc = ncsuc + 1
               if (ratio .ge. p5 .or. ncsuc .gt. 1) &
                  delta = dmax1(delta,pnorm/p5)
               if (dabs(ratio-one) .le. p1) delta = pnorm/p5
  240       continue
 
            !test for successful iteration.
 
            if (ratio .lt. p0001) go to 260
 
            !successful iteration. update x, fvec, and their norms.
 
            do 250 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
               fvec(j) = wa4(j)
  250          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  260       continue
 
            !determine the progress of the iteration.
 
            nslow1 = nslow1 + 1
            if (actred .ge. p001) nslow1 = 0
            if (jeval) nslow2 = nslow2 + 1
            if (actred .ge. p1) nslow2 = 0
 
            !test for convergence.
 
            if (delta .le. xtol*xnorm .or. fnorm .eq. zero) info = 1
            if (info .ne. 0) go to 300
 
            !tests for termination and stringent tolerances.
 
            if (nfev .ge. maxfev) info = 2
            if (p1*dmax1(p1*delta,pnorm) .le. epsmch*xnorm) info = 3
            if (nslow2 .eq. 5) info = 4
            if (nslow1 .eq. 10) info = 5
            if (info .ne. 0) go to 300
 
            !criterion for recalculating jacobian approximation
            !by forward differences.
 
            if (ncfail .eq. 2) go to 290
 
            !calculate the rank one modification to the jacobian
            !and update qtf if necessary.
 
            do 280 j = 1, n
               sum = zero
               do 270 i = 1, n
                  sum = sum + fjac(i,j)*wa4(i)
  270             continue
               wa2(j) = (sum - wa3(j))/pnorm
               wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
               if (ratio .ge. p0001) qtf(j) = sum
  280          continue
 
            !compute the qr factorization of the updated jacobian.
 
            call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
            call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)
            call r1mpyq(1,n,qtf,1,wa2,wa3)
 
            !end of the inner loop.
 
            jeval = .false.
            go to 180
  290    continue
 
         !end of the outer loop.
 
         go to 30
  300 continue
 
      !termination, either normal or user imposed.
 
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (nprint .gt. 0) call fcn(n,x,fvec,iflag,params)
      return
 
      !last card of subroutine hybrd.

      end
