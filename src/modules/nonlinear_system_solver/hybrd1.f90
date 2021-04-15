      subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa,params)
      integer n,info,lwa
      double precision tol
      double precision x(n),fvec(n),wa(lwa),params(14)
      external fcn
      !     **********
      
      !     subroutine hybrd1
      
      !     the purpose of hybrd1 is to find a zero of a system of
      !     n nonlinear functions in n variables by a modification
      !     of the powell hybrid method. this is done by using the
      !     more general nonlinear equation solver hybrd. the user
      !     must provide a subroutine which calculates the functions.
      !     the jacobian is then calculated by a forward-difference
      !     approximation.
      
      !     the subroutine statement is
      
      !       subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
      
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
      !         the user wants to terminate execution of hybrd1.
      !         in this case set iflag to a negative integer.
      
      !       n is a positive integer input variable set to the number
      !         of functions and variables.
      
      !       x is an array of length n. on input x must contain
      !         an initial estimate of the solution vector. on output x
      !         contains the final estimate of the solution vector.
      
      !       fvec is an output array of length n which contains
      !         the functions evaluated at the output x.
      
      !       tol is a nonnegative input variable. termination occurs
      !         when the algorithm estimates that the relative error
      !         between x and the solution is at most tol.
      
      !       info is an integer output variable. if the user has
      !         terminated execution, info is set to the (negative)
      !         value of iflag. see description of fcn. otherwise,
      !         info is set as follows.
      
      !         info = 0   improper input parameters.
      
      !         info = 1   algorithm estimates that the relative error
      !                    between x and the solution is at most tol.
      
      !         info = 2   number of calls to fcn has reached or exceeded
      !                    200*(n+1).
      
      !         info = 3   tol is too small. no further improvement in
      !                    the approximate solution x is possible.
      
      !         info = 4   iteration is not making good progress.
      
      !       wa is a work array of length lwa.
      
      !       lwa is a positive integer input variable not less than
      !         (n*(3*n+13))/2.
      
      !     subprograms called
      
      !       user-supplied ...... fcn
      
      !       minpack-supplied ... hybrd
      
      !     argonne national laboratory. minpack project. march 1980.
      !     burton s. garbow, kenneth e. hillstrom, jorge j. more
      
      !     **********
      integer index,j,lr,maxfev,ml,mode,mu,nfev,nprint
      double precision epsfcn,factor,one,xtol,zero
      data factor,one,zero /1.0d2,1.0d0,0.0d0/
      info = 0
 
      !check the input parameters for errors.
 
      if (n .le. 0 .or. tol .lt. zero .or. lwa .lt. (n*(3*n + 13))/2) &
         go to 20
 
      !call hybrd.
 
      maxfev = 200*(n + 1)
      xtol = tol
      ml = n - 1
      mu = n - 1
      epsfcn = zero
      mode = 2
      do 10 j = 1, n
         wa(j) = one
   10    continue
      nprint = 0
      lr = (n*(n + 1))/2
      index = 6*n + lr
      call hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,wa(1),mode, &
                 factor,nprint,info,nfev,wa(index+1),n,wa(6*n+1),lr, &
                 wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1),params)
      if (info .eq. 5) info = 4
   20 continue
      return
 
      !last card of subroutine hybrd1.
 
      end
