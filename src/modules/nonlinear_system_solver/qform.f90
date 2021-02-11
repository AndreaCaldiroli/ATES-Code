      subroutine qform(m,n,q,ldq,wa)
      integer m,n,ldq
      double precision q(ldq,m),wa(m)
      !     **********
      
      !     subroutine qform
      
      !     this subroutine proceeds from the computed qr factorization of
      !     an m by n matrix a to accumulate the m by m orthogonal matrix
      !     q from its factored form.
      
      !     the subroutine statement is
      
      !       subroutine qform(m,n,q,ldq,wa)
      
      !     where
      
      !       m is a positive integer input variable set to the number
      !         of rows of a and the order of q.
      
      !       n is a positive integer input variable set to the number
      !         of columns of a.
      
      !       q is an m by m array. on input the full lower trapezoid in
      !         the first min(m,n) columns of q contains the factored form.
      !         on output q has been accumulated into a square matrix.
      
      !       ldq is a positive integer input variable not less than m
      !         which specifies the leading dimension of the array q.
      
      !       wa is a work array of length m.
      
      !     subprograms called
      
      !       fortran-supplied ... min0
      
      !     argonne national laboratory. minpack project. march 1980.
      !     burton s. garbow, kenneth e. hillstrom, jorge j. more
      
      !     **********
      integer i,j,jm1,k,l,minmn,np1
      double precision one,sum,temp,zero
      data one,zero /1.0d0,0.0d0/
 
      !zero out upper triangle of q in the first min(m,n) columns.
 
      minmn = min0(m,n)
      if (minmn .lt. 2) go to 30
      do 20 j = 2, minmn
         jm1 = j - 1
         do 10 i = 1, jm1
            q(i,j) = zero
   10       continue
   20    continue
   30 continue
 
      !initialize remaining columns to those of the identity matrix.
 
      np1 = n + 1
      if (m .lt. np1) go to 60
      do 50 j = np1, m
         do 40 i = 1, m
            q(i,j) = zero
   40       continue
         q(j,j) = one
   50    continue
   60 continue
 
      !accumulate q from its factored form.
 
      do 120 l = 1, minmn
         k = minmn - l + 1
         do 70 i = k, m
            wa(i) = q(i,k)
            q(i,k) = zero
   70       continue
         q(k,k) = one
         if (wa(k) .eq. zero) go to 110
         do 100 j = k, m
            sum = zero
            do 80 i = k, m
               sum = sum + q(i,j)*wa(i)
   80          continue
            temp = sum/wa(k)
            do 90 i = k, m
               q(i,j) = q(i,j) - temp*wa(i)
   90          continue
  100       continue
  110    continue
  120    continue
      return
 
      !last card of subroutine qform.

      end
