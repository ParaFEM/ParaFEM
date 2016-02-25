c  ***********************************************************************
c  *                                                                     *
c  *                            Subroutine dspgss                        *
c  *                                                                     *
c  ***********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute the solution to the system of equations Ax = b through
c           Gaussian Elimination employing scaling and partial pivoting.
c
c  This routine estimates the solution to the set of simultaneous equations
c
c     [A]{x} = {b}
c
c  by reducing [A] to an upper triangular form through forward elimination,
c  then by back-substitution to obtain {x}. The solution vector {x} is
c  written over the input vector {b} (in fact {x} never appears). The
c  algorithm employs partial pivoting in which the pivots are chosen to
c  be the maximum scaled element in each sub-diagonal. Rows are not
c  actually swapped, only the row indices are swapped as recorded by
c  the integer vector `indx'. The solution in {b} is re-ordered before
c  returning so that b(1) corresponds to x(1), b(2) to x(2) etc.
c  Arguments to the routine are as follows;
c
c     A    real array of size at least n x n containing on input the
c          equation coefficients. On ouput, A will contain the LU
c          decomposition of the row permuted version of A. (input/output)
c
c    ia    leading dimension of A exactly as specified in the calling routine.
c          (input)
c
c     n    integer giving the number of equations in the problem (ie, the
c          number of rows in A). (input)
c
c     b    real vector of length at least n which on input contains the
c          RHS vector. On output, b will contain the solution. (input/output)
c
c     x    real vector of length at least n which on output will contain
c          the solution to [A]{x} = {b}. (output)
c
c  indx    integer vector of length at least n which on output will contain
c          the row indices of the row-permuted version of [A]. That is, if
c          indx(2) = 3, then row 3 was used in place of row 2 during the
c          Gaussian elimination (the rows are not actually swapped, just
c          the row indices). (output)
c
c  ierr    integer flag that is set to 0 if all goes well. ierr is set to -1
c          if the matrix is found to be algorithmically singular.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine dspgss( A, ia, n, b, x, indx, ierr )
      implicit real*8 (a-h,o-z)
      dimension A(ia,*), b(*), x(*)
      integer indx(*)
      data zero/0.d0/, one/1.d0/
      abs(xx) = dabs(xx)
c					assume all will not go well
      ierr = -1
c					initialize
      do 20 i = 1, n
c						set indx(i) = 1, 2, ..., n
         indx(i) = i
c						find scale factors  ===> x
         smax = abs( A(i,1) )
         do 10 j = 2, n
            s = abs( A(i,j) )
            if( s .gt. smax ) smax = s
  10     continue
         if( smax .eq. zero ) return
         x(i) = one/smax
  20  continue
c					forward reduce A and b
      do 60 i = 1, n - 1
         ii = indx(i)
c						find maximum scaled pivot
         k = i
         smax = x(ii)*abs( A(ii,i) )
         do 30 j = i+1, n
            jj = indx(j)
            s  = x(jj)*abs( A(jj,i) )
            if( s .gt. smax ) then
               smax = s
               k = j
            endif
  30     continue
c						swap rows (ie swap indices)
         if( k .ne. i ) then
            j       = indx(i)
            indx(i) = indx(k)
            indx(k) = j
            ii      = indx(i)
         endif
c						check pivotal element
         if( A(ii,i) .eq. zero ) return
c						and eliminate (subtract rows)
         do 50 j = i+1, n
            jj  = indx(j)
            A(jj,i)  = A(jj,i)/A(ii,i)
            do 40 k = i+1, n
               A(jj,k) = A(jj,k) - A(jj,i)*A(ii,k)
  40        continue
            b(jj) = b(jj) - A(jj,i)*b(ii)
  50     continue
  60  continue
c					check last diagonal element
      nn = indx(n)
      if( A(nn,n) .eq. zero ) return
c					now backsubstitute for the solution
c					(use x properly now)
      x(n) = b(nn)/A(nn,n)
      do 80 i = n-1, 1, -1
         ii = indx(i)
         s = b(ii) - A(ii,n)*x(n)
         do 70 j = i+1, n-1
            s  = s - A(ii,j)*x(j)
  70     continue
         x(i) = s/A(ii,i)
  80  continue
c					all went well, set flag and return
      ierr = 0
      return
      end
