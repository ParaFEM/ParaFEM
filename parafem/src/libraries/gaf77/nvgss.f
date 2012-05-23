c  ***********************************************************************
c  *                                                                     *
c  *                            Subroutine nvgss                         *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute the solution to the system of equations [A]{x} = {b}
c           through naive Gaussian Elimination.
c
c  This routine estimates the solution to the set of simultaneous equations
c
c     [A]{x} = {b}
c
c  by reducing [A] to an upper triangular form through forward elimination,
c  then by back-substitution to obtain {x}. The solution vector {x} is
c  written over the input vector {b} (in fact {x} never appears).
c  Arguments to the routine are as follows;
c
c     A    real array of size at least n x n containing on input the
c          equation coefficients. On ouput, A will contain its own LU
c          decomposition. (input/output)
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
c  ierr    integer flag that is set to 0 if all goes well. ierr is set to -1
c          if the matrix is found to be algorithmically singular.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine nvgss( A, ia, n, b, ierr )
      dimension A(ia,*), b(*)
      data zero/0.0/
c					assume all will not go well
      ierr = -1
c					forward reduce A and b
      do 20 k = 1, n - 1
         if( A(k,k) .eq. zero ) return
c					first, compute L(i,k) and reduce b
         do 10 i = k + 1, n
            A(i,k) = A(i,k)/A(k,k)
            b(i)   = b(i) - A(i,k)*b(k)
  10     continue
c					second, reduce the rest of A
         do 20 j = k + 1, n
         do 20 i = k + 1, n
            A(i,j) = A(i,j) - A(i,k)*A(k,j)
  20  continue
c					check the last diagonal element
      if( A(n,n) .eq. zero ) return
c					now backsubstitute for the solution
      b(n) = b(n)/A(n,n)
      do 40 i = n-1, 1, -1
         s = b(i) - A(i,n)*b(n)
         do 30 j = i+1, n-1
            s = s - A(i,j)*b(j)
  30     continue
         b(i) = s/A(i,i)
  40  continue
c					all went well, set flag
      ierr = 0
      return
      end
