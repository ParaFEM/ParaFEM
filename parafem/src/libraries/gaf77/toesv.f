c  **********************************************************************
c  *                                                                    *
c  *                         Subroutine toesv                           *
c  *                                                                    *
c  **********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to solve the matrix equation [A]{x} = {b} for [A] symmetric,
c           positive definite and Toeplitz.
c
c  This routine follows the algorithm of Levinson and computes a solution
c  to the Toeplitz system of equations in 4*n**2 FLOPS. Arguments to the
c  routine are as follows;
c
c     r    real vector of length at least n which contains the value of the
c          elements of each diagonal of the matrix [A] (all equal) starting
c          from the main diagonal. (input)
c
c     b    real vector of length at least n which contains the RHS vector.
c          (input)
c
c     x    real vector of length at least n which will contain the solution
c          vector. (output)
c
c     y    temporary vector of length at least n used for workspace.
c
c     n    number of equations in the system. (input)
c
c  ierr    integer flag = 0 if all goes well, = -1 if a singularity is
c          encountered.
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine toesv( r, b, x, y, n, ierr )
      dimension r(*), b(*), x(*), y(*)
      data zero/0.0/, one/1.0/
c					assume all will not go well
      ierr = -1
      if( r(1) .eq. zero ) return
c					begin recursive solution
      dd   = one/r(1)
      y(1) = -r(2)*dd
      x(1) = b(1)*dd
      bb   = one
      aa   = y(1)
      do 60 k = 1, n - 2
         bb = (one - aa*aa)*bb
         if( bb .eq. zero ) return
         db = dd/bb
         u  = b(k+1) - r(2)*x(k)
         do 20 i = 2, k
            u = u - r(i+1)*x(k-i+1)
  20     continue
         u = db*u
         do 30 i = 1, k
            x(i) = x(i) + u*y(k-i+1)
  30     continue
         x(k+1) = u

         aa = r(k+2) + r(k+1)*y(1)
         do 40 i = 2, k
            aa = aa + r(k-i+2)*y(i)
  40     continue
         aa = -db*aa
         kb2 = k/2
         do 50 i = 1, kb2
            l = k - i + 1
            v = y(i)
            y(i) = y(i) + aa*y(l)
            y(l) = y(l) + aa*v
  50     continue
         if( 2*kb2 .ne. k ) y(1+kb2) = y(1+kb2)*(one + aa)
         y(k+1) = aa

  60  continue

      bb = (one - aa*aa)*bb
      if( bb .eq. zero ) return
      u  = b(n) - r(2)*x(n-1)
      do 70 i = 2, n-1
         u = u - r(i+1)*x(n-i)
  70  continue
      u = dd*u/bb
      do 80 i = 1, n-1
         x(i) = x(i) + u*y(n-i)
  80  continue
      x(n) = u

      ierr = 0
      return
      end
