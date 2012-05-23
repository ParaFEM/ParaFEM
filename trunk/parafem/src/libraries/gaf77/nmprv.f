c  ************************************************************************
c  *                                                                      *
c  *                           Subroutine nmprv                           *
c  *                                                                      *
c  ************************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE to improve a solution to [A]{x} = {b} by a residual correction
c	   technique. (naive version)
c
c  This routine accepts a solution {x} and attempts to improve its accuracy
c  by computing the residual {r} = {b} - [A]{x}, solving [A]{z} = {r} for
c  {z} and then correcting {x} according to
c
c      {x} = {x} + {z}
c
c  Unless {x} has no significant digits, this should result in some additional
c  significant digits in the solution. Note, provision has been made in this
c  program to iterate the above scheme more than once (see `nit' below).
c  Arguments to this routine are as follows;
c
c     A    real array of size at least n x n which contains the elements of
c          the matrix [A]. (input)
c
c  ALUD    real array of size at least n x n which contains the LU
c          decomposition of [A]. This array is passed directly to the
c          solver `nsolv' and is not interpreted here. (input)
c
c    ia    integer giving the column dimension of A and ALUD exactly as
c          specified in the calling routine. (input)
c
c     n    integer giving the size of the matrix [A]. (input)
c
c     X    real vector of length at least n which on input contains the
c          original solution to [A]{x} = {b}. On output, X will contain
c          the improved solution. (input/output)
c
c     B    real vector of length at least n which contains the RHS elements
c          of {b}. (input)
c
c     R    temporary real vector of length at least n which is used to store
c          the residual. (output)
c
c   nit    integer denoting the number of times the improvement algorithm
c          is to be repeated (there is usually not much point in doing it
c          more than once). (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      subroutine nmprv( A, ALUD, ia, n, X, B, R, nit )
      dimension A(ia,*), ALUD(ia,*), X(*), B(*), R(*)
      real*8 dsum, dble

      do 40 k = 1, nit
c					compute residual (in extended prec.)
         do 20 i = 1, n
            dsum = dble(B(i)) - dble(A(i,1))*dble(X(1))
            do 10 j = 2, n
               dsum = dsum - dble(A(i,j))*dble(X(j))
  10        continue
            R(i) = dsum
  20     continue
c					compute correction
         call nsolv( ALUD, ia, n, R )
c					correct {x}
         do 30 i = 1, n
            X(i) = X(i) + R(i)
  30     continue

  40  continue

      return
      end
