c  ********************************************************************
c  *                                                                  *
c  *                         Subroutine dqrhde                        *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to compute the QR decomposition of a matrix A using Householder
c            reflections
c
c  This routine accepts an arbitrary n x m matrix [A] and reduces it to
c  upper triangular form using Householder transformations. The n x n
c  orthogonal transformation matrix [Q] which rotates the upper triangular
c  matrix [R] (n x m) back into [A] is stored in factorized form in the
c  lower triangle of [A]. [Q] and [R] are defined by the relationship
c
c     [A] = [Q][R]     <== or ==>    [Q^t][A] = [R]
c
c  The Householder matrices [H] are defined by
c
c                   V V^t
c     [H] = [I] - 2-------
c                   V^t V
c
c  where V is a Householder vector specifically chosen to reduce the
c  sub-diagonal elements of [A] to zero a column at a time. Thus
c
c    [H_m] ... [H_1] [A] = [R]
c
c  and so we see that [Q^t] is just the product of m Householder transformation
c  matrices.
c  Arguments to this routine are as follows;
c
c     A    real array of size at least n x m which on input contains the
c          elements of the array to be factorized. On output, A will
c          contain the elements of [R] in its upper triangle and the
c          essential components of each Householder vector V in its
c          columns below the diagonal (this is the factorized form of Q).
c          (input/output)
c
c    ia    leading dimension of A exactly as specified in the calling
c          routine. (input)
c
c     n    column dimension of the matrix [A]. (input)
c
c     m    row dimension of the matrix [A]. (input)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine dqrhde( A, ia, n, m )
      implicit real*8 (a-h,o-z)
      dimension A(ia,*)
      data one/1.d0/, two/2.d0/, small/1.d-40/
      sign( x, y ) = dsign( x, y )

c				for each column of [A] ...
      do 60 k = 1, m
c					calculate V and store in A
         s  = A(k,k)*A(k,k)
         do 10 i = k+1, n
            s = s + A(i,k)*A(i,k)
  10     continue
         s = sqrt(s)
c					skip the transform if column is zero
         if( s .lt. small ) go to 60
         t = sign( s, A(k,k) )
         r = one/(t + A(k,k))
         A(k,k) = -t
c					calculate V and 2/(V^t*V)
         s = one
         do 20 i = k+1, n
            A(i,k) = r*A(i,k)
            s = s + A(i,k)*A(i,k)
  20     continue
         s = -two/s
c					update k'th column of [A]
         do 50 j = k+1, m
            r = A(k,j)
            do 30 i = k+1, n
               r = r + A(i,j)*A(i,k)
  30        continue
            w = s*r
            A(k,j) = A(k,j) + w
            do 40 i = k+1, n
               A(i,j) = A(i,j) + w*A(i,k)
  40        continue
  50    continue
  60  continue

      return
      end
