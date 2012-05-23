c  ********************************************************************
c  *                                                                  *
c  *                         Subroutine qrgde                         *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 1.11
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   to compute the QR decomposition of a matrix A using Given's
c            rotations
c
c  This routine accepts an arbitrary n x m matrix [A] and reduces it to
c  upper triangular form using Givens transformations. The n x n
c  orthogonal transformation matrix [Q] which rotates the upper triangular
c  matrix [R] (n x m) back into [A] is produced along with [R].
c  [Q] and [R] are defined by the relationship
c
c     [A] = [Q][R]     <== or ==>    [Q^t][A] = [R]
c
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
c     Q    real array of size at least n x n which on output will contain
c          the elements of the orthogonal transformation matrix required to
c          rotate [A] into [R], ie [Q^t][A] = [R]. (output)
c
c    iq    leading dimension of Q exactly as specified in the calling routine.
c          (input)
c
c     n    column dimension of the matrix [A]. (input)
c
c     m    row dimension of the matrix [A]. (input)
c
c  REVISION HISTORY:
c  1.1	eliminated unused variable `two' (Dec 5/96)
c  1.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine qrgde( A, ia, Q, iq, n, m )
      dimension A(ia,*), Q(iq,*)
      data zero/0.0/, one/1.0/

c				initialize [Q]
      do 20 j = 1, n
         do 10 i = 1, j-1
            Q(i,j) = zero
            Q(j,i) = zero
  10     continue
         Q(j,j) = one
  20  continue
c				for each column of [A] ...
      do 50 k = 1, m
c					for each sub-diagonal element
         do 50 i = n, k+1, -1
c						skip transform if already 0
            if( A(i,k) .ne. zero ) then
c						find s & c using (i-1) and (i)
               if( abs(A(i,k)) .gt. abs(A(i-1,k)) ) then
                  t = A(i-1,k)/A(i,k)
                  s = one/sqrt(one + t*t)
                  c = t*s
               else
                  t = A(i,k)/A(i-1,k)
                  c = one/sqrt(one + t*t)
                  s = t*c
               endif
c						rotate affected rows of [A]
               do 30 j = k, m
                  tmp      = A(i-1,j)
                  A(i-1,j) = c*A(i-1,j) + s*A(i,j)
                  A(i,j)   = -s*tmp + c*A(i,j)
  30           continue
c						rotate affected columns of [Q]
               do 40 j = 1, n
                  tmp      = Q(j,i-1)
                  Q(j,i-1) = c*Q(j,i-1) + s*Q(j,i)
                  Q(j,i)   = -s*tmp + c*Q(j,i)
  40           continue
            endif
  50  continue

      return
      end
