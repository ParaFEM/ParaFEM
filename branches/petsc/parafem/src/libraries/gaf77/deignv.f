c  *****************************************************************
c  *                                                               *
c  *                     Subroutine Deignv                         *
c  *                                                               *
c  *****************************************************************
c  Double Precision Version 1.11
c  Written by Gordon A. Fenton, Princeton, Jan. 20, 1988.
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  computes the inverse of a matrix given its eigenvalues and
c           vectors.
c
c  Finds the inverse of a real symmetric matrix A of dimension N x N
c  given its eigenvalues and eigenvectors. Both the eigenvalues and
c  eigenvectors are destroyed in this routine and the inverse of A is
c  overwritten on its lower triangle. Eigenvalues and vectors are found
c  by calling EIGVL.
c  Variables are described as follows;
c
c      A    N x N real symmetric matrix which on output contains its
c           inverse in its lower triangle
c
c      EIG  vector of length N which on input contains the eigenvalues
c           of A. On output EIG contains its reciprocals.
c
c      Q    N x N real orthogonal matrix containing the eigenvectors of
c           A in the order dictated by the order of EIG.
c
c      N    dimension of A, EIG and Q
c
c      IERR if the inversion is successfull then ierr = 0. However if one
c           of the eigenvalues is zero then the matrix is singular and no
c           inversion is possible in which case ierr = -1.
c
c  REVISION HISTORY:
c  1.1	properly implemented ierr (checking eigenvalues). (Dec 5/96)
c  1.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine deignv( A, Q, eig, N, ierr )
      implicit real*8 (a-h,o-z)
      dimension A(n,*), Q(n,*), eig(*)
      data zero/0.d0/, one/1.d0/
c						assume failure
      ierr = -1

      s = zero
      do 10 k = 1,n
         if( eig(k) .eq. zero ) return
         eig(k) = one/eig(k)
         y = Q(1,k)
         Q(1,k) = eig(k) * Q(1,k)
         s = s + y * Q(1,k)
  10  continue

      A(1,1) = s

      do 30 i = 2,n
         s = zero
         do 20 k = 1,n
            s = s + Q(i,k) * Q(1,k)
  20     continue
         A(i,1) = s
  30  continue

      do 70 j = 2,n
         s = zero
         do 40 k = 1,n
            y = Q(j,k)
            Q(j,k) = eig(k) * Q(j,k)
            s = s + y * Q(j,k)
  40     continue
c
         A(j,j) = s
c
         do 60 i = j+1, n
            s = zero
            do 50 k = 1,n
               s = s + Q(i,k) * Q(j,k)
  50        continue
            A(i,j) = s
  60     continue
  70  continue
c
      ierr = 0
      return
      end
