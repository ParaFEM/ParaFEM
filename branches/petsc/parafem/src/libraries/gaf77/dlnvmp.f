c  ********************************************************************
c  *                                                                  *
c  *                        Subroutine dlnvmp                         *
c  *                                                                  *
c  ********************************************************************
c  Double Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  compute the inverse of a matrix using its LU decomposition with
c           residual improvement
c
c  This routine takes a matrix [A] and computes its inverse [A-inv] by
c  first finding its LU decomposition and then solving for each column
c  of [A-inv] sequentially. Each column is improved by one pass through
c  a residual solution algorithm. Arguments are as follows;
c
c    A    real array of size at least n x n which on input contains the
c         elements of the matrix to be inverted. (input)
c
c  ALU    real array of size at least n x n which on output will contain
c         the LU decomposition of [A] or of [P][A] (see LPP), where [P]
c         is a permutation matrix. (output)
c
c   AI    real array of size at least n x n which will contain the elements
c         of [A-inv]. (output)
c
c   ia    integer giving the leading dimension of A, ALU, and AI as specified
c         in the calling routine. (input)
c
c  LPP    logical flag which is true if partial pivoting is to be used in
c         the LU factorization stage. (input)
c
c indx    integer vector of length at least n used to store the partial
c         pivot row indices (ignored if LPP is false). (output)
c
c b,r,z   temporary real vectors of length at least n.
c
c    n    integer giving the order of the matrix [A]. (input)
c
c ierr    integer error flag which is set to 0 if all goes well and to
c         -1 if the matrix is found to be singular. (output)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c---------------------------------------------------------------------------
      subroutine dlnvmp( A, ALU, AI, ia, LPP, indx, b, r, z, n, ierr )
      implicit real*8 (a-h,o-z)
      dimension A(ia,*), ALU(ia,*), AI(ia,*), b(*), r(*), z(*)
      integer indx(*)
      logical LPP
      data zero/0.d0/, one/1.d0/

      ierr = -1
c					copy A into ALU
      do 10 j = 1, n
      do 10 i = 1, n
         ALU(i,j) = A(i,j)
  10  continue
c					use partial pivoting?
      if( LPP ) then
c						yes, use partial pivoting
         det = dspplu( ALU, ia, n, indx, b )
         if( det .eq. zero ) return
         do 20 i = 1, n
            b(i) = zero
  20     continue
         do 30 j = 1, n
            b(j) = one
            call dpsolv( ALU, ia, n, indx, b, AI(1,j) )
            call dimprv( A, ALU, ia, n, indx, AI(1,j), b, r, z, 1 )
            b(j) = zero
  30     continue
      else
c						no, use naive Gauss Elim.
         det = dgselu( ALU, ia, n )
         if( det .eq. zero ) return
         do 40 j = 1, n
            AI(j,j) = one
            b(j)    = zero
            do 40 i = j+1, n
               AI(i,j) = zero
               AI(j,i) = zero
  40     continue
         do 50 j = 1, n
            b(j) = one
            call dnsolv( ALU, ia, n, AI(1,j) )
            call dnmprv( A, ALU, ia, n, AI(1,j), b, r, 1 )
            b(j) = zero
  50     continue
      endif

      ierr = 0
      return
      end
