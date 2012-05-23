c  ********************************************************************
c  *                                                                  *
c  *                        Subroutine luinv                          *
c  *                                                                  *
c  ********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  compute the inverse of a matrix using its LU decomposition
c
c  This routine takes a matrix [A] and computes its inverse [A-inv] by
c  first finding its LU decomposition and then solving for each column
c  of [A-inv] sequentially. Arguments are as follows;
c
c    A    real array of size at least n x n which on input contains the
c         elements of the matrix to be inverted. On output A will contain
c         its own LU decomposition. (input/output)
c
c   AI    real array of size at least n x n which will contain the elements
c         of [A-inv]. (output)
c
c   ia    integer giving the leading dimension of A and AI as specified
c         in the calling routine. (input)
c
c  LPP    logical flag which is true if partial pivoting is to be used in
c         the LU factorization stage. (input)
c
c indx    integer vector of length at least n used to store the partial
c         pivot row indices (ignored if LPP is false). (output)
c
c    b    temporary real vector of length at least n.
c
c    n    integer giving the order of the matrix [A]. (input)
c
c ierr    integer error flag which is set to 0 if all goes well and to
c         -1 if the matrix is found to be singular. (output)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c---------------------------------------------------------------------------
      subroutine luinv( A, AI, ia, LPP, indx, b, n, ierr )
      dimension A(ia,*), AI(ia,*), b(*)
      integer indx(*)
      logical LPP
      data zero/0.0/, one/1.0/
c					use partial pivoting?
      ierr = -1
      if( LPP ) then
c						yes ...
         det = spplu( A, ia, n, indx, b )
         if( det .eq. zero ) return
         do 10 i = 1, n
            b(i) = zero
  10     continue
         do 20 j = 1, n
            b(j) = one
            call psolv( A, ia, n, indx, b, AI(1,j) )
            b(j) = zero
  20     continue
      else
c						no, use naive Gauss Elim.
         det = gselu( A, ia, n )
         if( det .eq. zero ) return
         do 30 j = 1, n
            AI(j,j) = one
            do 30 i = j+1, n
               AI(i,j) = zero
               AI(j,i) = zero
  30     continue
         do 40 j = 1, n
            call nsolv( A, ia, n, AI(1,j) )
  40     continue
      endif

      ierr = 0
      return
      end
