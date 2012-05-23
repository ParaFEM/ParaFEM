c  ***********************************************************************
c  *                                                                     *
c  *                             Function cpplu                          *
c  *                                                                     *
c  ***********************************************************************
c  Single Precision/Complex Version 1.01
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  to compute the complex determinant of a complex square matrix
c           A through its LU decomposition employing scaling and partial
c           pivoting. The LU decomposition of the row-permuted version of
c           [A] is stored directly in A.
c
c  This function reduces a matrix A into its (Doolittle) LU decomposition
c  `in place',
c
c      [L][U] = [P][A]    ===>   [A] <== [L\U]
c
c  Scaling and partial-pivoting are performed during the decomposition
c  to improve accuracy. The matrix [P] is the resulting permutation matrix
c  associated with the row interchanges. As a by-product of the decomposition,
c  the determinant is returned in the function name. If the matrix [A] is
c  found to be algorithmically singular, then the determinant is set to
c  zero and the function exited.
c  The matrix [L] is assumed to be unit lower triangular (ie., having 1's
c  on the diagonal) -- thus only the subdiagonal terms of [L] are stored
c  in [A]. Arguments to the routine are as follows;
c
c     A    complex array of size at least n x n containing on input the complex
c          matrix coefficients. On ouput, A will contain the LU
c          decomposition of the row permuted version of A. (input/output)
c
c    ia    leading dimension of A exactly as specified in the calling routine.
c          (input)
c
c     n    integer giving the size of the matrix A. (input)
c
c  indx    integer vector of length at least n which on output will contain
c          the indices of the rows in the row-permuted version of [A]. Note
c          that row interchanges are not actually performed, only the row
c          indices are swapped. Thus indx(2) = 3 implies that row 3 of [A]
c          would occupy the 2nd row of [P][A]. (output)
c
c scale    temporary real vector of length at least n which is used to store
c          the scaling factors used for each row. (output)
c
c  REVISION HISTORY:
c  1.01	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c-----------------------------------------------------------------------------
      complex function cpplu( A, ia, n, indx, scale )
      complex A(ia,*)
      dimension scale(*)
      integer indx(*)
      data zero/0.0/, one/1.0/
c					initialize cpplu
      cpplu   = cmplx( zero, zero )

      do 20 i = 1, n
c						set indx(i) = 1, 2, ..., n
         indx(i) = i
c						find scale factors
         smax = cabs( A(i,1) )
         do 10 j = 2, n
            s = cabs( A(i,j) )
            if( s .gt. smax ) smax = s
  10     continue
         if( smax .eq. zero ) return
         scale(i) = one/smax
  20  continue
c					forward reduce A
      sn = float(n)
      cpplu   = cmplx( one, zero )
      do 60 i = 1, n - 1
         ii = indx(i)
c						find maximum pivot element
         k = i
         smax = scale(ii)*cabs( A(ii,i) )
         do 30 j = i+1, n
            jj = indx(j)
            s  = scale(jj)*cabs( A(jj,i) )
            if( s .gt. smax ) then
               smax = s
               k = j
            endif
  30     continue
c						check for alg. singularity
         if( (smax + sn) .eq. sn ) then
            cpplu = cmplx( zero, zero )
            return
         endif
c						swap rows (ie swap indices)
         if( k .ne. i ) then
            j       = indx(i)
            indx(i) = indx(k)
            indx(k) = j
            ii      = indx(i)
            cpplu  = -cpplu
         endif
c						update determinant
         cpplu = cpplu*A(ii,i)
c						and eliminate (subtract rows)
         do 50 j = i+1, n
            jj  = indx(j)
            A(jj,i)  = A(jj,i)/A(ii,i)
            do 40 k = i+1, n
               A(jj,k) = A(jj,k) - A(jj,i)*A(ii,k)
  40        continue
  50     continue
  60  continue
c					include last diagonal element
      nn = indx(n)
      cpplu = cpplu*A(nn,n)

      return
      end
