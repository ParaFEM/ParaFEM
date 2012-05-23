c  ***********************************************************************
c  *                                                                     *
c  *                         Subroutine  dgjinv                          *
c  *                                                                     *
c  ***********************************************************************
c  Double Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Nov 1, 1996
c
c  PURPOSE  to compute the inverse of a square matrix A through Gauss-Jordan
c           transformations employing row scaling and partial pivoting.
c
c  A square matrix [A] is taken and reduced through Gauss-Jordan
c  transformations to the identity matrix [I]. The columns of the Gauss-Jordan
c  tranforms are the columns of [A-inv] and are stored in place in the original
c  matrix [A]. Since partial pivoting is performed, the inverse obtained is
c  actually the inverse of [P][A], where [P] is a row permutation matrix. Thus
c  the true inverse of [A] is given by
c
c           [A-inv] = [[PA]-inv][P]
c
c  that is, by permuting the columns of [PA]-inv. This final permutation is
c  optionally performed by this routine to yield [A-inv] (see LSWAP).
c  Arguments to the routine are as follows;
c
c     A    real array of size at least n x n containing on input the
c          matrix coefficients. On ouput, A will contain its inverse
c          or the inverse of a row permuted version of A. (input/output)
c
c    ia    leading dimension of A exactly as specified in the calling routine.
c          (input)
c
c     n    integer giving the size of the matrix A. (input)
c
c  indx    integer vector of length at least n which on output will contain
c          the indices of the row changes preformed in the row-permuted
c          version of [A]. Note that row interchanges are not actually
c          performed, only the row indices are swapped. (output)
c
c   idx    integer vector of length at least 2*n used to record the
c          individual row swap indices in the order performed. This
c          addition storage is required in order to unswap the final
c          inverse (not needed if LSWAP is false). (output)
c
c scale    temporary real vector of length at least n which is used to store
c          the scaling factors used for each row. (output)
c
c LSWAP    logical flag which is true if the columns of [PA]-inv are to be
c          swapped prior to returning to give [A-inv] directly. (input)
c
c  ierr    integer error flag which is zero if all goes well. If the matrix
c          is singular, ierr is set to -1.
c
c  REVISION HISTORY:
c  1.1	corrected unswapping algorithm, need additional idx storage (Nov 1/96)
c-----------------------------------------------------------------------------
      subroutine dgjinv( A, ia, n, indx, idx, scale, LSWAP, ierr )
      implicit real*8 (a-h,o-z)
      dimension A(ia,*), scale(*)
      integer indx(*), idx(2,*)
      logical LSWAP
      data zero/0.d0/, one/1.d0/
      abs(x) = dabs(x)
c					initialize
      ierr = -1
      do 20 k = 1, n
c						set indx(i) = 1, 2, ..., n
         indx(k) = k
c						find scale factors
         smax = abs( A(k,1) )
         do 10 j = 2, n
            s = abs( A(k,j) )
            if( s .gt. smax ) smax = s
  10     continue
         if( smax .eq. zero ) return
         scale(k) = one/smax
  20  continue
c					forward reduce A
      do 70 k = 1, n
         kk = indx(k)
c						find maximum pivot element
            j = k
            smax = scale(kk)*abs( A(kk,k) )
            do 30 i = k+1, n
               ii = indx(i)
               s  = scale(ii)*abs( A(ii,k) )
               if( s .gt. smax ) then
                  smax = s
                  j = i
               endif
  30        continue
c						swap rows (ie swap indices)
            if( j .ne. k ) then
               indx(k)  = indx(j)
               indx(j)  = kk
               kk       = indx(k)
               idx(1,k) = indx(j)
               idx(2,k) = indx(k)
            endif
c						check pivotal element
         if( A(kk,k) .eq. zero ) return
c						set up Gauss-Jordan Transform
         s = one/A(kk,k)
         t = (one - A(kk,k))*s
         do 40 i = 1, n
            ii = indx(i)
            A(ii,k) = -A(ii,k)*s
  40     continue
         A(kk,k) = t
c						update [PA]-inv
         do 50 j = 1, k-1
            akj = A(kk,j)
            do 50 i = 1, n
               ii = indx(i)
               A(ii,j) = A(ii,j) + A(ii,k)*akj
  50     continue
c						update [PA]
         do 60 j = k+1, n
            akj = A(kk,j)
            do 60 i = 1, n
               ii = indx(i)
               A(ii,j) = A(ii,j) + A(ii,k)*akj
  60     continue
         A(kk,k) = s
  70  continue
c					swap columns for true [A-inv]
      if( LSWAP ) then
         do 90 i = 1, n
            ii = idx(1,i)
            jj = idx(2,i)
            if( ii .ne. jj ) then
               do 80 j = 1, n
                  s       = A(j,ii)
                  A(j,ii) = A(j,jj)
                  A(j,jj) = s
  80           continue
            endif
  90     continue
c					now swap rows for true [A-inv]
         do 110 i = n, 1, -1
            ii = idx(1,i)
            jj = idx(2,i)
            if( ii .ne. jj ) then
               do 100 j = 1, n
                  s       = A(ii,j)
                  A(ii,j) = A(jj,j)
                  A(jj,j) = s
 100           continue
            endif
 110     continue
      endif

      ierr = 0
      return
      end
