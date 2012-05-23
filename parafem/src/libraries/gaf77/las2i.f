c  *********************************************************************
c  *                                                                   *
c  *                         subroutine las2i                          *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 3.41
c  Written by Gordon A. Fenton, TUNS, Aug. 26, 1992
c  Latest Update: Jun 9, 1999
c
c  PURPOSE   initializes parameters for LAS2G
c
c  This routine sets up the matrices required by LAS2G to construct
c  realizations of the random field. The covariances between local averages
c  at each subdivision stage are computed in double precision for accuracy,
c  but the final construction matrices are stored in single precision and
c  returned to the calling routine via the argument list. The general
c  recursive field construction follows the relationship
c
c         {Z^j} = [A^T]{Z^(j-1)} + [C]{U}
c
c  where {Z^j} is a vector of length 4 representing the values assigned to
c  the 2 x 2 cell subdivision and {Z^(j-1)} are the parent cell values in
c  some neighbourhood. The following figure illustrates the cell subdivision,
c  the neighbourhood, and the numbering scheme for an interior cell (special
c  subsets of the neighbourhood are used for the corners and sides):
c
c                ----------------------------------------
c                |            |            |            |
c                |            |            |            |
c                |     7      |     8      |      9     |
c                |            |            |            |
c                |            |            |            |
c                |------------|------------|------------|
c                |            | 3   |    4 |            |
c                |            |     |      |            |
c                |     4      |-----5------|      6     |
c                |            |     |      |            |
c                |            | 1   |    2 |            |
c                |------------|------------|------------|
c                |            |            |            |
c                |            |            |            |
c                |     1      |     2      |     3      |
c                |            |            |            |
c                |            |            |            |
c                ----------------------------------------
c
c  We see that for the upper left corner, the parent cell neighbourhood used
c  consists of just cells {4,5,7,8} and similarly for the other corners and
c  sides.
c
c  The first stage of the simulation involves the direct generation of a
c  k1 x k2 cell array, where k1 and k2 are integers which satisfying the
c  decomposition N1 = k1*2**m, N2 = k2*2**m for a common factor 2**m. The
c  integers k1 and k2 are chosen to be as large as possible while requiring
c  the product k1*k2 to be less than or equal to MXK. Note that the direct
c  simulation involves the inversion of a MXK x MXK matrix (at the upper
c  limit) and so MXK should not be overly large.
c  This formulation is somewhat less restrictive than simply requiring
c  N1 and N2 to be powers of 2. Also N1 and N2 do not have to be equal.
c  However N1 and N2 cannot still be chosen arbitrarily, for example the
c  set (N1,N2) = (336,256) results in k1 = 21, k2 = 16, m = 4 which is
c  not acceptable here (since k1*k2 > MXK for MXK = 256), while the set
c  (N1,N2) = (160,256) is acceptable since k1 = 10, k2 = 16, m = 4. In general
c  it may be easier to choose k1, k2, and m before specifying N1 and N2. In
c  the event that an unacceptable (k1,k2,m) combination is selected, IERR
c  is set to -1 and control returned to the calling routine.
c  The maximum value of m is set by the calling routine in the argument MXM.
c
c  Arguments to this routine are as follows;
c
c dvarfn  external real*8 function which returns the variance of the random
c         process averaged over a given area. dvarfn is referenced as follows
c
c                var = dvarfn( V1, V2 )
c
c         where (V1,V2) are the side dimensions of the rectangular averaging
c         domain. Any other parameters to the function may be passed by
c         common block from the calling routine. Note that the variance of
c         the random process averaged over the area (V1 x V2) is the product
c         of the point variance and the traditionally defined "variance"
c         function, as discussed by Vanmarcke (pg 186).
c
c  N1,N2  number of cells to discretize the field in the x and y directions
c         respectively (corresponding to the first and second indices of Z
c         respectively). Both N1 and N2 must have the form N1 = k1 * 2**m
c         and N2 = k2 * 2**m where m is common to both and k1 and k2 are
c         positive integers satisfying k1*k2 <= MXK. Generally k1 and k2
c         are chosen to be as large as possible and still satisfy the above
c         requirements so the the first stage involves directly simulating
c         a k1 x k2 field by inversion of a covariance matrix.
c         A potential example is (N1,N2) = (160,256) which gives k1 = 5,
c         k2 = 8, and m = 5. Note that in general N1 and N2 cannot be chosen
c         arbitrarily - it is usually best to choose m first then
c         k1 and k2 so as to satisfy or exceed the problem requirements. Note
c         that because of the requirements on k1*k2, N1 cannot be more than
c         MXK times as big (or smaller) than N2. (input)
c
c  XL,YL  physical dimensions of the process. (input)
c
c  KSEED  integer seed to be used for the pseudo-random number generator.
c         If KSEED = 0, then a random seed will be used (based on the
c         process ID of the current program -- see iseed.f).
c         On output, KSEED is set to the value of the actual seed used.
c
c    MXM  integer giving the largest value that M can take. An error is
c         generated if the process size is such that M > MXM. (input)
c
c    C0   real vector containing the upper triangular values of the Cholesky
c         decomposition of the covariance matrix for the initial stage (0) of
c         k1 x k2 cells. (output)
c
c    CT   real vector containing the upper triangular values of the Cholesky
c         decomposition of the covariance matrix for the k1 or k2 = 1 special
c         case. (output)
c
c    CC   real vector containing the upper triangular values of the Cholesky
c         decomposition of the covariance matrix for the corner cell 2 x 2
c         subdivisions. (output)
c
c    CS   real vector containing the upper triangular values of the Cholesky
c         decomposition of the covariance matrix for the side cell 2 x 2
c         subdivisions. (output)
c
c    CI   real vector containing the upper triangular values of the Cholesky
c         decomposition of the covariance matrix for the interior cell 2 x 2
c         subdivisions. (output)
c
c    AT   real array containing the best linear estimation coefficients for
c         the k1 or k2 = 1 special case. (output)
c
c    AC   real array containing the best linear estimation coefficients for
c         the corner cell subdivisions. (output)
c
c    AS   real array containing the best linear estimation coefficients for
c         the side cell subdivisions. (output)
c
c    AI   real array containing the best linear estimation coefficients for
c         the interior cell subdivisions. (output)
c
c     M   the number of 2 x 2 subdivisions to perform. It is an error for
c         M to be greater than MXM. (output)
c
c k1,k2   integers giving the size of the initial field (see C0). It is an
c         error for the product k1*k2 to exceed MXK. (output)
c
c    kk   integer giving the size of the initial field covariance matrix
c         (kk = k1*k2). (output)
c
c  iout   unit number to which error and warning messages are logged. (input)
c
c   tol   maximum relative error allowed in the Cholesky decomposition of
c         covariance matrices before a warning message is issued. (input)
c---------------------------------------------------------------------------
c  PARAMETERS:
c   MXK   represents the maximum number of cells (k1 x k2) in the initial
c         field, k1*k2 <= MXK. If the value of MXK is changed here, it must
c         also be changed in LAS2G.
c
c  Requires:
c    1) from libGAFsim:	ISEED, DCVIT2, DCVMT2, DCHOL2, CORN2D, SIDE2D,
c			INTR2D, DCVAA2, DCVAB2, DSIFA, DSISL, DAXPY, DSWAP,
c			IDAMAX, DDOT
c    2) user defined external variance function (see DVARFN)
c
c  REVISION HISTORY:
c  3.41	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c----------------------------------------------------------------------------
      subroutine las2i( dvarfn, N1, N2, XL, YL, kseed, MXM,
     >                  C0, CT, CC, CS, CI, AT, AC, AS, AI, M,
     >                  k1, k2, kk, iout, tol )
      parameter ( MXK = 256 )
      real C0(*), CT(6,*), CC(6,4,*), CS(6,4,*), CI(6,*)
      real AT(3,3,*), AC(4,3,4,*), AS(6,3,4,*), AI(9,3,*)
      real*8 R0(MXK*MXK)
      real*8 R(9,9,2), B(4,4), S(9,4)
      real*8 T1, T2, dvarfn, dble
      logical lformR
      integer mc(4,4), ms(6,4), mi(9)
      external dvarfn
      data mc/5,6,8,9, 4,5,7,8, 2,3,5,6, 1,2,4,5/
      data ms/4,5,6,7,8,9, 2,3,5,6,8,9, 1,2,4,5,7,8, 1,2,3,4,5,6/
      data mi/1,2,3,4,5,6,7,8,9/

   1  format(a,a,a)
   2  format(a,e13.4)
   3  format(a,i4,a,i4,a,i4,a)
c						decompose N1 and N2
      k1 = N1
      k2 = N2
      do 10 m = 0, MXM
         kk = k1*k2
         if( kk .le. MXK ) go to 20
         j1 = k1/2
         j2 = k2/2
         if( 2*j1 .ne. k1 .or. 2*j2 .ne. k2 ) go to 15
         k1 = j1
         k2 = j2
  10  continue
  15  write(iout,1)'Error: unable to determine an acceptable combination
     > of k1, k2 and m'
      write(iout,1)'       such that k1*2**m = N1, k2*2**m = N2.'
      write(iout,3)'       k1 = ',k1,', k2 = ',k2,', m = ',m
      write(iout,3)'       (k1*k2 must be less than ',MXK,' and m must b
     >e less than ',MXM,')'
      write(iout,1)'       Try changing N1 and N2.'
      stop
c						initialize internal generator
  20  kseed = iseed( kseed )
c						form initial covariance matrix
      T1 = dble(XL)/dble(k1)
      T2 = dble(YL)/dble(k2)
      call dcvit2( dvarfn, R0, kk, R, 9, k1, k2, T1, T2 )

c						and compute its cholesky decomp
      call dchol2( R0, kk, kk, rerr )
      if( rerr .gt. tol ) then
         write(iout,1)'Warning: Cholesky decomposition of stage 0 covari
     >ance matrix'
         write(iout,2)'         has maximum relative error of ',rerr
      endif
c						save in real*4 for LAS2G
      L = 0
      do 30 j = 1, kk
      do 30 i = 1, j
         L = L + 1
         C0(L) = R0(i+(j-1)*kk)
  30  continue
c						setup for subsequent stages
      in = 2
      io = 1
      nm = 1
      if( (k1 .eq. 1 .or. k2 .eq. 1) .and. M .gt. 0 ) then
c							special k1,k2 = 1 case
         T1 = 0.5d0*T1
         T2 = 0.5d0*T2
c							get basic cov matrices
         lformR = (1 .lt. M)
         call dcvmt2(dvarfn,R(1,1,in),9,B,4,S,9,T1,T2,lformR)
         i2 = 5
         if( k1 .eq. 1 ) then
            i1 = 2
            i3 = 8
         else
            i1 = 4
            i3 = 6
         endif
         call thin1d(R(1,1,io),9,B,4,S,9,AT,3,AT(1,1,2),3,CT,CT(1,2),
     >               i1,i2,i3,3,iout,tol)
         in = 1
         io = 2
         nm = 2
      endif

      do 40 k = nm, M
         T1 = 0.5d0*T1
         T2 = 0.5d0*T2
c							get basic cov matrices
         lformR = (k .lt. M)
         call dcvmt2(dvarfn,R(1,1,in),9,B,4,S,9,T1,T2,lformR)

c							corner parameters

         call corn2d(R(1,1,io),9,B,4,S,9,CC(1,1,k),6,3,AC(1,1,1,k),
     >               mc,iout,tol)
c							side parameters

         call side2d(R(1,1,io),9,B,4,S,9,CS(1,1,k),6,3,AS(1,1,1,k),
     >               ms,iout,tol)
c							interior parameters

         call intr2d(R(1,1,io),9,B,4,S,9,CI(1,k),3,AI(1,1,k),mi,
     >               iout,tol)
c							swap old/new indices
         ii = in
         in = io
         io = ii
  40  continue

      return
      end
