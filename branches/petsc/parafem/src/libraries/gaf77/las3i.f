c  *********************************************************************
c  *                                                                   *
c  *                         subroutine las3i                          *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 3.2
c  Written by Gordon A. Fenton, TUNS, Aug. 26, 1992
c  Latest Update: Feb. 22, 1994
c
c  PURPOSE   initializes parameters for LAS3G
c
c  This routine sets up the matrices required by LAS3G to construct
c  realizations of the random field. The covariances between local averages
c  at each subdivision stage are computed in double precision for accuracy,
c  but the final construction matrices are stored in single precision and
c  returned to the calling routine via the argument list. The general
c  recursive field construction follows the relationship
c
c         {Z^j} = [A^T]{Z^(j-1)} + [C]{U}
c
c  where {Z^j} is a vector of length 8 representing the values assigned to
c  the 2 x 2 x 2 cell subdivision and {Z^(j-1)} are the parent cell values in
c  some neighbourhood. The following figure illustrates the cell subdivision,
c  the neighbourhood, and the numbering scheme for an interior cell (special
c  subsets of the neighbourhood are used for the corners, edges, and sides).
c  Only one plane is shown for simplicity, however the top and side views are
c  similar, although the numbering will change (there being 27 parent elements
c  and 8 subdivided elements).
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
c  The first stage of the simulation involves the direct generation of a
c  k1 x k2 x k3 cell array, where k1, k2, and k3 are integers which satisfy
c  the decomposition N1 = k1*2**m, N2 = k2*2**m, and N3 = k3*2**m for a
c  common factor 2**m. The integers k1, k2, and k3 are chosen to be as large
c  as possible while requiring the product k1*k2*k3 to be less than or equal
c  to MXK. Note that the direct simulation involves the Cholesky Decomp. of a
c  MXK x MXK matrix (at the upper limit) and so MXK should not be overly large.
c  This formulation is somewhat less restrictive than simply requiring
c  N1, N2, and N3 to be powers of 2. Also N1, N2, and N3 do not have to be
c  equal. However N1, N2, and N3 cannot still be chosen arbitrarily, for
c  example the set (N1,N2,N3) = (144,256,256) results in k1 = 9, k2 = 16,
c  k3 = 16, m = 4 which is not acceptable here (since k1*k2*k3 > MXK for
c  MXK = 1000), while the set (N1,N2) = (160,256,256) is acceptable since
c  k1 = 5, k2 = 8, k3 = 8, m = 5. In general it may be easier to choose k1,
c  k2, k3, and m before specifying N1, N2, and N3. In the event that an
c  unacceptable (k1,k2,k3,m) combination is selected, IERR is set to -1 and
c  control returned to the calling routine.
c  The maximum value of m is set by the calling routine in the argument MXM.
c
c  Arguments to this routine are as follows;
c
c  dvfn   external real*8 function which returns the variance of the random
c         process averaged over a given volume. dvfn is referenced as follows
c
c                var = dvfn( V1, V2, V3 )
c
c         where (V1,V2,V3) are the side dimensions of the rectangular averaging
c         domain. Any other parameters to the function may be passed by
c         common block from the calling routine. Note that the variance of
c         the random process averaged over the volume (V1 x V2 x V3) is the
c         product of the point variance and the traditionally defined
c         "variance" function, as discussed by Vanmarcke (pg 186).
c
c  N1     number of cells to discretize the field in the x, y, and z directions
c  N2     respectively (corresponding to the first, second, and third indices
c  N3     of Z respectively). N1, N2, and N3 must have the form
c            N_i = k_i * 2**m
c         where m is common to N1, N2, and N3 and k1, k2, and k3 are
c         positive integers satisfying k1*k2*k3 <= MXK. Generally k1, k2 and k3
c         are chosen to be as large as possible and still satisfy the above
c         requirements so the the first stage involves directly simulating
c         a k1 x k2 x k3 field by inversion of a covariance matrix. (input)
c
c  XL     physical dimensions of the field in the x-direction. (input)
c  YL     physical dimensions of the field in the y-direction. (input)
c  ZL     physical dimensions of the field in the z-direction. (input)
c
c  KSEED  integer seed to be used for the pseudo-random number generator.
c         If KSEED = 0, then a random seed will be used (based on the
c         clock time when this routine is called for the first time).
c         On output, KSEED is set to the value of the actual seed used.
c
c  MXM    integer giving the largest value that M can take. An error is
c         generated if the process size is such that M > MXM. (input)
c
c  C0     real vector containing the lower triangular values of the Cholesky
c         decomposition of the covariance matrix for the initial stage of
c         k1 x k2 x k3 cells. C0 must be at least k1*k2*k3 x k1*k2*k3 in
c         length. (output)
c
c  CC     real vector containing the lower triangular values of the Cholesky
c         decomposition of the covariance matrix for the corner cell 2 x 2 x 2
c         subdivisions. (output)
c
c  CE     real vector containing the lower triangular values of the Cholesky
c         decomposition of the covariance matrix for the edge cell 2 x 2 x 2
c         subdivisions. (output)
c
c  CS     real vector containing the lower triangular values of the Cholesky
c         decomposition of the covariance matrix for the side cell 2 x 2 x 2
c         subdivisions. (output)
c
c  CI     real vector containing the lower triangular values of the Cholesky
c         decomposition of the covariance matrix for the interior cell
c         2 x 2 x 2 subdivisions. (output)
c
c  AC     real array containing the best linear estimation coefficients for
c         the corner cell subdivisions. (output)
c
c  AE     real array containing the best linear estimation coefficients for
c         the edge cell subdivisions. (output)
c
c  AS     real array containing the best linear estimation coefficients for
c         the side cell subdivisions. (output)
c
c  AI     real array containing the best linear estimation coefficients for
c         the interior cell subdivisions. (output)
c
c  M      the number of 2 x 2 x 2 subdivisions to perform. It is an error for
c         M to be greater than MXM. (output)
c
c  k1     integers giving the size of the initial field (see C0). It is an
c  k2     error for the product k1*k2*k3 to exceed MXK.
c  k3     (output)
c
c  kk     integer giving the size of the initial field covariance matrix
c         (kk = k1*k2*k3). (output)
c
c  iout   unit number to which error and warning messages are logged. (input)
c
c   tol   maximum relative error allowed in the Cholesky decomposition of
c         covariance matrices before a warning message is issued. (input)
c---------------------------------------------------------------------------
c  PARAMETERS:
c   MXK   represents the maximum number of cells (k1 x k2 x k3) in the initial
c         field, k1*k2*k3 <= MXK. If the value of MXK is changed here, it must
c         also be changed in LAS2G.
c
c  Requires:
c    1) from libGAFsim:	ISEED, DCVIT3, DCVMT3, DCHOL2, CORN3D, EDGE3D,
c			SIDE3D, INTR3D, DCVAA3, DCVAB3, DSIFA, DSISL, DAXPY,
c			DSWAP, IDAMAX, DDOT
c    4) user defined external variance function (see DVFN)
c----------------------------------------------------------------------------
      subroutine las3i( dvfn, N1, N2, N3, XL, YL, ZL, kseed, MXM,
     >                  C0, CC, CE, CS, CI, AC, AE, AS, AI,
     >                  ATC, ATS, ATI, CTC, CTS, CTI,
     >                  M, k1, k2, k3, kk, iout, tol )
      parameter ( MXK = 512 )
      real C0(*), CC(28,8,*), CE(28,12,*), CS(28,6,*), CI(28,*)
      real AC(8,7,8,*), AE(12,7,12,*), AS(18,7,6,*), AI(27,7,*)
      real ATC(4,7,*), ATS(6,7,*), ATI(9,*)
      real CTC(28,*), CTS(28,*), CTI(*)
      real*8 R0(MXK*MXK)
      real*8 R(27,27,2), B(8,8), S(27,8)
      real*8 T1, T2, T3, dvfn, dble
      logical lformR, lk1, lk2, lk3
      integer mc(4,4,3), ms(6,4,3), mi(9,3)
      external dvfn

      data mc/14,17,23,26,   11,14,20,23,    5, 8,14,17,   2, 5,11,14,
     >        14,15,23,24,   13,14,22,23,    5, 6,14,15,   4, 5,13,14,
     >        14,15,17,18,   13,14,16,17,   11,12,14,15,  10,11,13,14/
      data ms/11,14,17,20,23,26,   5, 8,14,17,23,26,
     >         2, 5,11,14,20,23,   2, 5, 8,11,14,17,
     >        13,14,15,22,23,24,   5, 6,14,15,23,24,
     >         4, 5,13,14,22,23,   4, 5, 6,13,14,15,
     >        13,14,15,16,17,18,  11,12,14,15,17,18,
     >        10,11,13,14,16,17,  10,11,12,13,14,15/
      data mi/ 2, 5, 8,11,14,17,20,23,26,
     >         4, 5, 6,13,14,15,22,23,24,
     >        10,11,12,13,14,15,16,17,18/

   1  format(a,a,a)
   2  format(a,e13.4)
   3  format(a,i4,a,i4,a,i4,a,i4,a)
c						decompose N1, N2, and N3
      k1 = N1
      k2 = N2
      k3 = N3
      do 10 m = 0, MXM
         kk = k1*k2*k3
         if( kk .le. MXK ) go to 20
         j1 = k1/2
         j2 = k2/2
         j3 = k3/2
         if( 2*j1 .ne. k1 .or. 2*j2 .ne. k2 .or. 2*j3 .ne. k3 ) go to 20
         k1 = j1
         k2 = j2
         k3 = j3
  10  continue
  15  write(iout,1)'Error: unable to determine an acceptable combination
     > of k1, k2, k3 and m'
      write(iout,1)'       such that k1*2**m = N1, k2*2**m = N2 and k3*2
     >**m = N3.'
      write(iout,3)'       k1 = ',k1,', k2 = ',k2,', k3 = ',k3,
     >             ', m = ',m
      write(iout,3)'       (k1*k2*k3 must be less than ',MXK,' and m mus
     >t be less than ',MXM,')'
      write(iout,1)'       Try changing N1, N2, and/or N3.'
      stop

  20  if(  kk .gt. MXK ) go to 15
      ifirst = 0
c						initialize internal generator
      kseed = iseed( kseed )
c						form initial covariance matrix
      T1 = dble(XL)/dble(k1)
      T2 = dble(YL)/dble(k2)
      T3 = dble(ZL)/dble(k3)
      call dcvit3( dvfn, R0, kk, R, 27, k1, k2, k3, T1, T2, T3 )

c						and compute its cholesky decomp
      call dchol2(R0,kk,kk,rerr)
      if( rerr .gt. tol ) then
         write(iout,1)'Warning: Cholesky decomposition of stage 0 covari
     >ance matrix'
         write(iout,2)'         has maximum relative error of ',rerr
      endif
c						save in real*4 for LAS2G
      L = 0
      do 30 j = 1, kk
         jj = (j-1)*kk
      do 30 i = 1, j
         L = L + 1
         C0(L) = R0(i+jj)
  30  continue
c						setup for subsequent stages
      in = 2
      io = 1
      nm = 1
c						(special k1,k2,k3 = 1 case)
      lk1 = (k1 .eq. 1)
      lk2 = (k2 .eq. 1)
      lk3 = (k3 .eq. 1)
      if( (lk1 .or. lk2 .or. lk3) .and. (m .gt. 0) ) then
         T1 = 0.5d0*T1
         T2 = 0.5d0*T2
         T3 = 0.5d0*T3
         lformR = (1 .lt. m)
         call dcvmt3(dvfn,R(1,1,in),27,B,8,S,27,T1,T2,T3,lformR)

         if((lk1.and.lk2) .or. (lk1.and.lk3) .or. (lk2.and.lk3)) then
            i2 = 14
            if( .not. lk1 ) then
               i1 = 13
               i3 = 15
            elseif( .not. lk2 ) then
               i1 = 11
               i3 = 17
            elseif( .not. lk3 ) then
               i1 = 5
               i3 = 23
            else
c						this should not happen...
               write(iout,1)'Error: not all k1, k2, k3 can be 1.'
               stop
            endif

            call thin1d(R(1,1,io),27,B,8,S,27,ATS,6,ATI,9,CTS,CTI,
     >                  i1,i2,i3,7,iout,tol)
         else
            if( lk1 ) k = 1
            if( lk2 ) k = 2
            if( lk3 ) k = 3
            call corn2d(R(1,1,io),27,B,8,S,27,CTC,28,7,ATC,
     >                  mc(1,1,k),iout,tol)

            call side2d(R(1,1,io),27,B,8,S,27,CTS,28,7,ATS,
     >                  ms(1,1,k),iout,tol)

            call intr2d(R(1,1,io),27,B,8,S,27,CTI,7,ATI,
     >                  mi(1,k),iout,tol)
         endif
         in = 1
         io = 2
         nm = 2
      endif

      do 40 k = nm, M
         T1 = 0.5d0*T1
         T2 = 0.5d0*T2
         T3 = 0.5d0*T3
c							get basic cov matrices
         lformR = (k .lt. M)

         call dcvmt3(dvfn,R(1,1,in),27,B,8,S,27,T1,T2,T3,lformR)

c							corner parameters

         call corn3d(R(1,1,io),27,B,8,S,27,CC(1,1,k),AC(1,1,1,k),
     >               iout,tol)
c							edge parameters

         call edge3d(R(1,1,io),27,B,8,S,27,CE(1,1,k),AE(1,1,1,k),
     >               iout,tol)
c							side parameters

         call side3d(R(1,1,io),27,B,8,S,27,CS(1,1,k),AS(1,1,1,k),
     >               iout,tol)
c							interior parameters

         call intr3d(R(1,1,io),27,B,8,S,27,CI(1,k),AI(1,1,k),
     >               iout,tol)
c							swap old/new indices
         ii = in
         in = io
         io = ii
  40  continue

      return
      end

