c  *********************************************************************
c  *                                                                   *
c  *                        subroutine las2g                           *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 3.42
c  Written by Gordon A. Fenton, TUNS, July 19, 1992
c  Latest Update: Jun 9, 1999
c
c  PURPOSE  produces a 2-D quadrant symmetric stationary Gaussian random field
c           using the Local Average Subdivision algorithm.
c
c  This routine creates a zero-mean realization of a 2-D random process given
c  its variance function (as defined by E.H. Vanmarcke in "Random Fields:
c  Analysis and Synthesis", MIT Press, 1984). Each discrete value generated
c  herein represents the local average of a realization of the process over
c  the area Dx x Dy, where (Dx,Dy) is the grid spacing of the desired field.
c  The construction of the realization proceeds recursively as follows;
c
c    1) generate a low resolution field of size k1 x k2. If (N1 x N2) is
c       the desired field resolution, then k1 and k2 are determined so
c       that N1 = k1*2**m, N2 = k2*2**m, where 2**m is a common factor and
c       k1*k2 <= MXK. Generally k1 and k2 are maximized where possible.
c       This is refered to subsequently as the Stage 0 generation.
c
c    2) subdivide the domain m times by dividing each cell into 4 equal
c       parts (2 x 2). In each subdivision, new random values are generated
c       for each new cell. The parent cells of the previous stage are used
c       to obtain best linear estimates of the mean of each new cell so that
c       the spatial correlation is approximated. Also upwards averaging is
c       preserved (ie the average of the 4 new values is the same as the
c       parent cell value). Only parent cells in a neighbourhood of 3 x 3
c       are considered (thus the approximation to the spatial correlation).
c
c  The linear estimation of the mean is accomplished by using the covariance
c  between local averages over each cell, consistent with the goal of
c  producing a local average field. Note that this conditioning process
c  implies that the construction of cells near the edge of the boundary will
c  require the use of values which are, strictly speaking, outside the
c  boundary of the field in question. This is handled by using special
c  reduced neighbourhoods along the boundaries (equivalent to saying that
c  what goes on beyond the boundary has no effect on the process within the
c  boundary).
c
c  Note that this routine sets up a number of parameters on the
c  first call and thus the time required to produce the first realization
c  is substantially greater than on subsequent calls (see INIT). For
c  more information on local average processes, see
c
c   1) G.A. Fenton, "Simulation and Analysis of Random Fields", Ph.D. thesis,
c      Dept. of Civil Eng. and Op. Research, Princeton University, Princeton,
c      New Jersey, 1990.
c   2) G.A. Fenton and E.H. Vanmarcke "Simulation of Random Fields
c      via Local Average Subdivision", ASCE Journal of Engineering Mechanics,
c      Vol. 116, No. 8, August 1990.
c   3) E. H. Vanmarcke, Random Fields: Analysis and Synthesis, MIT Press,
c      1984
c   4) G.A. Fenton, Error evaluation of three random field generators,
c      ASCE Journal of Engineering Mechanics (to appear), 1993.
c
c  Arguments to this routine are as follows;
c
c      Z  real array of size at least N1 x (5*N2/4) which on output will
c         contain the realization of the 2-D process in the first N1 x N2
c         locations. When dimensioned as Z(N1,5*N2/4) in the calling routine
c         and indexed as Z(i,j), Z(1,1) is the lower left cell, Z(2,1) is
c         the next cell in the X direction (to the right), etc.,
c         while Z(1,2) is the next cell in the Y direction (upwards), etc.
c         The extra space is required for workspace (specifically to store
c         the previous stage in the subdivision). (output)
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
c dvarfn  external real*8 function which returns the variance of the random
c         process averaged over a given area. dvarfn is referenced as follows
c
c                var = dvarfn( V1, V2 )
c
c         where (V1,V2) are the side dimensions of the rectangular averaging
c         domain. Any other parameters to the function must be passed by
c         common block from the calling routine. Note that the variance of
c         the random process averaged over the area (V1 x V2) is the product
c         of the point variance and the traditionally defined "variance"
c         function, as discussed by Vanmarcke (pg 186).
c
c  KSEED  integer seed to be used to initialize the pseudo-random number
c         generator. The generator is only initialized on the first call
c         to this routine or when abs(INIT) = 1.
c         If KSEED = 0, then a random seed will be used (based on the
c         process ID of the current program invocation - see iseed.f)
c         On output, KSEED is set to the value of the actual seed used.
c
c   INIT  integer flag which must be 1 when parameters of the process
c         are to be calculated or recalculated. If multiple realizations
c         of the same process are desired, then subsequent calls should
c         use INIT equal to 0 and M less than or equal to the value used 
c         initially. Note that if INIT = -1 is used, then only the
c         initialization is performed.
c
c   IOUT  unit number to which error and warning messages are to be logged.
c         (input)
c  -------------------------------------------------------------------------
c  PARAMETERS:
c   MXM   represents the maximum value that m can have in the decomposition
c         N = k * 2**m.
c
c   MXK   represents the maximum number of cells (k1 x k2) in the initial
c         field, k1*k2 <= MXK. If the value of MXK is changed here, it must
c         also be changed in LAS2I.
c
c   NGS   represents the maximum number of random Gaussian variates that can
c         be produced on each call to VNORM. NGS should be 3*[9*2**(MXM-1)],
c         but not less than MXK.
c
c  Notes: 
c    1) Simulation timing is available through common block LASTYM where
c       TI is the time required to initialize the parameters of the
c       process and TS is the cumulative simulation time. These timings
c       are obtained through the function SECOND which returns elapsed
c       user time in seconds since the start of the program.
c    2) In 2 and higher dimensions, the LAS method yeilds a slight pattern
c       in the point variance field. This is due to the neighborhood
c       approximations. This error lessens for larger scales of fluctuation
c       and smaller numbers of subdivisions.
c    3) The parameter `tol' is the maximum relative error allowed on the
c       Cholesky decomposition of covariance matrices before a warning
c       message is emitted. The relative error is estimated by computing
c       the lower-right most element of L*L' and comparing to the original
c       element of A (where L*L' = A is the decomposition). This give some
c       measure of the roundoff errors accumulated in the computation of L.
c
c  Requires:
c    1) from libGAFsim:	ISEED, LAS2I, DCVIT2, DCVMT2, DCHOL2, CORN2D, SIDE2D,
c			INTR2D, DCVAA2, DCVAB2, DSIFA, DSISL, DAXPY, DSWAP,
c			IDAMAX, DDOT, VNORM, RANDF, SECOND
c    2) user defined external variance function (see DVARFN)
c
c  REVISION HISTORY:
c  3.41	export saved parameters via /las2gb/ (previously just static)
c  3.42	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c ==========================================================================
      subroutine las2g( Z, N1, N2, XL, YL, dvarfn, kseed, init, iout )
      parameter( MXM = 9, MXK = 256, NGS = 6912 )
      real Z(*)
      real C0((MXK*(MXK + 1))/2), U(NGS)
      real CT(6,2), CC(6,4,MXM), CS(6,4,MXM), CI(6,MXM)
      real AT(3,3,2), AC(4,3,4,MXM), AS(6,3,4,MXM), AI(9,3,MXM)
      external dvarfn
      common/las2gb/C0,CT,CC,CS,CI,AT,AC,AS,AI,M,K1,K2,KK,NN
      common/LASTYM/ ti, ts
      data ifirst/1/, zero/0.0/, four/4.0/
      data tol/1.e-3/
c-------------------------------------- initialize LAS generator -------------

      if( ifirst .eq. 1 .or. iabs(init) .eq. 1 ) then
c						start timer
         ti = second()
         call las2i( dvarfn, N1, N2, XL, YL, kseed, MXM,
     >               C0, CT, CC, CS, CI, AT, AC, AS, AI, M,
     >               k1, k2, kk, iout, tol )
         NN = N1*N2
         ifirst = 0
c						all done, set timers
         ts   = zero
         ti   = second() - ti
         if( init .eq. -1 ) return
      endif
c-------------------------------------- generate realization -----------------
      tt = second()
      if( mod(M,2) .eq. 0 ) then
         iz = 0
         jz = NN
      else
         iz = NN
         jz = 0
      endif
c					generate stage 0 field
      call vnorm( U, kk )
      L = 1
      do 20 i = 1, kk
         Z(iz+i) = C0(L)*U(1)
         do 10 j = 2, i
            Z(iz+i) = Z(iz+i) + C0(L+j-1)*U(j)
  10     continue
         L = L + i
  20  continue
c					generate stage 1, 2, ..., M sub-fields
      jx = k1
      jy = k2
      nm = 1
      if( (k1.eq.1 .or. k2.eq.1) .and. M .gt. 0 ) then
         jq = max0(k1,k2)
         call vnorm( U, 3*jq )
         it = jz
         jz = iz
         iz = it
         ix = 2*k1
         iy = 2
         if( k1 .eq. 1 ) iy = 4
         i0 = iz + 1
         i1 = i0 + ix
         j0 = jz + 1
         Z(i0)   = AT(1,1,1)*Z(j0) + AT(2,1,1)*Z(j0+1)
     >             + CT(1,1)*U(1)
         Z(i0+1) = AT(1,2,1)*Z(j0) + AT(2,2,1)*Z(j0+1)
     >             + CT(2,1)*U(1) + CT(3,1)*U(2)
         Z(i1)   = AT(1,3,1)*Z(j0) + AT(2,3,1)*Z(j0+1)
     >             + CT(4,1)*U(1) + CT(5,1)*U(2) + CT(6,1)*U(3)
         Z(i1+1) = four*Z(j0) - Z(i0) - Z(i0+1) - Z(i1)
         L = 4
         do 30 js = 1, jq-2
            i0 = i0 + iy
            i1 = i1 + iy
            Z(i0)   = AT(1,1,2)*Z(j0  ) + AT(2,1,2)*Z(j0+1)
     >              + AT(3,1,2)*Z(j0+2)
     >              + CT(1,2)*U(L)
            Z(i0+1) = AT(1,2,2)*Z(j0  ) + AT(2,2,2)*Z(j0+1)
     >              + AT(3,2,2)*Z(j0+2)
     >              + CT(2,2)*U(L) + CT(3,2)*U(L+1)
            Z(i1)   = AT(1,3,2)*Z(j0  ) + AT(2,3,2)*Z(j0+1)
     >              + AT(3,3,2)*Z(j0+2)
     >              + CT(4,2)*U(L) + CT(5,2)*U(L+1) + CT(6,2)*U(L+2)
            Z(i1+1) = four*Z(j0+1) - Z(i0) - Z(i0+1) - Z(i1)
            j0 = j0 + 1
            L  = L  + 3
  30     continue
         i0 = i0 + iy
         i1 = i1 + iy
         Z(i0)   = AT(2,1,1)*Z(j0) + AT(1,1,1)*Z(j0+1)
     >             + CT(1,1)*U(1)
         Z(i0+1) = AT(2,2,1)*Z(j0) + AT(1,2,1)*Z(j0+1)
     >             + CT(2,1)*U(1) + CT(3,1)*U(2)
         Z(i1)   = AT(2,3,1)*Z(j0) + AT(1,3,1)*Z(j0+1)
     >             + CT(4,1)*U(1) + CT(5,1)*U(2) + CT(6,1)*U(3)
         Z(i1+1) = four*Z(j0) - Z(i0) - Z(i0+1) - Z(i1)
         jx = 2*k1
         jy = 2*k2
         nm = 2
      endif

      do 80 i = nm, M
c						swap current and prev fields
         it = jz
         jz = iz
         iz = it
c						new field dimensions
         ix = 2*jx
         iy = 2*jy
c						pointers into Z
         j0 = jz + 1
         j1 = j0 + jx
         i0 = iz + 1
         i1 = i0 + ix
         call vnorm( U, 3*jx )
c								corner #1
         Z(i0)   = AC(1,1,1,i)*Z(j0) + AC(2,1,1,i)*Z(j0+1)
     >           + AC(3,1,1,i)*Z(j1) + AC(4,1,1,i)*Z(j1+1)
     >           + CC(1,  1,i)*U(1)
         Z(i0+1) = AC(1,2,1,i)*Z(j0) + AC(2,2,1,i)*Z(j0+1)
     >           + AC(3,2,1,i)*Z(j1) + AC(4,2,1,i)*Z(j1+1)
     >           + CC(2,  1,i)*U(1) +  CC(3,  1,i)*U(2)
         Z(i1)   = AC(1,3,1,i)*Z(j0) + AC(2,3,1,i)*Z(j0+1)
     >           + AC(3,3,1,i)*Z(j1) + AC(4,3,1,i)*Z(j1+1)
     >           + CC(4,  1,i)*U(1) +  CC(5,  1,i)*U(2)
     >           + CC(6,  1,i)*U(3)
         Z(i1+1) = four*Z(j0) - Z(i0) - Z(i0+1) - Z(i1)
c								side #1
         L = 4
         do 40 js = 1, jx-2
            i0 = i0 + 2
            i1 = i1 + 2
            Z(i0)   = AS(1,1,1,i)*Z(j0  ) + AS(2,1,1,i)*Z(j0+1)
     >              + AS(3,1,1,i)*Z(j0+2) + AS(4,1,1,i)*Z(j1  )
     >              + AS(5,1,1,i)*Z(j1+1) + AS(6,1,1,i)*Z(j1+2)
     >              + CS(1,  1,i)*U(L)
            Z(i0+1) = AS(1,2,1,i)*Z(j0  ) + AS(2,2,1,i)*Z(j0+1)
     >              + AS(3,2,1,i)*Z(j0+2) + AS(4,2,1,i)*Z(j1  )
     >              + AS(5,2,1,i)*Z(j1+1) + AS(6,2,1,i)*Z(j1+2)
     >              + CS(2,  1,i)*U(L)    + CS(3,  1,i)*U(L+1)
            Z(i1)   = AS(1,3,1,i)*Z(j0  ) + AS(2,3,1,i)*Z(j0+1)
     >              + AS(3,3,1,i)*Z(j0+2) + AS(4,3,1,i)*Z(j1  )
     >              + AS(5,3,1,i)*Z(j1+1) + AS(6,3,1,i)*Z(j1+2)
     >              + CS(4,  1,i)*U(L)    + CS(5,  1,i)*U(L+1)
     >              + CS(6,  1,i)*U(L+2)
            Z(i1+1) = four*Z(j0+1) - Z(i0) - Z(i0+1) - Z(i1)
            j0 = j0 + 1
            j1 = j1 + 1
            L  = L  + 3
  40     continue
c								corner #2
         i0 = i0 + 2
         i1 = i1 + 2
         Z(i0)   = AC(1,1,2,i)*Z(j0) + AC(2,1,2,i)*Z(j0+1)
     >           + AC(3,1,2,i)*Z(j1) + AC(4,1,2,i)*Z(j1+1)
     >           + CC(1,  2,i)*U(L)
         Z(i0+1) = AC(1,2,2,i)*Z(j0) + AC(2,2,2,i)*Z(j0+1)
     >           + AC(3,2,2,i)*Z(j1) + AC(4,2,2,i)*Z(j1+1)
     >           + CC(2,  2,i)*U(L)  + CC(3,  2,i)*U(L+1)
         Z(i1)   = AC(1,3,2,i)*Z(j0) + AC(2,3,2,i)*Z(j0+1)
     >           + AC(3,3,2,i)*Z(j1) + AC(4,3,2,i)*Z(j1+1)
     >           + CC(4,  2,i)*U(L)  + CC(5,  2,i)*U(L+1)
     >           + CC(6,  2,i)*U(L+2)
         Z(i1+1) = four*Z(j0+1) - Z(i0) - Z(i0+1) - Z(i1)

         j0 = jz + 1
         do 60 ks = 1, jy-2
            j1 = j0 + jx
            j2 = j1 + jx
            i0 = i1 + 2
            i1 = i0 + ix
            call vnorm( U, 3*jx )
c								side #2
            Z(i0)   = AS(1,1,2,i)*Z(j0) + AS(2,1,2,i)*Z(j0+1)
     >              + AS(3,1,2,i)*Z(j1) + AS(4,1,2,i)*Z(j1+1)
     >              + AS(5,1,2,i)*Z(j2) + AS(6,1,2,i)*Z(j2+1)
     >              + CS(1,  2,i)*U(1)
            Z(i0+1) = AS(1,2,2,i)*Z(j0) + AS(2,2,2,i)*Z(j0+1)
     >              + AS(3,2,2,i)*Z(j1) + AS(4,2,2,i)*Z(j1+1)
     >              + AS(5,2,2,i)*Z(j2) + AS(6,2,2,i)*Z(j2+1)
     >              + CS(2,  2,i)*U(1)  + CS(3,  2,i)*U(2)
            Z(i1)   = AS(1,3,2,i)*Z(j0) + AS(2,3,2,i)*Z(j0+1)
     >              + AS(3,3,2,i)*Z(j1) + AS(4,3,2,i)*Z(j1+1)
     >              + AS(5,3,2,i)*Z(j2) + AS(6,3,2,i)*Z(j2+1)
     >              + CS(4,  2,i)*U(1)  + CS(5,  2,i)*U(2)
     >              + CS(6,  2,i)*U(3)
            Z(i1+1) = four*Z(j1) - Z(i0) - Z(i0+1) - Z(i1)
c								interior
            L = 4
            do 50 js = 1, jx-2
               i0 = i0 + 2
               it = i0 + 1
               i1 = i1 + 2
               Z(i0)=AI(1,1,i)*Z(j0)+AI(2,1,i)*Z(j0+1)+AI(3,1,i)*Z(j0+2)
     >              +AI(4,1,i)*Z(j1)+AI(5,1,i)*Z(j1+1)+AI(6,1,i)*Z(j1+2)
     >              +AI(7,1,i)*Z(j2)+AI(8,1,i)*Z(j2+1)+AI(9,1,i)*Z(j2+2)
     >              +CI(1,  i)*U(L)
               Z(it)=AI(1,2,i)*Z(j0)+AI(2,2,i)*Z(j0+1)+AI(3,2,i)*Z(j0+2)
     >              +AI(4,2,i)*Z(j1)+AI(5,2,i)*Z(j1+1)+AI(6,2,i)*Z(j1+2)
     >              +AI(7,2,i)*Z(j2)+AI(8,2,i)*Z(j2+1)+AI(9,2,i)*Z(j2+2)
     >              +CI(2,  i)*U(L) +CI(3,  i)*U(L+1)
               Z(i1)=AI(1,3,i)*Z(j0)+AI(2,3,i)*Z(j0+1)+AI(3,3,i)*Z(j0+2)
     >              +AI(4,3,i)*Z(j1)+AI(5,3,i)*Z(j1+1)+AI(6,3,i)*Z(j1+2)
     >              +AI(7,3,i)*Z(j2)+AI(8,3,i)*Z(j2+1)+AI(9,3,i)*Z(j2+2)
     >              +CI(4,  i)*U(L) +CI(5,  i)*U(L+1) +CI(6,  i)*U(L+2)
               Z(i1+1) = four*Z(j1+1) - Z(i0) - Z(it) - Z(i1)
               j0 = j0 + 1
               j1 = j1 + 1
               j2 = j2 + 1
               L  = L  + 3
  50        continue
c								side #3
            i0 = i0 + 2
            i1 = i1 + 2
            Z(i0)   = AS(1,1,3,i)*Z(j0) + AS(2,1,3,i)*Z(j0+1)
     >              + AS(3,1,3,i)*Z(j1) + AS(4,1,3,i)*Z(j1+1)
     >              + AS(5,1,3,i)*Z(j2) + AS(6,1,3,i)*Z(j2+1)
     >              + CS(1,  3,i)*U(L)
            Z(i0+1) = AS(1,2,3,i)*Z(j0) + AS(2,2,3,i)*Z(j0+1)
     >              + AS(3,2,3,i)*Z(j1) + AS(4,2,3,i)*Z(j1+1)
     >              + AS(5,2,3,i)*Z(j2) + AS(6,2,3,i)*Z(j2+1)
     >              + CS(2,  3,i)*U(L)  + CS(3,  3,i)*U(L+1)
            Z(i1)   = AS(1,3,3,i)*Z(j0) + AS(2,3,3,i)*Z(j0+1)
     >              + AS(3,3,3,i)*Z(j1) + AS(4,3,3,i)*Z(j1+1)
     >              + AS(5,3,3,i)*Z(j2) + AS(6,3,3,i)*Z(j2+1)
     >              + CS(4,  3,i)*U(L)  + CS(5,  3,i)*U(L+1)
     >              + CS(6,  3,i)*U(L+2)
            Z(i1+1) = four*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)
            j0 = j0 + 2
  60     continue

         j1 = j0 + jx
         i0 = i1 + 2
         i1 = i0 + ix
         call vnorm( U, 3*jx )
c								corner #3
         Z(i0)   = AC(1,1,3,i)*Z(j0) + AC(2,1,3,i)*Z(j0+1)
     >           + AC(3,1,3,i)*Z(j1) + AC(4,1,3,i)*Z(j1+1)
     >           + CC(1,  3,i)*U(1)
         Z(i0+1) = AC(1,2,3,i)*Z(j0) + AC(2,2,3,i)*Z(j0+1)
     >           + AC(3,2,3,i)*Z(j1) + AC(4,2,3,i)*Z(j1+1)
     >           + CC(2,  3,i)*U(1)  + CC(3,  3,i)*U(2)
         Z(i1)   = AC(1,3,3,i)*Z(j0) + AC(2,3,3,i)*Z(j0+1)
     >           + AC(3,3,3,i)*Z(j1) + AC(4,3,3,i)*Z(j1+1)
     >           + CC(4,  3,i)*U(1)  + CC(5,  3,i)*U(2)
     >           + CC(6,  3,i)*U(3)
         Z(i1+1) = four*Z(j1) - Z(i0) - Z(i0+1) - Z(i1)
c								side #4
         L = 4
         do 70 js = 1, jx-2
            i0 = i0 + 2
            i1 = i1 + 2
            Z(i0)   = AS(1,1,4,i)*Z(j0  ) + AS(2,1,4,i)*Z(j0+1)
     >              + AS(3,1,4,i)*Z(j0+2) + AS(4,1,4,i)*Z(j1  )
     >              + AS(5,1,4,i)*Z(j1+1) + AS(6,1,4,i)*Z(j1+2)
     >              + CS(1,  4,i)*U(L)
            Z(i0+1) = AS(1,2,4,i)*Z(j0  ) + AS(2,2,4,i)*Z(j0+1)
     >              + AS(3,2,4,i)*Z(j0+2) + AS(4,2,4,i)*Z(j1  )
     >              + AS(5,2,4,i)*Z(j1+1) + AS(6,2,4,i)*Z(j1+2)
     >              + CS(2,  4,i)*U(L)    + CS(3,  4,i)*U(L+1)
            Z(i1)   = AS(1,3,4,i)*Z(j0  ) + AS(2,3,4,i)*Z(j0+1)
     >              + AS(3,3,4,i)*Z(j0+2) + AS(4,3,4,i)*Z(j1  )
     >              + AS(5,3,4,i)*Z(j1+1) + AS(6,3,4,i)*Z(j1+2)
     >              + CS(4,  4,i)*U(L)    + CS(5,  4,i)*U(L+1)
     >              + CS(6,  4,i)*U(L+2)
            Z(i1+1) = four*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)
            j0 = j0 + 1
            j1 = j1 + 1
            L  = L  + 3
  70     continue
c								corner #4
         i0 = i0 + 2
         i1 = i1 + 2
         Z(i0)   = AC(1,1,4,i)*Z(j0) + AC(2,1,4,i)*Z(j0+1)
     >           + AC(3,1,4,i)*Z(j1) + AC(4,1,4,i)*Z(j1+1)
     >           + CC(1,  4,i)*U(L)
         Z(i0+1) = AC(1,2,4,i)*Z(j0) + AC(2,2,4,i)*Z(j0+1)
     >           + AC(3,2,4,i)*Z(j1) + AC(4,2,4,i)*Z(j1+1)
     >           + CC(2,  4,i)*U(L)  + CC(3,  4,i)*U(L+1)
         Z(i1)   = AC(1,3,4,i)*Z(j0) + AC(2,3,4,i)*Z(j0+1)
     >           + AC(3,3,4,i)*Z(j1) + AC(4,3,4,i)*Z(j1+1)
     >           + CC(4,  4,i)*U(L)  + CC(5,  4,i)*U(L+1)
     >           + CC(6,  4,i)*U(L+2)
         Z(i1+1) = four*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)

         jx = ix
         jy = iy
  80  continue
c						all done, compute elapsed time
      ts = ts + (second() - tt)

      return
      end
