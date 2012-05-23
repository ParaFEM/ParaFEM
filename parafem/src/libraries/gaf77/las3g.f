c  *******************************************************************
c  *                                                                 *
c  *                       Subroutine las3g                          *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 3.42
c  Written by Gordon A. Fenton, Jan. 1990.
c  Latest Update: Jun 21, 2000
c
c  PURPOSE  produces a 3-D quadrant symmetric homogeneous Gaussian random field
c           using the Local Average Subdivision algorithm.
c
c  This routine creates a zero-mean realization of a 3-D random process given
c  its variance function (as defined by E.H. Vanmarcke in "Random Fields:
c  Analysis and Synthesis", MIT Press, 1984). Each discrete value generated
c  herein represents the local average of a realization of the process over
c  the volume Dx x Dy x Dz, where (Dx,Dy,Dz) is the grid spacing of the
c  desired field.
c  The construction of the realization proceeds recursively as follows;
c
c    1) generate a low resolution field of size k1 x k2 x k3. If (N1 x N2 x N3)
c       is the desired field resolution, then k1, k2, and k3 are determined so
c       that N1 = k1*2**m, N2 = k2*2**m, and N3 = k3*2**m, where 2**m is a
c       common factor and k1*k2*k3 <= MXK. Generally k1, k2, and k3 are
c       maximized where possible since this improves the covariance structure
c       of anisotropic fields. However, the direct generation of the k1 x k2 x
c       k3 field involves the Cholesky decomposition of a (k1*k2*k3) x
c       (k1*k2*k3) matrix, and so MXK should not be too large. This is
c       refered to subsequently as the Stage 0 generation.
c
c    2) subdivide the domain m times by dividing each cell into 8 equal
c       parts (2 x 2 x 2). In each subdivision, new random values are generated
c       for each new cell. The parent cells of the previous stage are used
c       to obtain best linear estimates of the mean of each new cell so that
c       the spatial correlation is approximated. Also upwards averaging is
c       preserved (ie the average of the 8 new values is the same as the
c       parent cell value). Only parent cells in a neighbourhood of 3 x 3 x 3
c       are considered (thus the approximation to the spatial correlation).
c
c  The linear estimation of the mean is accomplished by using the covariance
c  between local averages over each cell, consistent with the goal of
c  producing a local average field. Note that this conditioning process
c  implies that the construction of cells near the edge of the boundary will
c  require the use of values which are, strictly speaking, outside the
c  boundary of the field in question. This is handled by using specially
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
c   Z     real array of size at least N1*N2*N3*(9/8) which on output will
c         contain the realization of the 3-D process in the first N1 x N2 x N3
c         locations. When dimensioned as Z(N1,N2,9*N3/8) in the calling routine
c         and indexed as Z(i,j,k), Z(1,1,1) is the lower left front cell,
c         Z(2,1,1) is the next cell in the X direction (to the right), etc.,
c         while Z(1,2,1) is the next cell in the Y direction (upwards), etc.
c         The extra space is required for workspace (specifically to store
c         the previous stage in the subdivision). (output)
c
c   N1    number of cells to discretize the field in the x, y, and z directions
c   N2    respectively (corresponding to the first, second, and third indices
c   N3    of Z respectively). Each N1, N2 and N3 must have the form
c         N_i = k_i * 2**m where m is common to all and k1, k2, and k3 are
c         positive integers satisfying k1*k2*k3 <= MXK. Generally k1 and k2
c         are chosen to be as large as possible and still satisfy the above
c         requirements so the the first stage involves directly simulating
c         a k1 x k2 x k3 field by Cholesky decomposition of a covariance
c         matrix. A possible example would be (N1,N2,N3) = (24,80,72) which
c         gives k1 = 3, k2 = 10, k3 = 9 (note k1*k2*k3 = 270 < 1000, OK) and
c         m = 3. Note that in general N1, N2, and N3 cannot be chosen
c         arbitrarily - it is usually best to choose m first then k1, k2, and
c         k3 so as to satisfy or exceed the problem requirements. (input)
c
c   XL    physical dimensions of the field in the x, y, and z directions.
c   YL    (input)
c   ZL
c
c   dvfn  external real*8 function which returns the variance of the random
c         process averaged over a given area. dvfn is referenced as follows
c
c                var = dvfn( V1, V2, V3 )
c
c         where (V1,V2,V3) are the side dimensions of the rectangular averaging
c         domain. Any other parameters to the function must be passed by
c         common block from the calling routine. Note that the variance of
c         the random process averaged over the volume (V1 x V2 x V3) is the
c         product of the point variance and the traditionally defined
c         "variance" function, as discussed by Vanmarcke (pg 186). This
c         function must return the final variance.
c
c   KSEED integer seed to be used for the pseudo-random number generator.
c         If KSEED = 0, then a random seed will be used (based on the
c         clock time when this routine is called for the first time).
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
c         N = k * 2**m. This is the maximum number of subdivisions. (Note
c         that MXM = 6 allows a field of size 512 x 512 x 512.)
c
c   MXK   represents the maximum number of cells (k1 x k2 x k3) in the initial
c         field, k1*k2*k3 <= MXK. If the value of MXK is changed here, it must
c         also be changed in LAS3I.
c
c   NGS   represents the maximum number of random Gaussian variates that can
c         be produced on each call to VNORM. NGS should be at least
c         7*N1/2, but not less than MXK. Since N1 is not known ahead of
c         time, I've taken the worst case NGS = 1000*2**5 = 32000
c         (corresponding to k1=MXK, k2=k3=1, M = 6)
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
c    1) from libGAFsim:	ISEED, LAS3I, DCVIT3, DCVMT3, DCHOL2, CORN3D, EDGE3D,
c			SIDE3D, INTR3D, DCVAA3, DCVAB3, DSIFA, DSISL, DAXPY,
c			DSWAP, IDAMAX, DDOT, DOT3D, VNORM, RANDF, SECOND
c    2) user defined external variance function (see DVFN)
c
c  REVISION HISTORY:
c  3.4	accommodated M = 0 case where k1, k2, and/or k3 = 1
c  3.41	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c  3.42	saved k1, k2, k3, kk, NN, and M from initialization (Jun 21/00)
c ==========================================================================
      subroutine las3g( Z, N1, N2, N3, XL, YL, ZL, dvfn, kseed, init,
     >                  iout )
      parameter( MXM = 6, MXK = 512, NGS = 32000 )
      real Z(*)
      real C0(MXK*(MXK + 1)/2), U(NGS)
      real CC(28,8,MXM), CE(28,12,MXM), CS(28,6,MXM), CI(28,MXM)
      real AC(8,7,8,MXM), AE(12,7,12,MXM), AS(18,7,6,MXM), AI(27,7,MXM)
      real ATC(4,7,4), ATS(6,7,4), ATI(9,7)
      real CTC(28,4), CTS(28,4), CTI(28)
      save C0, CC, CS, CE, CI, AC, AE, AS, AI
      save ATC, ATS, ATI, CTC, CTS, CTI
      save M, k1, k2, k3, kk, NN
      external dvfn, dot1, dot2, dot3, dot4, dot5, dot6, dot7
      common/LASTYM/ ti, ts
      data ifirst/1/, zero/0.0/, eight/8.0/
      data tol/1.e-3/
c-----------------------------------------------------------------------------
c-------------------------------------- initialize LAS generator -------------

      if( ifirst .eq. 1 .or. iabs(init) .eq. 1 ) then
c						start timer
         ti = second()
         call las3i( dvfn, N1, N2, N3, XL, YL, ZL, kseed, MXM,
     >               C0, CC, CE, CS, CI, AC, AE, AS, AI,
     >               ATC, ATS, ATI, CTC, CTS, CTI,
     >               M, k1, k2, k3, kk, iout, tol )
         NN = N1*N2*N3
         ifirst = 0
c						all done, set timers
         ts   = zero
         ti   = second() - ti
         if( init .eq. -1 ) return
      endif
c-----------------------------------------------------------------------------
c-------------------------------------- generate realization -----------------
      tt = second()
      if( mod(M,2) .eq. 0 ) then
         iq = 0
         jq = NN
      else
         iq = NN
         jq = 0
      endif
c					generate stage 0 field
      call vnorm( U, kk )
      L = 1
      do 20 i = 1, kk
         Z(iq+i) = C0(L)*U(1)
         do 10 j = 2, i
            Z(iq+i) = Z(iq+i) + C0(L+j-1)*U(j)
  10     continue
         L = L + i
  20  continue
c					generate stage 1, 2, ..., M sub-fields
      if( M .eq. 0 ) return
      if( (k1 .eq. 1) .or. (k2 .eq. 1) .or. (k3 .eq. 1) ) then
         ms = 2
         it = jq
         jq = iq
         iq = it
         call plan3d(Z,iq,jq,k1,k2,k3,ATC,ATS,ATI,CTC,CTS,CTI,U,iout)
         jx = 2*k1
         jy = 2*k2
         jz = 2*k3
      else
         jx = k1
         jy = k2
         jz = k3
         ms = 1
      endif

      do 160 i = ms, M
c						swap current and prev fields
         it = jq
         jq = iq
         iq = it
c						new field dimensions
         ix  = 2*jx
         iy  = 2*jy
         iz  = 2*jz
         ixy = ix*iy
         jxy = jx*jy
c						pointers into parent field
         j0 = jq + 1
         j1 = j0 + jx
         j3 = j0 + jxy
         j4 = j3 + jx
c						pointers into new field
         i0 = iq + 1
         i1 = i0 + ix
         i2 = i0 + ixy
         i3 = i2 + ix
         call vnorm( U, 7*jx )
c								corner #1
         Z(i0)   = dot8c(Z,AC(1,1,1,i),CC( 1,1,i),U(1),dot1,j0,j1,j3,j4)
         Z(i0+1) = dot8c(Z,AC(1,2,1,i),CC( 2,1,i),U(1),dot2,j0,j1,j3,j4)
         Z(i1)   = dot8c(Z,AC(1,3,1,i),CC( 4,1,i),U(1),dot3,j0,j1,j3,j4)
         Z(i1+1) = dot8c(Z,AC(1,4,1,i),CC( 7,1,i),U(1),dot4,j0,j1,j3,j4)
         Z(i2)   = dot8c(Z,AC(1,5,1,i),CC(11,1,i),U(1),dot5,j0,j1,j3,j4)
         Z(i2+1) = dot8c(Z,AC(1,6,1,i),CC(16,1,i),U(1),dot6,j0,j1,j3,j4)
         Z(i3)   = dot8c(Z,AC(1,7,1,i),CC(22,1,i),U(1),dot7,j0,j1,j3,j4)
         Z(i3+1) = eight*Z(j0) - Z(i0) - Z(i0+1) - Z(i1)
     >               - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
c								edge #1
         L = 8
         do 30 js = 1, jx-2
          i0 = i0 + 2
          i1 = i1 + 2
          i2 = i2 + 2
          i3 = i3 + 2
          Z(i0)  =dot12h(Z,AE(1,1,1,i),CE( 1,1,i),U(L),dot1,j0,j1,j3,j4)
          Z(i0+1)=dot12h(Z,AE(1,2,1,i),CE( 2,1,i),U(L),dot2,j0,j1,j3,j4)
          Z(i1)  =dot12h(Z,AE(1,3,1,i),CE( 4,1,i),U(L),dot3,j0,j1,j3,j4)
          Z(i1+1)=dot12h(Z,AE(1,4,1,i),CE( 7,1,i),U(L),dot4,j0,j1,j3,j4)
          Z(i2)  =dot12h(Z,AE(1,5,1,i),CE(11,1,i),U(L),dot5,j0,j1,j3,j4)
          Z(i2+1)=dot12h(Z,AE(1,6,1,i),CE(16,1,i),U(L),dot6,j0,j1,j3,j4)
          Z(i3)  =dot12h(Z,AE(1,7,1,i),CE(22,1,i),U(L),dot7,j0,j1,j3,j4)
          Z(i3+1)=eight*Z(j0+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
          j0 = j0 + 1
          j1 = j1 + 1
          j3 = j3 + 1
          j4 = j4 + 1
          L  = L  + 7
  30     continue
c								corner #2
         i0 = i0 + 2
         i1 = i1 + 2
         i2 = i2 + 2
         i3 = i3 + 2
         Z(i0)   = dot8c(Z,AC(1,1,2,i),CC( 1,2,i),U(L),dot1,j0,j1,j3,j4)
         Z(i0+1) = dot8c(Z,AC(1,2,2,i),CC( 2,2,i),U(L),dot2,j0,j1,j3,j4)
         Z(i1)   = dot8c(Z,AC(1,3,2,i),CC( 4,2,i),U(L),dot3,j0,j1,j3,j4)
         Z(i1+1) = dot8c(Z,AC(1,4,2,i),CC( 7,2,i),U(L),dot4,j0,j1,j3,j4)
         Z(i2)   = dot8c(Z,AC(1,5,2,i),CC(11,2,i),U(L),dot5,j0,j1,j3,j4)
         Z(i2+1) = dot8c(Z,AC(1,6,2,i),CC(16,2,i),U(L),dot6,j0,j1,j3,j4)
         Z(i3)   = dot8c(Z,AC(1,7,2,i),CC(22,2,i),U(L),dot7,j0,j1,j3,j4)
         Z(i3+1) = eight*Z(j0+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                 - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

         j0 = jq + 1
         do 50 ks = 1, jy-2
            j1 = j0 + jx
            j2 = j1 + jx
            j3 = j0 + jxy
            j4 = j3 + jx
            j5 = j4 + jx
            i0 = i1 + 2
            i1 = i0 + ix
            i2 = i0 + ixy
            i3 = i2 + ix
            call vnorm( U, 7*jx )
c								edge #2
            Z(i0)   = dot12v(Z,AE(1,1,2,i),CE( 1,2,i),U(1),dot1,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i0+1) = dot12v(Z,AE(1,2,2,i),CE( 2,2,i),U(1),dot2,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i1)   = dot12v(Z,AE(1,3,2,i),CE( 4,2,i),U(1),dot3,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i1+1) = dot12v(Z,AE(1,4,2,i),CE( 7,2,i),U(1),dot4,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i2)   = dot12v(Z,AE(1,5,2,i),CE(11,2,i),U(1),dot5,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i2+1) = dot12v(Z,AE(1,6,2,i),CE(16,2,i),U(1),dot6,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i3)   = dot12v(Z,AE(1,7,2,i),CE(22,2,i),U(1),dot7,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i3+1) = eight*Z(j1  ) - Z(i0) - Z(i0+1) - Z(i1)
     >                    - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
c								side #1
            L = 8
            do 40 js = 1, jx-2
               i0 = i0 + 2
               i1 = i1 + 2
               i2 = i2 + 2
               i3 = i3 + 2
               Z(i0  ) = dot18a(Z,AS(1,1,1,i),CS( 1,1,i),U(L),dot1,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i0+1) = dot18a(Z,AS(1,2,1,i),CS( 2,1,i),U(L),dot2,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i1  ) = dot18a(Z,AS(1,3,1,i),CS( 4,1,i),U(L),dot3,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i1+1) = dot18a(Z,AS(1,4,1,i),CS( 7,1,i),U(L),dot4,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i2  ) = dot18a(Z,AS(1,5,1,i),CS(11,1,i),U(L),dot5,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i2+1) = dot18a(Z,AS(1,6,1,i),CS(16,1,i),U(L),dot6,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i3  ) = dot18a(Z,AS(1,7,1,i),CS(22,1,i),U(L),dot7,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i3+1) = eight*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                       - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

               j0 = j0 + 1
               j1 = j1 + 1
               j2 = j2 + 1
               j3 = j3 + 1
               j4 = j4 + 1
               j5 = j5 + 1
               L  = L  + 7
  40        continue
c								edge #3
            i0 = i0 + 2
            i1 = i1 + 2
            i2 = i2 + 2
            i3 = i3 + 2
            Z(i0)   = dot12v(Z,AE(1,1,3,i),CE( 1,3,i),U(L),dot1,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i0+1) = dot12v(Z,AE(1,2,3,i),CE( 2,3,i),U(L),dot2,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i1)   = dot12v(Z,AE(1,3,3,i),CE( 4,3,i),U(L),dot3,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i1+1) = dot12v(Z,AE(1,4,3,i),CE( 7,3,i),U(L),dot4,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i2)   = dot12v(Z,AE(1,5,3,i),CE(11,3,i),U(L),dot5,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i2+1) = dot12v(Z,AE(1,6,3,i),CE(16,3,i),U(L),dot6,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i3)   = dot12v(Z,AE(1,7,3,i),CE(22,3,i),U(L),dot7,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i3+1) = eight*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                    - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

            j0 = j0 + 2
  50     continue

         j1 = j0 + jx
         j3 = j0 + jxy
         j4 = j3 + jx
         i0 = i1 + 2
         i1 = i0 + ix
         i2 = i0 + ixy
         i3 = i2 + ix
         call vnorm( U, 7*jx )
c								corner #3
         Z(i0)   = dot8c(Z,AC(1,1,3,i),CC( 1,3,i),U(1),dot1,j0,j1,j3,j4)
         Z(i0+1) = dot8c(Z,AC(1,2,3,i),CC( 2,3,i),U(1),dot2,j0,j1,j3,j4)
         Z(i1)   = dot8c(Z,AC(1,3,3,i),CC( 4,3,i),U(1),dot3,j0,j1,j3,j4)
         Z(i1+1) = dot8c(Z,AC(1,4,3,i),CC( 7,3,i),U(1),dot4,j0,j1,j3,j4)
         Z(i2)   = dot8c(Z,AC(1,5,3,i),CC(11,3,i),U(1),dot5,j0,j1,j3,j4)
         Z(i2+1) = dot8c(Z,AC(1,6,3,i),CC(16,3,i),U(1),dot6,j0,j1,j3,j4)
         Z(i3)   = dot8c(Z,AC(1,7,3,i),CC(22,3,i),U(1),dot7,j0,j1,j3,j4)
         Z(i3+1) = eight*Z(j1) - Z(i0) - Z(i0+1) - Z(i1)
     >               - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
c								edge #4
         L = 8
         do 60 js = 1, jx-2
          i0 = i0 + 2
          i1 = i1 + 2
          i2 = i2 + 2
          i3 = i3 + 2
          Z(i0)  =dot12h(Z,AE(1,1,4,i),CE( 1,4,i),U(L),dot1,j0,j1,j3,j4)
          Z(i0+1)=dot12h(Z,AE(1,2,4,i),CE( 2,4,i),U(L),dot2,j0,j1,j3,j4)
          Z(i1)  =dot12h(Z,AE(1,3,4,i),CE( 4,4,i),U(L),dot3,j0,j1,j3,j4)
          Z(i1+1)=dot12h(Z,AE(1,4,4,i),CE( 7,4,i),U(L),dot4,j0,j1,j3,j4)
          Z(i2)  =dot12h(Z,AE(1,5,4,i),CE(11,4,i),U(L),dot5,j0,j1,j3,j4)
          Z(i2+1)=dot12h(Z,AE(1,6,4,i),CE(16,4,i),U(L),dot6,j0,j1,j3,j4)
          Z(i3)  =dot12h(Z,AE(1,7,4,i),CE(22,4,i),U(L),dot7,j0,j1,j3,j4)
          Z(i3+1)=eight*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
          j0 = j0 + 1
          j1 = j1 + 1
          j3 = j3 + 1
          j4 = j4 + 1
          L  = L  + 7
  60     continue
c								corner #4
         i0 = i0 + 2
         i1 = i1 + 2
         i2 = i2 + 2
         i3 = i3 + 2
         Z(i0)   = dot8c(Z,AC(1,1,4,i),CC( 1,4,i),U(L),dot1,j0,j1,j3,j4)
         Z(i0+1) = dot8c(Z,AC(1,2,4,i),CC( 2,4,i),U(L),dot2,j0,j1,j3,j4)
         Z(i1)   = dot8c(Z,AC(1,3,4,i),CC( 4,4,i),U(L),dot3,j0,j1,j3,j4)
         Z(i1+1) = dot8c(Z,AC(1,4,4,i),CC( 7,4,i),U(L),dot4,j0,j1,j3,j4)
         Z(i2)   = dot8c(Z,AC(1,5,4,i),CC(11,4,i),U(L),dot5,j0,j1,j3,j4)
         Z(i2+1) = dot8c(Z,AC(1,6,4,i),CC(16,4,i),U(L),dot6,j0,j1,j3,j4)
         Z(i3)   = dot8c(Z,AC(1,7,4,i),CC(22,4,i),U(L),dot7,j0,j1,j3,j4)
         Z(i3+1) = eight*Z(j1+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                 - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

         j00 = jq + 1
         do 110 ls = 1, jz-2
            j0 = j00
            j1 = j0 + jx
            j3 = j0 + jxy
            j4 = j3 + jx
            j6 = j3 + jxy
            j7 = j6 + jx
            i0 = i3 + 2
            i1 = i0 + ix
            i2 = i0 + ixy
            i3 = i2 + ix
            call vnorm( U, 7*jx )
c								edge #5
            Z(i0)   = dot12v(Z,AE(1,1,5,i),CE( 1,5,i),U(1),dot1,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i0+1) = dot12v(Z,AE(1,2,5,i),CE( 2,5,i),U(1),dot2,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i1)   = dot12v(Z,AE(1,3,5,i),CE( 4,5,i),U(1),dot3,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i1+1) = dot12v(Z,AE(1,4,5,i),CE( 7,5,i),U(1),dot4,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i2)   = dot12v(Z,AE(1,5,5,i),CE(11,5,i),U(1),dot5,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i2+1) = dot12v(Z,AE(1,6,5,i),CE(16,5,i),U(1),dot6,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i3)   = dot12v(Z,AE(1,7,5,i),CE(22,5,i),U(1),dot7,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i3+1) = eight*Z(j3  ) - Z(i0) - Z(i0+1) - Z(i1)
     >                    - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
c								side #2
            L = 8
            do 70 js = 1, jx-2
               i0 = i0 + 2
               i1 = i1 + 2
               i2 = i2 + 2
               i3 = i3 + 2
               Z(i0  ) = dot18a(Z,AS(1,1,2,i),CS( 1,2,i),U(L),dot1,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i0+1) = dot18a(Z,AS(1,2,2,i),CS( 2,2,i),U(L),dot2,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i1  ) = dot18a(Z,AS(1,3,2,i),CS( 4,2,i),U(L),dot3,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i1+1) = dot18a(Z,AS(1,4,2,i),CS( 7,2,i),U(L),dot4,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i2  ) = dot18a(Z,AS(1,5,2,i),CS(11,2,i),U(L),dot5,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i2+1) = dot18a(Z,AS(1,6,2,i),CS(16,2,i),U(L),dot6,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i3  ) = dot18a(Z,AS(1,7,2,i),CS(22,2,i),U(L),dot7,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i3+1) = eight*Z(j3+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                       - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

               j0 = j0 + 1
               j1 = j1 + 1
               j3 = j3 + 1
               j4 = j4 + 1
               j6 = j6 + 1
               j7 = j7 + 1
               L  = L  + 7
  70        continue
c								edge #6
            i0 = i0 + 2
            i1 = i1 + 2
            i2 = i2 + 2
            i3 = i3 + 2
            Z(i0)   = dot12v(Z,AE(1,1,6,i),CE( 1,6,i),U(L),dot1,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i0+1) = dot12v(Z,AE(1,2,6,i),CE( 2,6,i),U(L),dot2,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i1)   = dot12v(Z,AE(1,3,6,i),CE( 4,6,i),U(L),dot3,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i1+1) = dot12v(Z,AE(1,4,6,i),CE( 7,6,i),U(L),dot4,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i2)   = dot12v(Z,AE(1,5,6,i),CE(11,6,i),U(L),dot5,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i2+1) = dot12v(Z,AE(1,6,6,i),CE(16,6,i),U(L),dot6,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i3)   = dot12v(Z,AE(1,7,6,i),CE(22,6,i),U(L),dot7,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i3+1) = eight*Z(j3+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                    - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

            j0 = j00
            do 90 ks = 1, jy-2
               j1 = j0 + jx
               j2 = j1 + jx
               j3 = j0 + jxy
               j4 = j3 + jx
               j5 = j4 + jx
               j6 = j3 + jxy
               j7 = j6 + jx
               j8 = j7 + jx
               i0 = i1 + 2
               i1 = i0 + ix
               i2 = i0 + ixy
               i3 = i2 + ix
               call vnorm( U, 7*jx )
c								side #3
               Z(i0  ) = dot18s(Z,AS(1,1,3,i),CS( 1,3,i),U(1),dot1,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i0+1) = dot18s(Z,AS(1,2,3,i),CS( 2,3,i),U(1),dot2,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i1  ) = dot18s(Z,AS(1,3,3,i),CS( 4,3,i),U(1),dot3,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i1+1) = dot18s(Z,AS(1,4,3,i),CS( 7,3,i),U(1),dot4,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i2  ) = dot18s(Z,AS(1,5,3,i),CS(11,3,i),U(1),dot5,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i2+1) = dot18s(Z,AS(1,6,3,i),CS(16,3,i),U(1),dot6,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i3  ) = dot18s(Z,AS(1,7,3,i),CS(22,3,i),U(1),dot7,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i3+1) = eight*Z(j4  ) - Z(i0) - Z(i0+1) - Z(i1)
     >                       - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

               L = 8
               do 80 js = 1, jx-2
                  i0 = i0 + 2
                  i1 = i1 + 2
                  i2 = i2 + 2
                  i3 = i3 + 2
c								interior
                  Z(i0  ) = dot27i(Z,AI(1,1,i),CI( 1,i),U(L),dot1,
     >                             j0,j1,j2,j3,j4,j5,j6,j7,j8)
                  Z(i0+1) = dot27i(Z,AI(1,2,i),CI( 2,i),U(L),dot2,
     >                             j0,j1,j2,j3,j4,j5,j6,j7,j8)
                  Z(i1  ) = dot27i(Z,AI(1,3,i),CI( 4,i),U(L),dot3,
     >                             j0,j1,j2,j3,j4,j5,j6,j7,j8)
                  Z(i1+1) = dot27i(Z,AI(1,4,i),CI( 7,i),U(L),dot4,
     >                             j0,j1,j2,j3,j4,j5,j6,j7,j8)
                  Z(i2  ) = dot27i(Z,AI(1,5,i),CI(11,i),U(L),dot5,
     >                             j0,j1,j2,j3,j4,j5,j6,j7,j8)
                  Z(i2+1) = dot27i(Z,AI(1,6,i),CI(16,i),U(L),dot6,
     >                             j0,j1,j2,j3,j4,j5,j6,j7,j8)
                  Z(i3  ) = dot27i(Z,AI(1,7,i),CI(22,i),U(L),dot7,
     >                             j0,j1,j2,j3,j4,j5,j6,j7,j8)
                  Z(i3+1) = eight*Z(j4+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                          - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
                  j0 = j0 + 1
                  j1 = j1 + 1
                  j2 = j2 + 1
                  j3 = j3 + 1
                  j4 = j4 + 1
                  j5 = j5 + 1
                  j6 = j6 + 1
                  j7 = j7 + 1
                  j8 = j8 + 1
                  L  = L + 7
  80           continue
c								side #4
               i0 = i0 + 2
               i1 = i1 + 2
               i2 = i2 + 2
               i3 = i3 + 2
               Z(i0  ) = dot18s(Z,AS(1,1,4,i),CS( 1,4,i),U(L),dot1,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i0+1) = dot18s(Z,AS(1,2,4,i),CS( 2,4,i),U(L),dot2,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i1  ) = dot18s(Z,AS(1,3,4,i),CS( 4,4,i),U(L),dot3,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i1+1) = dot18s(Z,AS(1,4,4,i),CS( 7,4,i),U(L),dot4,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i2  ) = dot18s(Z,AS(1,5,4,i),CS(11,4,i),U(L),dot5,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i2+1) = dot18s(Z,AS(1,6,4,i),CS(16,4,i),U(L),dot6,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i3  ) = dot18s(Z,AS(1,7,4,i),CS(22,4,i),U(L),dot7,
     >                          j0,j1,j2,j3,j4,j5,j6,j7,j8)
               Z(i3+1) = eight*Z(j4+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                       - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
               j0 = j0 + 2
  90        continue
            j1 = j0 + jx
            j3 = j0 + jxy
            j4 = j3 + jx
            j6 = j3 + jxy
            j7 = j6 + jx
            i0 = i1 + 2
            i1 = i0 + ix
            i2 = i0 + ixy
            i3 = i2 + ix
            call vnorm( U, 7*jx )
c								edge #7
            Z(i0)   = dot12v(Z,AE(1,1,7,i),CE( 1,7,i),U(1),dot1,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i0+1) = dot12v(Z,AE(1,2,7,i),CE( 2,7,i),U(1),dot2,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i1)   = dot12v(Z,AE(1,3,7,i),CE( 4,7,i),U(1),dot3,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i1+1) = dot12v(Z,AE(1,4,7,i),CE( 7,7,i),U(1),dot4,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i2)   = dot12v(Z,AE(1,5,7,i),CE(11,7,i),U(1),dot5,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i2+1) = dot12v(Z,AE(1,6,7,i),CE(16,7,i),U(1),dot6,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i3)   = dot12v(Z,AE(1,7,7,i),CE(22,7,i),U(1),dot7,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i3+1) = eight*Z(j4  ) - Z(i0) - Z(i0+1) - Z(i1)
     >                    - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
c								side #5
            L = 8
            do 100 js = 1, jx-2
               i0 = i0 + 2
               i1 = i1 + 2
               i2 = i2 + 2
               i3 = i3 + 2
               Z(i0  ) = dot18a(Z,AS(1,1,5,i),CS( 1,5,i),U(L),dot1,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i0+1) = dot18a(Z,AS(1,2,5,i),CS( 2,5,i),U(L),dot2,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i1  ) = dot18a(Z,AS(1,3,5,i),CS( 4,5,i),U(L),dot3,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i1+1) = dot18a(Z,AS(1,4,5,i),CS( 7,5,i),U(L),dot4,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i2  ) = dot18a(Z,AS(1,5,5,i),CS(11,5,i),U(L),dot5,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i2+1) = dot18a(Z,AS(1,6,5,i),CS(16,5,i),U(L),dot6,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i3  ) = dot18a(Z,AS(1,7,5,i),CS(22,5,i),U(L),dot7,
     >                          j0,j1,j3,j4,j6,j7)
               Z(i3+1) = eight*Z(j4+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                       - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

               j0 = j0 + 1
               j1 = j1 + 1
               j3 = j3 + 1
               j4 = j4 + 1
               j6 = j6 + 1
               j7 = j7 + 1
               L  = L  + 7
 100        continue
c								edge #8
            i0 = i0 + 2
            i1 = i1 + 2
            i2 = i2 + 2
            i3 = i3 + 2
            Z(i0)   = dot12v(Z,AE(1,1,8,i),CE( 1,8,i),U(L),dot1,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i0+1) = dot12v(Z,AE(1,2,8,i),CE( 2,8,i),U(L),dot2,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i1)   = dot12v(Z,AE(1,3,8,i),CE( 4,8,i),U(L),dot3,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i1+1) = dot12v(Z,AE(1,4,8,i),CE( 7,8,i),U(L),dot4,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i2)   = dot12v(Z,AE(1,5,8,i),CE(11,8,i),U(L),dot5,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i2+1) = dot12v(Z,AE(1,6,8,i),CE(16,8,i),U(L),dot6,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i3)   = dot12v(Z,AE(1,7,8,i),CE(22,8,i),U(L),dot7,
     >                       j0,j1,j3,j4,j6,j7)
            Z(i3+1) = eight*Z(j4+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                    - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

            j00 = j00 + jxy
 110     continue
         j0 = j00
         j1 = j0 + jx
         j3 = j0 + jxy
         j4 = j3 + jx
         i0 = i3 + 2
         i1 = i0 + ix
         i2 = i0 + ixy
         i3 = i2 + ix
         call vnorm( U, 7*jx )
c								corner #5
         Z(i0)   = dot8c(Z,AC(1,1,5,i),CC( 1,5,i),U(1),dot1,j0,j1,j3,j4)
         Z(i0+1) = dot8c(Z,AC(1,2,5,i),CC( 2,5,i),U(1),dot2,j0,j1,j3,j4)
         Z(i1)   = dot8c(Z,AC(1,3,5,i),CC( 4,5,i),U(1),dot3,j0,j1,j3,j4)
         Z(i1+1) = dot8c(Z,AC(1,4,5,i),CC( 7,5,i),U(1),dot4,j0,j1,j3,j4)
         Z(i2)   = dot8c(Z,AC(1,5,5,i),CC(11,5,i),U(1),dot5,j0,j1,j3,j4)
         Z(i2+1) = dot8c(Z,AC(1,6,5,i),CC(16,5,i),U(1),dot6,j0,j1,j3,j4)
         Z(i3)   = dot8c(Z,AC(1,7,5,i),CC(22,5,i),U(1),dot7,j0,j1,j3,j4)
         Z(i3+1) = eight*Z(j3) - Z(i0) - Z(i0+1) - Z(i1)
     >               - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
c								edge #9
         L = 8
         do 120 js = 1, jx-2
          i0 = i0 + 2
          i1 = i1 + 2
          i2 = i2 + 2
          i3 = i3 + 2
          Z(i0)  =dot12h(Z,AE(1,1,9,i),CE( 1,9,i),U(L),dot1,j0,j1,j3,j4)
          Z(i0+1)=dot12h(Z,AE(1,2,9,i),CE( 2,9,i),U(L),dot2,j0,j1,j3,j4)
          Z(i1)  =dot12h(Z,AE(1,3,9,i),CE( 4,9,i),U(L),dot3,j0,j1,j3,j4)
          Z(i1+1)=dot12h(Z,AE(1,4,9,i),CE( 7,9,i),U(L),dot4,j0,j1,j3,j4)
          Z(i2)  =dot12h(Z,AE(1,5,9,i),CE(11,9,i),U(L),dot5,j0,j1,j3,j4)
          Z(i2+1)=dot12h(Z,AE(1,6,9,i),CE(16,9,i),U(L),dot6,j0,j1,j3,j4)
          Z(i3)  =dot12h(Z,AE(1,7,9,i),CE(22,9,i),U(L),dot7,j0,j1,j3,j4)
          Z(i3+1)=eight*Z(j3+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
          j0 = j0 + 1
          j1 = j1 + 1
          j3 = j3 + 1
          j4 = j4 + 1
          L  = L  + 7
 120     continue
c								corner #6
         i0 = i0 + 2
         i1 = i1 + 2
         i2 = i2 + 2
         i3 = i3 + 2
         Z(i0)   = dot8c(Z,AC(1,1,6,i),CC( 1,6,i),U(L),dot1,j0,j1,j3,j4)
         Z(i0+1) = dot8c(Z,AC(1,2,6,i),CC( 2,6,i),U(L),dot2,j0,j1,j3,j4)
         Z(i1)   = dot8c(Z,AC(1,3,6,i),CC( 4,6,i),U(L),dot3,j0,j1,j3,j4)
         Z(i1+1) = dot8c(Z,AC(1,4,6,i),CC( 7,6,i),U(L),dot4,j0,j1,j3,j4)
         Z(i2)   = dot8c(Z,AC(1,5,6,i),CC(11,6,i),U(L),dot5,j0,j1,j3,j4)
         Z(i2+1) = dot8c(Z,AC(1,6,6,i),CC(16,6,i),U(L),dot6,j0,j1,j3,j4)
         Z(i3)   = dot8c(Z,AC(1,7,6,i),CC(22,6,i),U(L),dot7,j0,j1,j3,j4)
         Z(i3+1) = eight*Z(j3+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                 - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

         j0 = j00
         do 140 ks = 1, jy-2
            j1 = j0 + jx
            j2 = j1 + jx
            j3 = j0 + jxy
            j4 = j3 + jx
            j5 = j4 + jx
            i0 = i1 + 2
            i1 = i0 + ix
            i2 = i0 + ixy
            i3 = i2 + ix
            call vnorm( U, 7*jx )
c								edge #10
            Z(i0)   = dot12v(Z,AE(1,1,10,i),CE( 1,10,i),U(1),dot1,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i0+1) = dot12v(Z,AE(1,2,10,i),CE( 2,10,i),U(1),dot2,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i1)   = dot12v(Z,AE(1,3,10,i),CE( 4,10,i),U(1),dot3,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i1+1) = dot12v(Z,AE(1,4,10,i),CE( 7,10,i),U(1),dot4,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i2)   = dot12v(Z,AE(1,5,10,i),CE(11,10,i),U(1),dot5,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i2+1) = dot12v(Z,AE(1,6,10,i),CE(16,10,i),U(1),dot6,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i3)   = dot12v(Z,AE(1,7,10,i),CE(22,10,i),U(1),dot7,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i3+1) = eight*Z(j4  ) - Z(i0) - Z(i0+1) - Z(i1)
     >                    - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
c								side #6
            L = 8
            do 130 js = 1, jx-2
               i0 = i0 + 2
               i1 = i1 + 2
               i2 = i2 + 2
               i3 = i3 + 2
               Z(i0  ) = dot18a(Z,AS(1,1,6,i),CS( 1,6,i),U(L),dot1,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i0+1) = dot18a(Z,AS(1,2,6,i),CS( 2,6,i),U(L),dot2,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i1  ) = dot18a(Z,AS(1,3,6,i),CS( 4,6,i),U(L),dot3,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i1+1) = dot18a(Z,AS(1,4,6,i),CS( 7,6,i),U(L),dot4,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i2  ) = dot18a(Z,AS(1,5,6,i),CS(11,6,i),U(L),dot5,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i2+1) = dot18a(Z,AS(1,6,6,i),CS(16,6,i),U(L),dot6,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i3  ) = dot18a(Z,AS(1,7,6,i),CS(22,6,i),U(L),dot7,
     >                          j0,j1,j2,j3,j4,j5)
               Z(i3+1) = eight*Z(j4+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                       - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

               j0 = j0 + 1
               j1 = j1 + 1
               j2 = j2 + 1
               j3 = j3 + 1
               j4 = j4 + 1
               j5 = j5 + 1
               L  = L  + 7
 130        continue
c								edge #11
            i0 = i0 + 2
            i1 = i1 + 2
            i2 = i2 + 2
            i3 = i3 + 2
            Z(i0)   = dot12v(Z,AE(1,1,11,i),CE( 1,11,i),U(L),dot1,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i0+1) = dot12v(Z,AE(1,2,11,i),CE( 2,11,i),U(L),dot2,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i1)   = dot12v(Z,AE(1,3,11,i),CE( 4,11,i),U(L),dot3,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i1+1) = dot12v(Z,AE(1,4,11,i),CE( 7,11,i),U(L),dot4,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i2)   = dot12v(Z,AE(1,5,11,i),CE(11,11,i),U(L),dot5,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i2+1) = dot12v(Z,AE(1,6,11,i),CE(16,11,i),U(L),dot6,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i3)   = dot12v(Z,AE(1,7,11,i),CE(22,11,i),U(L),dot7,
     >                       j0,j1,j2,j3,j4,j5)
            Z(i3+1) = eight*Z(j4+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                    - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

            j0 = j0 + 2
 140     continue

         j1 = j0 + jx
         j3 = j0 + jxy
         j4 = j3 + jx
         i0 = i1 + 2
         i1 = i0 + ix
         i2 = i0 + ixy
         i3 = i2 + ix
         call vnorm( U, 7*jx )
c								corner #7
         Z(i0)   = dot8c(Z,AC(1,1,7,i),CC( 1,7,i),U(1),dot1,j0,j1,j3,j4)
         Z(i0+1) = dot8c(Z,AC(1,2,7,i),CC( 2,7,i),U(1),dot2,j0,j1,j3,j4)
         Z(i1)   = dot8c(Z,AC(1,3,7,i),CC( 4,7,i),U(1),dot3,j0,j1,j3,j4)
         Z(i1+1) = dot8c(Z,AC(1,4,7,i),CC( 7,7,i),U(1),dot4,j0,j1,j3,j4)
         Z(i2)   = dot8c(Z,AC(1,5,7,i),CC(11,7,i),U(1),dot5,j0,j1,j3,j4)
         Z(i2+1) = dot8c(Z,AC(1,6,7,i),CC(16,7,i),U(1),dot6,j0,j1,j3,j4)
         Z(i3)   = dot8c(Z,AC(1,7,7,i),CC(22,7,i),U(1),dot7,j0,j1,j3,j4)
         Z(i3+1) = eight*Z(j4) - Z(i0) - Z(i0+1) - Z(i1)
     >               - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
c								edge #12
         L = 8
         do 150 js = 1, jx-2
          i0 = i0 + 2
          i1 = i1 + 2
          i2 = i2 + 2
          i3 = i3 + 2
          Z(i0)  =dot12h(Z,AE(1,1,12,i),CE( 1,12,i),U(L),dot1,
     >                   j0,j1,j3,j4)
          Z(i0+1)=dot12h(Z,AE(1,2,12,i),CE( 2,12,i),U(L),dot2,
     >                   j0,j1,j3,j4)
          Z(i1)  =dot12h(Z,AE(1,3,12,i),CE( 4,12,i),U(L),dot3,
     >                   j0,j1,j3,j4)
          Z(i1+1)=dot12h(Z,AE(1,4,12,i),CE( 7,12,i),U(L),dot4,
     >                   j0,j1,j3,j4)
          Z(i2)  =dot12h(Z,AE(1,5,12,i),CE(11,12,i),U(L),dot5,
     >                   j0,j1,j3,j4)
          Z(i2+1)=dot12h(Z,AE(1,6,12,i),CE(16,12,i),U(L),dot6,
     >                   j0,j1,j3,j4)
          Z(i3)  =dot12h(Z,AE(1,7,12,i),CE(22,12,i),U(L),dot7,
     >                   j0,j1,j3,j4)
          Z(i3+1)=eight*Z(j4+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)
          j0 = j0 + 1
          j1 = j1 + 1
          j3 = j3 + 1
          j4 = j4 + 1
          L  = L  + 7
 150     continue
c								corner #8
         i0 = i0 + 2
         i1 = i1 + 2
         i2 = i2 + 2
         i3 = i3 + 2
         Z(i0)   = dot8c(Z,AC(1,1,8,i),CC( 1,8,i),U(L),dot1,j0,j1,j3,j4)
         Z(i0+1) = dot8c(Z,AC(1,2,8,i),CC( 2,8,i),U(L),dot2,j0,j1,j3,j4)
         Z(i1)   = dot8c(Z,AC(1,3,8,i),CC( 4,8,i),U(L),dot3,j0,j1,j3,j4)
         Z(i1+1) = dot8c(Z,AC(1,4,8,i),CC( 7,8,i),U(L),dot4,j0,j1,j3,j4)
         Z(i2)   = dot8c(Z,AC(1,5,8,i),CC(11,8,i),U(L),dot5,j0,j1,j3,j4)
         Z(i2+1) = dot8c(Z,AC(1,6,8,i),CC(16,8,i),U(L),dot6,j0,j1,j3,j4)
         Z(i3)   = dot8c(Z,AC(1,7,8,i),CC(22,8,i),U(L),dot7,j0,j1,j3,j4)
         Z(i3+1) = eight*Z(j4+1) - Z(i0) - Z(i0+1) - Z(i1)
     >                 - Z(i1+1) - Z(i2) - Z(i2+1) - Z(i3)

         jx = ix
         jy = iy
         jz = iz
 160  continue
c						all done, compute elapsed time
      ts = ts + (second() - tt)

      return
      end
