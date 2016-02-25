c  *******************************************************************
c  *                                                                 *
c  *                       Subroutine las1g                          *
c  *                                                                 *
c  *******************************************************************
c  Single Precision Version 3.34
c  Written by Gordon A. Fenton, TUNS, Feb. 4, 1993
c  Latest Update: Jun 3, 2006
c
c  PURPOSE  generates a realization of a 1-D Gaussian stationary random
c           local average process
c
c  This routine creates a realization of a zero-mean 1-D random process given
c  its variance function (see DVARFN below). Each discrete
c  value generated represents the local average of a realization of
c  the process over the cell length XL/N, where XL is the physical length
c  of the process and N is the number of cells desired. N must be such that
c  N = k1*2**m, for some integer k1 in the range [1,256] and for some integer
c  m in the range [0,16]. Thus the largest value of N is 2^(24) = 16,777,216.
c  The values of k1 and m are computed internally. The construction of the
c  realization proceeds in a recursive fashion as follows;
c    1) generate the first k1 cells directly,
c    2) subdivide each of the k1 cells into two equal parts and generate
c       random values for each part preserving the cell average and closely
c       approximating the covariance structure (the latter is exact within
c       the cell),
c    3) subdivide the 2*k1 cells obtained in step (2) into two equal parts
c       and again generate random values for each part,
c    4) continue the subdivision process until the domain is divided into
c       N equal cells.
c
c  Note that this routine sets up a number of parameters on the first
c  call and thus the time required to produce the first realization
c  is substantially greater than on subsequent calls (see INIT). For
c  more details on this Local Average Subdivision algorithm, see
c
c    1) Fenton, G.A., and Vanmarcke, E.H., "Simulation of Random Fields
c          via Local Average Subdivision", ASCE Journal of Engineering
c          Mechanics, 116(8), 1733-1749, 1990.
c    2) Fenton, G.A., Simulation and Analysis of Random Fields, Ph.D. Thesis,
c          Dept. of Civil Engineering and Operations Research, Princeton
c          University, Princeton, NJ, 1990.
c    3) Vanmarcke, E.H., Random Fields: Analysis and Synthesis, MIT Press,
c          Boston, Massachusetts, 1984.
c    4) Fenton, G.A., Error evaluation of three random field generators,
c          ASCE Journal of Engineering Mechanics (to appear), 1993.
c
c  Arguments to this routine are as follows;
c
c      Z  real vector of length at least (3*N/2) which on output will contain
c         the desired realization in the first N locations. The remainder is
c         used for workspace. (output)
c
c      N  desired number of cells discretizing the process. N must be such
c         that it can be written N = k1*2**m, k1 and m are integers with
c         1 <= k1 <= MXK and 0 <= m <= MXM (see parameters below). (input)
c
c     XL  physical length of the process. (input)
c
c DVARFN  external real*8 function which returns the covariance of the random
c         process between two points separated by a distance V1.
c         DVARFN is referenced as follows
c
c                var = dvarfn( V1 )
c
c         where (V1) is the separation distance. Any other
c         parameters to the function must be passed by common block from the
c         calling routine. The current version of LAS uses the sign on `var'
c         to tell dvarfn to return either the covariance, or the variance of
c         a local average over distance V1. The latter is returned if var < 0,
c         although this feature is no longer used by LAS.
c
c  NBH    integer giving the desired size of the neighborhood to include
c         in the determination of the best linear estimate of the mean value
c         of a sub-cell. At present only two neighborhood sizes are supported:
c         NBH = 3 or 5. Note that the code will run significantly slower using
c         NBH = 5, but it yields improved results for some types of processes
c         (such as damped oscillatory noise). (input)
c
c  KSEED  integer seed to be used to initialize the pseudo-random number
c         generator (which is assumed to be performed by calling rand(kseed)).
c         If KSEED = 0, then a random seed will be used (based on the
c         clock time when this routine is called for the first time).
c         On output, KSEED is set to the value of the actual seed used.
c         (input/output)
c
c  INIT   integer flag which must be set to +/- 1 when parameters of the
c         process are to be calculated or recalculated. If multiple
c         realizations of the same process are desired, then subsequent
c         calls should use INIT not equal to 1 and N1 the same as used
c         initially. If INIT = +1, then both the parameters are calculated
c         and a random realization is returned. If INIT = -1, then just
c         the initial parameters are calculated. (input)
c
c  IOUT   unit number to which error and warning messages are to be logged.
c         (input)
c  ---------------------------------------------------------------------------
c
c  PARAMETERS:
c
c    MXK  the maximum value of k1. Currently 256.
c
c    MXM  represents the maximum number of subdivisions (ie the maximum value
c         of m) that the routine can carry out. Currently 16. Note that this
c         means that N cannot exceed MXK*2**MXM = 16*2**16 = 1,048,576.
c
c    NGS  this routine will generate NGS independent Gaussian random variates
c         at a time (this is because some machines slow down drastically when
c         making multiple calls to an external routine - that is we attempt to
c         minimize the number of calls to VNORM herein).
c
c  Notes:
c    1) Simulation timing is available through common block LASTYM where
c       TI is the time required to set up the parameters of the process
c       and TS is the cumulative simulation time. These timings are
c       obtained through the function SECOND which returns elapsed user
c       time in seconds since the start of the program.
c    2) all variables that start with "d" are double precision - this is
c       used to perform the covariance calculations in extended precision.
c    3) The parameter `tol' is the maximum relative error allowed on the
c       Cholesky decomposition of covariance matrices before a warning
c       message is emitted. The relative error is estimated by computing
c       the lower-right most element of L*L' and comparing to the original
c       element of A (where L*L' = A is the decomposition). This give some
c       measure of the roundoff errors accumulated in the computation of L.
c
c  Requires:
c    1) from libGAFsim:	ISEED, DCVIT1, DCHOL2, DSIFA, DSISL, DAXPY, DSWAP,
c			IDAMAX, DDOT, LAS1I, VNORM, RANDF, SECOND
c    2) external user defined variance function (see DVARFN).
c
c  REVISION HISTORY:
c  3.31	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c  3.32	modified writeup above for new dvarfn return value. (Apr 17/01)
c  3.33	increased MXK from 16 to 256. (Aug 10/05)
c  3.34	just initialize if init = -1. (Jun 3/06)
c--------------------------------------------------------------------
      subroutine las1g( z, n, xl, dvarfn, nbh, kseed, init, iout )
      implicit real*8 (d)
      parameter( MXM = 16, MXK = 256, NGS = 4096 )
      real z(*), A(9,MXM), C(3,MXM), C0(MXK,MXK), g(NGS), ge(4*MXM)
      external dvarfn
      save A, C, C0, m, k1
      common/LASTYM/ ti, ts
      data ifirst/1/, two/2.0/
      data tol/1.e-3/
c------------------------------ initialization ---------------------------

      if( ifirst .eq. 1 .or. iabs(init) .eq. 1 ) then
         ti = second()
         ifirst = 0
c					compute parameters for this routine
         call las1i( dvarfn, n, xl, MXM, MXK, nbh, A, C, C0,
     >                m, k1, iout, tol )
c					initialize the random number generator
         kseed = iseed( kseed )
c					set timers
         ts = 0.
         ti = second() - ti
         if( init .eq. -1 ) return
      endif
c
c------------------------------ create the realization ---------------------
c
      tt = second()
c					create first k1 cells directly
      if( mod(m,2) .eq. 0 ) then
         iz = 0
         jz = n
      else
         iz = n
         jz = 0
      endif
c					generate stage 0 field
      call vnorm( g, k1 )
      do 20 i = 1, k1
         Z(iz+i) = C0(1,i)*g(1)
         do 10 j = 2, i
            Z(iz+i) = Z(iz+i) + C0(j,i)*g(j)
  10     continue
  20  continue
c					generate stage 1, 2, ... M fields
      jx = k1
      if ( nbh .eq. 3 ) then
c					for a neighborhood of 3 (NBH = 3)...
         call vnorm( ge, 2*M )
         do 40 i = 1, M
            ii = 2*i
            it = jz
            jz = iz
            iz = it
            i0 = iz + 1
            j0 = jz + 1
            Z(i0)   = A(1,i)*Z(j0) + A(2,i)*Z(j0+1) + C(1,i)*ge(ii-1)
            Z(i0+1) = two*Z(j0) - Z(i0)
            do 30 jj = 1, jx-2, NGS
               kk = min0( jx-jj-1, NGS )
               call vnorm( g, kk )
               do 30 j = 1, kk
                  i0      = i0 + 2
                  Z(i0)   = A(3,i)*(Z(j0)-Z(j0+2))+Z(j0+1) + C(2,i)*g(j)
                  Z(i0+1) = two*Z(j0+1) - Z(i0)
                  j0      = j0 + 1
  30        continue
            i0      = i0 + 2
            Z(i0+1) = A(1,i)*Z(j0+1) + A(2,i)*Z(j0) + C(1,i)*ge(ii)
            Z(i0)   = two*Z(j0+1) - Z(i0+1)
            jx      = 2*jx
  40     continue
      else
c					for a neighborhood of 5 (NBH = 5)...
         call vnorm( ge, 4*M )
         do 60 i = 1, M
            ii = 4*i
            it = jz
            jz = iz
            iz = it
            i0 = iz + 1
            j0 = jz + 1
            Z(i0)   = A(1,i)*Z(j0) + A(2,i)*Z(j0+1) + A(3,i)*Z(j0+2)
     >                             + C(1,i)*ge(ii-3)
            Z(i0+1) = two*Z(j0) - Z(i0)
            i0 = i0 + 2
            Z(i0)   = A(4,i)*Z(j0) + A(5,i)*Z(j0+1) + A(6,i)*Z(j0+2)
     >                             + A(7,i)*Z(j0+3) + C(2,i)*ge(ii-2)
            Z(i0+1) = two*Z(j0+1) - Z(i0)
            do 50 jj = 1, jx-4, NGS
               kk = min0( jx-jj-3, NGS )
               call vnorm( g, kk )
               do 50 j = 1, kk
                  i0      = i0 + 2
                  Z(i0)   = A(8,i)*(Z(j0  ) - Z(j0+4)) + Z(j0+2)
     >                    + A(9,i)*(Z(j0+1) - Z(j0+3)) + C(3,i)*g(j)
                  Z(i0+1) = two*Z(j0+2) - Z(i0)
                  j0      = j0 + 1
  50        continue
            i0      = i0 + 2
            Z(i0+1) = A(4,i)*Z(j0+3) + A(5,i)*Z(j0+2) + A(6,i)*Z(j0+1)
     >                               + A(7,i)*Z(j0) + C(2,i)*ge(ii-1)
            Z(i0)   = two*Z(j0+2) - Z(i0+1)
            i0      = i0 + 2
            Z(i0+1) = A(1,i)*Z(j0+3) + A(2,i)*Z(j0+2) + A(3,i)*Z(j0+1)
     >                               + C(1,i)*ge(ii)
            Z(i0)   = two*Z(j0+3) - Z(i0+1)
            jx      = 2*jx
  60     continue
      endif
c					all done, update timer
      ts = ts + (second() - tt)

      return
      end
