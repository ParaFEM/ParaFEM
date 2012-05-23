c  *********************************************************************
c  *                                                                   *
c  *                       subroutine dcvit                            *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.31
c  Written by Gordon A. Fenton, TUNS, Feb. 11, 1993
c  Latest Update: Apr 17, 2001
c
c  PURPOSE  computes the covariance vector between a sequence of
c           local averages in 1-D. Used by LAS1G.
c
c  This routine computes the covariances between a sequence of adjacent
c  equal-sized local averages of a 1-D random function. As well, the
c  covariances between the half cell 1' and the full cells 1, 2, ..., NBH
c  are computed. That is, for local average cells arranged as follows,
c
c    |--T--|
c     ------------------------------------
c    |  1  |  2  |  3  |  4  | ... |   k  |
c     ------------------------------------
c
c     ------------------------------------
c    |  1  |  2  |1'|2'|  4  | ... |   k  | (subdivided cells 1' and 2')
c     ------------------------------------
c
c  and where Z_i represents the local average of cell i, this routine
c  finds the vector {R} such that
c
c    R_j = E[ Z_1 * Z_j ]  ( = Cov[Z_1,Z_j] for Z zero-mean)
c
c  and the vector {S} such that
c
c    S_j = E[ Z_1' * Z_j]  ( = Cov[Z_1',Z_j] for Z zero-mean);
c
c  The actual location of the subdivided cell is given by the index (NBH+1)/2
c  (ie, NBH = 3 implies that the cells 1' and 2' occur in the full cell
c  number 2). Thus NBH must be odd.
c
c  Arguments to this routine are as follows;
c
c    vfn      external real*8 function which returns the covariance of the
c             random process between two points separated by a distance T.
c             VFN is referenced as follows
c
c                C = vfn(T)
c
c             where (T) is the separation distance. Any other parameters
c             to the function must be passed by common block from the
c             calling routine. The current version of LAS uses the sign on
c             `var' to tell vfn to return either the covariance, or the
c             variance of a local average over distance V1. The latter is
c             returned if var < 0, although this feature is no longer used
c             by LAS. VFN must be an even function, ie vfn(-T) = vfn(T).
c
c    R        real vector of length at least k which on output will
c             contain the covariances discussed above. (output)
c
c    S        real vector of length at least NBH which on output will contain
c             the covariances between the half cell 1' and the full cells
c             1, 2, ..., NBH. (output)
c
c    k        the number of adjacent cells considered in the calculation of
c             R. (input)
c
c    NBH      the cell which is subdivided in the calculation of {S} is
c             (NBH+1)/2. NBH must be odd. (input)
c
c    T        the cell dimension. (input)
c
c  REVISION HISTORY:
c  1.11	replaced dummy dimensions with a (*) for GNU's compiler (Jun 9/99)
c  1.2	now including a Gaussian Quadrature integration option as an
c	alternative to the variance function approach to evaluate
c	covariances between local averages. (Jun 16/00)
c  1.3	revised above docs to reflect elimination of lvarfn (Mar 27/01)
c  1.31	modified writeup above for new vfn return value. (Apr 17/01)
c---------------------------------------------------------------------------
      subroutine dcvit1( vfn, R, S, k, NBH, T )
      implicit real*8 (a-h,o-z)
      dimension R(*), S(*)
      external vfn
      data half/0.5d0/

c					compute cell 1 - cell i covariances
      do 10 i = 1, k
         C1   = dble(i-1)
         R(i) = dcvaa1(vfn,T,C1)
  10  continue
c					compute cell 1' - cell i covariances
      n1 = -1 - 2*NBH
      t2 = half*T
      do 20 i = 1, NBH
         C1   = half*dble(n1+4*i)
         S(i) = dcvab1(vfn,t2,C1)
  20  continue

      return
      end
