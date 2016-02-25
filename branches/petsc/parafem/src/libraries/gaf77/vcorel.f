c  *********************************************************************
c  *                                                                   *
c  *                         subroutine vcorel                         *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.01
c  Written by Gordon A. Fenton, TUNS, Thu Jul 18 11:57:34 1996
c  Latest Update: Aug 12, 1997
c
c  PURPOSE  simulates a set of correlated Gaussian random variables given its
c           covariance matrix
c
c  This routine produces a realization of a set of correlated, normally
c  distributed random variables (stationary) given the covariance matrix
c  associated with the set of RV's. This is accomplished using a Cholesky
c  decomposition of the covariance matrix.
c  NOTES:
c    1)	the pseudo-random number generator used by vnorm must be initialized
c	before calling this routine.
c    2)	although this is nominally a single precision routine, the array C
c	is double precision to minimize numerical round-off problems, which
c	covariance matrices are notorious for, in its Cholesky decomposition.
c
c  Arguments to this routine are as follows;
c
c	X	real vector of length at least n, which on output will contain
c		the correlated realization. (output)
c
c   xmean	the process mean (assumed stationary). (input)
c
c	C	real*8 array of size at least n x n containing on input
c		the process covariance matrix. On the first call to this
c		routine, the Cholesky decomposition of C is computed,
c		replacing it's upper triangle. On this and subsequent
c		calls, C should not be changed, since this decomposition
c		is used to produce the realizations. If C is recomputed
c		for another problem, this routine can be explicitly instructed
c		to recompute the decomposition by setting init = 1.
c		(input/output)
c
c      ic	the leading dimension of the array C as prescribed in the
c		calling routine. (input)
c
c	n	the number of elements in the vector X and the size of the
c		array C. (input)
c
c    init	integer flag which should be set to 1 if the Cholesky
c		decomposition of the matrix C is to be performed on this
c		call. Subsequent calls (for subsequent realizations of the
c		same process) should have init not equal to 1. (input)
c
c    ierr	integer error flag which is set to 0 if all goes well.
c		Otherwise ierr is set to 1 if the relative error in the
c		Cholesky decomposition (see dchol2) exceeds 20%. (output)
c
c  REVISION HISTORY:
c  1.01	eliminated unused variable 'zero' (Aug 12/97)
c-------------------------------------------------------------------------
      subroutine vcorel( X, xmean, C, ic, n, init, ierr )
      real X(*)
      real*8 C(ic,*), dble, t
      data emax/0.2/, ifirst/1/
c					assume all will go well
      ierr = 0
c					compute Cholesky Decomposition?
      if( (ifirst .eq. 1) .or. (init .eq. 1) ) then
         ifirst = 0
         call dchol2( C, ic, n, rerr )
         if( rerr .gt. emax ) ierr = 1
      endif
c					get independent standard normal vars
      call vnorm( X, n )
c					form into correlated variates
      do 20 j = n, 1, -1
         t = dble(xmean)
         do 10 i = 1, j
            t = t + C(i,j)*dble(X(i))
  10     continue
         X(j) = sngl(t)
  20  continue
c					all done
      return
      end
