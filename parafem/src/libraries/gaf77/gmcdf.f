c  *********************************************************************
c  *                                                                   *
c  *                       subroutine gmcdf                            *
c  *                                                                   *
c  *********************************************************************
c  Single Precision Version 1.1
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: Dec 5, 1996
c
c  PURPOSE  computes bounds on the cumulative distribution function of an
c           n-variate jointly normal random vector using GRCDF.
c
c  This routine computes an upper and lower bound on the cumulative
c  distribution function of a jointly standard normal n-variate. The bounds
c  are obtained by approximating the correlation matrix R by a special
c  multiplicative form, ie R_{ij} ~ g_i*g_j + (1 - g_i*g_i)*d_{ij},
c  where d_{ij} is the Chronecker delta and g_i is some constant between 0 and
c  1. For such an approximation, the multi-dimensional integral reduces to
c  a one-dimensional integral. However, the vector {g} will only uniquely
c  represent a general correlation matrix of order 3 or less, higher order
c  matrices can only be approximated. Two approximations are made, one in
c  which the elements of {g} are chosen using a maximum ratio technique,
c  the other using a minimum ratio technique (see Madsen etal, "Methods of
c  Structural Safety", Prentice-Hall, pg 62) - these yield the upper and
c  lower bounds on the multi-dimensional integral respectively. See also
c  GRCDF.
c  NOTE: at this time, the elements of the correlation matrix R must all be
c        positive to ensure real results.
c
c  Arguments to this routine are as follows;
c
c    B     real vector of length at least N containing the upper bounds of the
c          integration. B is re-ordered so that smallest elements appear first.
c          The lower bound of the integration is -infinity (ie this is a
c          distribution function). (input/output)
c
c    N     integer giving the number of elements in the random vector for
c          which the cumulative distribution is desired. (input)
c
c    R     real array of size at least N x N containing the elements of the
c          correlation matrix defining the n-variate process. R is re-ordered
c          to remain consistent with the re-ordering of B. (input/output)
c
c    ir    integer giving the leading dimension of the array R as specified
c          in the calling routine. (input)
c
c    G     real vector of length at least N to contain the vector used
c          to approximate R. (output)
c
c    bl    real value giving the lower bound on the integral. (output)
c
c    bu    real value giving the upper bound on the integral. (output)
c
c    ierr  integer error flag taking possible values;
c           =  0  if all goes well
c           = -1  if the integral cannot be calculated by GRCDF (usually means
c                 that a cubic spline could not be fitted - one or more data
c                 points are identical?)
c           =  1  if, in the computation of the upper bound, one of the
c                 elements of {G} was found to be greater than 1.0. The
c                 lower bound is still valid, but the upper bound is set
c                 to zero.
c           =  2  if, in the computation of the lower bound, one of the
c                 elements of {G} was found to be greater than 1.0. The
c                 upper bound is still valid, but the lower bound is set
c                 to zero.
c           = -2  if elements of {G} were found to be greater than 1.0 for
c                 both the upper and lower bounds. No result is possible.
c
c  REVISION HISTORY:
c  1.1	eliminated unused argument `T' (Dec 5/96)
c-----------------------------------------------------------------------------
      subroutine gmcdf( B, N, R, ir, G, bl, bu, ierr )
      dimension B(*), R(ir,*), G(*)
      data zero/0.0/, half/0.50/, one/1.0/, two/2.0/
c							if N < 4, unique...
      ierr = 0
      if( N .le. 1 ) then
         bl = phi( B(1) )
         bu = bl
      elseif( N .eq. 2 ) then
         G(1) = half
         G(2) = two*R(1,2)
         bl   = grcdf( B, G, N, ierr )
         bu   = bl
      elseif( N .eq. 3 ) then
         G(1) = sqrt(R(1,2)*R(1,3)/R(2,3))
         G(2) = sqrt(R(1,2)*R(2,3)/R(1,3))
         G(3) = sqrt(R(1,3)*R(2,3)/R(1,2))
         bl   = grcdf( B, G, N, ierr )
         bu   = bl
      else
c							set defaults
         bl = zero
         bu = zero
c							re-order B and R
         do 40 i = 1, N
            bmin = abs(B(i))
            k = i
            do 10 j = i+1, N
               if( abs(B(j)) .lt. bmin ) then
                  bmin = abs(B(j))
                  k    = j
               endif
  10        continue
c							swap B elements
            if( k .ne. i ) then
               tmp  = B(i)
               B(i) = B(k)
               B(k) = tmp
c							swap R rows
               do 20 j = 1, N
                  tmp    = R(i,j)
                  R(i,j) = R(k,j)
                  R(k,j) = tmp
  20           continue
c							swap R columns
               do 30 j = 1, N
                  tmp    = R(j,i)
                  R(j,i) = R(j,k)
                  R(j,k) = tmp
  30           continue
            endif
  40     continue
c							find elements of {G}
         G(1) = sqrt(R(1,2)*R(1,3)/R(2,3))
         G(2) = sqrt(R(1,2)*R(2,3)/R(1,3))
         G(3) = sqrt(R(1,3)*R(2,3)/R(1,2))
c							for upper bound...
         do 60 i = 4, N
            G(i) = R(i,1)/G(1)
            do 50 j = 2, i-1
               tmp = R(i,j)/G(j)
               if( tmp .gt. G(i) ) G(i) = tmp
  50        continue
c							this is an error
            if( abs(G(i)) .gt. one ) then
               ierr = 1
               go to 70
            endif
  60     continue
c							calculate upper bound
         bu = grcdf( B, G, N, ierr )
         if( ierr .ne. 0 ) bu = zero
c							for lower bound...
  70     do 90 i = 4, N
            G(i) = R(i,1)/G(1)
            do 80 j = 2, i-1
               tmp = R(i,j)/G(j)
               if( tmp .lt. G(i) ) G(i) = tmp
  80        continue
c							this is an error
            if( abs(G(i)) .gt. one ) then
               if( ierr .ne. 0 ) then
                  ierr = -2
               else
                  ierr = 2
               endif
               return
            endif
  90     continue
c							calculate lower bound
         bl = grcdf( B, G, N, ierr )
         if( ierr .ne. 0 ) bl = zero
      endif

      return
      end
