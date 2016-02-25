c  *********************************************************************
c  *                                                                   *
c  *                       subroutine dgmcdf                           *
c  *                                                                   *
c  *********************************************************************
c  Double Precision Version 1.2
c  Written by Gordon A. Fenton, TUNS, 1991
c  Latest Update: May 18, 1997
c
c  PURPOSE  computes bounds on the cumulative distribution function of an
c           n-variate jointly normal random vector.
c
c  This routine computes an upper and lower bound on the cumulative
c  distribution function of a jointly standard normal n-variate. The bounds
c  are obtained by approximating the correlation matrix Q by the special
c  multiplicative form, Q_{ij} ~ g_i*g_j + (1 - g_i*g_i)*d_{ij},
c  where d_{ij} is the Chronecker delta and g_i is some constant between 0 and
c  1. For such an approximation, the multi-dimensional integral reduces to
c  a one-dimensional integral. However, the vector {g} will only uniquely
c  represent a general correlation matrix of order 3 or less, higher order
c  matrices can only be approximated. Two approximations are made, one in
c  which the elements of {g} are chosen using a maximum ratio technique,
c  the other using a minimum ratio technique (see Madsen etal, "Methods of
c  Structural Safety", Prentice-Hall, pg 62) - these yield the upper and
c  lower bounds on the multi-dimensional integral respectively.
c  NOTES:
c  1) at this time, the elements of the correlation matrix Q must all be
c     positive to ensure real results (ie to ensure this routine works).
c  2) the single precision version of this routine employs a different
c     architecture. You may want to run both and compare the results...
c
c  Arguments to this routine are as follows;
c
c    B     real vector of length at least N containing the upper bounds of the
c          integration. B is re-ordered so that smallest (in absolute value)
c          elements appear first. Elements of B in excess of 16 (see parameter
c          BLRG) are eliminated and the marginal distribution is employed
c          (this reduces the dimension of the problem), in which case the
c          order of Q is also reduced by eliminating the appropriate rows
c          and columns and N is reduced. The lower bound of the integration
c          is -infinity (since this is a distribution function). (input/output)
c
c    N     integer giving the number of elements in the random vector for
c          which the cumulative distribution is desired. Note that N is
c          modified to denote the actual number of elements employed in
c          the integration. (input/output)
c
c    Q     real array of size at least N x N containing the elements of the
c          correlation matrix defining the n-variate process. Q is re-ordered
c          and reduced to remain consistent with the re-ordering/reduction of
c          B. (input/output)
c
c    iq    integer giving the leading dimension of the array Q as specified
c          in the calling routine. (input)
c
c    G     real vector of length at least N to contain the vector used
c          to approximate Q. (output)
c
c    T     temporary real vector of length at least N used for workspace.
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
c  1.1	changed name of dsplnt to dspint (Dec 5/96)
c  1.2	brought contents of dgblck.f directly into this routine; some compilers
c	don't like externally defined common blocks. (May 18/97)
c-----------------------------------------------------------------------------
      subroutine dgmcdf( B, N, Q, iq, G, T, bl, bu, ierr )
      implicit real*8 (a-h,o-z)
      parameter (NDIS = 257, BLRG = 16.d0)
      dimension B(*), Q(iq,*), G(*), T(*)
      dimension vx(NDIS), F(NDIS), TT(4*(NDIS-1))
c						real*8 conversion stuff
      data zero/0.d0/, half/0.5d0/, one/1.d0/
      data rt2pi/2.5066282746310005024d0/
c						[-10,10] discretization
      data (vx(i),i=1,35)/
     >-10.00d0, -9.50d0, -9.00d0, -8.50d0, -8.00d0, -7.50d0, -7.00d0,
     > -6.50d0, -6.00d0, -5.90d0, -5.80d0, -5.70d0, -5.60d0, -5.50d0,
     > -5.40d0, -5.30d0, -5.20d0, -5.10d0, -5.00d0, -4.90d0, -4.80d0,
     > -4.70d0, -4.60d0, -4.50d0, -4.40d0, -4.30d0, -4.20d0, -4.10d0,
     > -4.00d0, -3.90d0, -3.80d0, -3.70d0, -3.60d0, -3.50d0, -3.40d0/
      data (vx(i),i=36,70)/
     > -3.30d0, -3.20d0, -3.10d0, -3.00d0, -2.95d0, -2.90d0, -2.85d0,
     > -2.80d0, -2.75d0, -2.70d0, -2.65d0, -2.60d0, -2.55d0, -2.50d0,
     > -2.45d0, -2.40d0, -2.35d0, -2.30d0, -2.25d0, -2.20d0, -2.15d0,
     > -2.10d0, -2.05d0, -2.00d0, -1.95d0, -1.90d0, -1.85d0, -1.80d0,
     > -1.75d0, -1.70d0, -1.65d0, -1.60d0, -1.55d0, -1.50d0, -1.45d0/
      data (vx(i),i=71,105)/
     > -1.40d0, -1.35d0, -1.30d0, -1.25d0, -1.20d0, -1.15d0, -1.10d0,
     > -1.05d0, -1.00d0, -0.98d0, -0.96d0, -0.94d0, -0.92d0, -0.90d0,
     > -0.88d0, -0.86d0, -0.84d0, -0.82d0, -0.80d0, -0.78d0, -0.76d0,
     > -0.74d0, -0.72d0, -0.70d0, -0.68d0, -0.66d0, -0.64d0, -0.62d0,
     > -0.60d0, -0.58d0, -0.56d0, -0.54d0, -0.52d0, -0.50d0, -0.48d0/
      data (vx(i),i=106,140)/
     > -0.46d0, -0.44d0, -0.42d0, -0.40d0, -0.38d0, -0.36d0, -0.34d0,
     > -0.32d0, -0.30d0, -0.28d0, -0.26d0, -0.24d0, -0.22d0, -0.20d0,
     > -0.18d0, -0.16d0, -0.14d0, -0.12d0, -0.10d0, -0.08d0, -0.06d0,
     > -0.04d0, -0.02d0,  0.00d0,  0.02d0,  0.04d0,  0.06d0,  0.08d0,
     >  0.10d0,  0.12d0,  0.14d0,  0.16d0,  0.18d0,  0.20d0,  0.22d0/
      data (vx(i),i=141,175)/
     >  0.24d0,  0.26d0,  0.28d0,  0.30d0,  0.32d0,  0.34d0,  0.36d0,
     >  0.38d0,  0.40d0,  0.42d0,  0.44d0,  0.46d0,  0.48d0,  0.50d0,
     >  0.52d0,  0.54d0,  0.56d0,  0.58d0,  0.60d0,  0.62d0,  0.64d0,
     >  0.66d0,  0.68d0,  0.70d0,  0.72d0,  0.74d0,  0.76d0,  0.78d0,
     >  0.80d0,  0.82d0,  0.84d0,  0.86d0,  0.88d0,  0.90d0,  0.92d0/
      data (vx(i),i=176,210)/
     >  0.94d0,  0.96d0,  0.98d0,  1.00d0,  1.05d0,  1.10d0,  1.15d0,
     >  1.20d0,  1.25d0,  1.30d0,  1.35d0,  1.40d0,  1.45d0,  1.50d0,
     >  1.55d0,  1.60d0,  1.65d0,  1.70d0,  1.75d0,  1.80d0,  1.85d0,
     >  1.90d0,  1.95d0,  2.00d0,  2.05d0,  2.10d0,  2.15d0,  2.20d0,
     >  2.25d0,  2.30d0,  2.35d0,  2.40d0,  2.45d0,  2.50d0,  2.55d0/
      data (vx(i),i=211,245)/
     >  2.60d0,  2.65d0,  2.70d0,  2.75d0,  2.80d0,  2.85d0,  2.90d0,
     >  2.95d0,  3.00d0,  3.10d0,  3.20d0,  3.30d0,  3.40d0,  3.50d0,
     >  3.60d0,  3.70d0,  3.80d0,  3.90d0,  4.00d0,  4.10d0,  4.20d0,
     >  4.30d0,  4.40d0,  4.50d0,  4.60d0,  4.70d0,  4.80d0,  4.90d0,
     >  5.00d0,  5.10d0,  5.20d0,  5.30d0,  5.40d0,  5.50d0,  5.60d0/
      data (vx(i),i=246,257)/
     >  5.70d0,  5.80d0,  5.90d0,  6.00d0,  6.50d0,  7.00d0,  7.50d0,
     >  8.00d0,  8.50d0,  9.00d0,  9.50d0, 10.00d0/

c						conversion statement functions
      phi(x)  = dphi(x)
      sqrt(x) = dsqrt(x)
      abs(x)  = dabs(x)
      exp(x)  = dexp(x)
c---------------------------------------- start executable statements --------
      ierr = 0
c							re-order B and Q
      do 40 i = 1, N
         bmin = abs(B(i))
         k = i
         do 10 j = i+1, N
            if( abs(B(j)) .lt. bmin ) then
               bmin = abs(B(j))
               k    = j
            endif
  10     continue
c							swap B elements
         if( k .ne. i ) then
            tmp  = B(i)
            B(i) = B(k)
            B(k) = tmp
c							swap rows of Q
            do 20 j = 1, N
               tmp    = Q(i,j)
               Q(i,j) = Q(k,j)
               Q(k,j) = tmp
  20        continue
c							swap columns of Q
            do 30 j = 1, N
               tmp    = Q(j,i)
               Q(j,i) = Q(j,k)
               Q(j,k) = tmp
  30        continue
         endif
  40  continue
c							eliminate B large
      NN = N
      do 50 i = 1, N
         if( abs(B(i)) .gt. BLRG ) then
            if( B(i) .lt. zero ) then
               bl = zero
               bu = zero
               return
            endif
            if( NN .eq. N ) NN = i
         endif
  50  continue
      N = NN
c							if N < 4, unique...
      if( N .le. 0 ) then
         bl = one
         bu = one
         return
      elseif( N .eq. 1 ) then
         bl = phi( B(1) )
         bu = bl
         return
      elseif( N .eq. 2 ) then
         G(1) = sqrt(Q(2,1))
         rs   = sqrt(one - Q(2,1))
         T(1) = B(1)/rs
         T(2) = B(2)/rs
         G(1) = G(1)/rs
         G(2) = G(1)
         do 60 i = 1, NDIS
            rr   = exp(-half*vx(i)*vx(i))/rt2pi
            F(i) = rr*phi(T(1)-G(1)*vx(i))
     >               *phi(T(2)-G(2)*vx(i))
  60     continue
         bl = dspint( q, F, TT, NDIS, .true., zero, zero, ierr )
         bu = bl
         return
      elseif( N .eq. 3 ) then
         G(1) = sqrt(Q(2,1)*Q(3,1)/Q(3,2))
         G(2) = sqrt(Q(2,1)*Q(3,2)/Q(3,1))
         G(3) = sqrt(Q(3,1)*Q(3,2)/Q(2,1))
         do 70 i = 1, 3
            rr   = sqrt(one - G(i)*G(i))
            T(i) = B(i)/rr
            G(i) = G(i)/rr
  70     continue
         do 80 i = 1, NDIS
            rr   = exp(-half*vx(i)*vx(i))/rt2pi
            F(i) = rr*phi(T(1)-G(1)*vx(i))*phi(T(2)-G(2)*vx(i))
     >               *phi(T(3)-G(3)*vx(i))
  80     continue
         bl = dspint( q, F, TT, NDIS, .true., zero, zero, ierr )
         bu = bl
         return
      else
c							find elements of {G}
         G1 = sqrt(Q(2,1)*Q(3,1)/Q(3,2))
         G2 = sqrt(Q(2,1)*Q(3,2)/Q(3,1))
         G3 = sqrt(Q(3,1)*Q(3,2)/Q(2,1))
c							for upper bound...
         G(1) = G1
         G(2) = G2
         G(3) = G3
         do 100 i = 4, N
            G(i) = Q(i,1)/G(1)
            do 90 j = 2, i-1
               tmp = Q(i,j)/G(j)
               if( tmp .gt. G(i) ) G(i) = tmp
  90        continue
            if( abs(G(i)) .gt. one ) then
               ierr = 1
               go to 130
            endif
 100     continue
         do 110 i = 1, N
            rr   = sqrt(one - G(i)*G(i))
            T(i) = B(i)/rr
            G(i) = G(i)/rr
 110     continue
         do 120 i = 1, NDIS
            F(i) = exp(-half*vx(i)*vx(i))/rt2pi
            do 120 j = 1, N
               F(i) = F(i)*phi(T(j)-G(j)*vx(i))
 120     continue
         bu = dspint( q, F, TT, NDIS, .true., zero, zero, ierr )
         if( ierr .ne. 0 ) bu = zero
c							for lower bound...
 130     G(1) = G1
         G(2) = G2
         G(3) = G3
         do 150 i = 4, N
            G(i) = Q(i,1)/G(1)
            do 140 j = 2, i-1
               tmp = Q(i,j)/G(j)
               if( tmp .lt. G(i) ) G(i) = tmp
 140        continue
            if( abs(G(i)) .gt. one ) then
               if( ierr .ne. 0 ) then
                  ierr = -2
               else
                  ierr = 2
               endif
               return
            endif
 150     continue
         do 160 i = 1, N
            rr   = sqrt(one - G(i)*G(i))
            T(i) = B(i)/rr
            G(i) = G(i)/rr
 160     continue
         do 170 i = 1, NDIS
            F(i) = exp(-half*vx(i)*vx(i))/rt2pi
            do 170 j = 1, N
               F(i) = F(i)*phi(T(j)-G(j)*vx(i))
 170     continue
         bl = dspint( q, F, TT, NDIS, .true., zero, zero, ierr )
         if( ierr .ne. 0 ) bl = zero
      endif

      return
      end
