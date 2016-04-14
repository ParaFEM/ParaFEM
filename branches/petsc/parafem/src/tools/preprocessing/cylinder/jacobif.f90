!$Id: jacobif.f90 31 2016-01-24 21:47:27Z mexas $
! Copyright Anton Shterenlikht, The University of Bristol, UK

! module with Jacobi functions, etc.
module jacobif
implicit none
private
public :: jefr, sq2ci

contains

!***************************************************************

subroutine jefr( x , m , tol, iter, sn, cn, dn )

! Jacobian elliptic functions of a real argument
! Using Arithmetic-Geometric Mean
! See 22.20(ii) Arithmetic-Geometric Mean:
! http://dlmf.nist.gov/22.20
!
! x - argument
! m - parameter, m=k**2
! tol - required tolerance. The convergence is very good,
!       so epsilon() is adequate
! iter - how many iterations was required
! sn, cn, dn - 3 Jacobian elliptic functions

double precision, intent(in) :: x , m, tol
integer, intent(out) :: iter
double precision, intent(out) :: sn, cn, dn

double precision, parameter :: zero = 0.0d0, one = 1.0d0, &
 half = 0.5d0
integer, parameter :: lim = 20

double precision :: a(0:lim) , b(0:lim) , c(0:lim), f(0:lim)
integer :: i

! sanity check. Must have 0 .ge. m .ge. 1
if ( m .lt. zero .or. m .gt. one ) then
  stop "ERROR: jefr: m is out of limits [0..1]"
end if

! starting point
i = 0
a(0) = one 
b(0) = sqrt( one - m )
! c(0) = zero ! could be anything, not used in calcucations.

! iterations
do
  i = i + 1
  a(i) = half * ( a(i-1) + b(i-1) )
  b(i) =    sqrt( a(i-1) * b(i-1) )
  c(i) = half * ( a(i-1) - b(i-1) )
!  write (*,*) "i", i, a(i) , b(i) , c(i)
  if ( abs( c(i) ) .lt. tol ) exit
end do

! return the number of iterations taken
iter = i

! calculate all f
f(iter) = 2**iter * a(iter) * x

do i = iter , 1 , -1
  f(i-1) = half * ( f(i) + asin( c(i) / a(i) * sin( f(i) ) ) )
!  write (*,*) "i", i-1 , f(i-1)
end do

sn = sin( f(0) )
cn = cos( f(0) )

! Both dn forms are mathematically identical. The second
! form avoids division by zero, but might be less numerically
! accurate.
!dn = cn / cos( f(1) - f(0) )
dn = sqrt( one - m * sn**2 )

end subroutine jefr

!***************************************************************

subroutine sq2ci( xsq, ysq, xci, yci )

! adapted from
! http://arxiv.org/ftp/arxiv/papers/1011/1011.3189.pdf

double precision, intent(in) :: xsq , ysq
double precision, intent(out) :: xci , yci

double precision, parameter :: zero = 0.0d0 , half = 0.5d0 ,   &
  one = 1.0d0 , sqrt2 = sqrt( 2.0d0 ), sqrt22 = half * sqrt2 , &
  m = half , ke = 1.854074677301372d0

double precision :: eps, tol, s, c, d, s1, c1, d1, delta, &
  xpr, ypr
integer :: iter

! machine epsilon
eps = epsilon( 1.0d0 )
tol = eps

xpr = half * ke * ( xsq - ysq ) + ke
ypr = half * ke * ( xsq + ysq )

if ( abs(ypr) .lt. tol ) then
  call jefr( xpr , m , tol , iter , s , c , d )
  xci = c
  yci = zero
else
  call jefr( xpr , m     , tol , iter , s  , c  , d  )
  call jefr( ypr , one-m , tol , iter , s1 , c1 , d1 )
  delta = c1**2 + m * s**2 * s1**2
    xci =   c * c1 / delta
    yci = - s * d * s1 * d1 / delta
end if

end subroutine sq2ci

end module jacobif
