PROGRAM DSYMMV
!
! Short test program for prototyping potential Xeon Phi code
! Matrix-vector multiply where the matrix is stored as a vector
! 
! Author: Lee Margetts
!

IMPLICIT NONE

INTEGER::l,n,m,i,j,k,nod,ndim
REAL,ALLOCATABLE:: a(:,:),av(:),x(:),b(:),c(:),d(:),work(:),coord(:,:)
INTEGER,ALLOCATABLE::num(:)

PRINT *, ""; PRINT *, "Enter a value for n"; READ *, n
PRINT *, ""; PRINT *, "n = ", n

!-------------------------------------------------------------------------
! 1. Initialization
!-------------------------------------------------------------------------

nod=4; ndim=3 

m=0; m=n+n-1

ALLOCATE(a(n,n),av(m),x(n),b(n),c(n),d(n),work(m),num(nod),coord(nod,ndim))

a=0.0; av=0.0; x=0.0; b=0.0; c=0.0; d=0.0; work=0.0; num=0; coord=0.0

!-------------------------------------------------------------------------
! 2. Create geometry
!-------------------------------------------------------------------------

num = (/1,2,3,4/)

PRINT *, ""; PRINT *, "num = ", num

!-------------------------------------------------------------------------
! 3. Create symmetrix matrix
!-------------------------------------------------------------------------

DO i=1,n
  DO j=1,n
    a(i,j) = j+i-1
  END DO
END DO

DO i=1,n
  x(i) = i
END DO

b= MATMUL(a,x)

av(1:n) = a(1,:)
l=n-1
DO i=2,n
  j=l+i
  av(j)=a(i,n)
END DO

DO i=1,n
  c(i) = DOT_PRODUCT(av(i:n),x)
END DO

!-------------------------------------------------------------------------
! 4. Another way
!-------------------------------------------------------------------------

DO i=1,n
  DO j=1,n
    k=i+j-1
    d(i) = d(i) + av(k)*x(j)
  END DO
END DO

PRINT *, "d =", d

!-------------------------------------------------------------------------
! 5. Write out matrix a, vector av, vector x 
!-------------------------------------------------------------------------

PRINT *, ""; PRINT *, "*** Values stored ***"
PRINT *, ""; PRINT *, "Matrix 'a'"
DO i=1,n
  PRINT *, a(i,:)
END DO

PRINT *, ""; PRINT *, "Matrix 'a' condensed to vector 'av'"
PRINT *, av

PRINT *, ""; PRINT *, "Vector 'x'"
PRINT *, x

PRINT *, ""; PRINT *, "*** Compare output ***"
PRINT *, ""
PRINT *, "b= MATMUL(a,x) "
PRINT *, ""; PRINT *, "ans= ", b

PRINT *, ""
PRINT *, "DO i=1,n"
PRINT *, "  b(i) = DOT_PRODUCT(av(i:n),x)"
PRINT *, "END DO"
PRINT *, ""; PRINT *, "ans= ", c 
           ! note c(i) is used in the code
           ! b(i) is used in the text output above 

PRINT *, ""
PRINT *, "DO i=1,n"
PRINT *, "  DO j=1,n"
PRINT *, "    k=i+j-1"
PRINT *, "    b(i) = b(i) + av(k)*x(j)"
PRINT *, "  END DO"
PRINT *, "END DO"
PRINT *, ""; PRINT *, "ans= ", d
        ! note d(i) is used in the code
        ! b(i) is used in the text output above 

PRINT *, ""; PRINT *, "*** Storage ***"; PRINT *, ""
PRINT *, "Size of matrix 'a'                ", n*n
PRINT *, "Size when stored as vector 'av'   ", m
PRINT *, "% Storage saved using vector form ", (1-(REAL(m)/(REAL(n)*REAL(n))))*100 

END PROGRAM DSYMMV