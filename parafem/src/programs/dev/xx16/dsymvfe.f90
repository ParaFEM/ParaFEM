PROGRAM DSYMVFE
!
! Short test program for prototyping potential Xeon Phi code
! Matrix-vector multiply where the matrix is stored as a vector
! 
! Author: Lee Margetts
!

 USE library
 IMPLICIT NONE

 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER::l,n,m,i,j,k,nod,ndim,nip,nst=6,nodof=3,ndof
 REAL(iwp),PARAMETER::pt5=0.5_iwp,pt1=0.1_iwp,zero=0.0_iwp
 REAL(iwp)::det,e,v
 REAL(iwp),ALLOCATABLE:: a(:,:),av(:),x(:),b(:),c(:),d(:),work(:),        &
   coord(:,:),bee(:,:),dee(:,:),der(:,:),deriv(:,:),fun(:),jac(:,:),      &
   km(:,:),points(:,:),weights(:)
 INTEGER,ALLOCATABLE::num(:)
!CHARACTER(LEN=15)::element='tetrahedron'
 CHARACTER(LEN=15)::element='hexahedron'
!CHARACTER(LEN=15)::cube='perfect'
 CHARACTER(LEN=15)::cube='arbitrary'
 LOGICAL::reduced=.false.
!LOGICAL::reduced=.true.

!-------------------------------------------------------------------------
! x. Initialization
!-------------------------------------------------------------------------

 nip=8
!nip=1 

 IF(element=='tetrahedron') THEN
   nod=4; m=0; m=10 ! m = nod+2 for tetrahedra
 ELSE IF(element=='hexahedron') THEN
   IF(cube=='perfect') THEN
     IF(nip==1) THEN
       nod=8; m=5  ! guess
     ELSE IF(nip==8) THEN
       IF(reduced) THEN
         nod=8; m=6
       ELSE
         nod=8; m=9
       END IF
     END IF
   ELSE
     nod=8; m=28 ! guess
   END IF
 END IF

 n=nod ! probably not needed
 
 ndim=3; ndof=nod*nodof
 e=100.0; v=0.3

 ALLOCATE(a(n,n),av(m),x(n),b(n),c(n),d(n),work(m),num(nod),             &
   coord(nod,ndim),points(nip,ndim),fun(nod),jac(ndim,ndim),             &
   der(ndim,nod),deriv(ndim,nod),bee(nst,ndof),km(ndof,ndof),                           &
   weights(nip),dee(nst,nst))

 a=0.0; av=0.0; x=0.0; b=0.0; c=0.0; d=0.0; work=0.0; num=0; coord=0.0

!-------------------------------------------------------------------------
! x. Create geometry
!-------------------------------------------------------------------------

 IF(element=='tetrahedron') THEN
   num = (/1,2,3,4/)
   coord(1,:) = (/0.1,0.11,0.12/)
   coord(2,:) = (/1.2,0.13,0.14/)
   coord(3,:) = (/1.3,1.4,0.15/)
   coord(4,:) = (/0.6,0.7,1.5/)
 ELSE IF(element=='hexahedron'.and.nod==8) THEN
   IF(cube=='perfect') THEN
     PRINT *, "Created geometry for perfect cube"
     coord(1,:) = (/0.0,0.0,0.0/); coord(2,:) = (/0.0,0.0,1.0/)
     coord(3,:) = (/1.0,0.0,1.0/); coord(4,:) = (/1.0,0.0,0.0/)
     coord(5,:) = (/0.0,1.0,0.0/); coord(6,:) = (/0.0,1.0,1.0/)
     coord(7,:) = (/1.0,1.0,1.0/); coord(8,:) = (/1.0,1.0,0.0/)
   ELSE  ! general cube
     PRINT *, "Created geometry for general cube"
     coord(1,:) = (/0.1,0.2,0.3/); coord(2,:) = (/0.4,0.5,1.1/)
     coord(3,:) = (/1.2,0.6,1.3/); coord(4,:) = (/1.4,0.5,0.4/)
     coord(5,:) = (/0.3,1.7,0.2/); coord(6,:) = (/0.1,1.6,1.5/)
     coord(7,:) = (/1.4,1.3,1.2/); coord(8,:) = (/1.2,1.3,0.1/)
   END IF
 END IF

 DO i=1,nod
   PRINT *, "Node ", i, "Coord =", coord(i,:)
 END DO

!-------------------------------------------------------------------------
! x. Create stiffness matrix 'a'
!-------------------------------------------------------------------------

 CALL sample(element,points,weights)
 CALL deemat(dee,e,v)
 km=zero

 int_pts_1: DO i=1,nip
   CALL shape_fun(fun,points,i)
   CALL shape_der(der,points,i)
   jac=MATMUL(der,coord)
   det=determinant(jac)
   CALL invert(jac)
   deriv=MATMUL(jac,der)
   CALL beemat(bee,deriv)
   km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
 END DO int_pts_1

 a=km
 
!-------------------------------------------------------------------------
! x. Create LHS vector 'x'
!-------------------------------------------------------------------------

 DO i=1,n
   x(i) = i
 END DO

!-------------------------------------------------------------------------
! x. b=ax using MATMUL
!-------------------------------------------------------------------------

 b= MATMUL(a,x)

!-------------------------------------------------------------------------
! x. Compress 'a' to vector form 'av'
!-------------------------------------------------------------------------

 IF(element=='tetrahedron') THEN
   av(1) = a(1,1);  av(2) = a(1,2);  av(3) = a(1,3)
   av(4) = a(1,4);  av(5) = a(2,2);  av(6) = a(2,3)
   av(7) = a(2,4);  av(8) = a(3,3);  av(9) = a(3,4); av(10) = a(4,4)
 ELSE IF(element=='hexahedron'.and.nod==8) THEN
   IF(cube=='perfect') THEN
     IF(nip==1) THEN
       av(1)  = a(1,1);  av(2)  = a(1,2);  av(3)  = a(1,4);  av(4) = a(1,6)
       av(5)  = a(2,8)
     ELSE IF(nip==8) THEN
       IF(reduced) THEN
         av(1)  = a(1,1);  av(2)  = a(1,2);  av(3)  = a(1,4);  av(4) = a(1,6)
         av(5)  = a(1,7);  av(6)  = a(3,6)
       ELSE
         av(1)  = a(1,1);  av(2)  = a(1,2);  av(3)  = a(1,4);  av(4) = a(1,5)
         av(5)  = a(1,6);  av(6)  = a(1,7);  av(7)  = a(1,8);  av(8) = a(2,8)
         av(9)  = a(3,6)
       END IF
     END IF
   ELSE
     av(1)  = a(1,1);  av(2)  = a(1,2);  av(3)  = a(1,3); av(4)  = a(1,4)
     av(5)  = a(1,5);  av(6)  = a(1,6);  av(7)  = a(1,7); av(8)  = a(1,8)
     av(9)  = a(2,7);  av(10) = a(2,8);  av(11) = a(3,3); av(12) = a(3,4)
     av(13) = a(3,6);  av(14) = a(3,7);  av(15) = a(3,8); av(16) = a(4,4)
     av(17) = a(4,5);  av(18) = a(4,6);  av(19) = a(4,7); av(20) = a(4,8)   
     av(21) = a(5,7);  av(22) = a(5,8);  av(23) = a(6,6); av(24) = a(6,7)
     av(25) = a(6,8);  av(26) = a(7,7);  av(27) = a(7,8); av(28) = a(8,8)
  END IF
 END IF

!-------------------------------------------------------------------------
! x. Another way
!-------------------------------------------------------------------------

 IF(element=='tetrahedron') THEN
   d(1) = av(1)*x(1) + av(2)*x(2) + av(3)*x(3) + av(4)*x(4)
   d(2) = av(2)*x(1) + av(5)*x(2) + av(6)*x(3) + av(7)*x(4)
   d(3) = av(3)*x(1) + av(6)*x(2) + av(8)*x(3) + av(9)*x(4)
   d(4) = av(4)*x(1) + av(7)*x(2) + av(9)*x(3) + av(10)*x(4)
 ELSE IF(element=='hexahedron'.and.nod==8) THEN
   IF(cube=='perfect') THEN
     IF(nip==1) THEN
       d(1) = av(1)*x(1)  + av(2)*x(2)  + av(2)*x(3)  + av(3)*x(4)  +    &
              av(2)*x(5)  + av(4)*x(6)  - av(3)*x(7)  - av(4)*x(8)
       d(2) = av(2)*x(1)  + av(1)*x(2)  + av(2)*x(3)  + av(2)*x(4)  +    &
              av(3)*x(5)  + av(4)*x(6)  + av(4)*x(7)  + av(5)*x(8)
       d(3) = av(2)*x(1)  + av(2)*x(2)  + av(1)*x(3)  - av(4)*x(4)  -    &
              av(4)*x(5)  - av(5)*x(6)  - av(2)*x(7)  - av(4)*x(8)
       d(4) = av(3)*x(1)  + av(2)*x(2)  - av(4)*x(3)  + av(1)*x(4)  +    &
              av(2)*x(5)  - av(2)*x(6)  - av(5)*x(7)  - av(4)*x(8)
       d(5) = av(2)*x(1)  + av(3)*x(2)  - av(4)*x(3)  + av(2)*x(4)  +    &
              av(1)*x(5)  - av(2)*x(6)  + av(4)*x(7)  + av(3)*x(8)
       d(6) = av(4)*x(1)  + av(4)*x(2)  - av(5)*x(3)  - av(2)*x(4)  -    &
              av(2)*x(5)  + av(1)*x(6)  - av(4)*x(7)  - av(2)*x(8)
       d(7) =-av(3)*x(1)  + av(4)*x(2)  - av(2)*x(3)  - av(5)*x(4)  +    &
              av(4)*x(5)  - av(4)*x(6)  + av(1)*x(7)  - av(2)*x(8)
       d(8) =-av(4)*x(1)  + av(5)*x(2)  - av(4)*x(3)  - av(4)*x(4)  +    &
              av(3)*x(5)  - av(2)*x(6)  - av(2)*x(7)  + av(1)*x(8)
     ELSE IF(nip==8) THEN
       IF(.not.reduced) THEN
       d(1) = av(1)*x(1)  + av(2)*x(2)  + av(2)*x(3)  + av(3)*x(4)  +    &
              av(4)*x(5)  + av(5)*x(6)  + av(6)*x(7)  + av(7)*x(8)
       d(2) = av(2)*x(1)  + av(1)*x(2)  + av(2)*x(3)  + av(4)*x(4)  +    &
              av(3)*x(5)  + av(5)*x(6)  - av(7)*x(7)  + av(8)*x(8)
       d(3) = av(2)*x(1)  + av(2)*x(2)  + av(1)*x(3)  - av(5)*x(4)  -    &
              av(5)*x(5)  + av(9)*x(6)  - av(2)*x(7)  + av(7)*x(8)
       d(4) = av(3)*x(1)  + av(4)*x(2)  - av(5)*x(3)  + av(1)*x(4)  +    &
              av(2)*x(5)  - av(2)*x(6)  + av(9)*x(7)  - av(5)*x(8)
       d(5) = av(4)*x(1)  + av(3)*x(2)  - av(5)*x(3)  + av(2)*x(4)  +    &
              av(1)*x(5)  - av(2)*x(6)  + av(5)*x(7)  + av(3)*x(8)
       d(6) = av(5)*x(1)  + av(5)*x(2)  + av(9)*x(3)  - av(2)*x(4)  -    &
              av(2)*x(5)  + av(1)*x(6)  - av(5)*x(7)  - av(4)*x(8)
       d(7) = av(6)*x(1)  - av(7)*x(2)  - av(2)*x(3)  + av(9)*x(4)  +    &
              av(5)*x(5)  - av(5)*x(6)  + av(1)*x(7)  - av(2)*x(8)
       d(8) = av(7)*x(1)  + av(8)*x(2)  + av(7)*x(3)  - av(5)*x(4)  +    &
              av(3)*x(5)  - av(4)*x(6)  - av(2)*x(7)  + av(1)*x(8)
       ELSE
       PRINT *, "Using reduced av vector"
       d(1) = av(1)*x(1)  + av(2)*x(2)  + av(2)*x(3)  + av(3)*x(4)  +    &
              av(2)*x(5)*pt5  + av(4)*x(6) + av(5)*x(7) - av(4)*x(8)*pt5
       d(2) = av(2)*x(1)  + av(1)*x(2)  + av(2)*x(3)  + av(2)*x(4)*pt5 + &
              av(3)*x(5)  + av(4)*x(6)  + av(4)*x(7)*pt5 - av(3)*x(8)*pt1
       d(3) = av(2)*x(1)  + av(2)*x(2)  + av(1)*x(3)  - av(4)*x(4)  -    &
              av(4)*x(5)  + av(6)*x(6)  - av(2)*x(7)  - av(4)*x(8)*pt5
       d(4) = av(3)*x(1)  + av(2)*x(2)*pt5  - av(4)*x(3)  + av(1)*x(4) + &
              av(2)*x(5)  - av(2)*x(6)  + av(6)*x(7)  - av(4)*x(8)
       d(5) = av(2)*x(1)*pt5 + av(3)*x(2)  - av(4)*x(3)  + av(2)*x(4)  +    &
              av(1)*x(5)  - av(2)*x(6)  + av(4)*x(7)  + av(3)*x(8)
       d(6) = av(4)*x(1)  + av(4)*x(2)  + av(6)*x(3)  - av(2)*x(4)  -    &
              av(2)*x(5)  + av(1)*x(6)  - av(4)*x(7)  - av(2)*x(8)*pt5
       d(7) = av(5)*x(1)  + av(4)*x(2)*pt5 - av(2)*x(3)  + av(6)*x(4)  +    &
              av(4)*x(5)  - av(4)*x(6)  + av(1)*x(7)  - av(2)*x(8)
       d(8) = av(4)*x(1)*pt5  - av(3)*x(2)*pt1 - av(4)*x(3)*pt5 - av(4)*x(4)  +    &
              av(3)*x(5)  - av(2)*x(6)*pt5  - av(2)*x(7)  + av(1)*x(8)
       END IF
     ELSE
       PRINT *, "nip = ", nip, " not supported"
     END IF
   ELSE
     d(1) = av(1)*x(1)  + av(2)*x(2)  + av(3)*x(3)  + av(4)*x(4)  +      &
            av(5)*x(5)  + av(6)*x(6)  + av(7)*x(7)  + av(8)*x(8)
     d(2) = av(2)*x(1)  + av(1)*x(2)  + av(3)*x(3)  + av(5)*x(4)  +      &
            av(4)*x(5)  + av(6)*x(6)  + av(9)*x(7)  + av(10)*x(8)
     d(3) = av(3)*x(1)  + av(3)*x(2)  + av(11)*x(3) + av(12)*x(4) +      &
            av(12)*x(5) + av(13)*x(6) + av(14)*x(7) + av(15)*x(8)
     d(4) = av(4)*x(1)  + av(5)*x(2)  + av(12)*x(3) + av(16)*x(4) +      &
            av(17)*x(5) + av(18)*x(6) + av(19)*x(7) + av(20)*x(8)
     d(5) = av(5)*x(1)  + av(4)*x(2)  + av(12)*x(3) + av(17)*x(4) +      &
            av(16)*x(5) + av(18)*x(6) + av(21)*x(7) + av(22)*x(8)
     d(6) = av(6)*x(1)  + av(6)*x(2)  + av(13)*x(3) + av(18)*x(4) +      &
            av(18)*x(5) + av(23)*x(6) + av(24)*x(7) + av(25)*x(8)
     d(7) = av(7)*x(1)  + av(9)*x(2)  + av(14)*x(3) + av(19)*x(4) +      &
            av(21)*x(5) + av(24)*x(6) + av(26)*x(7) + av(27)*x(8)
     d(8) = av(8)*x(1)  + av(10)*x(2) + av(15)*x(3) + av(20)*x(4) +      &
            av(22)*x(5) + av(25)*x(6) + av(27)*x(7) + av(28)*x(8)
  END IF
 END IF
 
 PRINT *, "d =", d

!-------------------------------------------------------------------------
! x. Write out matrix a, vector av, vector x 
!-------------------------------------------------------------------------

PRINT *, ""; PRINT *, "Element = ", element, " with ", nod, " nodes" 
PRINT *, ""; PRINT *, "*** Values stored ***"
PRINT *, ""; PRINT *, "Matrix 'a'"
DO i=1,n
  WRITE(*,'(8F10.6)') a(i,:)
END DO

PRINT *, ""; PRINT *, "Matrix 'a' condensed to vector 'av'"
PRINT *, av

PRINT *, ""; PRINT *, "Vector 'x'"
PRINT *, x

PRINT *, ""; PRINT *, "*** Compare output ***"
PRINT *, ""
PRINT *, "Using full matrix b= MATMUL(a,x) "
PRINT *, ""; PRINT *, b

PRINT *, ""; PRINT *, "Using compressed form 'av'"; PRINT *, ""; PRINT *, d
        ! note d(i) is used in the code
        ! b(i) is used in the text output above 

PRINT *, ""; PRINT *, "*** Storage ***"; PRINT *, ""
PRINT *, "Size of matrix 'a'                ", n*n
PRINT *, "Size when stored as vector 'av'   ", m
PRINT *, "% Storage saved using vector form ", (1-(REAL(m)/(REAL(n)*REAL(n))))*100 

END PROGRAM DSYMVFE