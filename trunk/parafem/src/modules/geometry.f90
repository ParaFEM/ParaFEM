MODULE GEOMETRY

   !/****h* /geometry 
   !*  NAME
   !*    MODULE: geometry
   !*  SYNOPSIS
   !*    Usage: USE pcg
   !*  FUNCTION
   !*    Contains subroutines to create simple meshes for programs in 
   !*    Chapter 12 of "Programming the Finite Element Method"
   !*
   !*    Subroutine             Purpose
   !*
   !*    GEOMETRY_8BXZ          Forms the coordinates and nodal vector for
   !*                           boxes of 8-node brick elements counting x-z 
   !*                           planes in y-direction
   !*
   !*    GEOMETRY_20BXZ         Nodal vector and nodal coordinates for boxes
   !*                           of 20-node bricks counting x-z planes in the 
   !*                           y-direction
   !*
   !*    GEOMETRY_20BXZV        Nodal vector and nodal coordinates for boxes 
   !*                           of 20-node bricks counting x-z planes in the 
   !*                           y-direction. Variable mesh version.
   !*
   !*    CUBE_BC20              Boundary conditions for simple 20-node element
   !*                           cuboid  
   !*
   !*    CUBE_BC8               Boundary conditions for simple 8-node element
   !*                           cuboid  
   !*
   !*    BOX_BC8                Boundary conditions for simple 8-node element
   !*                           box
   !*
   !*    NS_LOADING             Simple lid displacement of a square 
   !*                           Navier-Stokes patch 
   !*
   !*    NS_CUBE_BC20           Boundary conditions for simple Navier-Stokes 
   !*                           20-node element cuboid
   !*  AUTHOR
   !*    L. Margetts
   !*    I.M. Smith
   !*    D.V. Griffiths
   !*  COPYRIGHT
   !*    2004-2010 University of Manchester
   !******
   !*  Place remarks that should not be included in the documentation here
   !*
   !*/

   USE precision

   IMPLICIT NONE

   CONTAINS

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE GEOMETRY_8BXZ(iel,nxe,nze,aa,bb,cc,coord,num)

  !/****f* geometry/geometry_8bxz
  !*  NAME
  !*    SUBROUTINE: geometry_8bxz
  !*  SYNOPSIS
  !*    Usage:      CALL geometry_8bxz(iel,nxe,nze,aa,bb,cc,coord,num)
  !*  FUNCTION
  !*    This subroutine forms the coordinates and nodal vector
  !*    for boxes of 8-node brick elements counting x-z planes in y-direction
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:    
  !*
  !*    iel                 : Integer
  !*                        : Element ID
  !*
  !*    nxe                 : Integer
  !*                        : Number of elements in the x-direction
  !*
  !*    nze                 : Integer
  !*                        : Number of elements in the z-direction
  !*
  !*    aa                  : Real
  !*                        : Element width
  !*
  !*    bb                  : Real
  !*                        : Element depth
  !*
  !*    cc                  : Real
  !*                        : Element height
  !*
  !*  OUTPUTS
  !*    The following arguments have the INTENT(OUT) attribute:    
  !*
  !*    num(:)              : Integer
  !*                        : Vector holding the nodal connectivities for a 
  !*                          single element
  !*
  !*    coord(:,:)          : Real
  !*                        : Array holding the nodal co-ordinates for a 
  !*                          single element
  !*
  !*  AUTHOR
  !*    I.M. Smith and D.V. Griffiths
  !******
  !*/    
 
   IMPLICIT NONE
 
   INTEGER,  INTENT(IN) :: iel,nxe,nze;integer,intent(out)::num(:)
   REAL(iwp),INTENT(IN) :: aa,bb,cc
   REAL(iwp),INTENT(OUT):: coord(:,:)
   INTEGER              :: ip,iq,is,iplane
   
   iq         = (iel-1)/(nxe*nze)+1 
   iplane     = iel -(iq-1)*nxe*nze
   is         = (iplane-1)/nxe+1
   ip         = iplane-(is-1)*nxe   
  
   num(1)     = (iq-1)*(nxe+1)*(nze+1)+is*(nxe+1)+ip 
   num(2)     = num(1)-nxe-1
   num(3)     = num(2)+1  
   num(4)     = num(1)+1 
   num(5)     = num(1)+(nxe+1)*(nze+1)
   num(6)     = num(5)-nxe-1 
   num(7)     = num(6)+1    
   num(8)     = num(5)+1
   
   coord(1,1) = (ip-1)*aa 
   coord(2,1) = (ip-1)*aa 
   coord(5,1) = (ip-1)*aa
   coord(6,1) = (ip-1)*aa 
   coord(3,1) = ip*aa 
   coord(4,1) = ip*aa
   coord(7,1) = ip*aa   
   coord(8,1) = ip*aa
   coord(1,2) = (iq-1)*bb  
   coord(2,2) = (iq-1)*bb 
   coord(3,2) = (iq-1)*bb
   coord(4,2) = (iq-1)*bb  
   coord(5,2) = iq*bb  
   coord(6,2) = iq*bb
   coord(7,2) = iq*bb 
   coord(8,2) = iq*bb 
   coord(1,3) = -is*cc
   coord(4,3) = -is*cc  
   coord(5,3) = -is*cc 
   coord(8,3) = -is*cc
  
   coord(2,3) = -(is-1)*cc 
   coord(3,3) = -(is-1)*cc  
   coord(6,3) = -(is-1)*cc
   coord(7,3) = -(is-1)*cc 
   
   RETURN
 
 END SUBROUTINE GEOMETRY_8BXZ  

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
 
  SUBROUTINE GEOMETRY_20BXZ(iel,nxe,nze,aa,bb,cc,coord,num)

    !/****f* geometry/geometry_20bxz
    !*  NAME
    !*    SUBROUTINE: geometry_20bxz
    !*  SYNOPSIS
    !*    Usage:      CALL geometry_20bxz(iel,nxe,nze,aa,bb,cc,coord,num)
    !*  FUNCTION
    !*    Nodal vector and nodal coordinates for boxes of 20-node bricks 
    !*    counting x-z planes in the y-direction.
    !*    Source: Smith and Griffiths Edition 3.
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:    
    !*
    !*    iel                 : Integer
    !*                        : Element ID
    !*
    !*    nxe                 : Integer
    !*                        : Number of elements in the x-direction
    !*
    !*    nze                 : Integer
    !*                        : Number of elements in the z-direction
    !*
    !*    width(:)            : Real
    !*                        : Vector of mesh spacings across the width
    !*
    !*    depth(:)            : Real
    !*                        : Vector of mesh spacings across the depth
    !*
    !*    height(:)           : Real
    !*                        : Vector of mesh spacings across the height
    !*
    !*  OUTPUTS
    !*    The following arguments have the INTENT(OUT) attribute:    
    !*
    !*    num(:)              : Integer
    !*                        : Vector holding the nodal connectivities for a 
    !*                          single element
    !*
    !*    coord(:,:)          : Real
    !*                        : Array holding the nodal co-ordinates for a 
    !*                          single element
    !*
    !*  AUTHOR
    !*    I.M. Smith and D.V. Griffiths
    !******
    !*/    

    IMPLICIT NONE
    
    INTEGER,  INTENT(IN)  :: iel,nxe,nze
    REAL(iwp),INTENT(IN)  :: aa,bb,cc
    REAL(iwp),INTENT(OUT) :: coord(:,:)
    INTEGER,  INTENT(OUT) :: num(:)
    INTEGER               :: fac1,fac2,ip,iq,is,iplane
    
    iq     = (iel-1)/(nxe*nze)+1
    iplane = iel-(iq-1)*nxe*nze
    is     = (iplane-1)/nxe+1 
    ip     = iplane-(is-1)*nxe
    fac1   = ((2*nxe+1)*(nze+1)+(2*nze+1)*(nxe+1))*(iq-1)
    fac2   = ((2*nxe+1)*(nze+1)+(2*nze+1)*(nxe+1))*iq

    num(1)  = fac1+(3*nxe+2)*is+2*ip-1
    num(2)  = fac1+(3*nxe+2)*is-nxe+ip-1
    num(3)  = num(1)-3*nxe-2
    num(4)  = num(3)+1
    num(5)  = num(4)+1
    num(6)  = num(2)+1
    num(7)  = num(1)+2
    num(8)  = num(1)+1
    num(9)  = fac2-(nxe+1)*(nze+1)+(nxe+1)*is+ip
    num(10) = num(9)-nxe-1
    num(11) = num(10)+1
    num(12) = num(9)+1
    num(13) = fac2+(3*nxe+2)*is+2*ip-1
    num(14) = fac2+(3*nxe+2)*is-nxe+ip-1
    num(15) = num(13)-3*nxe-2
    num(16) = num(15)+1
    num(17) = num(16)+1
    num(18) = num(14)+1
    num(19) = num(13)+2
    num(20) = num(13)+1 
    
    coord(1:3,1)    = (ip-1)*aa
    coord(9:10,1)   = (ip-1)*aa
    coord(13:15,1)  = (ip-1)*aa
    coord(5:7,1)    = ip*aa
    coord(11:12,1)  = ip*aa
    coord(17:19,1)  = ip*aa
    coord(4,1)      = .5*(coord(3,1)+coord(5,1))
    coord(8,1)      = .5*(coord(1,1)+coord(7,1))
    coord(16,1)     = .5*(coord(15,1)+coord(17,1))
    coord(20,1)     = .5*(coord(13,1)+coord(19,1))
    coord(1:8,2)    = (iq-1)*bb
    coord(13:20,2)  = iq*bb
    coord(9,2)      = .5*(coord(1,2)+coord(13,2))
    coord(10,2)     = .5*(coord(3,2)+coord(15,2))
    coord(11,2)     = .5*(coord(5,2)+coord(17,2))
    coord(12,2)     = .5*(coord(7,2)+coord(19,2))
    coord(1,3)      = -is*cc
    coord(7:9,3)    = -is*cc
    coord(12:13,3)  = -is*cc
    coord(19:20,3)  = -is*cc
    coord(3:5,3)    = -(is-1)*cc
    coord(10:11,3)  = -(is-1)*cc
    coord(15:17,3)  = -(is-1)*cc
    coord(2,3)      = .5*(coord(1,3)+coord(3,3))
    coord(6,3)      = .5*(coord(5,3)+coord(7,3))
    coord(14,3)     = .5*(coord(13,3)+coord(15,3))
    coord(18,3)     = .5*(coord(17,3)+coord(19,3))
   
    RETURN
    
  END SUBROUTINE GEOMETRY_20BXZ

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE GEOMETRY_20BXZV(iel,nxe,nze,width,depth,height,coord,num)

    !/****f* geometry/geometry_20bxzv
    !*  NAME
    !*    SUBROUTINE: geometry_20bxzv
    !*  SYNOPSIS
    !*    Usage:      CALL geometry_20bxzv(iel,nxe,nze,width,depth,height,   &
    !*                                    coord,num)
    !*  FUNCTION
    !*    Nodal vector and nodal coordinates for boxes of 20-node bricks 
    !*    counting x-z planes in the y-direction. Variable mesh version.
    !*    Written in the style of Smith and Griffiths Edition 3.
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:    
    !*
    !*    iel                 : Integer
    !*                        : Element ID
    !*
    !*    nxe                 : Integer
    !*                        : Number of elements in the x-direction
    !*
    !*    nze                 : Integer
    !*                        : Number of elements in the z-direction
    !*
    !*    width(:)            : Real
    !*                        : Vector of mesh spacings across the width
    !*
    !*    depth(:)            : Real
    !*                        : Vector of mesh spacings across the depth
    !*
    !*    height(:)           : Real
    !*                        : Vector of mesh spacings across the height
    !*
    !*  OUTPUTS
    !*    The following arguments have the INTENT(OUT) attribute:    
    !*
    !*    num(:)              : Integer
    !*                        : Vector holding the nodal connectivities for a 
    !*                          single element
    !*
    !*    coord(:,:)          : Real
    !*                        : Array holding the nodal co-ordinates for a 
    !*                          single element
    !*
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    19.04.2002
    !*  COPYRIGHT
    !*    (c) University of Manchester
    !******
    !*/    

    IMPLICIT NONE

    INTEGER,  INTENT(IN)  :: iel,nxe,nze
    REAL(iwp),INTENT(IN)  :: width(:),depth(:),height(:)
    REAL(iwp),INTENT(OUT) :: coord(:,:)
    INTEGER,  INTENT(OUT) :: num(:)
    INTEGER               :: fac1,fac2,ip,iq,is,iplane
    
    iq     = (iel-1)/(nxe*nze)+1
    iplane = iel-(iq-1)*nxe*nze
    is     = (iplane-1)/nxe+1 
    ip     = iplane-(is-1)*nxe
    fac1   = ((2*nxe+1)*(nze+1)+(2*nze+1)*(nxe+1))*(iq-1)
    fac2   = ((2*nxe+1)*(nze+1)+(2*nze+1)*(nxe+1))*iq

    !-----------------------------------------------------------------------
    ! 1. Find num(:)
    !-----------------------------------------------------------------------
    
    num(1)  = fac1+(3*nxe+2)*is+2*ip-1
    num(2)  = fac1+(3*nxe+2)*is-nxe+ip-1
    num(3)  = num(1)-3*nxe-2
    num(4)  = num(3)+1
    num(5)  = num(4)+1
    num(6)  = num(2)+1
    num(7)  = num(1)+2
    num(8)  = num(1)+1
    num(9)  = fac2-(nxe+1)*(nze+1)+(nxe+1)*is+ip
    num(10) = num(9)-nxe-1
    num(11) = num(10)+1
    num(12) = num(9)+1
    num(13) = fac2+(3*nxe+2)*is+2*ip-1
    num(14) = fac2+(3*nxe+2)*is-nxe+ip-1
    num(15) = num(13)-3*nxe-2
    num(16) = num(15)+1
    num(17) = num(16)+1
    num(18) = num(14)+1
    num(19) = num(13)+2
    num(20) = num(13)+1 

    !-----------------------------------------------------------------------
    ! 2. Compute the nodal co-ordinates and populate coord(:,:)
    !-----------------------------------------------------------------------
    
    coord(1:3,1)   = width(ip)
    coord(9:10,1)  = width(ip)
    coord(13:15,1) = width(ip)
    coord(5:7,1)   = width(ip+1)
    coord(11:12,1) = width(ip+1)
    coord(17:19,1) = width(ip+1)
    coord(4,1)     = .5*(coord(3,1)+coord(5,1))
    coord(8,1)     = .5*(coord(1,1)+coord(7,1))
    coord(16,1)    = .5*(coord(15,1)+coord(17,1))
    coord(20,1)    = .5*(coord(13,1)+coord(19,1))
    coord(1:8,2)   = depth(iq)
    coord(13:20,2) = depth(iq+1)
    coord(9,2)     = .5*(coord(1,2)+coord(13,2))
    coord(10,2)    = .5*(coord(3,2)+coord(15,2))
    coord(11,2)    = .5*(coord(5,2)+coord(17,2))
    coord(12,2)    = .5*(coord(7,2)+coord(19,2))
    coord(1,3)     = -height(is+1)
    coord(7:9,3)   = -height(is+1)
    coord(12:13,3) = -height(is+1)
    coord(19:20,3) = -height(is+1)
    coord(3:5,3)   = -height(is)
    coord(10:11,3) = -height(is)
    coord(15:17,3) = -height(is)
    coord(2,3)     = .5*(coord(1,3)+coord(3,3))
    coord(6,3)     = .5*(coord(5,3)+coord(7,3))
    coord(14,3)    = .5*(coord(13,3)+coord(15,3))
    coord(18,3)    = .5*(coord(17,3)+coord(19,3))
   
    RETURN
    
  END SUBROUTINE GEOMETRY_20BXZV

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE cube_bc20(rest,nxe,nye,nze)

   !/****f* geometry/cube_bc20
   !*  NAME
   !*    SUBROUTINE: geometry_20bxzv
   !*  SYNOPSIS
   !*    Usage:      CALL cube_bc20(rest,nxe,nye,nze)
   !*  FUNCTION
   !*    Boundary conditions for simple 20-node element cuboid  
   !*  INPUTS
   !*  OUTPUTS
   !*  AUTHOR
   !*    Smith & Griffiths, Edition 4
   !******
   !*/      ! 
  
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)  :: nxe,nye,nze
  INTEGER,INTENT(OUT) :: rest(:,:)
  INTEGER             :: i,j,k,l,m,n,face,face1,face2,count

  face1 = 3*nxe*nze+2*(nxe+nze)+1 
  face2 = (nxe+1)*(nze+1)
  face  = face1 + face2  
  l     = nze*(nxe+1) 
  m     = 3*nxe + 2
  n     = 3*nxe*nze + 2*nze 
  count = 0

  DO i = 0 , nye 
     DO j = i*face+1,i*face+face1
        k = j - i*face
         IF(i==0 .OR. i==nye) THEN
          IF((k+m-1)/m*m==(k+m-1) .AND. k<=n ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,0,1/) ; CYCLE 
          ELSE IF((k+nxe+1)/m*m==(k+nxe+1) .AND. k<=n ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,0,1/) ; CYCLE 
          ELSE IF((k+nxe)/m*m==(k+nxe) .AND. k<=n ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,0,1/) ; CYCLE
          ELSE IF(k/m*m==k .AND. k<=n) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,0,1/) ; CYCLE
          ELSE IF (k<=n) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/1,0,1/) ; CYCLE
          ELSE IF (k>n) THEN
            count = count + 1 ; rest(count,1) = j; rest(count,2:) = (/0,0,0/)   
          END IF
         ELSE
          IF((k+m-1)/m*m==(k+m-1) .AND. k<=n ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,1,1/) ; CYCLE 
          ELSE IF((k+nxe+1)/m*m==(k+nxe+1) .AND. k<=n ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,1,1/) ; CYCLE 
          ELSE IF((k+nxe)/m*m==(k+nxe) .AND. k<=n ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,1,1/) ; CYCLE
          ELSE IF(k/m*m==k .AND. k<=n) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,1,1/) ; CYCLE
          ELSE IF (k>n) THEN
            count = count + 1 ; rest(count,1) = j; rest(count,2:) = (/0,0,0/)   
          END IF
         END IF
     END DO
     IF(i<nye) THEN
       DO j= i*face + face1 + 1, (i+1)*face
          k = j - (face1 + i*face)
          IF((k+nxe)/(nxe+1)*(nxe+1)==(k+nxe) .AND. k<=l ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,1,1/); CYCLE
          ELSE IF(k/(nxe+1)*(nxe+1)==k .AND. k<=l) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,1,1/); CYCLE
          ELSE IF (k>l) THEN
            count = count + 1 ; rest(count,1) = j; rest(count,2:) = (/0,0,0/)
          END IF
       END DO
     END IF
  END DO

  END SUBROUTINE cube_bc20

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE cube_bc8(rest,nxe,nye,nze)

   !/****f* geometry/cube_bc8
   !*  NAME
   !*    SUBROUTINE: cube_bc8
   !*  SYNOPSIS
   !*    Usage:      CALL cube_bc8(rest,nxe,nye,nze)
   !*  FUNCTION
   !*    Boundary conditions for simple 8-node element cuboid  
   !*  INPUTS
   !*  OUTPUTS
   !*  AUTHOR
   !*    L. Margetts 
   !******
   !*/       
  
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)  :: nxe,nye,nze
  INTEGER,INTENT(OUT) :: rest(:,:)
  INTEGER             :: i,j,k,l,m,n,face,face1,face2,count

  face  = (nxe+1)*(nze+1)
  m     = 2*nxe + 2
  n     = (nxe+1)*nze
  count = 0

  DO i = 0 , nye 
     DO j = i*face+1,i*face+face
        k = j - i*face
         IF(i==0 .OR. i==nye) THEN
          IF((k+m-1)/m*m==(k+m-1) .AND. k<=n ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,0,1/) ; CYCLE 
          ELSE IF((k+nxe+1)/m*m==(k+nxe+1) .AND. k<=n ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,0,1/) ; CYCLE 
          ELSE IF((k+nxe)/m*m==(k+nxe) .AND. k<=n ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,0,1/) ; CYCLE
          ELSE IF(k/m*m==k .AND. k<=n) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,0,1/) ; CYCLE
          ELSE IF (k<=n) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/1,0,1/) ; CYCLE
          ELSE IF (k>n) THEN
            count = count + 1 ; rest(count,1) = j; rest(count,2:) = (/0,0,0/)   
          END IF
         ELSE
          IF((k+m-1)/m*m==(k+m-1) .AND. k<=n ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,1,1/) ; CYCLE 
          ELSE IF((k+nxe+1)/m*m==(k+nxe+1) .AND. k<=n ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,1,1/) ; CYCLE 
          ELSE IF((k+nxe)/m*m==(k+nxe) .AND. k<=n ) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,1,1/) ; CYCLE
          ELSE IF(k/m*m==k .AND. k<=n) THEN
            count=count+1; rest(count,1)= j; rest(count,2:)=(/0,1,1/) ; CYCLE
          ELSE IF (k>n) THEN
            count = count + 1 ; rest(count,1) = j; rest(count,2:) = (/0,0,0/)   
          END IF
         END IF
     END DO
  END DO

  END SUBROUTINE cube_bc8

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE box_bc8(rest,nxe,nye,nze)

   !/****f* geometry/box_bc8
   !*  NAME
   !*    SUBROUTINE: box_bc8
   !*  SYNOPSIS
   !*    Usage:      CALL box_bc8(rest,nxe,nye,nze)
   !*  FUNCTION
   !*    Boundary conditions for simple 8-node element box  
   !*  INPUTS
   !*  OUTPUTS
   !*  AUTHOR
   !*    Smith & Griffiths, Edition 4
   !******
   !*/   

  IMPLICIT NONE

  INTEGER,INTENT(IN)  :: nxe,nye,nze
  INTEGER,INTENT(OUT) :: rest(:,:)
  INTEGER             :: i,j,face,count
 
  face  = (nxe+1)*(nze+1) 
  count = 0
  
  DO i = 0 , nye - 1
    DO j = i*face+1,(i+1)*face
      IF (j/(nxe+1)*(nxe+1) ==j .OR. j < i*face + nxe + 1) THEN
        count         = count + 1
        rest(count,1) = j    
        rest(count,2) = 0
      END IF
    END DO
  END DO
 
  DO j = nye*face+1, (nye+1)*face
    count         = count + 1 
    rest(count,1) = j      
    rest(count,2) = 0
  END DO
 
  END SUBROUTINE box_bc8 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE ns_cube_bc20(rest,nxe,nye,nze)

  !/****f* geometry/ns_cube_bc20
  !*  NAME
  !*    SUBROUTINE: ns_cube_bc20
  !*  SYNOPSIS
  !*    Usage:      CALL ns_cube_bc20(rest,nxe,nye,nze)
  !*  FUNCTION
  !*    Boundary conditions for simple Navier-Stokes 20-node element cuboid  
  !*  INPUTS
  !*  OUTPUTS
  !*  AUTHOR
  !*    Smith & Griffiths, Edition 4
  !******
  !*/   
  
  IMPLICIT NONE
  INTEGER,INTENT(IN)  :: nxe,nye,nze
  INTEGER,INTENT(OUT) :: rest(:,:)
  INTEGER             :: i,j,k,l,m, count,node

  node = 0;  count = 0; k = 2*nxe+1; l = nxe+1
  DO m = 0 , nye
         IF(m==0 .OR. m==nye) THEN    ! front or back
          DO j=1,nze
             DO i=1,k
                IF((j==1.AND.i==1).OR.(j==1.AND.MOD(i,2)==0)) THEN 
                   node = node +1; count = count + 1
                   rest(count,1)=node; rest(count,2:)=(/1,0,0,0/)
                ELSE IF (j==1) THEN
                   node = node +1; count = count + 1
                   rest(count,1)=node; rest(count,2:)=(/1,1,0,0/)
                ELSE IF ((j>1.AND.mod(i,2)/=0)) THEN
                   node = node +1 ; count = count + 1
                   rest(count,1)=node; rest(count,2:)=(/0,1,0,0/)
                ELSE IF (mod(i,2)==0) THEN
                   node = node +1; count = count + 1
                   rest(count,1)=node; rest(count,2:)=(/0,0,0,0/)
                END IF
             END DO
             DO i=1,l
                node = node +1; count = count + 1
                rest(count,1)=node; rest(count,2:)=(/0,0,0,0/)
             END DO
          END DO
          DO i = 1,k
             IF(mod(i,2)/=0) THEN
                node = node +1; count = count + 1
                rest(count,1)=node; rest(count,2:)=(/0,1,0,0/)
             ELSE
                node = node + 1; count = count + 1
                rest(count,1)=node; rest(count,2:)=(/0,0,0,0/)
             END IF
          END DO
         ELSE
          DO j=1,nze        !  interiors
             DO i=1,k
                IF((j==1.AND.i==1).OR.(j==1.AND.mod(i,2)==0)) THEN 
                   node = node +1; count = count + 1
                   rest(count,1)=node; rest(count,2:)=(/1,0,0,0/)
                ELSE IF (j==1) THEN
                   node = node +1; count = count + 1
                   rest(count,1)=node; rest(count,2:)=(/1,1,0,0/)
                ELSE IF ((j>1.AND.i==1).OR.(j>1.AND.i==k)) THEN
                   node = node +1; count = count + 1
                   rest(count,1)=node; rest(count,2:)=(/0,1,0,0/)
                ELSE IF (j>1.AND.mod(i,2)/=0) THEN
                   node = node +1 
                ELSE IF (j>1.AND.mod(i,2)==0) THEN
                   node = node +1; count = count + 1
                   rest(count,1)=node; rest(count,2:)=(/1,0,1,1/)
                END IF
             END DO
             DO i=1,l
                IF(i==1.OR.i==l) THEN 
                   node = node +1; count = count + 1
                   rest(count,1)=node; rest(count,2:)=(/0,0,0,0/)
                ELSE
                   node = node +1; count = count + 1
                   rest(count,1)=node; rest(count,2:)=(/1,0,1,1/)
                END IF
             END DO
          END DO
          DO i = 1,k    ! base
             IF(mod(i,2)/=0) THEN
                node = node +1; count = count + 1
                rest(count,1)=node; rest(count,2:)=(/0,1,0,0/)
             ELSE
                node = node + 1; count = count + 1
                rest(count,1)=node; rest(count,2:)=(/0,0,0,0/)
             END IF
          END DO
         END IF
         IF(m<nye) THEN
          DO j=1,nze
             DO i=1,l
                IF(j==1) THEN
                 node = node + 1; count = count + 1
                 rest(count,1)=node; rest(count,2:)=(/1,0,0,0/)
                ELSE IF(i==1.OR.i==l) THEN
                 node = node + 1; count = count + 1
                 rest(count,1)=node; rest(count,2:)=(/0,0,0,0/)
                ELSE
                 node = node + 1; count = count + 1
                 rest(count,1)=node; rest(count,2:)=(/1,0,1,1/)
                END IF
             END DO
          END DO
          DO i=1,l
              node = node + 1; count = count + 1
              rest(count,1)=node; rest(count,2:)=(/0,0,0,0/)
          END DO
         END IF
  END DO
  
  END SUBROUTINE ns_cube_bc20

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

 SUBROUTINE ns_loading(no,nxe,nye,nze)

   !/****f* geometry/ns_loading
   !*  NAME
   !*    SUBROUTINE: ns_cube_bc20
   !*  SYNOPSIS
   !*    Usage:      CALL ns_loading(no,nxe,nye,nze)
   !*  FUNCTION
   !*    Simple lid displacement of a square  Navier-Stokes patch  
   !*  INPUTS
   !*  OUTPUTS
   !*  AUTHOR
   !*    Smith & Griffiths, Edition 4
   !******
   !*/  
   
   IMPLICIT NONE
   INTEGER,INTENT(IN)  :: nxe,nye,nze
   INTEGER,INTENT(OUT) :: no(:)
   INTEGER             :: f1,f2,f3,count,i,j
   
   f1 = 3*nxe+1 + (nxe+1)*nze
   f2 = 3*(nxe-1)*(nze-1)+(nxe+1)
   f3 = 10*nxe*nze -3*nxe - 5*nze + 4
   
   no(1) = 1
   count = 2
   
   DO  i = 1,nxe   ! front face
     no(count) = 3*i -1
     count     = count+1
     no(count) = 3*i
     count     = count + 1
   END DO
 
   DO j = 0 , nye-1    ! other face pairs
     DO i = 1 , nxe + 1   ! mid - planes
       no(count) = i + f1 + j*(f2 + f3)    
       count     = count + 1
     END DO
     no(count) = f1+(j+1)*f2+j*f3 +1
     count     = count + 1
     DO i = 1 , nxe 
       no(count) = 3*i-1 + f1 +(j+1) * f2 + j * f3   
       count     = count + 1
       no(count) = 3*i + f1 +(j+1) * f2 + j * f3   
       count     = count + 1
     END DO
   END DO

   RETURN
   
 END SUBROUTINE ns_loading

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

END MODULE GEOMETRY
