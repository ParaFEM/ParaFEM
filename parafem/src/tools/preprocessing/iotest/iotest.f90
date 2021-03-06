PROGRAM iotest
! Program to prototype, implement and test I/O subroutines
! Loads array need sorting out

 USE precision
 USE, INTRINSIC :: ISO_C_BINDING

 IMPLICIT NONE
 
 INTEGER   :: i, j, ndim, nod, nn, nels, nodof, nfe
 INTEGER(4):: parts=1
 INTEGER(4):: nn_c=8
 REAL(iwp) :: zero = 0.0_iwp
 
 INTEGER, ALLOCATABLE   :: g_num(:,:), nf(:,:), etype(:)
 REAL(iwp), ALLOCATABLE :: g_coord(:,:), loads(:,:), disp(:,:)

 REAL(c_float), ALLOCATABLE  :: g_coord_c(:,:)
 INTEGER(c_int), ALLOCATABLE :: g_num_c(:,:)
 
 LOGICAL   :: hardcoded=.true.

 CHARACTER(80) :: cbuffer

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 1. Hard coded data for one element
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  PRINT *, ""
  PRINT *, "Program IOTEST"
  PRINT *, ""
  PRINT *, "    - Used to test the implementation of data formats"
  PRINT *, "    - Currently writes ascii data into file named binary"
  PRINT *, "    - Under development"
  PRINT *, ""  

  ndim  = 3
  nod   = 8
  nn    = 8
  nels  = 1
  nodof = ndim + 1
  
  ALLOCATE (g_coord(ndim,nn))
  ALLOCATE (g_coord_c(ndim,nn))
  ALLOCATE (g_num(nod,nels))
  ALLOCATE (g_num_c(nod,nels))
  ALLOCATE (nf(nodof,nn))
  ALLOCATE (etype(nels))
  ALLOCATE (disp(ndim,nn))
  
  g_coord = zero
  g_num   = 0
  nf      = 0
  etype   = 0
  disp    = zero
  
  IF(hardcoded) THEN
  
    g_coord(:,1) = (/0.0,0.0,0.0/)
    g_coord(:,2) = (/0.0,0.0,1.0/)
    g_coord(:,3) = (/1.0,0.0,1.0/)
    g_coord(:,4) = (/1.0,0.0,0.0/)
    g_coord(:,5) = (/0.0,1.0,0.0/)
    g_coord(:,6) = (/0.0,1.0,1.0/)
    g_coord(:,7) = (/1.0,1.0,1.0/)
    g_coord(:,8) = (/1.0,1.0,0.0/)
 
    g_coord_c    = g_coord
 
    g_num(:,1) = (/1,2,3,4,5,6,7,8/)
    g_num_c    = g_num
  
    nf(:,1) = (/1,0,0,0/)
    nf(:,2) = (/2,1,1,1/)
    nf(:,3) = (/3,1,1,1/)
    nf(:,4) = (/4,1,1,0/)
    nf(:,5) = (/5,1,1,0/)
    nf(:,6) = (/6,1,1,1/)
    nf(:,7) = (/7,1,1,1/)
    nf(:,8) = (/8,1,1,0/)
  
    etype(1) = 1
  
    disp(:,1) = (/0.0,0.0,0.0/)
    disp(:,2) = (/0.0,0.0,-0.1/)
    disp(:,3) = (/0.0,0.0,-0.2/)
    disp(:,4) = (/0.0,0.0,0.0/)
    disp(:,5) = (/0.0,0.0,0.0/)
    disp(:,6) = (/0.0,0.0,-0.1/)
    disp(:,7) = (/0.0,0.0,-0.2/)
    disp(:,8) = (/0.0,0.0,0.0/)
    
  END IF
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 2. ASCII Ensight Gold Format
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! 2.1 Write ASCII Ensight Gold Case File
!------------------------------------------------------------------------------

  OPEN(20,FILE='ASCII.ensi.case')
  
  WRITE(20,'(A/A)')    "#", "# Post-processing file generated by program      &
                             &iotest.f90 "
  WRITE(20,'(A/A/A)')  "#","# Ensight Gold Format","#"
  WRITE(20,'(2A/A)')   "# Problem name: ","iotest","#"
  WRITE(20,'(A/A/A)')  "FORMAT","type:  ensight gold","GEOMETRY"
  WRITE(20,'(2A/A)')   "model: 1  ","ASCII.ensi.geo","VARIABLE"
  WRITE(20,'(2A)')     "scalar per element:  material      ",                 &
                       "ASCII.ensi.MATID"
  WRITE(20,'(2A)')     "scalar per node:     restraint     ",                 &
                       "ASCII.ensi.NDBND"
  WRITE(20,'(2A)')     "vector per node:     displacement  ",                 &
                       "ASCII.ensi.DISPL-******"
  WRITE(20,'(A/A)')    "TIME", "time set:   1"
  WRITE(20,'(A)')      "number of steps:   1"
  WRITE(20,'(A)')      "filename start number:   1"
  WRITE(20,'(A)')      "filename increment:   1"
  WRITE(20,'(A)')      "time values:"
  WRITE(20,'(A)')      "1.0"
  
  CLOSE(20)

!------------------------------------------------------------------------------
! 2.2 Write ASCII geometry file
!------------------------------------------------------------------------------
  
  OPEN(21,FILE='ASCII.ensi.geo')
  
  WRITE(21,'(/2A)')   "Problem name: ", "ASCII"
  WRITE(21,'(A/A/A)') "Geometry files","node id given","element id given"
  WRITE(21,'(A/A)')   "part","      1"
  WRITE(21,'(A)')     "Volume Mesh"
  WRITE(21,'(A)')     "coordinates"
  
  WRITE(21,'(I10)') nn
  DO j=1,ndim
    DO i=1,nn  
      WRITE(21,'(E12.5)') g_coord(j,i)
    END DO
  END DO
  
  WRITE(21,'(A/I10)') "hexa8",nels

  DO i = 1,nels
    WRITE(21,'(8I10)') g_num(1,i),g_num(4,i),g_num(8,i),g_num(5,i),           &
                       g_num(2,i),g_num(3,i),g_num(7,i),g_num(6,i)
  END DO

  CLOSE(21)
    
!------------------------------------------------------------------------------
! 2.3 Write ASCII Boundary Conditions
!------------------------------------------------------------------------------

  OPEN(22,FILE='ASCII.ensi.NDBND')
  WRITE(22,'(A)')     "Alya Ensight Gold --- Scalar per-node variable file"
  WRITE(22,'(A/A/A)') "part", "      1","coordinates"

  DO i=1,UBOUND(g_coord,2) 
    nfe=0
    IF(nf(1,i)==0) nfe=nfe+1
    IF(nf(2,i)==0) nfe=nfe+2
    IF(nf(3,i)==0) nfe=nfe+4
    WRITE(22,'(I2)') nfe
  END DO

  CLOSE(22)

!------------------------------------------------------------------------------
! 2.4 Write ASCII Loaded Nodes
!------------------------------------------------------------------------------

! OPEN(23,FILE='ASCII.ensi.NDLDS')
! WRITE(23,'(A)')     "Alya Ensight Gold --- Vector per-node variable file"
! WRITE(23,'(A/A/A)') "part", "      1","coordinates"

! DO j=1,UBOUND(nf,1)
!   DO i=1, UBOUND(nf,2)
!     WRITE(23,'(E12.5)') loads(nf(j,i))
!   END DO
! END DO

! CLOSE(23)

!------------------------------------------------------------------------------
! 2.5 Write ASCII file containing material IDs
!------------------------------------------------------------------------------
  
  OPEN(24,FILE='ASCII.ensi.MATID')

  WRITE(24,'(A)') "Alya Ensight Gold --- Scalar per-element variable file"
  WRITE(24,'(A/A)') "part", "      1"
  WRITE(24,'(A)') "hexa8"
  DO i=1,nels; WRITE(24,'(I10)') etype(i); END DO  

  CLOSE(24)
  
!------------------------------------------------------------------------------
! 2.6 Write ASCII Displacements
!------------------------------------------------------------------------------

  OPEN(25,file='ASCII.ensi.DISPL-000001',status='replace',action='write')
       
  WRITE(25,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
  WRITE(25,'(A/A/A)') "part", "     1","coordinates"  

  DO j=1,ndim
    DO i=1,nn
      WRITE(25,'(e12.4)') disp(j,i)
    END DO
  END DO

  CLOSE(25)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! 3. BINARY Ensight Gold Format
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! 3.1 Write BINARY Ensight Gold Case File
!------------------------------------------------------------------------------

  OPEN(26,FILE='BINARY.ensi.case')
  
  WRITE(26,'(A/A)')    "#", "# Post-processing file generated by program      &
                             &iotest.f90 "
  WRITE(26,'(A/A/A)')  "#","# Ensight Gold Format","#"
  WRITE(26,'(2A/A)')   "# Problem name: ","Binary IO Test","#"
  WRITE(26,'(A/A/A)')  "FORMAT","type:  ensight gold","GEOMETRY"
  WRITE(26,'(2A/A)')   "model: 1  ","BINARY.ensi.geo"!,"VARIABLE"
! WRITE(26,'(2A)')     "scalar per element:  material      ",                 &
!                      "BINARY.ensi.MATID"
! WRITE(26,'(2A)')     "scalar per node:     restraint     ",                 &
!                      "BINARY.ensi.NDBND"
! WRITE(26,'(2A)')     "vector per node:     displacement  ",                 &
!                      "BINARY.ensi.DISPL-******"
! WRITE(26,'(A/A)')    "TIME", "time set:   1"
! WRITE(26,'(A)')      "number of steps:   1"
! WRITE(26,'(A)')      "filename start number:   1"
! WRITE(26,'(A)')      "filename increment:   1"
! WRITE(26,'(A)')      "time values:"
! WRITE(26,'(A)')      "1.0"

  CLOSE(26)

!------------------------------------------------------------------------------
! 3.2 Write BINARY Geometry File
!------------------------------------------------------------------------------

  OPEN(27,FILE="BINARY.ensi.geo",STATUS="REPLACE", FORM="UNFORMATTED",        &
                ACTION="WRITE", &
                ACCESS="STREAM")
  
  cbuffer = "C Binary"                    ; WRITE(27) cbuffer
  cbuffer = "Problem name: IOTEST BINARY" ; WRITE(27) cbuffer
  cbuffer = "Geometry files"              ; WRITE(27) cbuffer
  cbuffer = "node id off"                 ; WRITE(27) cbuffer
  cbuffer = "element id off"              ; WRITE(27) cbuffer
  cbuffer = "part"                        ; WRITE(27) cbuffer
  WRITE(27) int(1,kind=c_int)
  cbuffer = "Volume"                      ; WRITE(27) cbuffer
  cbuffer = "coordinates"                 ; WRITE(27) cbuffer

  WRITE(27) nn_c
  DO j=1,ndim
    DO i=1,nn  
!     WRITE(27) real(g_coord_c(j,i),kind=4)
      WRITE(27) real(g_coord(j,i),kind=c_float)
!     WRITE(27) g_coord(j,i)
    END DO
  END DO
 
  cbuffer = "hexa8" ; WRITE(27) cbuffer 
  WRITE(27) int(1,kind=c_int)

  DO i = 1,nels
    WRITE(27) g_num_c(1,i),g_num_c(4,i),g_num_c(8,i),g_num_c(5,i),             &
              g_num_c(2,i),g_num_c(3,i),g_num_c(7,i),g_num_c(6,i)
  END DO


  CLOSE(27)
  
!------------------------------------------------------------------------------
! 3.3 Write BINARY Boundary Conditions
!------------------------------------------------------------------------------

  OPEN(28,FILE='BINARY.ensi.NDBND')
  WRITE(28,'(A)')     "Alya Ensight Gold --- Scalar per-node variable file"
  WRITE(28,'(A/A/A)') "part", "      1","coordinates"

  DO i=1,UBOUND(g_coord,2) 
    nfe=0
    IF(nf(1,i)==0) nfe=nfe+1
    IF(nf(2,i)==0) nfe=nfe+2
    IF(nf(3,i)==0) nfe=nfe+4
    WRITE(22,'(I2)') nfe
  END DO

  CLOSE(28)
  
!------------------------------------------------------------------------------
! 3.4 Write BINARY Loaded Nodes
!------------------------------------------------------------------------------

! OPEN(29,FILE='BINARY.ensi.NDLDS')
! WRITE(29,'(A)')     "Alya Ensight Gold --- Vector per-node variable file"
! WRITE(29,'(A/A/A)') "part", "      1","coordinates"

! DO j=1,UBOUND(nf,1)
!   DO i=1, UBOUND(nf,2)
!     WRITE(23,'(E12.5)') loads(nf(j,i))
!   END DO
! END DO

! CLOSE(29)

!------------------------------------------------------------------------------
! 2.5 Write BINARY file containing material IDs
!------------------------------------------------------------------------------
  
  OPEN(30,FILE='BINARY.ensi.MATID')

  WRITE(30,'(A)') "Alya Ensight Gold --- Scalar per-element variable file"
  WRITE(30,'(A/A)') "part", "      1"
  WRITE(30,'(A)') "hexa8"
  DO i=1,nels; WRITE(24,'(I10)') etype(i); END DO  

  CLOSE(30)

!------------------------------------------------------------------------------
! 2.6 Write ASCII Displacements
!------------------------------------------------------------------------------

  OPEN(31,file='BINARY.ensi.DISPL-000001',status='replace',action='write')
       
  WRITE(31,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
  WRITE(31,'(A/A/A)') "part", "     1","coordinates"  

  DO j=1,ndim
    DO i=1,nn
      WRITE(25,'(e12.4)') disp(j,i)
    END DO
  END DO

  CLOSE(31)
  
END PROGRAM iotest
