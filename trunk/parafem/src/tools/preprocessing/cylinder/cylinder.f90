!$Id: cylinder.f90 31 2016-01-24 21:47:27Z mexas $
!
! Copyright (c) 2015 Anton Shterenlikht, The University of Bristol, UK
!
! Mesh a circular cylinder.
! Takes 4 input arguments:
! 1 - radius
! 2 - length
! 3 - max uncertainty (error) in coordinates. The actual coord
!     will be coord = coord * ( 1 +- alpha ) where
!     alpha is a random number in the interval [ - input3 .. input3 ].
!     The intention is that this input is << 1.
!     For perfectly cylindrical models set it to 0.
! 4 - label to attach to all output files: nodes, elements, bc
! 5 - nel - number of element on a side of the square

program cylinder
use mesh
implicit none

!Variables for subroutine MESH_ENSI_CASE  
INTEGER          :: nlen,nstep,npri
REAL(8)          :: dtim
LOGICAL          :: solid

!Variables for MESH_ENSI_MATID_BIN, MESH_ENSI_NDLDS_BIN

double precision, ALLOCATABLE  :: loads(:)
double precision, ALLOCATABLE  :: etype_1(:)
INTEGER          :: nod_1
INTEGER          :: nels_1
INTEGER, ALLOCATABLE :: nf_1(:,:)


!Variables for subroutine ESH_ENSI_GEO_BIN
INTEGER, ALLOCATABLE           :: g_num_1(:,:)
double precision, ALLOCATABLE  :: g_coord_1(:,:)
CHARACTER(LEN=15)              :: element1 

! length of the end bits of the bars. The surface of
! the bar at each length is fixed in 12 plane
 double precision, parameter :: fixed_len = 20.0d0

! temp files
character( len=80 ) :: filebot="tmp_file_bot_nodes"

!integer, parameter :: nel=8, arrlen = (nel+1)**2

double precision :: rad, length, dist, total_len, cerr
integer :: nel, arglen, errstat, i, nel3, counter, springel1
character( len=80 ) :: value , nodefile, elfile, bcfile, springfile,   &
 filesuff
!type( node ) :: circ( arrlen )
type( node ), allocatable :: circ( : )

! read the values from command line

do i = 1, 5
   call get_command_argument( i, value, arglen, errstat )
   if ( errstat .gt. 0  ) then
     write (*,*) "Not enough arguments, must call as cylynder.dx&
                 & 1)radius 2)length 3)error 4)label 5)nel"
     stop
   end if
   if ( errstat .eq. -1 ) stop "ERROR: argument too long"
   if ( errstat .lt. -1 )                                              &
      stop "ERROR: unknown error, should never end up here"

   if ( i .eq. 1 ) then
     ! add a decimal point if it's not there
     if ( index( value, "." ) .eq. 0 ) value = trim(value)//"."
     read( value , "(es30.15)" ) rad
   end if

   if ( i .eq. 2 ) then
     ! add a decimal point if it's not there
     if ( index( value, "." ) .eq. 0 ) value = trim(value)//"."
     read( value , "(es30.15)" ) length
   end if

   if ( i .eq. 3 ) then
     ! add a decimal point if it's not there
     if ( index( value, "." ) .eq. 0 ) value = trim(value)//"."
     read( value , "(es30.15)" ) cerr
   end if

   if ( i .eq. 4 ) filesuff = value

   if ( i .eq. 5 ) then
     ! must be integerbrickel
     read( value , "(i3)" ) nel
   end if

end do

! set other values
  total_len = length + 2.0d0 * fixed_len
       nel3 = int( total_len ) ! FE size of 1mm
       dist = total_len / nel3
   nodefile = "nodes"//trim(filesuff)
     elfile = "elements"//trim(filesuff)
     bcfile = "bc"//trim(filesuff)
 springfile = "springel"//trim(filesuff)

allocate( circ( (nel+1)**2 ), stat=errstat )
if ( errstat .ne. 0 ) stop "ERROR: cylinder: allocate( circ )" 

call circle( nel , circ )
 !subroutine cyldmp( fname, n , rad, dist, array )
!call cyldmp( trim(nodefile), nel3, rad, dist, circ ) 
 !subroutine cyldmp_err( fname, nbot, nmid, ntop, rad, dist, err, array )

call cyldmp_err( trim( nodefile ), int( fixed_len), int( length ),     &
                 int( fixed_len ), rad, dist, cerr, circ, g_coord_1 )

ALLOCATE(g_num_1(8,nel3*nel*nel))

call brickel( trim(elfile), nel, nel, nel3, g_num_1 )


! Nodes in nset=top are moving along 3
! Nodes in nset=bot are constrained, or attached to spring1 elements

 !subroutine bound( fread, fwrite, nset, start, end, rad )
call bound( trim(nodefile), trim(filebot), "bot", 0.0d0, fixed_len,    &
            rad, counter )
call bound( trim(nodefile), "tmpzzz2", "top", length+fixed_len,        &
            total_len, rad, counter )
write( *,* ) "nodes in nset=top:", counter

! Do not delete bot nodes file, filebot, until spring
! elements have been created.

! The first spring element number
springel1 = nel**2 * nel3 + 1
call springel( trim(filebot), trim(springfile), springel1 )

! merge bot and top node sets into a single file and delete temp files
call bcmerge( trim(filebot), "tmpzzz2", trim(bcfile) )

!Calculating the length of the name file
nlen = len_trim(filesuff)

!Defining variables for MESH_ENSI_CASE and MESH_ENSI_GEO_BIN subroutines
nstep = 1
npri = 1
dtim = 1.0
solid = .true.
element1 = 'hexahedron'
nod_1 = 8

!Allocation of variables for subroutines 
ALLOCATE(etype_1(nel3*nel*nel))
ALLOCATE(nf_1(4,nels_1))
ALLOCATE(loads(nels_1))


!Material properties assignation
do i = 1, nel3*nel*nel
	etype_1(i) = 10000
end do


do i = 1, nels_1

	loads(i) = 1.0
	nf_1(1,i) = i
	nf_1(2,i) = 0
	nf_1(3,i) = 0
	nf_1(4,i) = 0

end do



!CALL ParaFEM subroutines to produce ensi.case files
CALL MESH_ENSI_CASE(filesuff,nlen,nstep,npri,dtim,solid)
CALL MESH_ENSI_GEO_BIN(filesuff,nlen,g_coord_1,g_num_1,element1)


CALL MESH_ENSI_MATID_BIN(filesuff,nlen,nod_1,element1,etype_1)
CALL MESH_ENSI_NDLDS_BIN(filesuff,nlen,nf_1,loads)
!CALL MESH_ENSI_NDBND_BIN(filesuff,nf_1,nlen,nod_1,solid)

CONTAINS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE MESH_ENSI_CASE(argv,nlen,nstep,npri,dtim,solid)

   !/****f* input/mesh_ensi_case
   !*  NAME
   !*    SUBROUTINE: mesh_ensi_case
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi_case(argv,nlen,nstep,npri,dtim,solid)
   !*  FUNCTION
   !*    This subroutine outputs the "case" file required for the C binary 
   !*    version of the Ensight gold format. Models in this format can be 
   !*    viewed in ParaView.
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*    npri             : print interval
   !*	 nstep            : number of time steps in analysis
   !*
   !*    Scalar reals
   !*    dtim             : time step
   !*
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*
   !*    Scalar logicals
   !*    solid            : type of analysis solid if .true. fluid if .false.
   !*
   !*  OUTPUTS
   !*  AUTHOR
   !*    I.M. Smith
   !*    D.V. Griffiths
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*
   !*/


    USE, INTRINSIC :: ISO_C_BINDING
    
    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER,   INTENT(IN)         :: nlen,nstep,npri
    INTEGER                       :: i,j,k,l,m,n,nfe,nod,nels,ndim,nn
    INTEGER                       :: prnwidth,remainder
    REAL(iwp), INTENT(IN)         :: dtim
    CHARACTER(LEN=15), INTENT(IN) :: argv
    CHARACTER(LEN=80)             :: cbuffer
    LOGICAL, INTENT(IN)           :: solid
    
  !-----------------------------------------------------------------------------
  ! 1. Write case file
  !-----------------------------------------------------------------------------
  
    OPEN(12,FILE=argv(1:nlen)//'.bin.ensi.case')
  
    WRITE(12,'(A/A)')    "#", "# Post-processing file generated by subroutine  &
                               &WRITE_ENSI in "
    WRITE(12,'(A,A,/A)') "#"," Smith, Griffiths and Margetts, 'Programming the &
                               &Finite Element Method',","# Wiley, 2013."        
    WRITE(12,'(A/A/A)')  "#","# Ensight Gold Format","#"
    WRITE(12,'(2A/A)')   "# Problem name: ",argv(1:nlen),"#"
    WRITE(12,'(A/A/A)')  "FORMAT","type:  ensight gold","GEOMETRY"
    WRITE(12,'(2A/A)')   "model: 1  ",argv(1:nlen)//'.bin.ensi.geo',"VARIABLE"
    WRITE(12,'(2A)')     "scalar per element:  material      ",                &
                          argv(1:nlen)//'.bin.ensi.MATID'
    IF(solid) THEN
      WRITE(12,'(2A)')   "scalar per node:     restraint     ",                &
                          argv(1:nlen)//'.bin.ensi.NDBND'
      WRITE(12,'(2A)')   "vector per node:     displacement  ",                &
                          argv(1:nlen)//'.bin.ensi.DISPL-******'
    ELSE
      WRITE(12,'(2A)')   "scalar per node:     pressure      ",                &
                          argv(1:nlen)//'.bin.ensi.PRESSURE-******'
    END IF
    WRITE(12,'(2A)')     "vector per node:     load          ",                &
                          argv(1:nlen)//'.bin.ensi.NDLDS'
    WRITE(12,'(A/A)')     "TIME","time set:     1"
    WRITE(12,'(A,I5)')    "number of steps:",nstep/npri
    WRITE(12,'(A,I5)')    "filename start number:",npri
    WRITE(12,'(A,I5)')    "filename increment:",npri
    WRITE(12,'(A)')       "time values:"
    prnwidth  = 5
    remainder = mod(nstep/npri,prnwidth)
    n         = ((nstep/npri) - remainder)/prnwidth
    IF(nstep/npri<=prnwidth) THEN
      DO i=1,nstep,npri
        IF(i==nstep) THEN
          WRITE(12,'(E12.5)') i*dtim
        ELSE
          WRITE(12,'(E12.5)',ADVANCE='no') i*dtim
        END IF
      END DO
    ELSE
      IF(remainder==0) THEN
        DO j=1,n
          m = ((j-1)*prnwidth)+1
          l = ((j-1)*prnwidth)+prnwidth
          WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
        END DO
      ELSE
  !     DO j=1,n-1
        DO j=1,n
          m = ((j-1)*prnwidth)+1
          l = ((j-1)*prnwidth)+prnwidth
          WRITE(12,'(5E12.5)') (k*dtim,k=m,l)
        END DO
        m = (n*prnwidth)+1
        l = (n*prnwidth)+remainder
        DO i=m,l
          IF(i==l) THEN
            WRITE(12,'(E12.5)') dtim*i
          ELSE
            WRITE(12,'(E12.5)',ADVANCE='no') dtim*i
          END IF
        END DO
      END IF
    END IF
   
    CLOSE(12)
    
    RETURN
  
  END SUBROUTINE MESH_ENSI_CASE


SUBROUTINE MESH_ENSI_GEO_BIN(argv,nlen,g_coord,g_num,element)

   !/****f* input/mesh_ensi_geo_bin
   !*  NAME
   !*    SUBROUTINE: mesh_ensi_geo_bin
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi_geo_bin(argv,nlen,g_coord,g_num,element)
   !*  FUNCTION
   !*    This subroutine outputs the "geo" file in the C binary version of the 
   !*    Ensight gold format. Models in this format can be viewed in ParaView.
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*	 element          : element type
   !*
   !*    Scalar logicals
   !*
   !*    Dynamic scalar arrays
   !*    g_num            : global element node numbers vector
   !* 
   !*    Dynamic real arrays
   !* 	 g_coord          : global nodal coordinates
   !*
   !*  OUTPUTS
   !*  AUTHOR
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !* 
   !*  Used in program p12meshgenbin
   !*/

    USE, INTRINSIC :: ISO_C_BINDING
    
    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER,   INTENT(IN)         :: nlen
    INTEGER,   INTENT(INOUT)      :: g_num(:,:)
    INTEGER                       :: i,j
    INTEGER                       :: nod,nels,ndim,nn
    REAL(iwp), INTENT(INOUT)      :: g_coord(:,:)
    CHARACTER(LEN=15), INTENT(IN) :: argv,element  
    CHARACTER(LEN=80)             :: cbuffer
    
  !----------------------------------------------------------------------------
  ! 1. Initialisation
  !----------------------------------------------------------------------------
  
    nn   = UBOUND(g_coord,2) ; ndim = UBOUND(g_coord,1)
    nels = UBOUND(g_num,2)   ; nod  = UBOUND(g_num,1)
  
  !----------------------------------------------------------------------------
  ! 2. Write geometry file
  !
  !    Only 8 node bricks tested
  !----------------------------------------------------------------------------
  
    OPEN(13,FILE=argv(1:nlen)//'.bin.ensi.geo',STATUS="REPLACE",              &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

    cbuffer = "C Binary"                     ; WRITE(13) cbuffer
    cbuffer = "Problem name: "//argv(1:nlen) ; WRITE(13) cbuffer
    cbuffer = "Geometry files"               ; WRITE(13) cbuffer
    cbuffer = "node id off"                  ; WRITE(13) cbuffer
    cbuffer = "element id off"               ; WRITE(13) cbuffer
    cbuffer = "part"                         ; WRITE(13) cbuffer
    WRITE(13) int(1,kind=c_int)
    IF(ndim==2) THEN 
       cbuffer = "2d-mesh"                   ; WRITE(13) cbuffer
    END IF
    IF(ndim==3) THEN
       cbuffer = "Volume"                    ; WRITE(13) cbuffer
    END IF
    cbuffer = "coordinates"                  ; WRITE(13) cbuffer
    
    WRITE(13) int(nn,kind=c_int)
    DO j=1,ndim
      DO i=1,nn  
        
        WRITE(13) real(g_coord(j,i),kind=c_float)
	
      END DO
    END DO


    IF(ndim==2) THEN ! ensight requires zeros for the z-ordinate
      DO i=1,nn
        WRITE(13,'(A)') " 0.00000E+00" ! needs fixing for binary
      END DO
    END IF
  
    SELECT CASE(element)
      CASE('triangle')
        SELECT CASE(nod)
!         CASE(3)
!           WRITE(13,'(A/I10)') "tria3", nels
!           DO i = 1,nels
!             WRITE(13,'(3I10)')g_num(3,i),g_num(2,i),g_num(1,i)
!           END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('quadrilateral')
        SELECT CASE(nod)
          CASE(4)
            WRITE(13,'(A/I10)') "quad4", nels
            DO i = 1,nels
              WRITE(13,'(4I10)')g_num(1,i),g_num(4,i),g_num(3,i),g_num(2,i)
            END DO
          CASE(8)
            WRITE(13,'(A/I10)') "quad8", nels
            DO i = 1,nels
              WRITE(13,'(8I10)')g_num(1,i),g_num(7,i),g_num(5,i),g_num(3,i),  &
                                g_num(8,i),g_num(6,i),g_num(4,i),g_num(2,i)
             

            END DO
          CASE DEFAULT
            WRITE(13,'(A)')   "# Element type not recognised"
        END SELECT
      CASE('hexahedron')
        SELECT CASE(nod)
          CASE(8)
            cbuffer = "hexa8"       ; WRITE(13) cbuffer
            WRITE(13) int(nels,kind=c_int)
            DO i = 1,nels
              WRITE(13) int(g_num(1,i),kind=c_int),int(g_num(4,i),kind=c_int),&
                        int(g_num(8,i),kind=c_int),int(g_num(5,i),kind=c_int),&
                        int(g_num(2,i),kind=c_int),int(g_num(3,i),kind=c_int),&
                        int(g_num(7,i),kind=c_int),int(g_num(6,i),kind=c_int)
		
            END DO 
          CASE(20)
            cbuffer = "hexa20"       ; WRITE(13) cbuffer
            WRITE(13) int(nels,kind=c_int)
            DO i = 1,nels
              WRITE(13)                                                       &
                int(g_num(1,i),kind=c_int), int(g_num(7,i),kind=c_int),       &
                int(g_num(19,i),kind=c_int),int(g_num(13,i),kind=c_int),      &
                int(g_num(3,i),kind=c_int),int(g_num(5,i),kind=c_int),        &
                int(g_num(17,i),kind=c_int),int(g_num(15,i),kind=c_int),      &
                int(g_num(8,i),kind=c_int),int(g_num(12,i),kind=c_int),       &
                int(g_num(20,i),kind=c_int),int(g_num(9,i),kind=c_int),       &
                int(g_num(4,i),kind=c_int),int(g_num(11,i),kind=c_int),       &
                int(g_num(16,i),kind=c_int),int(g_num(10,i),kind=c_int),      &
                int(g_num(2,i),kind=c_int),int(g_num(6,i),kind=c_int),        &
                int(g_num(18,i),kind=c_int),int(g_num(14,i),kind=c_int) 
	
            END DO		
          CASE DEFAULT
            cbuffer = "# Element type not recognised" ; WRITE(13) cbuffer
        END SELECT
      CASE('tetrahedron')
        SELECT CASE(nod)
          CASE(4)
            cbuffer = "tetra4" ; WRITE(13)
            WRITE(13) int(nels,kind=c_int)
            DO i = 1,nels
              WRITE(13) int(g_num(1,i),kind=c_int),int(g_num(3,i),kind=c_int),&
                        int(g_num(2,i),kind=c_int),int(g_num(4,i),kind=c_int)
            END DO
          CASE DEFAULT
            cbuffer = "# Element type not recognised" ; WRITE(13)
        END SELECT
      CASE DEFAULT
        cbuffer = "# Element type not recognised" ; WRITE(13)
    END SELECT
  
    CLOSE(13)

    RETURN
  
  END SUBROUTINE MESH_ENSI_GEO_BIN

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE MESH_ENSI_MATID_BIN(argv,nlen,nod,element,etype)

   !/****f* input/mesh_ensi_matid_bin
   !*  NAME
   !*    SUBROUTINE: mesh_ensi_matid_bin
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi_matid_bin(argv,nlen,nod,element,etype)
   !*  FUNCTION
   !*    This subroutine outputs material type for each element in the mesh
   !*    in the C binary version of the Ensight gold format. Models in this 
   !*    format can be viewed in ParaView.
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*    nod              : number of nodes in the element
   !*
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*	 element          : element type
   !*
   !*    Dynamic scalar arrays
   !*    etype            : element property type vector
   !* 
   !*  OUTPUTS
   !*  AUTHOR
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*
   !*  ParaView has a bug which prevents the 'Material' section of the 
   !*  ENSIGHT Gold to be read and therefore integer MATID values
   !*  http://www.paraview.org/Bug/view.php?id=15151
   !*  http://www3.ensight.com/EnSight10_Docs/UserManual.pdf pp.713
   !*  Workaround - use MATIDs as reals and convert into int for ParaFEM
   !*
   !*/

    USE, INTRINSIC :: ISO_C_BINDING
    
    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER,   INTENT(IN)         :: nlen,nod
    REAL(iwp), INTENT(IN)         :: etype(:)
    CHARACTER(LEN=15), INTENT(IN) :: argv,element  
    CHARACTER(LEN=80)             :: cbuffer
  
!------------------------------------------------------------------------------
! 1. Write file containing material IDs
!------------------------------------------------------------------------------
  
    OPEN(14,FILE=argv(1:nlen)//'.bin.ensi.MATID',STATUS="REPLACE",            &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

    cbuffer = "Alya Ensight Gold --- Scalar per-element variable file"
    WRITE(14) cbuffer
    cbuffer = "part" ;  WRITE(14)
    WRITE(14) cbuffer
    WRITE(14) int(1,kind=c_int)
    cbuffer = "coordinates" ; WRITE(14)
  
    SELECT CASE(element)
      CASE('triangle')
        SELECT CASE(nod) 
          CASE(3)
            cbuffer = "tria3" ; WRITE(14) cbuffer
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('quadrilateral')
        SELECT CASE(nod) 
          CASE(4)
            cbuffer = "quad4" ; WRITE(14) cbuffer
          CASE(8)
            cbuffer = "quad8" ; WRITE(14) cbuffer
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('hexahedron')
        SELECT CASE(nod) 
          CASE(8)
            cbuffer = "hexa8"   ; WRITE(14) cbuffer
          CASE(20)
            cbuffer = "hexa20"  ; WRITE(14) cbuffer
          CASE DEFAULT
            WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE('tetrahedron')
        SELECT CASE(nod)
          CASE(4)
            WRITE(14,'(A)') "tetra4"
          CASE DEFAULT
          WRITE(14,'(A)') "# Element type not recognised"
        END SELECT
      CASE DEFAULT
        WRITE(14,'(A)')   "# Element type not recognised"
    END SELECT
   
    WRITE(14) real(etype(:),kind=c_float) 
    !PRINT *, etype
    CLOSE(14)
  
    RETURN
  
  END SUBROUTINE MESH_ENSI_MATID_BIN
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE MESH_ENSI_NDBND_BIN(argv,nf,nlen,nod,solid)

   !/****f* input/mesh_ensi_ndbnd_bin
   !*  NAME
   !*    SUBROUTINE: mesh_ensi_ndbnd_bin
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi_ndbnd_bin(argv,nf,nlen,nod,solid)
   !*  FUNCTION
   !*    This subroutine outputs a file of restrained nodes in the C binary 
   !*    version of the Ensight gold format. Models in this format can be 
   !*    viewed in ParaView.
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*    nod              : number of nodes per element
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*
   !*    Scalar logicals
   !*    solid            : type of analysis solid if .true. fluid if .false.
   !*
   !*    Dynamic scalar arrays
   !*    nf               : nodal freedom matrix
   !*  OUTPUTS
   !*    File: <job_name>.bin.ensi.NDBND 
   !*  AUTHOR
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*
   !*/

    USE, INTRINSIC :: ISO_C_BINDING
    
    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER,   INTENT(IN)         :: nlen,nod
    INTEGER,   INTENT(IN)         :: nf(:,:)
    INTEGER                       :: i,nfe,nn,ndim
    CHARACTER(LEN=15), INTENT(IN) :: argv  
    CHARACTER(LEN=80)             :: cbuffer
    LOGICAL, INTENT(IN)           :: solid
    
  !------------------------------------------------------------------------------
  ! 1. Initialisation
  !------------------------------------------------------------------------------
  
    ndim = UBOUND(nf,1)-1  
    nn   = UBOUND(nf,2)
  
  !------------------------------------------------------------------------------
  ! 2. Write boundary conditions. Encoded using formula: 4z + 2y + 1x
  !
  !    110 = 1   010 = 2   100 = 3   011 = 4   101 = 5   001 = 6   000 = 7
  !-----------------------------------------------------------------------------
  
    IF(solid) THEN

      OPEN(15,FILE=argv(1:nlen)//'.bin.ensi.NDBND',STATUS="REPLACE",           &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

      cbuffer = "Alya Ensight Gold --- Scalar per-node variable file"
      WRITE(15)
      cbuffer = "part"         ; WRITE(15)
      WRITE(15) int(1,kind=c_int)
      cbuffer = "coordinates"  ; WRITE(15)

      IF(ndim==3) THEN
        DO i=1,nod 
          nfe=0
          IF(nf(1,i)==0) nfe=nfe+1
          IF(nf(2,i)==0) nfe=nfe+2
          IF(nf(3,i)==0) nfe=nfe+4
          WRITE(15) int(nfe,kind=c_int)
        END DO
      ELSE IF(ndim==2) THEN
        DO i=1,nn
          nfe=0
          IF(nf(1,i)==0) nfe=nfe+1
          IF(nf(2,i)==0) nfe=nfe+2
          WRITE(15) int(nfe,kind=c_int)
        END DO
      ELSE
        PRINT *, "Wrong number of dimensions in mesh_ensi"
      END IF   
    END IF
  
	
    CLOSE(15)
  
    RETURN
  
  END SUBROUTINE MESH_ENSI_NDBND_BIN

 SUBROUTINE MESH_ENSI_NDLDS_BIN(argv,nlen,nf,loads)

   !/****f* input/mesh_ensi_ndlds_bin
   !*  NAME
   !*    SUBROUTINE: mesh_ensi_ndlds_bin
   !*  SYNOPSIS
   !*    Usage:      CALL mesh_ensi_ndlds_bin(argv,nlen,nf,loads)
   !*  FUNCTION
   !*    This subroutine outputs a file of loads in the C binary version of the 
   !*    Ensight gold format. Models in this format can be viewed in ParaView.
   !*  INPUTS
   !*    Scalar integers
   !*    nlen             : number of characters in data file base name
   !*
   !*    Scalar characters
   !*    argv             : holds data file base name
   !*
   !*    Dynamic scalar arrays
   !*    nf               : nodal freedom matrix
   !* 
   !*    Dynamic real arrays
   !*	 oldlds           : initial loads vector
   !*
   !*  OUTPUTS
   !*  AUTHOR
   !*    L. Margetts
   !*  COPYRIGHT
   !*    (c) University of Manchester 2004-2014
   !******
   !*  Place remarks that should not be included in the documentation here.
   !*
   !*/

    USE, INTRINSIC :: ISO_C_BINDING
    
    IMPLICIT none
  
    INTEGER,PARAMETER             :: iwp=SELECTED_REAL_KIND(15)
    INTEGER,   INTENT(IN)         :: nlen
    INTEGER,   INTENT(IN)         :: nf(:,:)
    INTEGER                       :: i,j
    REAL(iwp), INTENT(IN)         :: loads(:)
    CHARACTER(LEN=15), INTENT(IN) :: argv
    CHARACTER(LEN=80)             :: cbuffer
    
  !-----------------------------------------------------------------------------
  ! 1. Write loaded nodes
  !-----------------------------------------------------------------------------
  
    OPEN(16,FILE=argv(1:nlen)//'.bin.ensi.NDLDS',STATUS="REPLACE",             &
                 FORM="UNFORMATTED", ACTION="WRITE", ACCESS="STREAM")

    cbuffer = "Alya Ensight Gold --- Vector per-node variable file"
    WRITE(16)
    cbuffer = "part"        ; WRITE(16)
    WRITE(16) int(1,kind=c_int)
    cbuffer = "coordinates" ; WRITE(16)

    DO j=1,UBOUND(nf,1)
      DO i=1, UBOUND(nf,2)
        WRITE(16) real(loads(nf(j,i)),kind=c_float)
      END DO
    END DO
    CLOSE(16)
  
    RETURN
  
  END SUBROUTINE MESH_ENSI_NDLDS_BIN



end program cylinder
