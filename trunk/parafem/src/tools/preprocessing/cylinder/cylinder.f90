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
use input

implicit none

!Variables for subroutine MESH_ENSI_CASE  
INTEGER          :: nlen,nstep,npri
REAL(8)          :: dtim
LOGICAL          :: solid

!Variables for MESH_ENSI_MATID_BIN, MESH_ENSI_NDLDS_BIN

double precision, ALLOCATABLE  :: loads(:)
double precision, ALLOCATABLE  :: etype_1(:)
INTEGER          :: nod_1
INTEGER          :: nn
INTEGER, ALLOCATABLE :: nf_1(:,:)

double precision, ALLOCATABLE  :: val(:,:)
INTEGER, ALLOCATABLE           :: node_1(:) 

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
                 int( fixed_len ), rad, dist, cerr, circ, g_coord_1, nn )

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
ALLOCATE(nf_1(4,nn))
ALLOCATE(loads(nn))

ALLOCATE(val(4,nn))


!Material properties assignation
do i = 1, nel3*nel*nel
	etype_1(i) = 10000
end do

!Automatic 3-2-1 boundary conditions 

do i = 1, nn


	nf_1(1,i) = i
	nf_1(2,i) = 1
	nf_1(3,i) = 1
	nf_1(4,i) = 1

end do

do i = 1, (nel+1)*(nel+1)

	nf_1(1,i) = i
	nf_1(2,i) = 1
	nf_1(3,i) = 1
	nf_1(4,i) = 0
        
end do

nf_1(2,1) = 0
nf_1(3,1) = 0
nf_1(2, (nel+1)+(nel+1)) = 0



!End automatic 3-2-1 boundary conditions 

ALLOCATE (node_1(1))

node_1(1) = 10 

do i = 1, nn

        val(1,i) = i
	val(2,i) = 0
	val(3,i) = 1
	val(4,i) = 0

end do

!Loads



!End loads

!CALL ParaFEM subroutines to produce ensi.case files
CALL MESH_ENSI_CASE(filesuff,nlen,nstep,npri,dtim,solid)
CALL MESH_ENSI_GEO_BIN(filesuff,nlen,g_coord_1,g_num_1,element1)
CALL MESH_ENSI_MATID_BIN(filesuff,nlen,nod_1,element1,etype_1)
CALL MESH_ENSI_NDBND_BIN(filesuff,nf_1,nlen,nod_1,solid)
CALL MESH_ENSI_NDLDS_BIN(filesuff,nlen,nn,val,node_1)

end program cylinder
