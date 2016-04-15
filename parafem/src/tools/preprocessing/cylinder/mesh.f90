!$Id: mesh.f90 31 2016-01-24 21:47:27Z mexas $
! Copyright Anton Shterenlikht, The University of Bristol, UK

! module with mesh routines

module mesh
implicit none
private
public :: node, brickel, circle, cyldmp, bound, bcmerge, springel,     &
  cyldmp_err

type node
  integer :: number
  double precision :: coord(3)
end type node

contains

!***************************************************************

subroutine brickel( fname , el1, el2 , el3, g_num)

! dump element definitions to file in abaqus standard
! Start from 1. Assume the model is brick, with el1, el2, el3 elements
! along 1,2,3.

character( len=* ), intent(in) :: fname
integer, intent(in) :: el1, el2, el3




integer :: funit, errstat, i, j, k , elnum, nd1, nd2, nd3, nd12,       &
 full_lay, full_row, fullt, node1, node2, node3, node4, node5,         &
 node6, node7, node8

!Variable for MESH_ENSI_GEO_BIN subroutine
!New variable nodestore - this variable is used to store the nodal order
!generated in this subroutine by node1, node2 to node 8
integer :: nodestore(8)
integer, intent(inout) :: g_num(8,el3)


open( newunit=funit , file=trim(fname) , status="replace" ,            &
      access="sequential", form="formatted" , iostat=errstat ) 
if ( errstat .ne. 0 ) then
  stop "ERROR: brickel: open"
end if

! number of nodes along 1, 2, 3 and on 12 plane
nd1 = el1 + 1
nd2 = el2 + 1
nd3 = el3 + 1
nd12 = nd1 * nd2

elnum = 1

do k = 1 , el3
do j = 1 , el2
do i = 1 , el1
  full_lay = (k-1) * nd12 ! nodes in full layers
  full_row = (j-1) * nd1  ! nodes in full rows
  fullt = full_lay + full_row

  node1 = fullt + i
  node2 = node1 + 1
  node4 = node1 + nd1
  node3 = node4 + 1 
  node5 = node1 + nd12
  node6 = node5 + 1
  node8 = node5 + nd1
  node7 = node8 + 1 
  
  nodestore(1) = node1
  nodestore(2) = node2
  nodestore(4) = node4
  nodestore(3) = node3
  nodestore(5) = node5
  nodestore(6) = node6
  nodestore(8) = node8
  nodestore(7) = node7

  !Storage of nodal order for MESH_ENSI_GEO_BIN subroutine
  g_num(:,elnum) = nodestore

  write ( funit , "(8(i7,',')i7)")  &
       elnum, node1, node2, node3, node4, node5, node6, node7, node8
  elnum = elnum+1

end do
end do
end do



close( funit, iostat=errstat )
if ( errstat .ne. 0 ) then
  stop "ERROR: brickel: close"
end if

end subroutine brickel

!***************************************************************

subroutine circle( nel , array )

! mesh a unit circle via a conformal mapping of a square

use jacobif

integer, intent(in) :: nel
type( node ), intent(inout), allocatable :: array( : )

double precision, parameter :: zero = 0.0d0, half = 0.5d0,             &
 one = 1.0d0, size = 2.0d0

integer :: i, j, nelp1, nn
double precision :: x, y, dx, dy, xci, yci

! check that array is allocated
if ( .not. allocated( array ) )                                        &
  stop "ERROR: circle/mesh: array not allocated"

   dx = size / nel
   dy = size / nel
nelp1 = nel + 1

nn = 1     ! node number
 y = - one

do j = 1, nelp1
  x = - one
  do i = 1, nelp1
    call sq2ci( x , y , xci , yci )
    array( (j-1) * nelp1 + i ) = node( nn , (/ xci , yci , zero /) )
    nn = nn + 1
     x = x + dx 
  end do
  y = y + dy
end do

end subroutine circle

!***************************************************************

subroutine cyldmp( fname, n , rad, dist, array )

! - replicate circle nodes "n" times along 3, increasing node numbers
!   accordingly.
!   The total number of layers of nodes along 3 is n+1 !
! - increase the radius to "rad" by multiplying both x and y coord.
!   by "rad". The unit circle is assumed on input.
! - Each layer is "dist" further than the previous along 3. 
! - dump to file "fname".
! The circle is passed in as "array".

character( len=* ), intent(in) :: fname
integer, intent(in) :: n
double precision, intent(in) :: rad, dist
type( node ), intent(in), allocatable :: array(:)

integer :: arrlen, errstat, nn, funit, i, j
double precision :: coord1, coord2, coord3

arrlen = size( array )

open( newunit=funit , file=trim(fname) , status="replace" , &
      access="sequential", form="formatted" , iostat=errstat ) 
if ( errstat .ne. 0 ) then
  stop "ERROR: cirep: open"
end if

do j = 0 , n
do i = 1 , arrlen
      nn = array(i)%number + j * arrlen
  coord1 = array(i)%coord(1) * rad
  coord2 = array(i)%coord(2) * rad
  coord3 = array(i)%coord(3) + j * dist
  write ( funit ,"(i6,3(',',es25.16))") nn, coord1, coord2, coord3
  
end do
end do

close( funit, iostat=errstat )
if ( errstat .ne. 0 ) then
  stop "ERROR: cirep: close"
end if

end subroutine cyldmp

!***************************************************************

subroutine cyldmp_err( fname, nbot, nmid, ntop, rad, dist, err, array, g_coord_1, nels_1 )

! - replicate circle nodes along 3. First "nbot" times a perfect
!   circle, then "nmid" times with a random error of +- "err" in
!   both x and y coordinates. Finally "top" times a perfect circle.
!   The node numbers are increasing accordingly.
!   The total number of layers of nodes along 3 is nbot+nmid+ntop+1 !
! - increase the radius to "rad" by multiplying both x and y coord.
!   by "rad". The unit circle is assumed on input.
! - Each layer is "dist" further than the previous along 3. 
! - dump to file "fname".
! The circle is passed in as "array".

character( len=* ), intent(in) :: fname
integer, intent(in) :: nbot, nmid, ntop
double precision, intent(in) :: rad, dist, err
type( node ), intent(in), allocatable :: array(:)

double precision, parameter :: one = 1.0d0, two = 2.0d0

integer :: arrlen, errstat, nn, funit, i, j
double precision :: coord1, coord2, coord3, rnderr(2)


!New Variables
integer, intent(out)   :: nels_1

double precision :: storagecoord(3)
double precision, Allocatable, intent(inout) :: g_coord_1(:,:)
INTEGER  :: findex
INTEGER  :: lengthg_coord
INTEGER :: addition

INTEGER :: formula1, formula2, formula3

storagecoord = 0
findex = 1

arrlen = size( array )

open( newunit=funit , file=trim(fname) , status="replace" , &
      access="sequential", form="formatted" , iostat=errstat ) 
if ( errstat .ne. 0 ) then
  stop "ERROR: cirep: open"
end if

!Calculation of the length for the array g_coord_1
lengthg_coord = (((nbot+1)*arrlen)+(((nbot+nmid+1)-(nbot+1))*arrlen)+(((nbot+nmid+ntop+1)-(nbot+nmid+1))*arrlen))

ALLOCATE(g_coord_1(3,lengthg_coord))

! bottom nodes with no error
do j = 0 , nbot

do i = 1 , arrlen
      nn = array(i)%number + j * arrlen
  coord1 = array(i)%coord(1) * rad
  coord2 = array(i)%coord(2) * rad
  coord3 = array(i)%coord(3) + j * dist
  write ( funit ,"(i6,3(',',es25.16))") nn, coord1, coord2, coord3

  storagecoord(1) = coord1
  storagecoord(2) = coord2
  storagecoord(3) = coord3

  g_coord_1(1,findex) = coord1
  g_coord_1(2,findex) = coord2
  g_coord_1(3,findex) = coord3

  findex = findex + 1

end do
end do
  
! middle nodes with error
do j = nbot+1 , nbot+nmid
do i = 1 , arrlen
      nn = array(i)%number + j * arrlen

  ! random error in coordnates
  call random_number( rnderr )   ! [  0 .. 1 )
  rnderr = rnderr * two - one    ! [ -1 .. 1 )
  rnderr = rnderr * err          ! [ -err .. err )

  coord1 = array(i)%coord(1) * rad * ( one + rnderr(1) )
  coord2 = array(i)%coord(2) * rad * ( one + rnderr(2) )
  coord3 = array(i)%coord(3) + j * dist
  write ( funit ,"(i6,3(',',es25.16))") nn, coord1, coord2, coord3

  storagecoord(1) = coord1
  storagecoord(2) = coord2
  storagecoord(3) = coord3

  g_coord_1(1,findex) = coord1
  g_coord_1(2,findex) = coord2
  g_coord_1(3,findex) = coord3

  findex = findex + 1  
  
end do
end do

! top nodes with no error
do j = nbot+nmid+1 , nbot+nmid+ntop
do i = 1 , arrlen
      nn = array(i)%number + j * arrlen
  coord1 = array(i)%coord(1) * rad
  coord2 = array(i)%coord(2) * rad
  coord3 = array(i)%coord(3) + j * dist
  write ( funit ,"(i6,3(',',es25.16))") nn, coord1, coord2, coord3

  storagecoord(1) = coord1
  storagecoord(2) = coord2
  storagecoord(3) = coord3

  g_coord_1(1,findex) = coord1
  g_coord_1(2,findex) = coord2
  g_coord_1(3,findex) = coord3
  
  findex = findex + 1
  
end do
end do
 
nels_1 = nn

close( funit, iostat=errstat )
if ( errstat .ne. 0 ) then
  stop "ERROR: cirep: close"
end if

end subroutine cyldmp_err

!PRINT *, g_coord_1

!***************************************************************

subroutine bound( fread, fwrite, nset, start, end, rad, counter )

! Read node coord. from "fread", select those which
! are on the surface of the cylinder of radius "rad"
! and between locations "start" and "end" along axis 3.
! Write those node labels into file "fwrite", and assign
! them to node set "nset".

character( len=* ), intent(in) :: fread, fwrite, nset
double precision, intent(in) :: start, end, rad
integer, intent(out) :: counter

double precision, parameter :: tol = 1.0d-5

integer :: label, fr, fw, errstat
double precision :: x1, x2, x3, calc_rad

open( newunit=fr, file=fread, status="old", &
      form="formatted", access="sequential", iostat=errstat )
if ( errstat .ne. 0 ) stop "ERROR: bound: open fread"
open( newunit=fw, file=fwrite, status="replace", &
      form="formatted", access="sequential", iostat=errstat )
if ( errstat .ne. 0 ) stop "ERROR: bound: open fwrite"

! write the command
write( fw, "(a,a)" ) "*nset, nset = ", nset

! write the node labels
counter = 0
do
  read( fr, *, iostat=errstat ) label, x1, x2, x3
  if ( errstat .ne. 0 ) exit
  calc_rad = sqrt( x1**2 + x2**2 )
  if ( abs( calc_rad - rad) .lt. tol .and. &
         x3 .ge. start .and. x3 .le. end ) then
    write( fw, * ) label
    counter = counter + 1
  end if
end do

close( fr, iostat=errstat )
if ( errstat .ne. 0 ) stop "ERROR: bound: close fread"
close( fw, iostat=errstat )
if ( errstat .ne. 0 ) stop "ERROR: bound: close fwrite"

end subroutine bound

!***************************************************************

subroutine bcmerge( file1, file2, outfile )

! Dump date from file1 and file2 into outfile and
! delete file1 and file2

character( len=* ), intent(in) :: file1, file2, outfile

integer, parameter :: nval = 8 ! values per line 
integer :: f1, f2, fout, errstat, label, counter, i, infile
character( len=80 ) :: line

open( newunit=fout, file=outfile, status="replace", &
      form="formatted", access="sequential", iostat=errstat )
if ( errstat .ne. 0 ) stop "ERROR: bound: open outfile"

open( newunit=f1, file=file1, status="old", &
      form="formatted", access="sequential", iostat=errstat )
if ( errstat .ne. 0 ) stop "ERROR: bound: open file1"
open( newunit=f2, file=file2, status="old", &
      form="formatted", access="sequential", iostat=errstat )
if ( errstat .ne. 0 ) stop "ERROR: bound: open file2"

filedo: do i = 1, 2
  counter = 0

  if ( i .eq. 1 ) then
    infile=f1
  else
    infile=f2
    write( fout, * )
  end if

  readdo: do
    if ( counter .eq. 0 ) then
      read( infile, '(a80)', iostat=errstat ) line
      if ( errstat .ne. 0 ) exit readdo
      write( fout,  '(a80)' ) line
    else
      read( infile, * , iostat=errstat ) label
      if ( errstat .ne. 0 ) exit readdo
 
      if ( counter .eq. 1 ) then
        ! first value in file
        write( fout, '(i10)', advance="no" ) label
      else if ( counter .gt. 1 .and. mod( counter-1, nval ) .eq. 0 ) then
        ! last value on the line
        write( fout , '(a)' ) ","
        write( fout , '(i10)', advance="no" ) label
      else
        ! normal value
        write( fout, '(a,i10)', advance="no" ) ",", label
      end if

    end if
    counter = counter + 1

  end do readdo

end do filedo

! close and delete the input files
close( f1, iostat=errstat, status="delete" )
 if ( errstat .ne. 0 ) stop "ERROR: bcmerge: close f1"
close( f2, iostat=errstat, status="delete" )
 if ( errstat .ne. 0 ) stop "ERROR: bcmerge: close f2"

! close and keep the outfile
close( fout, iostat=errstat )
 if ( errstat .ne. 0 ) stop "ERROR: bcmerge: close fout"

end subroutine bcmerge

!***************************************************************

subroutine springel( nfile, springfile, springel1 )

! Read node labels from "elfile" and write a spring element,
! spring1, for each node into "springfile". The first element
! number is "springel1"

character( len=* ), intent(in) :: nfile, springfile
integer, intent(in) :: springel1

integer :: fnode, fel, label, elno, errstat
character( len=80 ) :: header

! open node and element files
open( newunit=fnode, file=trim(nfile), status="old", &
      form="formatted", access="sequential", iostat=errstat )
if ( errstat .ne. 0 ) stop "ERROR: springel: open nfile"

open( newunit=fel, file=trim(springfile), status="replace", &
      form="formatted", access="sequential", iostat=errstat )
if ( errstat .ne. 0 ) stop "ERROR: springel: open springfile"

! read and write headers
read( fnode, * ) header

! read nodes and write elements
elno = springel1
do
  read( fnode, * , iostat=errstat ) label
  if ( errstat .ne. 0 ) exit
  write( fel, * ) elno, ",", label
  elno = elno + 1
end do

close( fnode, iostat=errstat )
 if ( errstat .ne. 0 ) stop "ERROR: springel: close nfile"
close( fel, iostat=errstat )
 if ( errstat .ne. 0 ) stop "ERROR: springel: close springfile"

end subroutine springel

end module mesh
