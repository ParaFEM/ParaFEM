!Program KUBCgen - takes parafem format *.d, and *.dat, and produces 6 test cases for testing Kinematic Uniform Boundary Conditions

!WARNING - there is a lot of repetitve code in here, mainly because I'm uncomfortable passing around dynamic arrays.
!          There's also room for decreasing memory overheads by running xcomp/ycomp etc seperately and allocating/deallocating arrays accordingly.
!          (search for "***TO IMPROVE MEMORY USAGE"...)
!          
!          Also, as of 30/8/11, this is written to use the currently undocumented *.nset and *.fix file formats.
!
!
!I will eventually write a small bash script to sort files into folders.  At the moment this software just dumps out 6*3 additional files, and requires the user to duplicate the original *.d and file six times

!Author : Peter L. Falkingham (peter.falkingham@manchester.ac.uk)

Program KUBCgen

IMPLICIT NONE

!Declarations
	TYPE :: node              !node info
		integer :: nodeID
		real :: x
		real :: y
		real :: z
		!integer :: flag
	End Type node
	
  CHARACTER(LEN=50) :: jobname
	CHARACTER(LEN=50) :: filename

!generic program variables
  INTEGER :: x, i, j, k, unique, temp
  INTEGER :: fixed_freedomsxc, fixed_freedomsyc, fixed_freedomszc
  INTEGER :: fixed_freedomsxs, fixed_freedomsys, fixed_freedomszs
  REAL :: motion, maxx, maxy, maxz, minx, miny, minz         !motion is the distance nodes will be moving.  this is hardcoded here to 10% (0.1)
  REAL :: xcorr, ycorr, zcorr  !these are for correcting the box to 0,0,0 
  INTEGER, DIMENSION(:), ALLOCATABLE :: concat_nset
  TYPE(node), DIMENSION(:), ALLOCATABLE :: sorted_nset
  TYPE(node), DIMENSION(:), ALLOCATABLE :: out_xcomp
  TYPE(node), DIMENSION(:), ALLOCATABLE :: out_ycomp
  TYPE(node), DIMENSION(:), ALLOCATABLE :: out_zcomp
  TYPE(node), DIMENSION(:), ALLOCATABLE :: out_xshear
  TYPE(node), DIMENSION(:), ALLOCATABLE :: out_yshear
  TYPE(node), DIMENSION(:), ALLOCATABLE :: out_zshear

!*.d values
	CHARACTER(LEN=50) :: dummy
	TYPE(node), DIMENSION(:), ALLOCATABLE :: node_list


!*.dat values:
	CHARACTER(LEN=50) :: element
	INTEGER :: mesh, nels, nn, nr, nip, nod, loaded_nodes, fixed_freedoms
	REAL :: e, v, tol, limit, partitioner
	
!*.nset values
  INTEGER :: nsets
  INTEGER, DIMENSION(6) :: nnset
  CHARACTER(LEN=50) :: nset_type
  INTEGER, DIMENSION(:), ALLOCATABLE :: xMax, xMin, yMax, yMin, zMax, zMin !
	
	
!Program proper	

 call GETARG(1,jobname)

!open and read *.dat

  filename = jobname(1:INDEX(jobname, " ")-1) // ".dat"

 OPEN(10, file=filename)
 
 READ(10,*) element
 READ(10,*) mesh
 READ(10,*) partitioner
 READ(10,*) nels, nn, nr, nip, nod, loaded_nodes, fixed_freedoms
 READ(10,*) e, v, tol, limit

 CLOSE(10)
 
!Allocate node array (type node, indlude 'done' flag - all set to '0')

 ALLOCATE(node_list(nn))

!open *.d and read nodes

 filename = jobname(1:INDEX(jobname, " ")-1) // ".d"
 
 OPEN(10, file=filename)
 
 READ(10,*)dummy
 READ(10,*)dummy
 
 READ(10,*)node_list
 
 CLOSE(10)
 
 maxx=maxval(node_list(:)%x)
 minx=minval(node_list(:)%x)
 maxy=maxval(node_list(:)%y)
 miny=minval(node_list(:)%y)
 maxz=maxval(node_list(:)%z)
 minz=minval(node_list(:)%z)
 
 if(minx.ne.0)then
   xcorr = -1*minx
 else 
   xcorr = 0
 end if
 
 if(miny.ne.0)then
  ycorr = -1*miny
 else
  ycorr = 0
 end if
 
 if(minz.ne.0)then
  zcorr = -1*minz
 else
  zcorr = 0
 end if
 
 
 
!now read nsets  *****MUST ADD IN READ DUMMY BEFORE AND AFTER NNSET*****
 filename = jobname(1:INDEX(jobname, " ")-1) // ".nset"

 OPEN(10,file=filename)
 
 READ(10,*)dummy
 if(dummy.ne."KUBC")then
  WRITE(*,*)"WARNING, KUBC FILE-TYPE NOT FOUND.  Error..."
 end if
 READ(10,*)nsets
 READ(10,*)dummy, nnset(1), dummy
 ALLOCATE(xmin(nnset(1)))
 READ(10,*)xmin
 READ(10,*)dummy, nnset(2), dummy
 ALLOCATE(xmax(nnset(2)))
 READ(10,*)xmax
 READ(10,*)dummy, nnset(3), dummy
 ALLOCATE(ymin(nnset(3)))
 READ(10,*)ymin
 READ(10,*)dummy, nnset(4), dummy
 ALLOCATE(ymax(nnset(4)))
 READ(10,*)ymax
 READ(10,*)dummy, nnset(5), dummy
 ALLOCATE(zmin(nnset(5)))
 READ(10,*)zmin
 READ(10,*)dummy, nnset(6), dummy
 ALLOCATE(zmax(nnset(6)))
 READ(10,*)zmax 
 CLOSE(10)

!concatenate those arrays and then search for duplicates, Then sort into unique_nset
!Concat
 ALLOCATE(concat_nset(nnset(1)+nnset(2)+nnset(3)+nnset(4)+nnset(5)+nnset(6)))
 concat_nset(1:nnset(1)) = xmin
 concat_nset(nnset(1)+1:nnset(1)+nnset(2)) = xmax
 concat_nset(nnset(1)+nnset(2)+1:nnset(1)+nnset(2)+nnset(3)) = ymin
 concat_nset(nnset(1)+nnset(2)+nnset(3)+1:nnset(1)+nnset(2)+nnset(3)+nnset(4))=ymax
 concat_nset(nnset(1)+nnset(2)+nnset(3)+nnset(4)+1:nnset(1)+nnset(2)+nnset(3)+nnset(4)+nnset(5))=zmin 
 concat_nset(nnset(1)+nnset(2)+nnset(3)+nnset(4)+nnset(5)+1:nnset(1)+nnset(2)+nnset(3)+nnset(4)+nnset(5)+nnset(6))=zmax 

  DEALLOCATE(xmin)
  DEALLOCATE(xmax)
  DEALLOCATE(ymin)
  DEALLOCATE(ymax)
  DEALLOCATE(zmin)
  DEALLOCATE(zmax)

!remove duplicates 
unique = size(concat_nset)
 do x=1,size(concat_nset)
  do i=x+1, size(concat_nset)
    if(concat_nset(i).eq.concat_nset(x))then
      concat_nset(x)=0
      unique = unique-1
    end if
  end do
 end do

!now sort
 do x=1,size(concat_nset)
  do i=x+1, size(concat_nset)
    if(concat_nset(x).gt.concat_nset(i))then
      temp = concat_nset(x)
      concat_nset(x) = concat_nset(i)
      concat_nset(i) = temp
    end if
  end do
 end do

 
 ALLOCATE(sorted_nset(unique))
 i = size(concat_nset)-unique+1
 do x=1,unique
  sorted_nset(x)%nodeid = concat_nset(i)
  i = i+1
 end do

 DEALLOCATE(concat_nset)

!We now have a sorted array containing only unique node ids on the boundary.  Next, get xyz coords from *.d node list
 do x=1, unique
  sorted_nset(x)%x = node_list(sorted_nset(x)%nodeid)%x
  sorted_nset(x)%y = node_list(sorted_nset(x)%nodeid)%y
  sorted_nset(x)%z = node_list(sorted_nset(x)%nodeid)%z
 end do


motion = 0.1
!Now produce an array that is infact displacements (xcomp)    ***TO IMPROVE MEMORY USAGE, MAKE ONE ARRAY, WRITE OUT FILES, THEN DEALLOCATE AND REMAKE
 ALLOCATE(out_xcomp(unique))
 out_xcomp = sorted_nset
 do x=1,unique
 out_xcomp(x)%x = (out_xcomp(x)%x+xcorr)*motion
 end do
 out_xcomp%y = 0
 out_xcomp%z = 0
 
!ycomp
 ALLOCATE(out_ycomp(unique))
 out_ycomp = sorted_nset
 out_ycomp%x = 0
 do x=1, unique
 out_ycomp(x)%y = (out_ycomp(x)%y+ycorr)*motion
 end do
 out_ycomp%z = 0

!zcomp
 ALLOCATE(out_zcomp(unique))
 out_zcomp = sorted_nset
 out_zcomp%x = 0
 out_zcomp%y = 0
 do x=1, unique
 out_zcomp(x)%z = (out_zcomp(x)%z+zcorr)*motion
 end do
 
!Xshear
 ALLOCATE(out_xshear(unique))
 out_xshear = sorted_nset
 do x=1, unique
 out_xshear(x)%x = (out_xshear(x)%x+xcorr)*motion
 out_xshear(x)%y = (out_xshear(x)%y+ycorr)*motion
 end do
 out_xshear%z = 0
 
!Yshear
 ALLOCATE(out_yshear(unique))
 out_yshear = sorted_nset
 do x=1, unique
 out_yshear(x)%x = (out_yshear(x)%x+xcorr)*motion
 out_yshear(x)%z = (out_yshear(x)%z+zcorr)*motion
 end do 
 out_yshear%y = 0
 
!Zshear
 ALLOCATE(out_zshear(unique))
 out_zshear = sorted_nset
 out_zshear%x = 0
 do x=1, unique
 out_zshear(x)%y = (out_zshear(x)%y+ycorr)*motion
 out_zshear(x)%z = (out_zshear(x)%z+zcorr)*motion
 end do
   

!now to output 6 bnd files, 6 fix files, and 6 dat files
!Xcompression
!Xcomp bnd
OPEN(11,file='xcomp.bnd')

do x=1, unique
  WRITE(11,'(I12)',advance='no')out_xcomp(x)%nodeid
  if(out_xcomp(x)%x.eq.0)then
    WRITE(11,'(A)',advance='no')' 0 '
  else
    WRITE(11,'(A)',advance='no')' 1 '
  end if
  if(out_xcomp(x)%y.eq.0)then
    WRITE(11,'(A)',advance='no')' 0 '
  else
    WRITE(11,'(A)',advance='no')' 1 '
  end if  
  if(out_xcomp(x)%z.eq.0)then
    WRITE(11,'(A)',advance='yes')' 0 '
  else
    WRITE(11,'(A)',advance='yes')' 1 '
  end if
end do
  CLOSE(11)  
  
!Ycomp bnd
OPEN(11,file='ycomp.bnd')
do x=1, unique
  WRITE(11,'(I12)',advance='no')out_ycomp(x)%nodeid
  if(out_ycomp(x)%x.eq.0)then
    WRITE(11,'(A)',advance='no')' 0 '
  else
    WRITE(11,'(A)',advance='no')' 1 '
  end if
  if(out_ycomp(x)%y.eq.0)then
    WRITE(11,'(A)',advance='no')' 0 '
  else
    WRITE(11,'(A)',advance='no')' 1 '
  end if  
  if(out_ycomp(x)%z.eq.0)then
    WRITE(11,'(A)',advance='yes')' 0 '
  else
    WRITE(11,'(A)',advance='yes')' 1 '
  end if
end do
  CLOSE(11)
!zcomp bnd

OPEN(11,file='zcomp.bnd')
do x=1, unique
  WRITE(11,'(I12)',advance='no')out_zcomp(x)%nodeid
  if(out_zcomp(x)%x.eq.0)then
    WRITE(11,'(A)',advance='no')' 0 '
  else
    WRITE(11,'(A)',advance='no')' 1 '
  end if
  if(out_zcomp(x)%y.eq.0)then
    WRITE(11,'(A)',advance='no')' 0 '
  else
    WRITE(11,'(A)',advance='no')' 1 '
  end if  
  if(out_zcomp(x)%z.eq.0)then
    WRITE(11,'(A)',advance='yes')' 0 '
  else
    WRITE(11,'(A)',advance='yes')' 1 '
  end if
end do
  CLOSE(11)
!xshear
OPEN(11,file='xshear.bnd')
do x=1, unique
  WRITE(11,'(I12)',advance='no')out_xshear(x)%nodeid
  if(out_xshear(x)%x.eq.0)then
    WRITE(11,'(A)',advance='no')' 0 '
  else
    WRITE(11,'(A)',advance='no')' 1 '
  end if
  if(out_xshear(x)%y.eq.0)then
    WRITE(11,'(A)',advance='no')' 0 '
  else
    WRITE(11,'(A)',advance='no')' 1 '
  end if  
  if(out_xshear(x)%z.eq.0)then
    WRITE(11,'(A)',advance='yes')' 0 '
  else
    WRITE(11,'(A)',advance='yes')' 1 '
  end if
end do
  CLOSE(11)
  
  !yshear
  OPEN(11,file='yshear.bnd')
do x=1, unique
  WRITE(11,'(I12)',advance='no')out_yshear(x)%nodeid
  if(out_yshear(x)%x.eq.0)then
    WRITE(11,'(A)',advance='no')' 0 '
  else
    WRITE(11,'(A)',advance='no')' 1 '
  end if
  if(out_yshear(x)%y.eq.0)then
    WRITE(11,'(A)',advance='no')' 0 '
  else
    WRITE(11,'(A)',advance='no')' 1 '
  end if  
  if(out_yshear(x)%z.eq.0)then
    WRITE(11,'(A)',advance='yes')' 0 '
  else
    WRITE(11,'(A)',advance='yes')' 1 '
  end if
end do
  CLOSE(11)
  
  
  !zshear
  OPEN(11,file='zshear.bnd')
do x=1, unique
  WRITE(11,'(I12)',advance='no')out_zshear(x)%nodeid
  if(out_zshear(x)%x.eq.0)then
    WRITE(11,'(A)',advance='no')' 0 '
  else
    WRITE(11,'(A)',advance='no')' 1 '
  end if
  if(out_zshear(x)%y.eq.0)then
    WRITE(11,'(A)',advance='no')' 0 '
  else
    WRITE(11,'(A)',advance='no')' 1 '
  end if  
  if(out_zshear(x)%z.eq.0)then
    WRITE(11,'(A)',advance='yes')' 0 '
  else
    WRITE(11,'(A)',advance='yes')' 1 '
  end if
end do
  CLOSE(11)
  
  
!Now output the Fix file
fixed_freedomsxc = 0
fixed_freedomsyc = 0
fixed_freedomszc = 0
fixed_freedomsxs = 0
fixed_freedomsys = 0
fixed_freedomszs = 0



OPEN(11,file='xcomp.fix')

do x=1, unique
  if(out_xcomp(x)%x+out_xcomp(x)%y+out_xcomp(x)%z.ne.0)then
  WRITE(11,'(I12)',advance='no')out_xcomp(x)%nodeid
  if(out_xcomp(x)%x.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 1 ', out_xcomp(x)%x
  end if
  if(out_xcomp(x)%y.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 2 ', out_xcomp(x)%y
  end if  
  if(out_xcomp(x)%z.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 3 ', out_xcomp(x)%z
  end if
  WRITE(11,*)' '
  fixed_freedomsxc = fixed_freedomsxc+1
  end if
end do
  CLOSE(11)  
  
  
OPEN(11,file='ycomp.fix')

do x=1, unique
  if(out_ycomp(x)%x+out_ycomp(x)%y+out_ycomp(x)%z.ne.0)then
  WRITE(11,'(I12)',advance='no')out_ycomp(x)%nodeid
  if(out_ycomp(x)%x.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 1 ', out_ycomp(x)%x
  end if
  if(out_ycomp(x)%y.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 2 ', out_ycomp(x)%y
  end if  
  if(out_ycomp(x)%z.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 3 ', out_ycomp(x)%z
  end if
  WRITE(11,*)' '
  fixed_freedomsyc = fixed_freedomsyc+1
  end if
end do
  CLOSE(11)  
  
  
  OPEN(11,file='zcomp.fix')

do x=1, unique
  if(out_zcomp(x)%x+out_zcomp(x)%y+out_zcomp(x)%z.ne.0)then
  WRITE(11,'(I12)',advance='no')out_zcomp(x)%nodeid
  if(out_zcomp(x)%x.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 1 ', out_zcomp(x)%x
  end if
  if(out_zcomp(x)%y.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 2 ', out_zcomp(x)%y
  end if  
  if(out_zcomp(x)%z.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 3 ', out_zcomp(x)%z
  end if
  WRITE(11,*)' '
  fixed_freedomszc = fixed_freedomszc+1
  end if
end do
  CLOSE(11)  
  
  
  OPEN(11,file='xshear.fix')

do x=1, unique
  if(out_xshear(x)%x+out_xshear(x)%y+out_xshear(x)%z.ne.0)then
  WRITE(11,'(I12)',advance='no')out_xshear(x)%nodeid
  if(out_xshear(x)%x.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 1 ', out_xshear(x)%x
  end if
  if(out_xshear(x)%y.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 2 ', out_xshear(x)%y
  end if  
  if(out_xshear(x)%z.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 3 ', out_xshear(x)%z
  end if
  WRITE(11,*)' '
    fixed_freedomsxs = fixed_freedomsxs+1
  end if
end do
  CLOSE(11)  
  
  
  
  OPEN(11,file='yshear.fix')

do x=1, unique
  if(out_yshear(x)%x+out_yshear(x)%y+out_yshear(x)%z.ne.0)then
  WRITE(11,'(I12)',advance='no')out_yshear(x)%nodeid
  if(out_yshear(x)%x.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 1 ', out_yshear(x)%x
  end if
  if(out_yshear(x)%y.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 2 ', out_yshear(x)%y
  end if  
  if(out_yshear(x)%z.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 3 ', out_yshear(x)%z
  end if
    WRITE(11,*)' '
        fixed_freedomsys = fixed_freedomsys+1
  end if
end do
  CLOSE(11)  
  
  
  
  OPEN(11,file='zshear.fix')

do x=1, unique
  if(out_zshear(x)%x+out_zshear(x)%y+out_zshear(x)%z.ne.0)then
  WRITE(11,'(I12)',advance='no')out_zshear(x)%nodeid
  if(out_zshear(x)%x.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 1 ', out_zshear(x)%x
  end if
  if(out_zshear(x)%y.ne.0)then
    WRITE(11,'(A, E16.8)',advance='no')' 2 ', out_zshear(x)%y
  end if  
  if(out_zshear(x)%z.ne.0)then
    WRITE(11,'(A, E16.8)',advance='yes')' 3 ', out_zshear(x)%z
  end if
    WRITE(11,*)' '
        fixed_freedomszs = fixed_freedomszs+1
  end if
end do
  CLOSE(11)  
  
  
  !now output the dat files - Xcomp
  
  OPEN(11,file='xcomp.dat')
  WRITE(11,*) element
  WRITE(11,*) mesh
  WRITE(11,*) partitioner
  WRITE(11,*) nels, nn, unique, nip, nod, '0', fixed_freedomsxc
  WRITE(11,*) e, v, tol, limit
  CLOSE(11)
  
  OPEN(11,file='ycomp.dat')
  WRITE(11,*) element
  WRITE(11,*) mesh
  WRITE(11,*) partitioner
  WRITE(11,*) nels, nn, unique, nip, nod, '0', fixed_freedomsyc
  WRITE(11,*) e, v, tol, limit
  CLOSE(11)
  
  OPEN(11,file='zcomp.dat')
  WRITE(11,*) element
  WRITE(11,*) mesh
  WRITE(11,*) partitioner
  WRITE(11,*) nels, nn, unique, nip, nod, '0', fixed_freedomszc
  WRITE(11,*) e, v, tol, limit
  CLOSE(11)    
  
  OPEN(11,file='xshear.dat')
  WRITE(11,*) element
  WRITE(11,*) mesh
  WRITE(11,*) partitioner
  WRITE(11,*) nels, nn, unique, nip, nod, '0', fixed_freedomsxs
  WRITE(11,*) e, v, tol, limit
  CLOSE(11)
    
  OPEN(11,file='yshear.dat')
  WRITE(11,*) element
  WRITE(11,*) mesh
  WRITE(11,*) partitioner
  WRITE(11,*) nels, nn, unique, nip, nod, '0', fixed_freedomsys
  WRITE(11,*) e, v, tol, limit
  CLOSE(11)
  
  OPEN(11,file='zshear.dat')
  WRITE(11,*) element
  WRITE(11,*) mesh
  WRITE(11,*) partitioner
  WRITE(11,*) nels, nn, unique, nip, nod, '0', fixed_freedomszs
  WRITE(11,*) e, v, tol, limit
  CLOSE(11)  
  
end program KUBCgen

