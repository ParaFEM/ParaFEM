!quick fortran program to convert palaeofem (smith and griffiths) output to abaqus style.
!This swaps node order thus:
!Mapping:
!S&G    1  2  3  4 5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20

!should be 1  7 19 13 3  5 17 15  8 12 20  9  4 11 16 10  2  6 18 14

!16/3/9 - version 2, for outputting useful data (maximum depth at each layer, maximum uplift at eaach layer, ability to select layer to view)

Program p2a

IMPLICIT NONE
 DOUBLE PRECISION, ALLOCATABLE :: nodes(:,:)
 INTEGER, ALLOCATABLE :: elements(:,:)
 CHARACTER(LEN=50) :: line
 CHARACTER(LEN=50) :: line2
 CHARACTER(LEN=50) :: line3
 CHARACTER(len=50) :: job_name
 CHARACTER(LEN=50) :: fname_dat
 CHARACTER(LEN=50) :: fname_d
 CHARACTER(LEN=50) :: new_d
 CHARACTER(LEN=50) :: buffer
 INTEGER :: numelements, numnodes, x, layer, layerdir, footmat
 
 call GETARG(1, job_name)

 call GETARG(2, buffer)
 READ(BUFFER,*) layer

 fname_dat = job_name(1:INDEX(job_name, " ")-1) // ".dat"
 fname_d = job_name(1:INDEX(job_name, " ")-1) // ".d"
 new_d = job_name(1:INDEX(job_name, " ")-1) // buffer 
 new_d = new_d(1:INDEX(new_d, " ")-1)  //"new.d"


 call GETARG(3, buffer)
 READ(BUFFER,*) layerdir

 OPEN(10,file=fname_dat)
 READ(10,*) numelements,numnodes
 CLOSE(10)

 ALLOCATE(nodes(4,numnodes))
 ALLOCATE(elements(25,numelements))


 OPEN(11,file=fname_d)
 READ(11,*) line,line2,nodes,line3,elements
 CLOSE(11)


10 format(F10.0,2x,a8,2x,a8,2x,a8)

 OPEN(12,file=new_d)
 WRITE(12,'(A18)')'*THREE_DIMENSIONAL'
 WRITE(12,'(A6)')'*NODES'
 do x=1,numnodes
   WRITE(12,*)x,nodes(2,x),nodes(3,x),nodes(4,x)
 end do
 WRITE(12,'(A9)')'*ELEMENTS'

20 format(i10,2x,i3,2x,i3,2x,i3,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,&
2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i10,2x,i3)

footmat = elements(25,numelements)

WRITE(*,*)'Writing elements of ', layer
if(layerdir.gt.0)then
WRITE(*,*)'and below'
else if(layerdir.lt.0)then
WRITE(*,*)'and above'
end if

if(layer.eq.0)then
 do x=1,numelements
  WRITE(12,*)elements(1,x), elements(2,x), elements(3,x), elements(4,x), &
             elements(5,x), elements(11,x), elements(23,x), elements(17,x), &
             elements(7,x), elements(9,x), elements(21,x), elements(19,x), &
             elements(12,x), elements(16,x), elements(24,x), elements(13,x), & 
             elements(8,x), elements(15,x), elements(20,x), elements(14,x), & 
             elements(6,x), elements(10,x), elements(22,x), elements(18,x), &
             elements(25,x)
 end do
else
   if(layerdir.eq.0)then

   do x=1,numelements
   if(elements(25,x).eq.layer)then
   WRITE(12,*)elements(1,x), elements(2,x), elements(3,x), elements(4,x), &
             elements(5,x), elements(11,x), elements(23,x), elements(17,x), &
             elements(7,x), elements(9,x), elements(21,x), elements(19,x), &
             elements(12,x), elements(16,x), elements(24,x), elements(13,x), & 
             elements(8,x), elements(15,x), elements(20,x), elements(14,x), & 
             elements(6,x), elements(10,x), elements(22,x), elements(18,x), &
             elements(25,x)
   end if
 end do
 
 else if(layerdir.lt.0)then
    do x=1,numelements
   if(elements(25,x).le.layer)then
    if(elements(25,x).ne.footmat)then
   WRITE(12,*)elements(1,x), elements(2,x), elements(3,x), elements(4,x), &
             elements(5,x), elements(11,x), elements(23,x), elements(17,x), &
             elements(7,x), elements(9,x), elements(21,x), elements(19,x), &
             elements(12,x), elements(16,x), elements(24,x), elements(13,x), & 
             elements(8,x), elements(15,x), elements(20,x), elements(14,x), & 
             elements(6,x), elements(10,x), elements(22,x), elements(18,x), &
             elements(25,x)
     end if
   end if
 end do
 
  else if(layerdir.gt.0)then
    do x=1,numelements
   if(elements(25,x).ge.layer)then
       if(elements(25,x).ne.footmat)then
   WRITE(12,*)elements(1,x), elements(2,x), elements(3,x), elements(4,x), &
             elements(5,x), elements(11,x), elements(23,x), elements(17,x), &
             elements(7,x), elements(9,x), elements(21,x), elements(19,x), &
             elements(12,x), elements(16,x), elements(24,x), elements(13,x), & 
             elements(8,x), elements(15,x), elements(20,x), elements(14,x), & 
             elements(6,x), elements(10,x), elements(22,x), elements(18,x), &
             elements(25,x)
             end if
   end if
 end do
 end if
end if


end program

