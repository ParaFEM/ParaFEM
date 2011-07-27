!PROGRAM MESHGENv6

!   ----------------------------------------------------------
!   MeshGen6.F90    22/3/11
!   Re-write of MeshGen for generating input data to PalaeoFEM
!	Generates a cuboid mesh with 2.5D foot/indenter on the
!	surface.
!   No error checking implemented yet
!	Changes: Comments standardised for ParaFEM
!   Update 22/3/11 : Now uses Abaqus/FEAR node ordering
!   ----------------------------------------------------------

  IMPLICIT NONE

!-----------------------!
!1.  Declarations:      !
!-----------------------!
  
  INTEGER :: denseResX, denseResZ
  INTEGER :: resolutionX, resolutionZ, layers, i, x, y, z, x1, x2, x3,&
             y1, y2, y3, z1, z2, z3, random, ce, le, n1, n1b, n2, n2b
  INTEGER :: num_nodes, num_elements, nodes_per_el, xdr2, zdr, stepStages, loadType, old_num_els, num_load_els
  INTEGER :: num_loads, rigid_nodes, c, old_num_nodes, lightRes, num_nodes_x, num_nodes_y, num_nodes_z, footcount, footlayers
  INTEGER :: allonum, allonod, biggestDenseRes, smooth, num_incs, shallow_layers
  DOUBLE PRECISION :: cu, e, v, footThickness, elsize, mass, layerThickness, vari, loadperel, twelth, third, scaleFact, dist, &
  						mult, dummy, scaleFactY
  REAL, ALLOCATABLE :: loadedEls(:,:)
  REAL, ALLOCATABLE :: increments(:)
  CHARACTER(LEN=50),ALLOCATABLE   ::  stepFiles(:)
  CHARACTER(LEN=10) :: p3
  CHARACTER(len=50) :: job_name
  CHARACTER(LEN=50) :: ldsname
  TYPE :: node
     integer :: id
     real :: x
     real :: y
     real :: z
     integer :: freedomx
     integer :: freedomy
     integer :: freedomz
     real :: loadx
     real :: loady
     real :: loadz
  END TYPE node

  TYPE :: element8			!This is a stub for 8-node hexahedral elements 
     integer :: id
     type (node) :: node1
     type (node) :: node2
     type (node) :: node3
     type (node) :: node4
     type (node) :: node5
     type (node) :: node6
     type (node) :: node7
     type (node) :: node8
  END TYPE element8

  TYPE :: element20
     integer :: id
     type (node) :: node1
     type (node) :: node2
     type (node) :: node3
     type (node) :: node4
     type (node) :: node5
     type (node) :: node6
     type (node) :: node7
     type (node) :: node8
     type (node) :: node9
     type (node) :: node10
     type (node) :: node11
     type (node) :: node12
     type (node) :: node13
     type (node) :: node14
     type (node) :: node15
     type (node) :: node16
     type (node) :: node17
     type (node) :: node18
     type (node) :: node19
     type (node) :: node20
  END TYPE element20

  TYPE(element8), DIMENSION(:,:,:), ALLOCATABLE :: el8			!Stub
  TYPE(element20), DIMENSION(:,:,:), ALLOCATABLE :: el20
  TYPE(node), DIMENSION(:,:,:), ALLOCATABLE :: nodes


 call GETARG(1, job_name)

!-----------------------!
!2. Read in input:		!
!-----------------------!
  OPEN(14,file=job_name)
  READ(14,*)  elSize, denseResX, denseResZ, lightRes, scaleFact, nodes_per_el, &
              layers, layerThickness, scaleFactY, shallow_layers, &
              footlayers, footThickness, &
              random, cu, e, v, vari
  READ(14,*)  loadType
  READ(14,*)stepStages, mass  
  ALLOCATE(stepFiles(stepStages))
  READ(14,*)stepFiles
  READ(14,*)smooth
  READ(14,*)num_incs
  ALLOCATE(increments(num_incs))
  READ(14,*)increments
  CLOSE(14)

!-----------------------------------------!
!2.1 Read in which elements are loaded	  !
!-----------------------------------------!
  OPEN(12,file=stepFiles(1)) !this will be dynamic later.
  READ(12,*)p3
  READ(12,*)denseResX,denseResZ,dummy
  ALLOCATE(loadedEls(denseResX*3,denseResZ))
  READ(12,*)loadedEls
  CLOSE(12)
  
  WRITE(*,*)'read denseRes X and Z from ppm: ', denseResX, denseResZ


  resolutionX = denseResX + lightRes*2
  resolutionZ = denseResZ + lightRes*2


!------------------------------------------------!
!3.	Generate mesh				 !
!------------------------------------------------!
  
  allonum = -1*(footlayers-1)
  allonod = (-2*footlayers)+1

  if(nodes_per_el.eq.20)then  
    num_nodes_x = 4*lightRes + 2*denseResX + 1
    num_nodes_y = 2*layers + 3
    num_nodes_z = 4*lightRes + 2*denseResZ + 1
  else 
    num_nodes_x = 2*lightRes + denseResX + 1
    num_nodes_y = layers + 2
    num_nodes_z = 2*lightRes + denseResZ + 1
  end if


  ALLOCATE(nodes(0:num_nodes_x,allonod:num_nodes_y, 0:num_nodes_z))

  if(nodes_per_el.eq.20)then
    ALLOCATE(el20(resolutionX, allonum:layers, resolutionZ))
  else
    ALLOCATE(el8(resolutionX, allonum:layers, resolutionZ))
  end if

  nodes = node(0,0,0,0,1,1,1,0,0,0) !just initialise 


  ![****Fill nodes****!]

  if(nodes_per_el.eq.20)then
    i = 1
    do y=1, layers*2+1, 2 
       do z=1, resolutionZ*2+1, 2
          do x=1, resolutionX*2+1
             nodes(x,y,z)%id = i
             nodes(x,y,z)%x = x-1
             nodes(x,y,z)%y = -1*(y-1)
             nodes(x,y,z)%z = z-1
             i=i+1
          end do
       end do
    end do
  !  [now create half way nodes, on horizontal...]
    do y=1, layers*2+1, 2
      do z=2, resolutionZ*2+1, 2
        do x=1, resolutionX*2+1, 2
           nodes(x,y,z)%id = i
           nodes(x,y,z)%x = x-1
           nodes(x,y,z)%y = -1*(y-1)
           nodes(x,y,z)%z = z-1
           i=i+1
        end do
      end do
    end do
  ![... and vertical.]
    do y=2, layers*2+1, 2
      do z=1, resolutionZ*2+1, 2
        do x=1, resolutionX*2+1, 2
           nodes(x,y,z)%id = i
           nodes(x,y,z)%x = x-1
           nodes(x,y,z)%y = -1*(y-1)
           nodes(x,y,z)%z = z-1
           i=i+1
        end do
      end do
    end do
  
  else 
    i=1
    do y=1, layers+1, 1
      do z=1, resolutionZ+1, 1
        do x=1, resolutionX+1, 1
          nodes(x,y,z)%id = i
          nodes(x,y,z)%x = x-1
          nodes(x,y,z)%y = -1*(y-1)
          nodes(x,y,z)%z = z-1
	  i = i+1
	end do
      end do
    end do
  end if

  num_nodes = i
  old_num_nodes = num_nodes	
  
  ![fix nodes on edges - 0 = fixed, 1 = free]
if(nodeS_per_el.eq.20)then
  nodes(1,:,:)%freedomx = 0
  nodes(resolutionX*2+1,:,:)%freedomx = 0
  nodes(:,layers*2+1,:)%freedomy = 0	!
  nodes(:,layers*2+1,:)%freedomx = 0    ! >Base elements
  nodes(:,layers*2+1,:)%freedomz = 0	!
  nodes(:,:,1)%freedomz = 0
  nodes(:,:,resolutionZ*2+1)%freedomz = 0
else
  nodes(1,:,:)%freedomx = 0
  nodes(resolutionX+1,:,:)%freedomx = 0
  nodes(:,layers+1,:)%freedomy = 0	!
  nodes(:,layers+1,:)%freedomx = 0      ! >Base elements
  nodes(:,layers+1,:)%freedomz = 0	!
  nodes(:,:,1)%freedomz = 0
  nodes(:,:,resolutionZ+1)%freedomz = 0
end if

  ![Fill Elements]

if(nodes_per_el.eq.20)then
  i=1
  do y=1, layers
     do z = 1, resolutionZ
        do x = 1, resolutionX
           x2 = x*2
           x1 = x2-1
           x3 = x2+1
           z2 = z*2
           z1 = z2-1
           z3 = z2+1
           y2 = y*2
           y1 = y2-1
           y3 = y2+1
           el20(x,y,z)%id = i
           el20(x,y,z)%node1 = nodes(x1,y3,z1)
           el20(x,y,z)%node2 = nodes(x1,y2,z1)
           el20(x,y,z)%node3 = nodes(x1,y1,z1)
           el20(x,y,z)%node4 = nodes(x1,y1,z2)
           el20(x,y,z)%node5 = nodes(x1,y1,z3)
           el20(x,y,z)%node6 = nodes(x1,y2,z3)
           el20(x,y,z)%node7 = nodes(x1,y3,z3)
           el20(x,y,z)%node8 = nodes(x1,y3,z2)
           el20(x,y,z)%node9 = nodes(x2,y3,z1)
           el20(x,y,z)%node10 = nodes(x2,y1,z1)
           el20(x,y,z)%node11 = nodes(x2,y1,z3)
           el20(x,y,z)%node12 = nodes(x2,y3,z3)
           el20(x,y,z)%node13 = nodes(x3,y3,z1)
           el20(x,y,z)%node14 = nodes(x3,y2,z1)
           el20(x,y,z)%node15 = nodes(x3,y1,z1)
           el20(x,y,z)%node16 = nodes(x3,y1,z2)
           el20(x,y,z)%node17 = nodes(x3,y1,z3)
           el20(x,y,z)%node18 = nodes(x3,y2,z3)
           el20(x,y,z)%node19 = nodes(x3,y3,z3)
           el20(x,y,z)%node20 = nodes(x3,y3,z2)
           i = i+1
        end do
     end do
  end do
else
  i=1
  do y=1, layers
     do z = 1, resolutionZ
        do x = 1, resolutionX
	   x2 = x+1
           x1 = x
	   z2 = z+1
           z1 = z
	   y2 = y+1
           y1 = y
           el8(x,y,z)%id = i
           el8(x,y,z)%node1 = nodes(x1,y2,z1)
           el8(x,y,z)%node2 = nodes(x1,y1,z1)
           el8(x,y,z)%node3 = nodes(x1,y1,z2)
           el8(x,y,z)%node4 = nodes(x1,y2,z2)
           el8(x,y,z)%node5 = nodes(x2,y2,z1)
           el8(x,y,z)%node6 = nodes(x2,y1,z1)
           el8(x,y,z)%node7 = nodes(x2,y1,z2)
           el8(x,y,z)%node8 = nodes(x2,y2,z2)
           i = i+1
        end do
     end do
  end do  
end if
  num_elements = i
  old_num_els = num_elements
  


!------------------------------------------------------!
!4. Create new elements (elements forming the indenter)!
!Only available for 20-node elements			!
!------------------------------------------------------!
if(footThickness.ne.0)then

if(nodes_per_el.eq.20)then
  
do footcount=1, footlayers*2, 2
    x2=1
    do x=2, denseResX*3,3
       do z=1, denseResZ
          if(loadedEls(x-1,z)+loadedEls(x,z)+loadedEls(x+1,z).ne.0)then  !this will have to become magnitude, as currently, 1+ -1 = 0
             xdr2 = x2+lightRes
             zdr = z + lightRes
             ce = footcount-1
             ce = ce/2
             ce = -1*ce
             le = ce+1
             n1 = -1*footcount
             n1b = footcount+1
             n2 = n1+1
             n2b = footcount
             el20(xdr2,ce,zdr)%id = num_elements
             num_elements = num_elements+1
             el20(xdr2,ce,zdr)%node1 = el20(xdr2,le,zdr)%node3
             el20(xdr2,ce,zdr)%node7 = el20(xdr2,le,zdr)%node5
             el20(xdr2,ce,zdr)%node8 = el20(xdr2,le,zdr)%node4
             el20(xdr2,ce,zdr)%node9 = el20(xdr2,le,zdr)%node10
             el20(xdr2,ce,zdr)%node12 = el20(xdr2,le,zdr)%node11
             el20(xdr2,ce,zdr)%node13 = el20(xdr2,le,zdr)%node15
             el20(xdr2,ce,zdr)%node19 = el20(xdr2,le,zdr)%node17
             el20(xdr2,ce,zdr)%node20 = el20(xdr2,le,zdr)%node16
             !create new nodes, layer 0

             if(nodes(xdr2*2-1, n1, zdr*2-1)%id.eq.0)then
                nodes(xdr2*2-1,n1,zdr*2-1)%id = num_nodes
                nodes(xdr2*2-1,n1,zdr*2-1)%x = nodes(xdr2*2-1,1,zdr*2-1)%x
                nodes(xdr2*2-1,n1,zdr*2-1)%y = n1b
                nodes(xdr2*2-1,n1,zdr*2-1)%z = nodes(xdr2*2-1,1,zdr*2-1)%z
                num_nodes = num_nodes+1
             end if
             el20(xdr2,ce,zdr)%node3 = nodes(xdr2*2-1,n1,zdr*2-1)				

             if(nodes(xdr2*2-1, n1, zdr*2)%id.eq.0)then
                nodes(xdr2*2-1,n1,zdr*2)%id = num_nodes
                nodes(xdr2*2-1,n1,zdr*2)%x = nodes(xdr2*2-1,1,zdr*2)%x
                nodes(xdr2*2-1,n1,zdr*2)%y = n1b
                nodes(xdr2*2-1,n1,zdr*2)%z = nodes(xdr2*2-1,1,zdr*2)%z
                num_nodes = num_nodes+1
             end if
             el20(xdr2,ce,zdr)%node4 = nodes(xdr2*2-1,n1,zdr*2)		

             if(nodes(xdr2*2-1,n1,  zdr*2+1)%id.eq.0)then
                nodes(xdr2*2-1,n1,zdr*2+1)%id = num_nodes
                nodes(xdr2*2-1,n1,zdr*2+1)%x = nodes(xdr2*2-1,1,zdr*2+1)%x
                nodes(xdr2*2-1,n1,zdr*2+1)%y = n1b
                nodes(xdr2*2-1,n1,zdr*2+1)%z = nodes(xdr2*2-1,1,zdr*2+1)%z
                num_nodes = num_nodes+1
             end if
             el20(xdr2,ce,zdr)%node5 = nodes(xdr2*2-1,n1,zdr*2+1)							

             if(nodes(xdr2*2,n1,  zdr*2-1)%id.eq.0)then
                nodes(xdr2*2,n1,zdr*2-1)%id = num_nodes
                nodes(xdr2*2,n1,zdr*2-1)%x = nodes(xdr2*2,1,zdr*2-1)%x
                nodes(xdr2*2,n1,zdr*2-1)%y = n1b
                nodes(xdr2*2,n1,zdr*2-1)%z = nodes(xdr2*2,1,zdr*2-1)%z
                num_nodes = num_nodes+1
             end if
             el20(xdr2,ce,zdr)%node10 = nodes(xdr2*2,n1,zdr*2-1)							

             if(nodes(xdr2*2,n1,  zdr*2+1)%id.eq.0)then
                nodes(xdr2*2,n1,zdr*2+1)%id = num_nodes
                nodes(xdr2*2,n1,zdr*2+1)%x = nodes(xdr2*2,1,zdr*2+1)%x
                nodes(xdr2*2,n1,zdr*2+1)%y = n1b
                nodes(xdr2*2,n1,zdr*2+1)%z = nodes(xdr2*2,1,zdr*2+1)%z
                num_nodes = num_nodes+1
             end if
             el20(xdr2,ce,zdr)%node11 = nodes(xdr2*2,n1,zdr*2+1)					

             if(nodes(xdr2*2+1,n1,  zdr*2-1)%id.eq.0)then
                nodes(xdr2*2+1,n1,zdr*2-1)%id = num_nodes
                nodes(xdr2*2+1,n1,zdr*2-1)%x = nodes(xdr2*2+1,1,zdr*2-1)%x
                nodes(xdr2*2+1,n1,zdr*2-1)%y = n1b
                nodes(xdr2*2+1,n1,zdr*2-1)%z = nodes(xdr2*2+1,1,zdr*2-1)%z
                num_nodes = num_nodes+1
             end if
             el20(xdr2,ce,zdr)%node15 = nodes(xdr2*2+1,n1,zdr*2-1)		

             if(nodes(xdr2*2+1,n1,  zdr*2)%id.eq.0)then
                nodes(xdr2*2+1,n1,zdr*2)%id = num_nodes
                nodes(xdr2*2+1,n1,zdr*2)%x = nodes(xdr2*2+1,1,zdr*2)%x
                nodes(xdr2*2+1,n1,zdr*2)%y = n1b
                nodes(xdr2*2+1,n1,zdr*2)%z = nodes(xdr2*2+1,1,zdr*2)%z
                num_nodes = num_nodes+1
             end if
             el20(xdr2,ce,zdr)%node16 = nodes(xdr2*2+1,n1,zdr*2)		

             if(nodes(xdr2*2+1,n1,  zdr*2+1)%id.eq.0)then
                nodes(xdr2*2+1,n1,zdr*2+1)%id = num_nodes
                nodes(xdr2*2+1,n1,zdr*2+1)%x = nodes(xdr2*2+1,1,zdr*2+1)%x
                nodes(xdr2*2+1,n1,zdr*2+1)%y = n1b
                nodes(xdr2*2+1,n1,zdr*2+1)%z = nodes(xdr2*2+1,1,zdr*2+1)%z
                num_nodes = num_nodes+1
             end if
             el20(xdr2,ce,zdr)%node17 = nodes(xdr2*2+1,n1,zdr*2+1)

             if(nodes(xdr2*2-1, n2, zdr*2-1)%id.eq.0)then
                nodes(xdr2*2-1,n2,zdr*2-1)%id = num_nodes
                nodes(xdr2*2-1,n2,zdr*2-1)%x = nodes(xdr2*2-1,1,zdr*2-1)%x
                nodes(xdr2*2-1,n2,zdr*2-1)%y = n2b
                nodes(xdr2*2-1,n2,zdr*2-1)%z = nodes(xdr2*2-1,1,zdr*2-1)%z
                num_nodes = num_nodes+1
             end if
             el20(xdr2,ce,zdr)%node2 = nodes(xdr2*2-1,n2,zdr*2-1)

             if(nodes(xdr2*2-1, n2, zdr*2+1)%id.eq.0)then
                nodes(xdr2*2-1,n2,zdr*2+1)%id = num_nodes
                nodes(xdr2*2-1,n2,zdr*2+1)%x = nodes(xdr2*2-1,1,zdr*2+1)%x
                nodes(xdr2*2-1,n2,zdr*2+1)%y = n2b
                nodes(xdr2*2-1,n2,zdr*2+1)%z = nodes(xdr2*2-1,1,zdr*2+1)%z
                num_nodes = num_nodes+1
             end if
             el20(xdr2,ce,zdr)%node6 = nodes(xdr2*2-1,n2,zdr*2+1)		

             if(nodes(xdr2*2+1, n2, zdr*2-1)%id.eq.0)then
                nodes(xdr2*2+1,n2,zdr*2-1)%id = num_nodes
                nodes(xdr2*2+1,n2,zdr*2-1)%x = nodes(xdr2*2+1,1,zdr*2-1)%x
                nodes(xdr2*2+1,n2,zdr*2-1)%y = n2b
                nodes(xdr2*2+1,n2,zdr*2-1)%z = nodes(xdr2*2+1,1,zdr*2-1)%z
                num_nodes = num_nodes+1
             end if
             el20(xdr2,ce,zdr)%node14 = nodes(xdr2*2+1,n2,zdr*2-1)							

             if(nodes(xdr2*2+1, n2, zdr*2+1)%id.eq.0)then
                nodes(xdr2*2+1,n2,zdr*2+1)%id = num_nodes
                nodes(xdr2*2+1,n2,zdr*2+1)%x = nodes(xdr2*2+1,1,zdr*2+1)%x
                nodes(xdr2*2+1,n2,zdr*2+1)%y = n2b
                nodes(xdr2*2+1,n2,zdr*2+1)%z = nodes(xdr2*2+1,1,zdr*2+1)%z
                num_nodes = num_nodes+1
             end if
             el20(xdr2,ce,zdr)%node18 = nodes(xdr2*2+1,n2,zdr*2+1)		
          end if
       end do
       x2 = x2+1
     end do
 
end do

else
  WRITE(*,*)'Generating an indenter is only available for 20-node elements &
 &at this time. No elements have been generated above the soil.'
end if     
   end if

  num_load_els = (num_elements - old_num_els)/footlayers

!-----------------------------------------------_-------!
!5.  Set all nodes to correct scale 			!
!-------------------------------------------------------!
if(nodes_per_el.eq.20)then
  do x=lightRes*2, 1, -2
  	dist = (nodes(x+3,1,1)%x-nodes(x+1,1,1)%x)*scaleFact
  	nodes(x,:,:)%x = nodes(x+1,:,:)%x-(0.5*dist)
  	nodes(x-1,:,:)%x = nodes(x+1,:,:)%x-dist
  end do
   
  do z=lightRes*2, 1, -2
  	dist = (nodes(1,1,z+3)%z-nodes(1,1,z+1)%z)*scaleFact
  	nodes(:,:,z)%z = nodes(:,:,z+1)%z-(0.5*dist)
  	nodes(:,:,z-1)%z = nodes(:,:,z+1)%z-dist
  end do
  
  do x=(denseResX+lightRes)*2+2, resolutionX*2+1, 2
  	dist = (nodes(x-1,1,1)%x-nodes(x-3,1,1)%x)*scaleFact
  	nodes(x,:,:)%x = nodes(x-1,:,:)%x+(0.5*dist)
  	nodes(x+1,:,:)%x = nodes(x-1,:,:)%x+dist
  end do
  
  do z=(lightRes+denseResZ)*2+2, resolutionZ*2+1, 2
  	dist = (nodes(1,1,z-1)%z-nodes(1,1,z-3)%z)*scaleFact
  	nodes(:,:,z)%z = nodes(:,:,z-1)%z+(0.5*dist)
  	nodes(:,:,z+1)%z = nodes(:,:,z-1)%z+dist
  end do
  
![y scaling]

  if(shallow_layers<0)then

    if(denseResX.ge.denseResZ)then
      biggestDenseRes = denseResX
    else
      biggestDenseRes = denseResZ
    end if

    do y=biggestDenseRes/2+2, layers*2+1, 2 !insert here to stop foot expanding
  	dist = (nodes(1,y-1,1)%y-nodes(1,y-3,1)%y)*scaleFactY	
  	nodes(:,y,:)%y = nodes(:,y-1,:)%y+(0.5*dist)
  	nodes(:,y+1,:)%y = nodes(:,y-1,:)%y+dist
    end do

  else

    do y=shallow_layers*2+2, layers*2+1, 2 !insert here to stop foot expanding
  	dist = (nodes(1,y-1,1)%y-nodes(1,y-3,1)%y)*scaleFactY	
  	nodes(:,y,:)%y = nodes(:,y-1,:)%y+(0.5*dist)
  	nodes(:,y+1,:)%y = nodes(:,y-1,:)%y+dist
    end do

  end if

else
!put scaling for 8-node here
end if

if(nodes_per_el.eq.20)then  
  mult = elSize/2
  nodes(:,:,:)%x = nodes(:,:,:)%x*mult
  nodes(:,:,:)%y = nodes(:,:,:)%y*layerthickness/2
  nodes(:,:,:)%z = nodes(:,:,:)%z*elSize/2

else 
  mult = elSize
  nodes(:,:,:)%x = nodes(:,:,:)%x*mult
  nodes(:,:,:)%y = nodes(:,:,:)%y*layerthickness
  nodes(:,:,:)%z = nodes(:,:,:)%z*elSize
end if


!---------------------------------------------------------------------------------!
!5. Attempt smoothing of corner els.  This is experimental! Use at your own risk! !
!	 (May cause numerical problems)
!	Only works with 20node elements.					   							  !
!---------------------------------------------------------------------------------!
if(nodes_per_el.eq.20)then
if(smooth.eq.1)then
WRITE(*,*)'Applying smoothing algorithm...  '
write(*,*)'***WARNING*** DO NOT USE ON LOW RESOLUTION MESHES'
y = 0
do z=lightRes, lightRes+denseResZ
  do x=lightRes, lightRes+denseResX
  	if(el20(x,y,z)%id.ne.0)then

  	  ![check upper left]
  	  if(el20(x-1,y,z-1)%id.eq.0)then
  	  	if(el20(x-1,y,z)%id.eq.0)then
  	  	  if(el20(x,y,z-1)%id.eq.0)then
  	  	    nodes(x*2-1,:,z*2-1)%x = nodes(x*2-1,:,z*2-1)%x+elSize/4
 			nodes(x*2-1,:,z*2-1)%z = nodes(x*2-1,:,z*2-1)%z+elSize/4

		  end if
		end if
		
		if(el20(x-1,y,z)%id.ne.0)then
  	  	  if(el20(x,y,z-1)%id.ne.0)then
  	  	    !nodes(x*2-1,:,z*2-1)%x = nodes(x*2-1,:,z*2-1)%x-elSize/4
 			!nodes(x*2-1,:,z*2-1)%z = nodes(x*2-1,:,z*2-1)%z-elSize/4

		  end if
		end if
		
 	  end if
 	  
 	  ![check upper right]
 	  if(el20(x+1,y,z-1)%id.eq.0)then
  	  	if(el20(x+1,y,z)%id.eq.0)then
  	  	  if(el20(x,y,z-1)%id.eq.0)then
  	  	    nodes(x*2+1,:,z*2-1)%x = nodes(x*2+1,:,z*2-1)%x-elSize/4
 			nodes(x*2+1,:,z*2-1)%z = nodes(x*2+1,:,z*2-1)%z+elSize/4

		  end if
		end if
		
		if(el20(x+1,y,z)%id.ne.0)then
  	  	  if(el20(x,y,z-1)%id.ne.0)then
  	  	    !nodes(x*2+1,:,z*2-1)%x = nodes(x*2+1,:,z*2-1)%x+elSize/4
 			!nodes(x*2+1,:,z*2-1)%z = nodes(x*2+1,:,z*2-1)%z-elSize/4

		  end if
		end if
		
 	  end if
 	  

  ![check lower left]
 	if(el20(x-1,y,z)%id.eq.0)then
 	  if(el20(x-1,y,z+1)%id.eq.0)then
  	   	  if(el20(x,y,z+1)%id.eq.0)then
  	  	    nodes(x*2-1,:,z*2+1)%x = nodes(x*2-1,:,z*2+1)%x+elSize/4
 			nodes(x*2-1,:,z*2+1)%z = nodes(x*2-1,:,z*2+1)%z-elSize/4

		  end if
		end if	
	  if(el20(x-1,y,z+1)%id.ne.0)then
	  	if(el20(x,y,z+1)%id.ne.0)then
  	  	    nodes(x*2-1,:,z*2+1)%x = nodes(x*2-1,:,z*2+1)%x-elSize/4
 			nodes(x*2-1,:,z*2+1)%z = nodes(x*2-1,:,z*2+1)%z-elSize/4

		end if	  	
		
 	  end if
  
 	 end if


   ![check lower right]
   if(el20(x+1,y,z+1)%id.eq.0)then
  	if(el20(x+1,y,z)%id.eq.0)then
  	  if(el20(x,y,z+1)%id.eq.0)then
  	    nodes(x*2+1,:,z*2+1)%x = nodes(x*2+1,:,z*2+1)%x-elSize/4
 		nodes(x*2+1,:,z*2+1)%z = nodes(x*2+1,:,z*2+1)%z-elSize/4
 	  end if
 	end if
   end if
   
  if(el20(x+1,y,z)%id.eq.0)then
  
    if(el20(x+1,y,z+1)%id.ne.0)then
    if(el20(x,y,z+1)%id.ne.0)then
        nodes(x*2+1,:,z*2+1)%x = nodes(x*2+1,:,z*2+1)%x+elSize/4
 		nodes(x*2+1,:,z*2+1)%z = nodes(x*2+1,:,z*2+1)%z-elSize/4
 	end if
 	end if
  end if
  
  if(el20(x+1,y,z+1)%id.eq.0)then
  
    if(el20(x+1,y,z)%id.ne.0)then
    if(el20(x,y,z+1)%id.ne.0)then
        nodes(x*2+1,:,z*2+1)%x = nodes(x*2+1,:,z*2+1)%x+elSize/4
 		nodes(x*2+1,:,z*2+1)%z = nodes(x*2+1,:,z*2+1)%z+elSize/4
 	end if
 	end if
  end if
  
  if(el20(x,y,z+1)%id.eq.0)then
  
    if(el20(x+1,y,z+1)%id.ne.0)then
    if(el20(x+1,y,z)%id.ne.0)then
        nodes(x*2+1,:,z*2+1)%x = nodes(x*2+1,:,z*2+1)%x-elSize/4
 		nodes(x*2+1,:,z*2+1)%z = nodes(x*2+1,:,z*2+1)%z+elSize/4
 	end if
 	end if
  end if
   

 	 
 	 
 end if
  end do
end do
end if
else
  if(smooth.eq.1)then
  WRITE(*,*)'only 20-node elements can use the scaling effect'
  end if
end if

! ----------------------------------------!
! end experimental section				  !
! ----------------------------------------!



!-------------------------------------------------------!
!6.   Write *.d file and *.bnd:							!
!     Writes to file the nodes,							!
!     elements, and boundary conditions used by ParaFEM !
!-------------------------------------------------------!

  open(10,file='dino.d')
  open(11,file='dino.bnd')
  write(10,'(A18)')'*THREE_DIMENSIONAL'
  write(10,'(A6)')'*NODES'
  rigid_nodes = 0

if(nodes_per_el.eq.20)then
 do y=1,layers*2+1, 2
    do z=1, resolutionZ*2+1, 2
       do x=1, resolutionX*2+1
          if(nodes(x,y,z)%id.ne.0)then
             write(10,*)nodes(x,y,z)%id, nodes(x,y,z)%x,nodes(x,y,z)%y,nodes(x,y,z)%z
             if(nodes(x,y,z)%freedomx+nodes(x,y,z)%freedomy+nodes(x,y,z)%freedomz.ne.3)then
                write(11,*) nodes(x,y,z)%id, nodes(x,y,z)%freedomx,nodes(x,y,z)%freedomy, &
                     nodes(x,y,z)%freedomz
                rigid_nodes = rigid_nodes+1
             end if
          end if
       end do
    end do
 end do
 do y=1, layers*2+1, 2
    do z=2, resolutionZ*2+1, 2
       do x=1, resolutionX*2+1, 2
          if(nodes(x,y,z)%id.ne.0)then
             write(10,*)nodes(x,y,z)%id, nodes(x,y,z)%x,nodes(x,y,z)%y,nodes(x,y,z)%z
             if(nodes(x,y,z)%freedomx+nodes(x,y,z)%freedomy+nodes(x,y,z)%freedomz.ne.3)then
                write(11,*) nodes(x,y,z)%id, nodes(x,y,z)%freedomx,nodes(x,y,z)%freedomy, &
                     nodes(x,y,z)%freedomz
                rigid_nodes = rigid_nodes+1
             end if
          end if
       end do
    end do
 end do

 do y=2,layers*2+1, 2
    do z=1, resolutionZ*2+1, 2
       do x=1, resolutionX*2+1
          if(nodes(x,y,z)%id.ne.0)then
             write(10,*)nodes(x,y,z)%id, nodes(x,y,z)%x,nodes(x,y,z)%y,nodes(x,y,z)%z
             if(nodes(x,y,z)%freedomx+nodes(x,y,z)%freedomy+nodes(x,y,z)%freedomz.ne.3)then
                write(11,*) nodes(x,y,z)%id, nodes(x,y,z)%freedomx,nodes(x,y,z)%freedomy, &
                     nodes(x,y,z)%freedomz
                rigid_nodes = rigid_nodes+1
             end if
          end if
       end do
    end do
 end do
 ![now output indenter nodes.]
 do c=old_num_nodes, num_nodes
 
 do y=allonod,0
    do z = lightRes*2, (lightRes+denseResZ)*2+1
       do x = lightRes*2, (lightRes+denseResX)*2+1
          if(nodes(x,y,z)%id.eq.c)then
             write(10,*)nodes(x,y,z)%id, nodes(x,y,z)%x,nodes(x,y,z)%y,nodes(x,y,z)%z
             if(nodes(x,y,z)%freedomx+nodes(x,y,z)%freedomy+nodes(x,y,z)%freedomz.ne.3)then
                write(11,*) nodes(x,y,z)%id, nodes(x,y,z)%freedomx,nodes(x,y,z)%freedomy, &
                     nodes(x,y,z)%freedomz
                rigid_nodes = rigid_nodes+1
             end if
          end if
       end do
    end do
 end do
 end do
else
 do y=1,layers+1
    do z=1, resolutionZ+1
       do x=1, resolutionX+1
          if(nodes(x,y,z)%id.ne.0)then
             write(10,*)nodes(x,y,z)%id, nodes(x,y,z)%x,nodes(x,y,z)%y,nodes(x,y,z)%z
             if(nodes(x,y,z)%freedomx+nodes(x,y,z)%freedomy+nodes(x,y,z)%freedomz.ne.3)then
                write(11,*) nodes(x,y,z)%id, nodes(x,y,z)%freedomx,nodes(x,y,z)%freedomy, &
                     nodes(x,y,z)%freedomz
                rigid_nodes = rigid_nodes+1
             end if
          end if
       end do
    end do
 end do

 ![now output indenter nodes.]
 do c=old_num_nodes, num_nodes
 
 do y=allonod,0
    do z = lightRes, (lightRes+denseResZ)+1
       do x = lightRes, (lightRes+denseResX)+1
          if(nodes(x,y,z)%id.eq.c)then
             write(10,*)nodes(x,y,z)%id, nodes(x,y,z)%x,nodes(x,y,z)%y,nodes(x,y,z)%z
             if(nodes(x,y,z)%freedomx+nodes(x,y,z)%freedomy+nodes(x,y,z)%freedomz.ne.3)then
                write(11,*) nodes(x,y,z)%id, nodes(x,y,z)%freedomx,nodes(x,y,z)%freedomy, &
                     nodes(x,y,z)%freedomz
                rigid_nodes = rigid_nodes+1
             end if
          end if
       end do
    end do
 end do
 end do
end if

 write(10,'(A9)')'*ELEMENTS'

if(nodes_per_el.eq.20)then
 do y=1,layers
    do z=1, resolutionZ
       do x=1, resolutionX
          if(el20(x,y,z)%id.ne.0)then
             write(10,*) el20(x,y,z)%id, ' 3 ', ' 20 ', ' 1 ', el20(x,y,z)%node1%id, &
                  el20(x,y,z)%node7%id, el20(x,y,z)%node19%id, el20(x,y,z)%node13%id, &
                  el20(x,y,z)%node3%id, el20(x,y,z)%node5%id, el20(x,y,z)%node17%id, &
                  el20(x,y,z)%node15%id, el20(x,y,z)%node8%id, el20(x,y,z)%node12%id, &
                  el20(x,y,z)%node20%id, el20(x,y,z)%node9%id, el20(x,y,z)%node4%id, &
                  el20(x,y,z)%node11%id, el20(x,y,z)%node16%id, el20(x,y,z)%node10%id, &
                  el20(x,y,z)%node2%id, el20(x,y,z)%node6%id, el20(x,y,z)%node18%id, &
                  el20(x,y,z)%node14%id, y
          end if
       end do
    end do
 end do
 ![output indenter elements]
 do ce=0, allonum, -1

 do x=1, resolutionX
    do z=1, resolutionZ
       if(el20(x,ce,z)%id.ne.0)then
          write(10,*) el20(x,ce,z)%id, ' 3 ', ' 20 ', ' 1 ', el20(x,ce,z)%node1%id, &
               el20(x,ce,z)%node7%id, el20(x,ce,z)%node19%id, el20(x,ce,z)%node13%id, &
               el20(x,ce,z)%node3%id, el20(x,ce,z)%node5%id, el20(x,ce,z)%node17%id, &
               el20(x,ce,z)%node15%id, el20(x,ce,z)%node8%id, el20(x,ce,z)%node12%id, &
               el20(x,ce,z)%node20%id, el20(x,ce,z)%node9%id, el20(x,ce,z)%node4%id, &
               el20(x,ce,z)%node11%id, el20(x,ce,z)%node16%id, el20(x,ce,z)%node10%id, &
               el20(x,ce,z)%node2%id, el20(x,ce,z)%node6%id, el20(x,ce,z)%node18%id, &
               el20(x,ce,z)%node14%id, y
       end if
    end do
 end do
 end do

else
 do y=1,layers
    do z=1, resolutionZ
       do x=1, resolutionX
          if(el8(x,y,z)%id.ne.0)then
             write(10,*) el8(x,y,z)%id, ' 3 ', ' 8 ', ' 1 ', el8(x,y,z)%node1%id, &
                  el8(x,y,z)%node2%id, el8(x,y,z)%node3%id, el8(x,y,z)%node4%id, &
                  el8(x,y,z)%node5%id, el8(x,y,z)%node6%id, el8(x,y,z)%node7%id, &
                  el8(x,y,z)%node8%id, y
          end if
       end do
    end do
 end do
 ![output indenter elements]
 do ce=0, allonum, -1

 do x=1, resolutionX
    do z=1, resolutionZ
       if(el8(x,ce,z)%id.ne.0)then
          write(10,*) el8(x,ce,z)%id, ' 3 ', ' 8 ', ' 1 ', el8(x,ce,z)%node1%id, &
               el8(x,ce,z)%node2%id, el8(x,ce,z)%node3%id, el8(x,ce,z)%node4%id, &
               el8(x,ce,z)%node5%id, el8(x,ce,z)%node6%id, el8(x,ce,z)%node7%id, &
               el8(x,ce,z)%node8%id, y
       end if
    end do
 end do
 end do
end if

  close(10)
  close(11)

!---------------------------------------------------------!
!7.  Write *.mat										  !
!     Writes material properties to ParaFEM format (*.mat)!
!---------------------------------------------------------!

  call writemat(footThickness, layers, cu, e, v)

if(loadType.eq.2)then
!-------------------------------------------------!
!8. Write *.lds (currently only vertical)	  !
!   Generates loading conditions in ParaFEM format!
!   Stubs included for non-vertical loading       !
!   Only Available for 20-node elements		  !
!-------------------------------------------------!
if(nodes_per_el.eq.20)then
  loadperel = mass/num_load_els
  third = (-1*0.333333333)*loadperel
  twelth = (0.08333333333)*loadperel

  do z=1, resolutionZ
     do x=1, resolutionX
        if(el20(x,allonum,z)%id.ne.0)then
           !nodes(x*2-1,-1,z*2-1)%loadx
           nodes(x*2-1,allonod,z*2-1)%loady	= nodes(x*2-1,allonod,z*2-1)%loady + twelth
           !nodes(x*2-1,-1,z*2-1)%loadz

           !nodes(x*2,-1,z*2-1)%loadx
           nodes(x*2,allonod,z*2-1)%loady	= nodes(x*2,allonod,z*2-1)%loady + third
           !nodes(x*2,-1,z*2-1)%loadz

           !nodes(x*2+1,-1,z*2-1)%loadx
           nodes(x*2+1,allonod,z*2-1)%loady	= nodes(x*2+1,allonod,z*2-1)%loady + twelth
           !nodes(x*2+1,-1,z*2-1)%loadz

           !nodes(x*2-1,-1,z*2)%loadx
           nodes(x*2-1,allonod,z*2)%loady	= nodes(x*2-1,allonod,z*2)%loady + third
           !nodes(x*2-1,-1,z*2)%loadz

           !nodes(x*2+1,-1,z*2)%loadx
           nodes(x*2+1,allonod,z*2)%loady	= nodes(x*2+1,allonod,z*2)%loady + third
           !nodes(x*2+1,-1,z*2)%loadz

           !nodes(x*2-1,-1,z*2+1)%loadx
           nodes(x*2-1,allonod,z*2+1)%loady	= nodes(x*2-1,allonod,z*2+1)%loady + twelth
           !nodes(x*2-1,-1,z*2+1)%loadz

           !nodes(x*2,-1,z*2+1)%loadx
           nodes(x*2,allonod,z*2+1)%loady	= nodes(x*2,allonod,z*2+1)%loady + third
           !nodes(x*2,-1,z*2+1)%loadz

           !nodes(x*2+1,-1,z*2+1)%loadx
           nodes(x*2+1,allonod,z*2+1)%loady	= nodes(x*2+1,allonod,z*2+1)%loady + twelth
           !nodes(x*2+1,-1,z*2+1)%loadz
        end if
     end do
  end do


  do c=1, num_incs

  WRITE(ldsname, '(A5,I1,A4)') "dino_", c, ".lds"
  
  open(13,file=ldsname)
  num_loads = 0
  do z=1, resolutionZ*2+1
     do x=1, resolutionX*2+1
        if(nodes(x,allonod,z)%id.ne.0)then
          if(nodes(x,allonod,z)%loady.ne.0)then
           WRITE(13,*)nodes(x,allonod,z)%id, increments(c)*nodes(x,allonod,z)%loadx, increments(c)*nodes(x,allonod,z)%loady, &
           				increments(c)*nodes(x,allonod,z)%loadz
           num_loads = num_loads+1
          end if
        end if
     end do
  end do
  close(13)
  end do

else
WRITE(*,*)'load conditions are only generated for 20-node element meshes currently. No load file has been created.'
end if

!-------------------------------------------------!
!8b. Write *.fix file (currently only vertical)	  !
!   Generates fixed conditions in ParaFEM format  !
!   Stubs included for non-vertical loading       !
!   Only Available for 20-node elements		  !
!-------------------------------------------------!
else

if(nodes_per_el.eq.20)then
  loadperel = mass


  do z=1, resolutionZ
     do x=1, resolutionX
        if(el20(x,allonum,z)%id.ne.0)then
           !nodes(x*2-1,-1,z*2-1)%loadx
           nodes(x*2-1,allonod,z*2-1)%loady	= loadperel
           !nodes(x*2-1,-1,z*2-1)%loadz

           !nodes(x*2,-1,z*2-1)%loadx
           nodes(x*2,allonod,z*2-1)%loady	= loadperel
           !nodes(x*2,-1,z*2-1)%loadz

           !nodes(x*2+1,-1,z*2-1)%loadx
           nodes(x*2+1,allonod,z*2-1)%loady	= loadperel
           !nodes(x*2+1,-1,z*2-1)%loadz

           !nodes(x*2-1,-1,z*2)%loadx
           nodes(x*2-1,allonod,z*2)%loady	= loadperel
           !nodes(x*2-1,-1,z*2)%loadz

           !nodes(x*2+1,-1,z*2)%loadx
           nodes(x*2+1,allonod,z*2)%loady	= loadperel
           !nodes(x*2+1,-1,z*2)%loadz

           !nodes(x*2-1,-1,z*2+1)%loadx
           nodes(x*2-1,allonod,z*2+1)%loady	= loadperel
           !nodes(x*2-1,-1,z*2+1)%loadz

           !nodes(x*2,-1,z*2+1)%loadx
           nodes(x*2,allonod,z*2+1)%loady	= loadperel
           !nodes(x*2,-1,z*2+1)%loadz

           !nodes(x*2+1,-1,z*2+1)%loadx
           nodes(x*2+1,allonod,z*2+1)%loady	= loadperel
           !nodes(x*2+1,-1,z*2+1)%loadz
        end if
     end do
  end do


  do c=1, num_incs

  WRITE(ldsname, '(A5,I1,A4)') "dino_", c, ".fix"
  
  open(13,file=ldsname)
  num_loads = 0
  do z=1, resolutionZ*2+1
     do x=1, resolutionX*2+1
        if(nodes(x,allonod,z)%id.ne.0)then
          if(nodes(x,allonod,z)%loady.ne.0)then
           WRITE(13,*)nodes(x,allonod,z)%id, increments(c)*nodes(x,allonod,z)%loadx, increments(c)*nodes(x,allonod,z)%loady, &
           				increments(c)*nodes(x,allonod,z)%loadz
           num_loads = num_loads+1
          end if
        end if
     end do
  end do
  close(13)
  end do

else
WRITE(*,*)'fixed conditions are only generated for 20-node element meshes currently. No fixed file has been created.'
end if
end if

!-----------------!
!9.  Write *.dat  !
!-----------------!
  open(14,file='dino.dat')
  WRITE(14,*) "'hexahedron'"
  WRITE(14,*) '2'
  WRITE(14,*) num_elements-1, num_nodes-1, rigid_nodes, '8'
  WRITE(14,*) '12 3 50 15000'
  WRITE(14,*) '1.e-6 1.e-7 -1.e-6 1.e-5'
  WRITE(14,*) layers, num_incs
if(loadType.eq.2)then
  do c = 1, num_incs
  WRITE(14,*) '0'
  WRITE(14,*) num_loads
  end do
else
  do c = 1, num_incs
  WRITE(14,*) num_loads
  WRITE(14,*) '0'
  end do
end if

!------------------!
!  end of program: !
!------------------!
END PROGRAM



!--------------------------!
!10. Subroutine write *.mat:!
!--------------------------!

subroutine writemat(footThickness, layers, cu,e,v)
	DOUBLE PRECISION :: footThickness, cu, e, v
	integer :: layers, i

	open(10,file='dino.mat')
	if(footThickness.ne.0)then
		layers = layers+1
	end if
	WRITE(10,*)layers
	do i=1, layers-1
		WRITE(10,*)i, cu, e, v
	end do
	WRITE(10,*)i, cu*1000, e*1000, '0.48'

end subroutine

