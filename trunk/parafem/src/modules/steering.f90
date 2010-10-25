MODULE STEERING

  !/****h* /steering
  !*  NAME
  !*    MODULE: steering
  !*  SYNOPSIS
  !*    Usage:      USE steering
  !*  FUNCTION
  !*    Contains subroutines that relate equations or nodes to elements
  !*
  !*    Subroutine             Purpose
  !*    
  !*    REARRANGE              Modifies REST array
  !*    FIND_G                 Finds g from node numbers and restraints "rest"
  !*    FIND_NO                Forms vector of loaded equations NO
  !*    ABAQUS2SG              Swaps node order from Abaqus to S&G convention
  !*
  !*  AUTHOR
  !*    F. Calvo
  !*    L. Margetts
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  USE precision
  USE mp_interface

  CONTAINS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE REARRANGE(rest)

    !/****f* steering/rearrange
    !*  NAME
    !*    SUBROUTINE: rearrange
    !*  SYNOPSIS
    !*    Usage:      CALL rearrange(rest)
    !*  FUNCTION
    !*    Modifies REST array, converting restrained flags into restrained
    !*    equation numbers.   
    !*  INPUTS
    !*    The following argument has the INTENT(INOUT) attribute:
    !*
    !*    rest                : Integer
    !*                        : Array of restrained nodes
    !*                        : Possible input values are 1 or 0
    !*                        : 1 - unrestrained
    !*                        : 0 - restrained
    !*  AUTHOR
    !*    I.M Smith
    !*  COPYRIGHT
    !*    (c) University of Manchester 2004-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/
    
    IMPLICIT NONE
    
    INTEGER, INTENT(INOUT) :: rest(:,:)
    INTEGER                :: nr, nodof, i, k, m

    nr    = ubound(rest,1)
    nodof = ubound(rest,2) - 1

    m = 0

    DO i = 1, nr
      DO k = 2, nodof+1
        IF (rest(i,k)/=0) THEN
          m = m + 1
          rest(i,k) = m
        END IF
      END DO
    END DO

    DO i = 1, nr
      k = rest(i,1)
      DO m = 2, nodof+1
        IF (rest(i,m)/=0) rest(i,m) = rest(i,m) + nodof*(k-i)
      END DO
    END DO

  END SUBROUTINE rearrange

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE find_g3(num,g,rest)

    !/****f* steering/find_g3
    !*  NAME
    !*    SUBROUTINE: find_g3
    !*  SYNOPSIS
    !*    Usage:      CALL find_g3(num,g,rest)
    !*  FUNCTION
    !*    Finds g from node numbers and restraints "rest"
    !*  INPUTS
    !*    The following argument has the INTENT(INOUT) attribute:
    !*
    !*    rest                : Integer
    !*                        : Array of restrained nodes
    !*                        : Possible input values are 1 or 0
    !*                        : 1 - unrestrained
    !*                        : 0 - restrained
    !*  AUTHOR
    !*    I.M Smith
    !*  COPYRIGHT
    !*    (c) University of Manchester 2004-2010
    !******
    !*  Seems to have a bug 
    !*  Smith and Griffiths "Programming the Finite Element Method", Edition 4
    !*
    !*/

    INTEGER, INTENT(IN)  :: rest(:,:),num(:)
    INTEGER, INTENT(OUT) :: g(:)
    INTEGER              :: i,j,l,nod,nodof,s1,s2,s3,only,first,last,half
    LOGICAL              :: found
    
    nod   = ubound(num,1) 
    nodof = ubound(rest,2) - 1

    DO i=1,nod 

      l     = num(i)
      first = 1 
      last  = ubound(rest,1)

      DO
        IF(first==last) EXIT ! silly array or converged
        half = (first + last)/2
        IF(l<=rest(half,1)) THEN
          last = half  ! discard second half
        ELSE
          first = half + 1  ! discard first half
        END IF
      END DO

      only  = first
      found = (l==rest(only,1))

      IF(found) THEN
        DO j = 1 , nodof
          g(nodof*i-nodof+j) = rest(only,j+1)  
        END DO
      ELSE
!       k = 1
!       DO
!         s1 = only - k
!         IF(sum(rest(s1,2:))/=0) EXIT
!         k = k + 1
!       END DO

        s1 = only -1
        s2 = maxval(rest(s1,2:))    
        s3 = l - rest(s1,1)
  
        DO j=1,nodof
          g(nodof*i-nodof+j) = s2+nodof*s3 - nodof + j
        END DO
 
      END IF

    END DO
    
  END SUBROUTINE find_g3

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
  subroutine find_g(num,g,rest)
  ! Finds g from node numbers and restraints "rest"
  ! Only for uncoupled problems
  implicit none
  integer,intent(in)::num(:),rest(:,:); integer,intent(out)::g(:)
  integer::inc,i,j,k,l,count,nr,nod,nodof
  nr=ubound(rest,1); nod=ubound(num,1);nodof=ubound(rest,2)-1
  g=1;inc=0
  do i=1,nod
   count=0;l=num(i)
   do j=1,nr
      if(l>rest(j,1)) then
         do k=2,nodof+1; if(rest(j,k)==0)count=count+1;end do
      end if
   end do
   do k=2,nodof+1
    inc=inc+1
    do j=1,nr
     if(l==rest(j,1).and.rest(j,k)==0) then
       g(inc)=0; count=count+1
     end if
    end do
    if(g(inc)/=0) g(inc) = l*nodof-count-(nodof+1-k)
   end do
  end do
  end subroutine find_g
 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE g_t_g_ns(nod,g_t,g)
  
    !/****f* structure_dof/g_t_g_ns
    !*  NAME
    !*    SUBROUTINE: g_t_g_ns
    !*  SYNOPSIS
    !*    Usage:      CALL g_t_g_ns(nod,g_t,g)
    !*  FUNCTION
    !*    Finds g from g_t (Navier - Stokes)
    !*  INPUTS
    !*  AUTHOR
    !*    Ian M. Smith
    !*  COPYRIGHT
    !*    (c) University of Manchester
    !******
    !*  
    !*/

    IMPLICIT NONE
    
    INTEGER,INTENT(IN)  :: nod,g_t(:)
    INTEGER,INTENT(OUT) :: g(:)
    INTEGER             :: i
    
    DO i=1,nod
      g(i)    = g_t(4*i-3)
      g(i+28) = g_t(4*i-1)
      g(i+48) = g_t(4*i) 
    END DO
    
    DO i=1,4 
      g(i+20) = g_t(8*i-6)
      g(i+24) = g_t(8*i+42)
    END DO
    
  END SUBROUTINE g_t_g_ns

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE FIND_NO(NODE,REST,SENSE,NO)

    !/****f* structure_dof/find_no
    !*  NAME
    !*    SUBROUTINE: find_no
    !*  SYNOPSIS
    !*    Usage:      CALL find_no(node,rest,sense,no)
    !*  FUNCTION
    !*    Forms vector of loaded equations NO from vector of loaded nodes NODE
    !*    using freedom information stored in REST.  Will work for serial and
    !*    parallel programs.
    !*  INPUTS
    !*    It is assumed that the node numbers in both NODE and REST are arranged
    !*    in ascending order.
    !*  AUTHOR
    !*    Lee Margetts
    !*  CREATION DATE
    !*    25.08.2009
    !*  COPYRIGHT
    !*    (c) University of Manchester 2009-2010
    !******
    !*  Needs fully testing, for all the cases.
    !*/

    IMPLICIT NONE
    INTEGER,INTENT(IN)    :: node(:)
    INTEGER,INTENT(IN)    :: rest(:,:)
    INTEGER,INTENT(IN)    :: sense(:)
    INTEGER,INTENT(OUT)   :: no(:)
    INTEGER               :: i,j,k          ! simple counters
    INTEGER               :: nrf            ! number of restrained freedoms
    INTEGER               :: method = 0     ! for case arguments
    INTEGER               :: fixed_freedoms ! number of fixed_freedoms
    INTEGER               :: nodof          ! number of degrees of freedom
    INTEGER               :: nr             ! number of restrained nodes
  
    
    fixed_freedoms = UBOUND(node,1)
    nodof          = UBOUND(rest,2) - 1
    nr             = UBOUND(rest,1)
    
    IF(node(fixed_freedoms) < rest(1,1)) method = 1
    IF(node(1) > rest(nr,1))             method = 2

    SELECT CASE (method)
    
    CASE(1)
    
      DO i = 1, fixed_freedoms
        no(i) = ((node(i) - 1) * nodof) + sense(i)
      END DO
    
    CASE(2)
      
      nrf = 0
      
      DO i = 1,nr
        DO j = 1,nodof
          IF(rest(i,j+1) == 0) nrf = nrf + 1
        END DO
      END DO
      
      DO i = 1,fixed_freedoms
        no(i) = ((node(i) - 1) * nodof) + sense(i) - nrf
      END DO
      
    CASE default
    
      DO i = 1,fixed_freedoms
     
        nrf = 0 
      
        DO j = 1,nr
          IF(rest(j,1) < node(i) ) THEN
            DO k = 1,nodof
              IF(rest(j,k+1) == 0) nrf = nrf + 1
            END DO
          ELSE IF (rest(j,1) == node(i)) THEN
            DO k = 1,sense(i) - 1
              IF(rest(j,k+1) == 0) nrf = nrf + 1
            END DO
          ELSE IF(rest(j,1) > node(i)) THEN
            EXIT
          END IF
        END DO
        
        no(i) = ((node(i) - 1) * nodof) + sense(i) - nrf
        
      END DO
      
    END SELECT
    
    PRINT*, "NO =", no
    
    RETURN
  
  END SUBROUTINE FIND_NO
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------  
  
    SUBROUTINE REINDEX_FIXED_NODES (ieq_start,no,no_local_temp,num_no,  &
                                    no_index_start,neq_pp)
  
      !/****f* structure_dof/reindex_fixed_nodes
      !*  NAME
      !*    SUBROUTINE: reindex_fixed_nodes
      !*  SYNOPSIS
      !*    Usage:      CALL reindex_fixed_nodes(ieq_start,no,no_local_temp,  &
      !*                                         num_no,no_index_start,neq_pp)
      !*  FUNCTION
      !*    Creates a local index of loaded equations. Will function in both
      !*    serial and MPI based programs.
      !*  INPUTS
      !*    The following arguments have the INTENT(IN) attribute:
      !*
      !*    ieq_start           : Integer
      !*                        : The starting equation number on the processor
      !*
      !*    neq_pp              : Integer
      !*                          Number of equations per processor
      !*
      !*    no(:)               : Integer
      !*                        : List of the nodes with some degree of freedom
      !*                          loaded.
      !*
      !*    The following arguments have the INTENT(INOUT) attribute:    
      !*
      !*    no_local_temp(:)    : Integer
      !*                        : A temporary array
      !*
      !*  OUTPUTS
      !*    The following arguments have the INTENT(OUT) attribute:
      !*
      !*    num_no              : Integer
      !*                        : The number of local fixed nodes
      !*
      !*    no_index_start      : Integer
      !*                        : The starting position in the array no_local_temp
      !*    
      !*  AUTHOR
      !*    Lee Margetts
      !*  CREATION DATE
      !*    24.04.2002
      !*  COPYRIGHT
      !*    (c) University of Manchester
      !******
      !*  Note: This looks more like REINDEX_FIXED_EQUATIONS. Rename?
      !*
      !*/
      
      IMPLICIT NONE
  
      INTEGER, INTENT(IN)     :: ieq_start, neq_pp, no(:)
      INTEGER, INTENT(INOUT)  :: no_local_temp(:)
      INTEGER, INTENT(OUT)    :: num_no, no_index_start
      INTEGER                 :: fixed_nodes, i
  
      fixed_nodes    = UBOUND(no,1)
      no_index_start = 0
      num_no         = 0
  
      DO i = 1,fixed_nodes
        IF (no(i) < ieq_start) THEN
          CYCLE
        ELSE IF (no(i) >= ieq_start + neq_pp) THEN
          EXIT
        ELSE IF (no_index_start==0) THEN
          no_index_start   = i
          no_local_temp(1) = no(i)
          num_no           = 1
        ELSE
          num_no                = num_no + 1
          no_local_temp(num_no) = no(i)
        END IF
      END DO
   
  END SUBROUTINE REINDEX_FIXED_NODES
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------  
  
  SUBROUTINE ABAQUS2SG(element,g_num)

    !/****f* steering/abaqus2sg
    !*  NAME
    !*    SUBROUTINE: abaqus2sg
    !*  SYNOPSIS
    !*    Usage:      CALL abaqus2sg(element,g_num)
    !*  FUNCTION
    !*    Converts the node steering array from the ABAQUS node ordering to 
    !*    the ordering in Smith and Griffiths "Programming the Finite Element
    !*    Method", 4th Edition. Works with both serial and parallel programs.
    !*
    !*    For 20-node hexahedra, the conversion is:
    !*
    !*    S&G   :   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
    !*    Abaqus:   1  7 19 13  3  5 17 15  8 12 20  9  4 11 16 10  2  6 18 14 
    !*  INPUTS
    !*
    !*  AUTHOR
    !*    L. Margetts
    !*  COPYRIGHT
    !*    (c) University of Manchester 2009-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/
  
    IMPLICIT NONE
    INTEGER, INTENT(INOUT)        :: g_num(:,:)
    INTEGER, ALLOCATABLE          :: temp(:)
    INTEGER                       :: iel
    INTEGER                       :: nod
    INTEGER                       :: nels
    CHARACTER(LEN=15), INTENT(IN) :: element

!------------------------------------------------------------------------------
! 1. Find dimensions of problem and allocate temp
!------------------------------------------------------------------------------

    nels = UBOUND(g_num,2)
    nod  = UBOUND(g_num,1)
    
!------------------------------------------------------------------------------
! 2. Identify which type of element is being dealt with and reorder the nodes
!------------------------------------------------------------------------------

    SELECT CASE (element)
    
    CASE('hexahedron')
    
      SELECT CASE (nod)
      
      CASE(20)
   
      ALLOCATE(temp(nod))
      temp = 0

      DO iel = 1,nels
        temp          = g_num(:,iel)
        g_num(1,iel)  = temp(4)
        g_num(2,iel)  = temp(12)
        g_num(3,iel)  = temp(1)
        g_num(4,iel)  = temp(9)
        g_num(5,iel)  = temp(2)
        g_num(6,iel)  = temp(10)
        g_num(7,iel)  = temp(3)
        g_num(8,iel)  = temp(11)
        g_num(9,iel)  = temp(20)
        g_num(10,iel) = temp(17)
        g_num(11,iel) = temp(18)
        g_num(12,iel) = temp(19)
        g_num(13,iel) = temp(8)
        g_num(14,iel) = temp(16)
        g_num(15,iel) = temp(5)
        g_num(16,iel) = temp(13)
        g_num(17,iel) = temp(6)
        g_num(18,iel) = temp(14)
        g_num(19,iel) = temp(7)
        g_num(20,iel) = temp(15)
      END DO
        
      DEALLOCATE(temp)
     
      CASE(8)
      
      ALLOCATE(temp(nod))
      temp = 0

      DO iel = 1,nels
        temp          = g_num(:,iel)
        g_num(1,iel)  = temp(1)
        g_num(2,iel)  = temp(5)
        g_num(3,iel)  = temp(6)
        g_num(4,iel)  = temp(2)
        g_num(5,iel)  = temp(4)
        g_num(6,iel)  = temp(8)
        g_num(7,iel)  = temp(7)
        g_num(8,iel)  = temp(3)
      END DO
     
      DEALLOCATE(temp)
      
      CASE default
      
        PRINT*
        PRINT*, "This element type not supported in ABAQUS2SG"
        PRINT*, "Program aborting"
        PRINT*

        CALL shutdown()
 
      END SELECT
    
    CASE('tetrahedron')
      
      SELECT CASE (nod)
      
      CASE(4)

      ALLOCATE(temp(nod))
      
      DO iel = 1,nels
        temp          = g_num(:,iel)
        g_num(1,iel)  = temp(1)
        g_num(2,iel)  = temp(3)
        g_num(3,iel)  = temp(4)
        g_num(4,iel)  = temp(2)
      END DO

      DEALLOCATE(temp)

      CASE default
        
        PRINT*
        PRINT*, "Wrong number of nodes for a tetrahedron"
        PRINT*, "Program aborting"
        PRINT*
        
        CALL shutdown()

      END SELECT
    
    CASE default
    
      PRINT*
      PRINT*, "Wrong type of element in subroutine ABAQUS2SG"
      PRINT*, "Program aborting"
      PRINT*
     
      CALL shutdown()
    
    END SELECT
    
    RETURN
    
    END SUBROUTINE ABAQUS2SG
    
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

END MODULE STEERING
