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

  SUBROUTINE find_g(num,g,rest)

    !/****f* steering/find_g
    !*  NAME
    !*    SUBROUTINE: find_g
    !*  SYNOPSIS
    !*    Usage:      CALL find_g(num,g,rest)
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
    !*  Was find_g3 in: 
    !*  Smith and Griffiths "Programming the Finite Element Method", Edition 4
    !*
    !*/

    INTEGER, INTENT(IN)  :: rest(:,:),num(:)
    INTEGER, INTENT(OUT) :: g(:)
    INTEGER              :: i,j,k,l,nod,nodof,s1,s2,s3,only,first,last,half
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
        k = 1
        DO
          s1 = only - k
          IF(sum(rest(s1,2:))/=0) EXIT
          k = k + 1
        END DO

        s2 = maxval(rest(s1,2:))    
        s3 = l - rest(s1,1)
  
        DO j=1,nodof
          g(nodof*i-nodof+j) = s2+nodof*s3 - k*nodof + j
        END DO
 
      END IF

    END DO
    
  END SUBROUTINE find_g
 
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
    
    ALLOCATE(temp(nod))

!------------------------------------------------------------------------------
! 2. Identify which type of element is being dealt with and reorder the nodes
!------------------------------------------------------------------------------

    SELECT CASE (element)
    
    CASE('hexahedron')
    
      SELECT CASE (nod)
      
      CASE(20)
   
!     This loop could be optimized by a compiler
   
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
        
      CASE default
      
        PRINT*, "Wrong number of nodes for a hexahedron"
        
      END SELECT
      
    CASE default
    
      PRINT*, "Wrong type of element in subroutine ABAQUS2SG"
    
    END SELECT
    
    DEALLOCATE(temp)
    
    RETURN
    
    END SUBROUTINE ABAQUS2SG
    
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

END MODULE STEERING
