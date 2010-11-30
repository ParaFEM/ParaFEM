MODULE LOADING

  !/****h* /loading
  !*  NAME
  !*    MODULE: loading
  !*  SYNOPSIS
  !*    Usage:      USE loading
  !*  FUNCTION
  !*    Contains subroutines used for loading. These subroutines are parallel 
  !*    and require MPI.
  !*    
  !*    Subroutine             Purpose
  !*
  !*    LOAD                   Creates the distributed applied loads vector
  !*  AUTHOR
  !*    L. Margetts
  !*    I.M. Smith
  !*    D.V. Griffiths
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/
  
  USE precision
  USE gather_scatter

  CONTAINS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE LOAD(g_g_pp,g_num_pp,load_node,load_value,fext_pp)

  !/****f* loading/load
  !*  NAME
  !*    SUBROUTINE: load
  !*  SYNOPSIS
  !*    Usage:      CALL load(g_g_pp,g_num_pp,load_node,load_value,fext_pp)
  !*  FUNCTION
  !*    Build the distributed vector containing the natural boundary
  !*    conditions (forces in mechanics, fluxes in thermal problems)
  !*  
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    g_g_pp(ntot,nels_pp)    : Integer
  !*                            : Array linking equations to elements
  !*
  !*    g_num_pp(nod,nels_pp)   : Integer
  !*                            : Elements connectivity
  !*
  !*    load_node(loaded_nodes) : Integer
  !*                            : Node numbers that have a load imposed
  !*                          
  !*    load_value(ndim,loaded_nodes)
  !*                            : Real
  !*                            : Value of the force at each loaded node
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*                          
  !*    fext_pp(neq_pp)         : Real
  !*                            : Natural boundary conditions vector
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    26.05.2008
  !*  COPYRIGHT
  !*    (c) University of Manchester 2008-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  IMPLICIT NONE

  INTEGER,   INTENT(IN)    :: g_num_pp(:,:)         ! node to element
  INTEGER,   INTENT(IN)    :: g_g_pp(:,:)           ! equation to element
  INTEGER,   INTENT(IN)    :: load_node(:)          ! loaded node IDs
  INTEGER                  :: i, j, k, l, iel       ! loop indices
  INTEGER                  :: nels_pp, ipos
  INTEGER                  :: loaded_nodes, ntot, nodof, nod
  INTEGER                  :: foundFreedoms
  REAL(iwp), INTENT(IN)    :: load_value(:,:)       ! loaded node values
  REAL(iwp), INTENT(INOUT) :: fext_pp(:)            ! force vector
  REAL(iwp), ALLOCATABLE   :: r_temp(:,:)           ! temporary array
  REAL(iwp)                :: zero = 0.0_iwp 
  LOGICAL                  :: found

  found = .false. 

!------------------------------------------------------------------------------
! 1. Find values of common variables from imported arrays
!------------------------------------------------------------------------------

  loaded_nodes = UBOUND(load_node,1)
  nod          = UBOUND(g_num_pp,1)
  nels_pp      = UBOUND(g_num_pp,2)
  nodof        = UBOUND(load_value,1)
  ntot         = nod*nodof

!------------------------------------------------------------------------------
! 2. Allocate and zero temporary array
!------------------------------------------------------------------------------
    
  ALLOCATE(r_temp(ntot,nels_pp))
  r_temp  = zero

!------------------------------------------------------------------------------
! 3. Populate temporary array with load values
!------------------------------------------------------------------------------
   
  DO i = 1, loaded_nodes
    DO iel = 1, nels_pp
      DO j = 1, nod
        IF (g_num_pp(j,iel) == load_node(i)) THEN  !if it's a loaded node
          DO k = 1, nodof
            ipos = (nodof * (j-1)) + k 
            IF (g_g_pp(ipos,iel) > 0) THEN
              r_temp(ipos,iel) = load_value(k,i)
            END IF
          END DO
        END IF
      END DO
    END DO
  END DO

!------------------------------------------------------------------------------
! 4. Populate external force vector with load values
!------------------------------------------------------------------------------
    
  fext_pp = zero
  CALL SCATTER_NOADD(r_temp,fext_pp)
  
  DEALLOCATE(r_temp)
  
  RETURN
  
  END SUBROUTINE LOAD
 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE loading_p121(no,val,nle,nxe,nze)

  !/****f* loading/loading_p121
  !*  NAME
  !*    SUBROUTINE loading_p121
  !*  FUNCTION
  !*    Returns loaded freedoms NO and values VAL from number of elements in
  !*    x and z and number of loaded elements along a side of a square
  !*  INPUTS
  !*    CALL loading_p121(no,val, nle,nxe,nze)
  !*    
  !*    NO         : Integer. 
  !*    VAL        : Real. 
  !*    NLE        : Integer. 
  !*    NXE        : Integer.
  !*    NZE        : Integer.
  !*  AUTHOR
  !*    Smith and Griffiths, 4th Edition
  !****
  !*/ 

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)    :: nxe,nze,nle
  INTEGER,INTENT(OUT)   :: no(:)
  REAL(iwp),INTENT(OUT) :: val(:)  
  INTEGER               :: f1,f2,f3,count,i,j
  
  f1 = 6*nxe*nze
  f2 = (3*nxe+1)*nze
  f3 = (9*nxe+2)*nze

  count = 1
  
  DO  i = 1,4*nle+1,2   ! front face
    no(count) = i
    IF(i==1.OR.i==4*nle+1) THEN
     val(count) = -1._iwp
    ELSE IF (mod(i+1,4)==0) THEN
     val(count) = 4._iwp
    ELSE
     val(count) = -2._iwp
    END IF
    count = count + 1
  END DO
 
  DO j = 0 , nle - 1    ! other face pairs
 
    DO i = 1 , 3*nle + 1 , 3
      no(count) = i + f1 + j*(f2 + f3) + 1
      IF(i==1.OR.i==3*nle+1) THEN
        val(count) = 4._iwp
      ELSE
        val(count) = 8._iwp
      END IF
      count = count + 1
    END DO
 
    DO i = 1 , 6*nle+1 , 3
      no(count) = i + f1 +(j+1) * f2 + j * f3  + 1
      IF(j/=nle-1) THEN
        IF(i==1.OR.i==6*nle+1) THEN
          val(count) = -2._iwp
        ELSE IF (mod(i,2)==0) THEN
          val(count) = 8._iwp
        ELSE
          val(count) = -4._iwp
        END IF
      ELSE
        IF(i==1.OR.i==6*nle+1) THEN
          val(count) = -1._iwp
        ELSE IF (mod(i,2)==0) THEN
          val(count) = 4._iwp
        ELSE
          val(count) = -2._iwp
        END IF
      END IF
      count = count + 1
    END DO
 
  END DO
 
  END SUBROUTINE loading_p121

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

 SUBROUTINE loading_p126(no,nxe,nye,nze)

   !/****f* loading/loading_p126
   !*  NAME
   !*    SUBROUTINE: loading_p126
   !*  SYNOPSIS
   !*    Usage:      CALL loading_p126(no,nxe,nye,nze)
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
   
 END SUBROUTINE loading_p126
 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE load_p121(node,val,nle,nxe,nze)

  !/****f* loading/load_p121
  !*  NAME
  !*    SUBROUTINE load_p121
  !*  FUNCTION
  !*    Returns loaded nodes NODE and values VAL from number of elements in
  !*    x and z and number of loaded elements along a side of a square
  !*  INPUTS
  !*    CALL load_p121(node,val, nle,nxe,nze)
  !*    
  !*    NODE       : Integer. 
  !*    VAL        : Real. 
  !*    NLE        : Integer. 
  !*    NXE        : Integer.
  !*    NZE        : Integer.
  !*  AUTHOR
  !*    Smith and Griffiths, 4th Edition
  !****
  !*/ 

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)    :: nxe,nze,nle
  INTEGER,INTENT(OUT)   :: node(:)
  REAL(iwp),INTENT(OUT) :: val(:)  
  INTEGER               :: f1,f2,f3,count,i,j
  
  f1 = ((2*nxe+1)*(nze+1)+(nxe+1)*nze)
  f2 = (nxe+1)*(nze+1)

  count = 1
  
  DO  i = 1,2*nle+1   ! front face
    node(count) = i
    IF(i==1.OR.i==2*nle+1) THEN
     val(count) = -1._iwp
    ELSE IF (mod(i,2)==0) THEN
     val(count) = 4._iwp
    ELSE
     val(count) = -2._iwp
    END IF
    count = count + 1
  END DO
 
  DO j = 0 , nle - 1    ! other face pairs
 
    DO i = 1 , nle + 1
      node(count) = i + f1 + j*(f1 + f2)
      IF(i==1.OR.i==nle+1) THEN
        val(count) = 4._iwp
      ELSE
        val(count) = 8._iwp
      END IF
      count = count + 1
    END DO
 
    DO i = 1 , 2*nle+1
      node(count) = i + f1 +(j+1) * f2 + j * f1
      IF(j/=nle-1) THEN
        IF(i==1.OR.i==2*nle+1) THEN
          val(count) = -1._iwp
        ELSE IF (mod(i,2)==0) THEN
          val(count) = 8._iwp
        ELSE
          val(count) = -4._iwp
        END IF
      ELSE
        IF(i==1.OR.i==2*nle+1) THEN
          val(count) = -1._iwp
        ELSE IF (mod(i,2)==0) THEN
          val(count) = 4._iwp
        ELSE
          val(count) = -2._iwp
        END IF
      END IF
      count = count + 1
    END DO
 
  END DO
 
  END SUBROUTINE load_p121
  
END MODULE LOADING
