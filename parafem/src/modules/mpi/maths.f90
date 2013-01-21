MODULE MATHS

  !/****h* /maths
  !*  NAME
  !*    MODULE: maths
  !*  SYNOPSIS
  !*    Usage:      USE maths
  !*  FUNCTION
  !*    Contains subroutines required for standard mathematical operations on
  !*    matrices and vectors. These subroutines are parallel and require MPI.
  !*    
  !*    Subroutine             Purpose
  !*
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

  USE mpi_wrapper
  USE precision
  USE mp_interface

  CONTAINS

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------

  FUNCTION MAX_INTEGER_P(n_pp) RESULT(n)

  !/****f* maths/max_integer_p
  !*  NAME
  !*    FUNCTION:   MAX_INTEGER_P
  !*  SYNOPSIS
  !*    Usage:      max_integer_p(n_pp)
  !*  FUNCTION
  !*    Maximum value of an integer across processors
  !*    Every process keeps the result
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    n_pp:          : Integer
  !*                   : Integer on each processor
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    n:             : Integer
  !*                   : Maximum value across processors
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2002-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n_pp
  INTEGER :: n, bufsize, ier

  bufsize = 1
  CALL MPI_ALLREDUCE(n_pp,n,bufsize,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

  END FUNCTION MAX_INTEGER_P

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  FUNCTION MAX_P(n_pp) RESULT(n)

  !/****f* maths/max_p
  !*  NAME
  !*    FUNCTION:   MAX_P
  !*  SYNOPSIS
  !*    Usage:      max_p(n_pp)
  !*  FUNCTION
  !*    Maximum value of an integer across processors
  !*    Every process keeps the result
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    n_pp:          : Integer
  !*                   : Integer on each processor
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    n:             : Integer
  !*                   : Maximum value across processors
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2002-2013
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  Same as MAX_INTEGER_P
  !*  Place holder subroutine for one that uses an interface to swap 
  !*  between REAL and INTEGER versions
  !*/

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n_pp
  INTEGER :: n, bufsize, ier

  bufsize = 1
  CALL MPI_ALLREDUCE(n_pp,n,bufsize,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

  END FUNCTION MAX_P

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  FUNCTION MAX_REAL_P(r_pp) RESULT(r)

  !/****f* maths/max_real_r
  !*  NAME
  !*    FUNCTION:   MAX_REAL_P
  !*  SYNOPSIS
  !*    Usage:      max_real_p(r_pp)
  !*  FUNCTION
  !*    Maximum value of a real across processors
  !*    Every process keeps the result
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    r_pp:          : Real
  !*                   : Real on each processor
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    r:             : Real
  !*                   : Maximum value across processors
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2002-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  IMPLICIT NONE

  REAL(iwp), INTENT(IN) :: r_pp
  REAL(iwp)             :: r
  INTEGER               :: bufsize, ier

  bufsize = 1
  CALL MPI_ALLREDUCE(r_pp,r,bufsize,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ier)

  END FUNCTION MAX_REAL_P

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  FUNCTION DOT_PRODUCT_P(vectora,vectorb)

  !/****f* maths/dot_product_p
  !*  NAME
  !*    FUNCTION:   DOT_PRODUCT_P
  !*  SYNOPSIS
  !*    Usage:      dot_product_p(r_pp)
  !*  FUNCTION
  !*    Parallel version of the serial fortran intrinsic function dot_product.
  !*    Each processor computes the local dot_product. A global sum is then
  !*    performed across processors with the result being returned to all 
  !*    processors.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    vectora(:)     : Real
  !*                   : Chunk of vector a on each processor
  !*
  !*    vectorb(:)     : Real
  !*                   : Chunk of vector b on each processor
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    dot_product_p  : Real
  !*                   : Dot product of the distributed vectors a and b
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2002-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  Parallel version of the serial intrinsic function DOT_PRODUCT
  !*
  !*/

  IMPLICIT NONE

  REAL(iwp), INTENT(IN) :: vectora(:), vectorb(:)
  INTEGER               :: bufsize, ier
  REAL(iwp)             :: local_product, dot_product_p

  local_product = DOT_PRODUCT(vectora,vectorb)

  bufsize = 1
  CALL MPI_ALLREDUCE(local_product,dot_product_p,bufsize,MPI_REAL8,MPI_SUM, &
                     MPI_COMM_WORLD,ier)

  END FUNCTION DOT_PRODUCT_P

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  FUNCTION NORM_P(vectora)

  !/****f* maths/norm_p
  !*  NAME
  !*    FUNCTION:   NORM_P
  !*  SYNOPSIS
  !*    Usage:      norm_p(vectora)
  !*  FUNCTION
  !*    Calculates L2-norm of a distrubited vector
  !*    Every process keeps the result
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    vectora(:)     : Real
  !*                   : Chunk of vector a on each processor
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    norm_p         : Real
  !*                   : Norm of the vector
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2002-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  IMPLICIT NONE

  REAL(iwp), INTENT(IN) :: vectora(:)
  INTEGER               :: bufsize, ier
  REAL(iwp)             :: local_product, global_product, norm_p

  local_product = DOT_PRODUCT(vectora,vectora)

  bufsize = 1
  CALL MPI_ALLREDUCE(local_product,global_product,bufsize,MPI_REAL8,MPI_SUM,&
                     MPI_COMM_WORLD,ier)

  norm_p = SQRT(global_product)
   
  END FUNCTION NORM_P

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  FUNCTION SUM_P(vectora)

  !/****f* maths/sum_p
  !*  NAME
  !*    FUNCTION:   SUM_P
  !*  SYNOPSIS
  !*    Usage:      sum_p(vectora)
  !*  FUNCTION
  !*    Parallel version of the serial fortran instrinsic SUM. The function
  !*    sums all the components of a distributed vector. The result of SUM_P
  !*    is returned to all processors. This function requires MPI.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    vectora(:)     : Real
  !*                   : Local part of vector stored on calling processor
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    sum_p          : Real
  !*                   : Sum of all the components of the vector
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2002-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  IMPLICIT NONE

  REAL(iwp), INTENT(IN) :: vectora(:)
  INTEGER               :: bufsize, ier
  REAL(iwp)             :: local_sum, global_sum , sum_p

  local_sum = SUM(vectora)

  bufsize = 1
  CALL MPI_ALLREDUCE(local_sum,global_sum,bufsize,MPI_REAL8,MPI_SUM,    &
                     MPI_COMM_WORLD,ier)

  sum_p = global_sum

  END FUNCTION SUM_P

  FUNCTION ISUM_P(vectora)

  !/****f* maths/isum_p
  !*  NAME
  !*    FUNCTION:   ISUM_P
  !*  SYNOPSIS
  !*    Usage:      isum_p(vectora)
  !*  FUNCTION
  !*    Parallel version of the serial fortran instrinsic SUM. The function
  !*    sums all the components of a distributed vector. The result of ISUM_P
  !*    is returned to all processors. This function requires MPI.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    vectora(:)     : Integer
  !*                   : Local part of vector stored on calling processor
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    isum_p         : Integer
  !*                   : Sum of all the components of the vector
  !*  AUTHOR
  !*    L. Margetts, Louise Lever
  !*  COPYRIGHT
  !*    (c) University of Manchester 2002-2012
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: vectora(:)
  INTEGER             :: bufsize, ier
  INTEGER             :: local_isum, global_isum , isum_p

  local_isum = SUM(vectora)

  bufsize = 1
  CALL MPI_ALLREDUCE(local_isum,global_isum,bufsize,MPI_INTEGER,MPI_SUM,    &
                     MPI_COMM_WORLD,ier)

  isum_p = global_isum

  END FUNCTION ISUM_P

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  FUNCTION MAXVAL_P(vectora)

  !/****f* maths/maxval_p
  !*  NAME
  !*    FUNCTION:   MAXVAL_P
  !*  SYNOPSIS
  !*    Usage:      maxval_p(vectora)
  !*  FUNCTION
  !*    Parallel version of the serial intrinsic fortran function MAXVAL.
  !*    Computes the maximum value of the components of a distributed vector, 
  !*    and then computes the maximum value across all vectors. The result
  !*    is returned to all processors.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    vectora(:)     : Real
  !*                   : Chunk of vector a on each processor
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    maxval_p       : Real
  !*                   : Maximum value of the components of the vector
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2002-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*  Parallel version of the serial intrinsic function MAXVAL
  !*
  !*/

  IMPLICIT NONE

  REAL(iwp), INTENT(IN) :: vectora(:)
  INTEGER :: bufsize, ier
  REAL(iwp) :: local_max, global_max , maxval_p

  local_max = MAXVAL(vectora)

  bufsize = 1
  CALL MPI_ALLREDUCE(local_max,global_max,bufsize,MPI_REAL8,MPI_MAX,    &
                     MPI_COMM_WORLD,ier)

  maxval_p = global_max

  END FUNCTION MAXVAL_P

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  REAL FUNCTION MAXABSVAL_P(vectora)

  !/****f* maths/maxabsval_p
  !*  NAME
  !*    FUNCTION:   MAXABSVAL_P
  !*  SYNOPSIS
  !*    Usage:      maxabsval_p(vectora)
  !*  FUNCTION
  !*    Computes the maximum absolute value of the components of a distributed
  !*    vector and then computes the maximum value across all vectors. 
  !*    The result is returned to all processors.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    vectora(:)     : Real
  !*                   : Chunk of vector a on each processor
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    maxabsval_p    : Real
  !*                   : Maximum absolute value of the components of the
  !*                     vector
  !*  AUTHOR
  !*    F. Calvo
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2008-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  IMPLICIT NONE

  REAL(iwp), INTENT(IN)  :: vectora(:)
  REAL(iwp)              :: val_pp, auxval, val
  INTEGER                :: bufsize, ier, i, n_pp

  n_pp        = UBOUND(vectora,1)
  val_pp      = 0.0_iwp
  maxabsval_p = 0.0_iwp

!------------------------------------------------------------------------------
! 1. Every processor finds its own local absolute maximum 
!------------------------------------------------------------------------------   
  
  DO i = 1,n_pp
    auxval = ABS(vectora(i))
    IF (auxval > val_pp) THEN
      val_pp = auxval
    END IF
  END DO

!------------------------------------------------------------------------------
! 2. Find the global absolute maximum across all processors
!------------------------------------------------------------------------------   

  bufsize = 1
  CALL MPI_ALLREDUCE(val_pp,val,bufsize,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ier)

  maxabsval_p = val

  RETURN

  END FUNCTION MAXABSVAL_P

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE INVERT(matrix)

    !/****f* maths/invert
    !*  NAME
    !*    SUBROUTINE: invert(matrix)
    !*  SYNOPSIS
    !*    Usage:      CALL invert(matrix)
    !*  FUNCTION
    !*    Invert a small square matrix onto itself.
    !*  AUTHOR
    !*    I.M. Smith
    !*    D.V. Griffiths
    !*  COPYRIGHT
    !*    (c) University of Manchester 2004-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    IMPLICIT NONE
 
    REAL(iwp),INTENT(INOUT) :: matrix(:,:)
    REAL(iwp)               :: det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
    INTEGER                 :: ndim,i,k
 
    ndim = UBOUND(matrix,1)
   
    IF(ndim==2)THEN
      det         =  matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
      j11         =  matrix(1,1)
      matrix(1,1) =  matrix(2,2)
      matrix(2,2) =  j11
      matrix(1,2) = -matrix(1,2)
      matrix(2,1) = -matrix(2,1)
      matrix      =  matrix/det
    ELSE IF(ndim==3)THEN
      det =  matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
      det =  det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
      det =  det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
      j11 =  matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
      j21 = -matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
      j31 =  matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
      j12 = -matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
      j22 =  matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
      j32 = -matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
      j13 =  matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
      j23 = -matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
      j33 =  matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
      matrix(1,1) = j11
      matrix(1,2) = j12
      matrix(1,3) = j13
      matrix(2,1) = j21
      matrix(2,2) = j22
      matrix(2,3) = j23
      matrix(3,1) = j31
      matrix(3,2) = j32
      matrix(3,3) = j33
      matrix      = matrix/det
    ELSE
      DO k=1,ndim
        con         = matrix(k,k)
        matrix(k,k) = 1.0_iwp
        matrix(k,:) = matrix(k,:)/con
        DO i=1,ndim
          IF(i/=k)THEN
            con         = matrix(i,k)
            matrix(i,k) = 0.0_iwp
            matrix(i,:) = matrix(i,:)-matrix(k,:)*con
          END IF
        END DO
      END DO
    END IF

    RETURN
  
  END SUBROUTINE INVERT

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  FUNCTION DETERMINANT (jac) RESULT(det)

    !/****f* maths/determinant
    !*  NAME
    !*    FUNCTION: determinant
    !*  SYNOPSIS
    !*    Usage:    determinant
    !*  FUNCTION
    !*    Returns the determinant of a 1x1 2x2 3x3 jacobian matrix
    !*    Source: "Programming the Finite Element Method"
    !*  INPUTS
    !*    The following argument has the INTENT(IN) attribute:
    !*
    !*    jac(:,:)       : Real
    !*                   : Matrix whose determinant is to be computed
    !*
    !*  AUTHOR
    !*    I.M. Smith
    !*    D.V. Griffiths
    !*  COPYRIGHT
    !*    (c) University of Manchester 2004-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/
    
    IMPLICIT NONE
    
    REAL(iwp)            :: det
    REAL(iwp),INTENT(IN) :: jac(:,:)
    INTEGER              :: it 
    
    it = UBOUND(jac,1)  
   
    SELECT CASE (it)
     
      CASE (1)
  
      det=1.0_iwp
     
      CASE (2)
     
      det=jac(1,1)*jac(2,2) - jac(1,2) * jac(2,1)
     
      CASE (3)
     
      det= jac(1,1)*(jac(2,2) * jac(3,3) -jac(3,2) * jac(2,3))
      det= det-jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
      det= det+jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
     
      CASE DEFAULT
     
      print*,' wrong dimension for jacobian matrix'
   
    END SELECT
  
    RETURN
  
  END FUNCTION DETERMINANT

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

  SUBROUTINE PRINCIPALSTRESS3D(sigma,principal)

    !/****f* maths/principalstress3D
    !*  NAME
    !*    SUBROUTINE: principalstress3D
    !*  SYNOPSIS
    !*    Usage:      CALL principalstress3D(sigma,principal)
    !*  FUNCTION
    !*    Compute principal stresses in 3D
    !*    
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    sigma(nst)            : Real
    !*                          : Stress tensor (Voigt notation)
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    principal(ndim)       : Real
    !*                          : Principal stresses
    !*  AUTHOR
    !*    Francisco Calvo
    !*    Lee Margetts
    !*  CREATION DATE
    !*    10.09.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2010
   !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    IMPLICIT NONE

    REAL(iwp), INTENT(IN)  :: sigma(:)
    REAL(iwp), INTENT(OUT) :: principal(:)
    INTEGER                :: i
    REAL(iwp)              :: bd(6), pi, pi23, thrd, tol
    REAL(iwp)              :: al, b1, b2, b3, c1, c2,c3

    pi        = 4.0_iwp * ATAN(1.0_iwp)  
    pi23      = (pi/3.0_iwp) * 2.0_iwp
    thrd      = 1.0_iwp/3.0_iwp
    tol       = 1.0e-12_iwp
    bd        = 0.0_iwp
    principal = 0.0_iwp

    b1 = (sigma(1) + sigma(2) + sigma(3))*thrd

    DO i = 1,6
      bd(i) = sigma(i)
    END DO

    DO i = 1,3
      bd(i) = bd(i) - b1
    END DO

    c1 = bd(4)*bd(4)
    c2 = bd(5)*bd(5)
    c3 = bd(6)*bd(6)
    b2 = 0.5_iwp*(bd(1)*bd(1)+bd(2)*bd(2)+bd(3)*bd(3))+ c1 + c2 + c3
    IF(b2.le.tol*b1*b1) THEN
      principal(1) = b1
      principal(2) = b1
      principal(3) = b1
    ELSE
      b3 = bd(1)*bd(2)*bd(3)+(bd(4)+bd(4))*bd(5)*bd(6) +        &
           bd(1)*(c1-c2) + bd(2)*(c1-c3)
      c1 = 2.0_iwp*SQRT(b2*thrd)
      c2 = 4.0_iwp*b3
      c3 = c1*c1*c1
      al = ATAN2(SQRT(ABS(c3*c3-c2*c2)),c2)*thrd
      principal(1) = b1 + c1*COS(al)
      principal(2) = b1 + c1*COS(al-pi23)
      principal(3) = b1 + c1*COS(al+pi23)
    END IF

  END SUBROUTINE PRINCIPALSTRESS3D

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE NODAL_PROJECTION(npes,nn,nels_pp,g_num_pp,nod,numvar,nodes_pp, &
                node_start,node_end,shape_integral,stress_integral,stressnodes)

    !/****f* stress/nodal_projection
    !*  NAME
    !*    SUBROUTINE: nodal_projection
    !*  SYNOPSIS
    !*    Usage:      CALL nodal_projection(npes,nn,nels_pp,g_num_pp,nod, &
    !*                            numvar,nodes_pp,node_start,node_end,    &
    !*                            shape_integral,stress_integral,stressnodes)
    !*  FUNCTION
    !*    Project a variable (i.e. stress) from Gauss points to nodes
    !*    
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    nels_pp               : Integer
    !*                          : Number of elements per processor
    !*
    !*    nn                    : Integer
    !*                          : Total number of nodes
    !*
    !*    nod                   : Integer
    !*                          : Number of nodes per element
    !*
    !*    node_end              : Integer
    !*                          : Last node stored in the process
    !*
    !*    node_start            : Integer
    !*                          : First node stored in the process
    !*
    !*    nodes_pp              : Integer
    !*                          : Number of nodes assigned to the process
    !*
    !*    npes                  : Integer
    !*                          : Number of processes
    !*
    !*    numvar                : Integer
    !*                          : Number of components of the variable
    !*                            (1-scalar, 3-vector, 6-tensor)
    !*
    !*    g_num_pp(nod,nels_pp) : Integer
    !*                          : Elements connectivity
    !*
    !*    shape_integral(nod,nels_pp) : Real
    !*                                : Integral of shape functions in the
    !*                                  element
    !*
    !*    stress_integral(nod*numvar,nels_pp) : Real
    !*                                        : Integral of the variable (i.e.
    !*                                          stress) in the element
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    stressnodes(nodes_pp*numvar) : Real
    !*                                 : Nodal value after projection
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    10.11.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2010
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    IMPLICIT NONE

    INTEGER,   INTENT(IN)  :: nels_pp, nn, nod, node_end, node_start, &
                              nodes_pp, npes, numvar, g_num_pp(:,:)
    REAL(iwp), INTENT(IN)  :: shape_integral(:,:), stress_integral(:,:)
    REAL(iwp), INTENT(OUT) :: stressnodes(:)
    INTEGER :: i, j, k, iel, startNode, endNode, idx1, idx2, idx3, bufsize, &
               size_scalar, size_vector, size_scalar_max, size_vector_max, ier
    REAL(iwp), ALLOCATABLE :: all_assembly_shape(:), local_assembly_shape(:), &
                              local_assembly_stress(:)
    REAL(iwp), PARAMETER   :: zero = 0.0_iwp

    !-----------------------------------------------------------------------
    ! 1. Allocate internal arrays
    !-----------------------------------------------------------------------

    size_scalar_max =  nodes_pp+1
    size_vector_max = (nodes_pp+1)*numvar
    ALLOCATE(all_assembly_shape(size_scalar_max))
    ALLOCATE(local_assembly_shape(size_scalar_max))
    ALLOCATE(local_assembly_stress(size_vector_max))


    all_assembly_shape = zero
    stressnodes        = zero

    !-----------------------------------------------------------------------
    ! 2. Every process (first loop) sends the range of nodes assigned to that
    !    process. The other processes check whether the elements contain nodes
    !    in that range for the contribution. The contributions are summed up
    !    (assembled) and the process collects the results
    !-----------------------------------------------------------------------

    bufsize = 1

    DO i = 1,npes
      startNode = node_start
      endNode   = node_end
      CALL MPI_BCAST(startNode,bufsize,MPI_INTEGER,i-1,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST(endNode,bufsize,MPI_INTEGER,i-1,MPI_COMM_WORLD,ier)

      size_scalar = nodes_pp
      size_vector = nodes_pp*numvar
      CALL MPI_BCAST(size_scalar,bufsize,MPI_INTEGER,i-1,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST(size_vector,bufsize,MPI_INTEGER,i-1,MPI_COMM_WORLD,ier)

      local_assembly_stress = 0.0_iwp
      local_assembly_shape  = 0.0_iwp

      DO iel = 1,nels_pp
        DO j = 1,nod
          IF (g_num_pp(j,iel)>=startNode .and. g_num_pp(j,iel)<=endNode) THEN
            idx1 = (g_num_pp(j,iel)-startNode)*numvar
            idx2 = (j-1)*numvar
            idx3 = g_num_pp(j,iel)-startNode
            DO k = 1,numvar
              local_assembly_stress(idx1+k) = local_assembly_stress(idx1+k) + &
                                              stress_integral(idx2+k,iel)
            END DO
            local_assembly_shape(idx3+1)  = local_assembly_shape(idx3+1)  + &
                                            shape_integral(j,iel)
          END IF
        END DO
      END DO

      CALL MPI_REDUCE(local_assembly_shape, all_assembly_shape, size_scalar, &
                      MPI_REAL8, MPI_SUM, i-1, MPI_COMM_WORLD, ier)
      CALL MPI_REDUCE(local_assembly_stress, stressnodes, size_vector ,  &
                      MPI_REAL8, MPI_SUM, i-1, MPI_COMM_WORLD, ier)
    END DO

    !-----------------------------------------------------------------------
    ! 3. Compute the projected nodal value (i. e. stress)
    !-----------------------------------------------------------------------

    DO i = 1,nodes_pp
      idx1 = (i-1)*numvar
      DO j = 1,numvar
        stressnodes(idx1+j) = stressnodes(idx1+j)/all_assembly_shape(i)
      END DO
    END DO

    !-----------------------------------------------------------------------
    ! 4. Deallocate internal arrays
    !-----------------------------------------------------------------------

    DEALLOCATE(all_assembly_shape,local_assembly_shape,local_assembly_stress)
      
  END SUBROUTINE NODAL_PROJECTION
  
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

 SUBROUTINE INVERT_MATRIX_3x3(a,a_1)

    !/****f* maths/invert_matrix_3x3
    !*  NAME
    !*    SUBROUTINE: invert_matrix_3x3
    !*  SYNOPSIS
    !*    Usage:      CALL invert_matrix_3x3(a,a_1)
    !*  FUNCTION
    !*    Invert a 3x3 matrix
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    a(3,3):        : Real
    !*                   : Matrix to be inverted
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    a_1(3,3):      : Real
    !*                   : Inverse of the matrix
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    05.06.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  A programming branch by Fran that needs to be removed from the main
    !*  ParaFEM distribution
    !*/

    IMPLICIT NONE
	  
    REAL(iwp), INTENT(IN)  :: a(3,3)
    REAL(iwp), INTENT(OUT) :: a_1(3,3)
    REAL(iwp) :: det

    CALL DETERMINANT_MATRIX_3x3(a,det)
    det = 1._iwp/det

    a_1(1,1) =  (a(2,2)*a(3,3) - a(3,2)*a(2,3))*det
    a_1(1,2) = -(a(1,2)*a(3,3) - a(3,2)*a(1,3))*det
    a_1(1,3) =  (a(1,2)*a(2,3) - a(2,2)*a(1,3))*det
    a_1(2,1) = -(a(2,1)*a(3,3) - a(3,1)*a(2,3))*det
    a_1(2,2) =  (a(1,1)*a(3,3) - a(3,1)*a(1,3))*det
    a_1(2,3) = -(a(1,1)*a(2,3) - a(2,1)*a(1,3))*det
    a_1(3,1) =  (a(2,1)*a(3,2) - a(3,1)*a(2,2))*det
    a_1(3,2) = -(a(1,1)*a(3,2) - a(3,1)*a(1,2))*det
    a_1(3,3) =  (a(1,1)*a(2,2) - a(2,1)*a(1,2))*det

  END SUBROUTINE INVERT_MATRIX_3x3

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

  SUBROUTINE DETERMINANT_MATRIX_3x3(a,det)

    !/****f* maths/determinant_matrix_3x3
    !*  NAME
    !*    SUBROUTINE: determinant_matrix_3x3
    !*  SYNOPSIS
    !*    Usage:      CALL determinant_matrix_3x3(a,det)
    !*  FUNCTION
    !*    Compute the determinant of a matrix 3x3
    !*  INPUTS
    !*    The following arguments have the INTENT(IN) attribute:
    !*
    !*    a(3,3):        : Real
    !*                   : Matrix whose determinant is to be computed
    !*
    !*    The following arguments have the INTENT(OUT) attribute:
    !*
    !*    det:           : Real
    !*                   : Determinant of the matrix
    !*  AUTHOR
    !*    Francisco Calvo
    !*  CREATION DATE
    !*    05.06.2007
    !*  COPYRIGHT
    !*    (c) University of Manchester 2007-2011
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  A programming branch from Fran that needs to be removed from the main
    !*  ParaFEM distribution
    !*/

    IMPLICIT NONE

    REAL(iwp), INTENT(IN)  :: a(3,3)
    REAL(iwp), INTENT(OUT) ::  det

    det = a(1,1)*( a(2,2)*a(3,3) - a(2,3)*a(3,2) ) + &
          a(2,1)*( a(3,2)*a(1,3) - a(1,2)*a(3,3) ) + &
          a(3,1)*( a(1,2)*a(2,3) - a(1,3)*a(2,2) )

  END SUBROUTINE DETERMINANT_MATRIX_3x3

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  FUNCTION reduce(neq_temp) RESULT(neq)
  ! Maximum neq on any processor
   
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: neq_temp
  INTEGER             :: neq

  CALL MPI_ALLREDUCE(neq_temp,neq,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

  END FUNCTION reduce

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE checon_par(loads,tol,converged,oldlds)

  !/****f* pcg/checon_par
  !*  NAME
  !*    SUBROUTINE: checon_par
  !*  SYNOPSIS
  !*    Usage:      CALL checon_par(loads,tol,converged,oldlds)
  !*  FUNCTION
  !*    Parallel version of checon
  !*    Sets converged to .false. if relative change in loads and
  !*    oldlds is greater than tol and updates oldlds
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    loads(:)                 : Real
  !*                             : Solution vector
  !*
  !*    tol                      : Real
  !*                             : Solution tolerance
  !*
  !*    The following argument has the INTENT(INOUT) attribute:
  !*
  !*    oldlds(:)                : Real
  !*                             : The old solution vector
  !*
  !*    The following argument has the INTENT(OUT) attribute:
  !*
  !*    converged                : Logical
  !*                             : Solution has converged if 
  !*                               converged = .true.
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*    iters                    : Integer
  !*                             : Number of PCG iterations
  !*
  !*  AUTHOR
  !*    I.M. Smith and D.V. Griffiths
  !*    Modified by L. Margetts
  !*  COPYRIGHT
  !*    2004-2010 University of Manchester
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/

  IMPLICIT NONE

  REAL(iwp)                :: maxloads_pp,maxdiff_pp,maxloads,maxdiff
  REAL(iwp), INTENT(IN)    :: loads(:),tol
  REAL(iwp), INTENT(INOUT) :: oldlds(:)
  LOGICAL, INTENT(OUT)     :: converged
  
  maxloads_pp = maxval(abs(loads))
  maxdiff_pp  = maxval(abs(loads-oldlds))

  CALL MPI_ALLREDUCE(maxloads_pp,maxloads,1,MPI_REAL8,MPI_MAX,                &
                     MPI_COMM_WORLD,ier)
  CALL MPI_ALLREDUCE(maxdiff_pp,maxdiff,1,MPI_REAL8,MPI_MAX,                  &
                     MPI_COMM_WORLD,ier)

  converged = .true.
  converged = ((maxdiff/maxloads)<=tol)
  oldlds    = loads

! PRINT *, "MAXDIFF/MAXLOADS = ", maxdiff/maxloads
! Need to modify this subroutine to output convergence history
 
  RETURN

  END SUBROUTINE checon_par

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
  
END MODULE MATHS
