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
 
    REAL(iwp),INTENT(IN OUT):: matrix(:,:)
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
  
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
  
END MODULE MATHS
