MODULE PCG

  !/****h* /pcg
  !*  NAME
  !*    MODULE: pcg
  !*  SYNOPSIS
  !*    Usage:      USE pcg
  !*  FUNCTION
  !*    Contains subroutines required by the preconditioned conjugate gradient
  !*    method. These subroutines are parallel and require MPI.
  !*    
  !*    Subroutine             Purpose
  !*
  !*    CHECON_PAR             Convergence test
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
  USE maths
  USE gather_scatter

  IMPLICIT NONE
  
  CONTAINS

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

END MODULE PCG
