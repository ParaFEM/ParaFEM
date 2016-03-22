MODULE choose_solvers
  
  !/****h* /choose_solvers
  !*  NAME
  !*    MODULE: choose_solvers
  !*  SYNOPSIS
  !*    Usage:      USE choose_solvers
  !*  FUNCTION
  !*    Contains data and subroutines that are used to choose between different
  !*    solver libraries, currently the solvers provided by ParaFEM and by
  !*    PETSc.
  !*    
  !*    Subroutine             Purpose
  !*
  !*    solvers_valid          Test the the solvers string is valid
  !*    solvers_list           Return a list of valid solvers
  !*  AUTHOR
  !*    Mark Filipiak
  !*  COPYRIGHT
  !*    2016 University of Edinburgh
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*  If you contribute to this module, add your author name.
  !*
  !*/
  
  IMPLICIT NONE
  
  PUBLIC
  
  ! Parameters
  CHARACTER(len=*),PARAMETER :: parafem_solvers   = "parafem"
  CHARACTER(len=*),PARAMETER :: petsc_solvers     = "petsc"
  
CONTAINS
  
  FUNCTION solvers_valid(solvers)

    !/****if* choose_solvers/solvers_valid
    !*  NAME
    !*    SUBROUTINE: solvers_valid
    !*  SYNOPSIS
    !*    Usage:      solvers_valid(solvers)
    !*  FUNCTION
    !*      Check if solvers matches one of the supported solvers
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*
    !*    solvers            : Character
    !*                         Name of the solvers
    !*  RESULT
    !*
    !*    solvers_valid      : Logical
    !*                         .true. if solvers matches one of the supported
    !*                         solvers, otherwise .false.
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    21.03.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 21.03.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    IMPLICIT NONE
    LOGICAL                      :: solvers_valid
    CHARACTER(len=*), INTENT(IN) :: solvers
    IF (     solvers == parafem_solvers                                        &
        .OR. solvers == petsc_solvers  ) THEN
      solvers_valid = .TRUE.
    ELSE
      solvers_valid = .FALSE.
    END IF
  END FUNCTION solvers_valid
  
  FUNCTION solvers_list()
    !/****if* choose_solvers/solvers_list
    !*  NAME
    !*    SUBROUTINE: solvers_list
    !*  SYNOPSIS
    !*    Usage:      solvers_list()
    !*  FUNCTION
    !*      Return list of supported solvers
    !*  ARGUMENTS
    !*    None
    !*  RESULT
    !*
    !*    solvers_list      : Character
    !*                        List of supported solvers
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    21.03.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 21.03.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    IMPLICIT NONE
    CHARACTER(:),allocatable :: solvers_list
    solvers_list = parafem_solvers                                             &
      // " or " // petsc_solvers
  END FUNCTION solvers_list
  
END MODULE choose_solvers
