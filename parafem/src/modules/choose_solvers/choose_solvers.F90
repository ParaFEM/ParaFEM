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
  !*    get_solvers            Gets the type of solvers to use
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
  
  USE global_variables, ONLY: numpe
  IMPLICIT NONE
  
  PUBLIC
  
  ! Parameters
  INTEGER,PARAMETER          :: choose_solvers_string_length = 1024
  CHARACTER(len=*),PARAMETER :: parafem_solvers   = "parafem"
  CHARACTER(len=*),PARAMETER :: petsc_solvers     = "petsc"
  
CONTAINS

  FUNCTION get_solvers()

    !/****if* choose_solvers/get_solvers
    !*  NAME
    !*    SUBROUTINE: get_solvers
    !*  SYNOPSIS
    !*    Usage:      get_solvers
    !*  FUNCTION
    !*      Gets the type of solvers (e.g. ParaFEM or PETSc) to use
    !*  ARGUMENTS
    !*    None
    !*  RESULT
    !*
    !*    solvers            : Character
    !*                         Name of the solvers
    !*  AUTHOR
    !*    Mark Filipiak
    !*  CREATION DATE
    !*    14.04.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 14.04.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*/

    ! get_solvers cannot be CHARACTER(:),ALLOCATABLE because it needs to be
    ! able to hold any length of string from the command argument.
    CHARACTER(len=choose_solvers_string_length) :: get_solvers

    INTEGER :: argc, status

    ! Default solvers are ParaFEM.
    get_solvers = parafem_solvers
    
    argc = command_argument_count()
    ! The second argument in the command line is the name of the solvers.
    IF (argc >= 2) THEN
      CALL GET_COMMAND_ARGUMENT(2,get_solvers,status=status)
    END IF

    IF (.NOT. solvers_valid(get_solvers)) THEN
      IF (numpe == 1) THEN
        WRITE(*,*) "Solvers can be " // solvers_list()
      END IF
    END IF
  END FUNCTION get_solvers
  
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
    !*  This needs to be kept consistent with solvers_list().
    !*/

    LOGICAL                      :: solvers_valid
    CHARACTER(len=*), INTENT(IN) :: solvers

    IF (     solvers == parafem_solvers                                        &
        .OR. solvers == petsc_solvers   ) THEN
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
    !*  This needs to be kept consistent with solvers_valid().
    !*/

    CHARACTER(:),ALLOCATABLE :: solvers_list

    solvers_list = parafem_solvers                                             &
      // " or " // petsc_solvers
  END FUNCTION solvers_list
  
END MODULE choose_solvers
