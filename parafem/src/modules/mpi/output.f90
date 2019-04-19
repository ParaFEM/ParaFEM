MODULE OUTPUT

  !/****h* /output
  !*  NAME
  !*    MODULE: output
  !*  SYNOPSIS
  !*    Usage:      USE output
  !*  FUNCTION
  !*    Contains subroutines that write out the results. These subroutines are 
  !*    parallel and require MPI.
  !*    
  !*    Subroutine                   Purpose
  !*    WRITE_P121                   Writes out basic program data and timing info
  !*    WRITE_P123                   Writes out basic program data and timing info
  !*    WRITE_P125                   Writes out basic program data and timing info 
  !*    WRITE_P126                   Writes out basic program data and timing info  
  !*    WRITE_P129                   Writes out basic program data and timing info
  !*    WRITE_P1210                   
  !*    WRITE_XX1                    Writes out basic program data and timing info
  !*    WRITE_XX12                   Writes out basic program data and timing info
  !*    WRITE_NODAL_VARIABLE         Writes out results computed at the nodes
  !*    WRITE_NODAL_VARIABLE2
  !*    WRITE_NODAL_VARIABLE_BINARY  Write the values of a nodal variable to a file
  !*    WRITE_X_PP                   Write the values of a nodal variable to a file
  !*    JOB_NAME_ERROR               Writes error message if job_name is missing
  !*    GETFILENUMBER                Returns the next file number
  !*    DISMSH_ENSI                  ASCII Ensight Gold output for Paraview
  !*    DISMSH_ENSI_PB               Fortran binary Ensight Gold output for Paraview
  !*	DISMSH_ENSI_PB2	             C binary Ensight Gold output for Paraview
  !*    DISMSH_ENSI_PB3              Write the values of a nodal variable to a 
  !*                                 binary file
  !*    DISMSH_ENSI_PB2_INT          Write the values of a nodal variable to a 
  !*                                 binary file
  !*
  !*  AUTHOR
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

  SUBROUTINE WRITE_P121(fixed_freedoms,iters,job_name,loaded_nodes,neq,nn,    &
                        npes,nr,numpe,timest,tload)

  !/****f* output/write_p121
  !*  NAME
  !*    SUBROUTINE: write_p121
  !*  SYNOPSIS
  !*    Usage:      CALL write_p121(fixed_freedoms,iters,job_name,loaded_nodes,&
  !*                                neq,nn,npes,nr,numpe,timest,tload)
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and 
  !*    some performance data
  !*  INPUTS
  !*    The following scalar integers have the INTENT(IN) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    iters                  : Number of PCG iterations taken to solve problem
  !*    loaded_nodes           : Number of loaded_nodes
  !*    neq                    : Total number of equations in the mesh
  !*    nn                     : Number of nodes in the mesh
  !*    npes                   : Number of processors used in the simulations
  !*    nr                     : Number of restrained nodes in the mesh
  !*    numpe                  : Processor number
  !*
  !*    The following scalar real has the INTENT(IN) attribute:
  !*
  !*    tload                  : Total applied load
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : Job name used to name output file
  !*
  !*    The following dynamic real array has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*    Based on Smith I.M. and Griffiths D.V. "Programming the Finite Element
  !*    Method", Edition 4, Wiley, 2004.
  !*  CREATION DATE
  !*    02.03.2010
  !*  COPYRIGHT
  !*    (c) University of Manchester 2010-2015
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: numpe,npes,nn,nr,neq,iters
  INTEGER, INTENT(IN)            :: fixed_freedoms,loaded_nodes
  REAL(iwp), INTENT(IN)          :: timest(:),tload

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  
  CHARACTER(LEN=50)              :: fname
  INTEGER                        :: i          ! loop counter
 
  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".res"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

!------------------------------------------------------------------------------
! 2. Write basic details about the problem
!------------------------------------------------------------------------------

    WRITE(11,'(/A)')   "BASIC JOB DATA                                  "     
 
    WRITE(11,'(A,I12)')    "Number of processors used                   ",npes 
    WRITE(11,'(A,I12)')    "Number of nodes in the mesh                 ",nn
    WRITE(11,'(A,I12)')    "Number of nodes that were restrained        ",nr
    WRITE(11,'(A,I12)')    "Number of equations solved                  ",neq
    WRITE(11,'(A,I12)')    "Number of PCG iterations                    ",iters
    IF(loaded_nodes > 0) THEN
      WRITE(11,'(A,I12)')    "Number of loaded nodes                      ",   &
                              loaded_nodes 
      WRITE(11,'(A,E12.4)')  "Total load applied                          ",   &
                              tload
    END IF
    IF(fixed_freedoms > 0) THEN
      WRITE(11,'(A,I12)')    "Number of fixed displacements               ",   &
                              fixed_freedoms 
    END IF

!------------------------------------------------------------------------------
! 3. Output timing data
!------------------------------------------------------------------------------

    WRITE(11,'(/3A)')   "PROGRAM SECTION EXECUTION TIMES                  ",  &
                        "SECONDS  ", "%TOTAL    "
    WRITE(11,'(A,F12.6,F8.2)') "Setup                                       ",&
                           timest(2)-timest(1),                               &
                           ((timest(2)-timest(1))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read element steering array                 ",&
                           timest(3)-timest(2),                               &
                           ((timest(3)-timest(2))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Convert Abaqus to S&G node ordering         ",&
                           timest(4)-timest(3),                               &
                           ((timest(4)-timest(3))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read nodal coordinates                      ",&
                           timest(5)-timest(4),                               &
                           ((timest(5)-timest(4))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read restrained nodes                       ",&
                           timest(6)-timest(5),                               &
                           ((timest(6)-timest(5))/(timest(14)-timest(1)))*100                             
    WRITE(11,'(A,F12.6,F8.2)') "Compute steering array and neq              ",&
                           timest(7)-timest(6),                               &
                          ((timest(7)-timest(6))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute interprocessor communication tables ",&
                           timest(8)-timest(7),                               &
                          ((timest(8)-timest(7))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Allocate neq_pp arrays                      ",&
                           timest(9)-timest(8),                               &
                          ((timest(9)-timest(8))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute element stiffness matrices          ",&
                            timest(10)-timest(9),                              &
                          ((timest(10)-timest(9))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Build the preconditioner                    ",&
                           timest(11)-timest(10),                               &
                          ((timest(11)-timest(10))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Get starting r                              ",&
                           timest(12)-timest(11),                               &
                          ((timest(12)-timest(11))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Solve equations                             ",&
                           timest(13)-timest(12),                               &
                           ((timest(13)-timest(12))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Output results                              ",&
                           timest(14)-timest(13),                              &
                          ((timest(14)-timest(13))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A/)')  "Total execution time                        ",  &
                          timest(14)-timest(1),"  100.00"
    CLOSE(11)
    
  END IF
  
  RETURN

  END SUBROUTINE WRITE_P121

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE WRITE_P123(fixed_freedoms,iters,job_name,loaded_freedoms,neq,nn, &
                        npes,nr,numpe,timest,tload)

  !/****f* output/write_p123
  !*  NAME
  !*    SUBROUTINE: write_p123
  !*  SYNOPSIS
  !*    Usage:      CALL write_p123(fixed_freedoms,iters,job_name,            &
  !*                                loaded_freedoms,neq,nn,npes,nr,numpe,     &
  !*                                timest,tload)
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and 
  !*    some performance data
  !*  INPUTS
  !*    The following scalar integers have the INTENT(IN) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    iters                  : Number of PCG iterations taken to solve problem
  !*    loaded_nodes           : Number of loaded_nodes
  !*    neq                    : Total number of equations in the mesh
  !*    nn                     : Number of nodes in the mesh
  !*    npes                   : Number of processors used in the simulations
  !*    nr                     : Number of restrained nodes in the mesh
  !*    numpe                  : Processor number
  !*
  !*    The following scalar real has the INTENT(IN) attribute:
  !*
  !*    tload                  : Total applied load
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : Job name used to name output file
  !*
  !*    The following dynamic real array has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*    Based on Smith I.M. and Griffiths D.V. "Programming the Finite Element
  !*    Method", Edition 4, Wiley, 2004.
  !*  CREATION DATE
  !*    13.03.2012
  !*  COPYRIGHT
  !*    (c) University of Manchester 2012
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: numpe,npes,nn,nr,neq,iters
  INTEGER, INTENT(IN)            :: fixed_freedoms,loaded_freedoms
  REAL(iwp), INTENT(IN)          :: timest(:),tload

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  
  CHARACTER(LEN=50)              :: fname
  INTEGER                        :: i          ! loop counter
 
  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".res"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

!------------------------------------------------------------------------------
! 2. Write basic details about the problem
!------------------------------------------------------------------------------

    WRITE(11,'(/A)')   "BASIC JOB DATA                                  "     
 
    WRITE(11,'(A,I12)')    "Number of processors used                   ",npes 
    WRITE(11,'(A,I12)')    "Number of nodes in the mesh                 ",nn
    WRITE(11,'(A,I12)')    "Number of nodes that were restrained        ",nr
    WRITE(11,'(A,I12)')    "Number of equations solved                  ",neq
    WRITE(11,'(A,I12)')    "Number of PCG iterations                    ",iters
    IF(loaded_freedoms > 0) THEN
      WRITE(11,'(A,I12)')    "Number of loaded freedoms                   ",   &
                              loaded_freedoms 
      WRITE(11,'(A,E12.4)')  "Total load applied                          ",   &
                              tload
    END IF
    IF(fixed_freedoms > 0) THEN
      WRITE(11,'(A,I12)')    "Number of fixed displacements               ",   &
                              fixed_freedoms 
    END IF

!------------------------------------------------------------------------------
! 3. Output timing data
!------------------------------------------------------------------------------

    WRITE(11,'(/3A)')   "PROGRAM SECTION EXECUTION TIMES                  ",  &
                        "SECONDS  ", "%TOTAL    "
    WRITE(11,'(A,F12.6,F8.2)') "Setup                                       ",&
                           timest(2)-timest(1),                               &
                           ((timest(2)-timest(1))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read element steering array                 ",&
                           timest(3)-timest(2),                               &
                           ((timest(3)-timest(2))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Convert Abaqus to S&G node ordering         ",&
                           timest(4)-timest(3),                               &
                           ((timest(4)-timest(3))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read nodal coordinates                      ",&
                           timest(5)-timest(4),                               &
                           ((timest(5)-timest(4))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read restrained nodes                       ",&
                           timest(6)-timest(5),                               &
                           ((timest(6)-timest(5))/(timest(14)-timest(1)))*100                             
    WRITE(11,'(A,F12.6,F8.2)') "Compute steering array and neq              ",&
                           timest(7)-timest(6),                               &
                          ((timest(7)-timest(6))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute interprocessor communication tables ",&
                           timest(8)-timest(7),                               &
                          ((timest(8)-timest(7))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Allocate neq_pp arrays                      ",&
                           timest(9)-timest(8),                               &
                          ((timest(9)-timest(8))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute element stiffness matrices          ",&
                            timest(10)-timest(9),                              &
                          ((timest(10)-timest(9))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Build the preconditioner                    ",&
                           timest(11)-timest(10),                               &
                          ((timest(11)-timest(10))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Get starting r                              ",&
                           timest(12)-timest(11),                               &
                          ((timest(12)-timest(11))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Solve equations                             ",&
                           timest(13)-timest(12),                               &
                           ((timest(13)-timest(12))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Output results                              ",&
                           timest(14)-timest(13),                              &
                          ((timest(14)-timest(13))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A/)')  "Total execution time                        ",  &
                          timest(14)-timest(1),"  100.00"
    CLOSE(11)
    
  END IF
  
  RETURN

  END SUBROUTINE WRITE_P123

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE WRITE_P125(fixed_freedoms,nstep,job_name,loaded_freedoms,neq,nn, &
                        npes,nr,numpe,timest)

  !/****f* output/write_p125
  !*  NAME
  !*    SUBROUTINE: write_p125
  !*  SYNOPSIS
  !*    Usage:      CALL write_p125(fixed_freedoms,iters,job_name,            &
  !*                                loaded_freedoms,neq,nn,npes,nr,numpe,     &
  !*                                timest)
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and 
  !*    some performance data
  !*  INPUTS
  !*    The following scalar integers have the INTENT(IN) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    nstep                  : Number of time steps carried out
  !*    loaded_nodes           : Number of loaded_nodes
  !*    neq                    : Total number of equations in the mesh
  !*    nn                     : Number of nodes in the mesh
  !*    npes                   : Number of processors used in the simulations
  !*    nr                     : Number of restrained nodes in the mesh
  !*    numpe                  : Processor number
  !*
  !*  ! The following scalar real has the INTENT(IN) attribute:
  !*  !
  !*  ! tload                  : Total applied load
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : Job name used to name output file
  !*
  !*    The following dynamic real array has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*    Based on Smith I.M. and Griffiths D.V. "Programming the Finite Element
  !*    Method", Edition 4, Wiley, 2004.
  !*  CREATION DATE
  !*    05.12.2012
  !*  COPYRIGHT
  !*    (c) University of Manchester 2012
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: numpe,npes,nn,nr,neq,nstep
  INTEGER, INTENT(IN)            :: fixed_freedoms,loaded_freedoms
  REAL(iwp), INTENT(IN)          :: timest(:)

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  
  CHARACTER(LEN=50)              :: fname
  INTEGER                        :: i          ! loop counter
 
  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".per"
    OPEN(12,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

!------------------------------------------------------------------------------
! 2. Write basic details about the problem
!------------------------------------------------------------------------------

    WRITE(12,'(/A)')   "BASIC JOB DATA                                  "     
 
    WRITE(12,'(A,I12)')    "Number of processors used                   ",npes 
    WRITE(12,'(A,I12)')    "Number of nodes in the mesh                 ",nn
    WRITE(12,'(A,I12)')    "Number of nodes that were restrained        ",nr
    WRITE(12,'(A,I12)')    "Number of equations solved                  ",neq
    WRITE(12,'(A,I12)')    "Number of time steps                        ",nstep
    IF(loaded_freedoms > 0) THEN
      WRITE(12,'(A,I12)')    "Number of loaded freedoms                   ",   &
                              loaded_freedoms 
!     WRITE(12,'(A,E12.4)')  "Total load applied                          ",   &
!                             tload
    END IF
    IF(fixed_freedoms > 0) THEN
      WRITE(12,'(A,I12)')    "Number of fixed displacements               ",   &
                              fixed_freedoms 
    END IF

!------------------------------------------------------------------------------
! 3. Output timing data
!------------------------------------------------------------------------------

    WRITE(12,'(/3A)')   "PROGRAM SECTION EXECUTION TIMES                  ",  &
                        "SECONDS  ", "%TOTAL    "
    WRITE(12,'(A,F12.6,F8.2)') "Setup                                       ",&
                           timest(2)-timest(1),                               &
                           ((timest(2)-timest(1))/(timest(12)-timest(1)))*100  
    WRITE(12,'(A,F12.6,F8.2)') "Read element steering array                 ",&
                           timest(3)-timest(2),                               &
                           ((timest(3)-timest(2))/(timest(12)-timest(1)))*100  
    WRITE(12,'(A,F12.6,F8.2)') "Convert Abaqus to S&G node ordering         ",&
                           timest(4)-timest(3),                               &
                           ((timest(4)-timest(3))/(timest(12)-timest(1)))*100  
    WRITE(12,'(A,F12.6,F8.2)') "Read nodal coordinates                      ",&
                           timest(5)-timest(4),                               &
                           ((timest(5)-timest(4))/(timest(12)-timest(1)))*100  
    WRITE(12,'(A,F12.6,F8.2)') "Read restrained nodes                       ",&
                           timest(6)-timest(5),                               &
                           ((timest(6)-timest(5))/(timest(12)-timest(1)))*100                             
    WRITE(12,'(A,F12.6,F8.2)') "Compute steering array and neq              ",&
                           timest(7)-timest(6),                               &
                          ((timest(7)-timest(6))/(timest(12)-timest(1)))*100  
    WRITE(12,'(A,F12.6,F8.2)') "Compute interprocessor communication tables ",&
                           timest(8)-timest(7),                               &
                          ((timest(8)-timest(7))/(timest(12)-timest(1)))*100  
    WRITE(12,'(A,F12.6,F8.2)') "Allocate neq_pp arrays                      ",&
                           timest(9)-timest(8),                               &
                          ((timest(9)-timest(8))/(timest(12)-timest(1)))*100  
    WRITE(12,'(A,F12.6,F8.2)') "Compute element stiffness matrices          ",&
                            timest(10)-timest(9),                             &
                          ((timest(10)-timest(9))/(timest(12)-timest(1)))*100  
    WRITE(12,'(A,F12.6,F8.2)') "Invert the mass matrix                      ",&
                           timest(11)-timest(10),                             &
                          ((timest(11)-timest(10))/(timest(12)-timest(1)))*100  
    WRITE(12,'(A,F12.6,F8.2)') "Time stepping recursion                     ",&
                           timest(12)-timest(11),                             &
                          ((timest(12)-timest(11))/(timest(12)-timest(1)))*100  
!    WRITE(12,'(A,F12.6,F8.2)')"Output results                              ",&
!                           timest(14)-timest(13),                            &
!                          ((timest(14)-timest(13))/(timest(14)-timest(1)))*100  
    WRITE(12,'(A,F12.6,A/)')   "Total execution time                        ",&
                          timest(12)-timest(1),"  100.00"
    CLOSE(12)
    
  END IF
  
  RETURN

  END SUBROUTINE WRITE_P125
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE WRITE_P126(bicg,bicg_pp,cj_tot,iters,job_name,neq,nn,npes,nr,    &
                        numpe,timest)

  !/****f* output/write_p126
  !*  NAME
  !*    SUBROUTINE: write_p126
  !*  SYNOPSIS
  !*    Usage:      CALL write_p126(bicg,bicg_pp,cj_tot,iters,job_name,neq,   &
  !*                                nn,npes,nr,numpe,timest)
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and 
  !*    some performance data
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    cj_tot                 : Integer
  !*                           : Total number of BiCGSTAB(l) iterations taken 
  !*                           : to solve problem
  !*
  !*    iters                  : Integer
  !*                           : Number of external iterations required
  !*
  !*    job_name               : Character
  !*                           : Job name used to name output file
  !*
  !*    neq                    : Integer
  !*                           : Total number of equations in the mesh
  !*
  !*    nn                     : Integer
  !*                           : Number of nodes in the mesh
  !*
  !*    npes                   : Integer
  !*                           : Number of processors used in the simulations
  !*
  !*    nr                     : Integer
  !*                           : Number of restrained nodes in the mesh
  !*
  !*    numpe                  : Integer
  !*                           : Processor number
  !*
  !*    timest(:)              : Real array
  !*                           : Holds timing information
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    26.07.2010
  !*  COPYRIGHT
  !*    (c) University of Manchester 2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: numpe,npes,nn,nr,neq,cj_tot,iters,bicg(:)
  REAL(iwp), INTENT(IN)          :: timest(:),bicg_pp(:)

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  
  CHARACTER(LEN=50)              :: fname
  INTEGER                        :: i          ! loop counter
 
  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".res"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

!------------------------------------------------------------------------------
! 2. Write basic details about the problem
!------------------------------------------------------------------------------

    WRITE(11,'(/A,I12)') "Number of processors used:                  ",npes 
    WRITE(11,'(A,I12)')  "Number of nodes in the mesh:                ",nn
    WRITE(11,'(A,I12)')  "Number of nodes that were restrained:       ",nr
    WRITE(11,'(A,I12/)') "Number of equations solved:                 ",neq

!------------------------------------------------------------------------------
! 3. Write output for the BiCGSTAB(l) iterations 
!------------------------------------------------------------------------------

    WRITE(11,'(A,I12/)') "Total BiCGSTAB(l) iterations:               ",cj_tot

    WRITE(11,'(2A)')  "                  iterations        norm"

    DO i = 1, iters
      WRITE(11,'(A,I12,I12,E12.4)') "    ",i,bicg(i),bicg_pp(i)
    END DO
    
!------------------------------------------------------------------------------
! 4. Output timing data
!------------------------------------------------------------------------------
 
    WRITE(11,'(/3A)')  "Program section execution times                   ",  &
                        "Seconds  ", "%Total    "
    WRITE(11,'(A,F12.6,F8.2)')"Setup                                        ",&
                          timest(2)-timest(1),                                &
                          ((timest(2)-timest(1))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)')"Compute coordinates and steering array       ",&
                          timest(3)-timest(2),                                &
                          ((timest(3)-timest(2))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)')"Compute interprocessor communication tables  ",&
                          timest(4)-timest(3),                                &
                          ((timest(4)-timest(3))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)')"Allocate neq_pp arrays                       ",&
                          timest(5)-timest(4),                                &
                          ((timest(5)-timest(4))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)')"Organize fixed nodes                         ",&
                          timest(6)-timest(5),                                &
                          ((timest(6)-timest(5))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)')"Element stiffness integration                ",&
                          timest(8),                                          &
                          ((timest(8))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)')"Build the preconditioner                     ",&
                          timest(9),                                          &
                          ((timest(9))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)')"Solve equations                              ",&
                          timest(10),                                         &
                          ((timest(10))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)')"Output results                               ",&
                          timest(12)-timest(11),                              &
                          ((timest(12)-timest(11))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A)') "Total execution time                         ",  &
                           timest(12)-timest(1),"  100.00"
    CLOSE(11)
    
  END IF
  
  RETURN
  END SUBROUTINE WRITE_P126
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE WRITE_P129(ans,job_name,neq,nn,npes,nr,numpe,timest,tload,       &
                        ttliters)

  !/****f* output/write_p129
  !*  NAME
  !*    SUBROUTINE: write_p129
  !*  SYNOPSIS
  !*    Usage:      CALL write_p129(ans,job_name,neq,nn,npes,nr,numpe,        &
  !*                                timest,tload,ttliters))
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and 
  !*    some performance data
  !*  INPUTS
  !*    The following scalar character argument has the INTENT(IN) attribute:
  !*
  !*    job_name               : Job name used to name output file
  !*
  !*    The following scalar integer arguments have the INTENT(IN) attribute:
  !*
  !*    neq                    : Number of equations
  !*    nn                     : Number of nodes in the mesh
  !*    npes                   : Number of processors used in the simulations
  !*    nr                     : Number of restrained nodes
  !*    numpe                  : Processor ID number
  !*
  !*    The following scalar real argument has the INTENT(IN) attribute:
  !*
  !*    tload                  : Total applied load
  !*
  !*    The following real array arguments have the INTENT(IN) attribute:
  !*
  !*    ans(:)                 : Holds a single displacement value at the 
  !*                           : equation number "it" for each time step
  !*    timest(:)              : Holds timing information
  !*    ttliters(:)            : Records the number of iterations required by
  !*                           : PCG for each time step
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    26.02.2010
  !*  COPYRIGHT
  !*    (c) University of Manchester 2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: numpe,npes,nn,nr,neq,ttliters(:)
  REAL(iwp), INTENT(IN)          :: ans(:,:),timest(:),tload

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  
  CHARACTER(LEN=50)              :: fname
  INTEGER                        :: i          ! loop counter
 
  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".res"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

!------------------------------------------------------------------------------
! 3. Write basic details about the problem
!------------------------------------------------------------------------------

    WRITE(11,'(A,I12)')    "Number of processors used:                  ",npes 
    WRITE(11,'(A,I12)')    "Number of nodes in the mesh:                ",nn
    WRITE(11,'(A,I12)')    "Number of nodes that were restrained:       ",nr
    WRITE(11,'(A,I12)')    "Number of equations solved:                 ",neq
    WRITE(11,'(A,E12.4)')  "Total load applied:                         ",tload

!------------------------------------------------------------------------------
! 4. Write output at the equation specified by the user, for each time step, 
!    and report the number of iterations required
!------------------------------------------------------------------------------

    WRITE(11,'(/A)')"        Time t  cos(omega*t)  Displacement    Iterations"
    DO i = 1,UBOUND(ans,2)
      WRITE(11,'(3E14.4,I14)') ans(1,i), ans(2,i), ans(3,i), ttliters(i)
    END DO

!------------------------------------------------------------------------------
! 5. Output timing data
!------------------------------------------------------------------------------

    WRITE(11,'(/3A)')   "Program section execution times                   ",  &
                        "Seconds  ", "%Total    "
    WRITE(11,'(A,F12.6,F8.2)') "Setup                                        ",&
                           timest(2)-timest(1),                                &
                          ((timest(2)-timest(1))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute steering array                       ",&
                           timest(3)-timest(2),                                &
                          ((timest(3)-timest(2))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute interprocessor communication tables  ",&
                           timest(4)-timest(3),                                &
                          ((timest(4)-timest(3))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Allocate neq_pp arrays                       ",&
                           timest(5)-timest(4),                                &
                          ((timest(5)-timest(4))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Element stiffness integration                ",&
                           timest(6)-timest(5),                                &
                          ((timest(6)-timest(5))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Build the diagonal preconditioner            ",&
                           timest(7)-timest(6),                                &
                          ((timest(7)-timest(6))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read the applied forces                      ",&
                           timest(8)-timest(7),                                &
                          ((timest(8)-timest(7))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Solve the equations                          ",&
                           timest(10),                                         &
                          ((timest(10))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Output the results                           ",&
                          timest(11),                                          &
                         ((timest(11))/(timest(12)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A)')  "Total execution time                         ",  &
                          timest(12)-timest(1),"  100.00"
    CLOSE(11)
    
  END IF
  
  RETURN
  END SUBROUTINE WRITE_P129

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE WRITE_P1210(job_name,meshgen,neq,nn,npes,nr,numpe,timest,tload)

  !/****f* output/write_p1210
  !*  NAME
  !*    SUBROUTINE: write_p1210
  !*  SYNOPSIS
  !*    Usage:      CALL write_p1210(job_name,meshgen,neq,nn,npes,nr,numpe,    &
  !*                                 timest,tload)
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and 
  !*    some performance data
  !*  INPUTS
  !*    The following scalar character argument has the INTENT(IN) attribute:
  !*
  !*    job_name               : Job name used to name output file
  !*
  !*    The following scalar integer arguments have the INTENT(IN) attribute:
  !*
  !*    meshgen                : Mesh numbering scheme
  !*                           : 1 = Smith and Griffiths numbering scheme
  !*                           : 2 = Abaqus numbering scheme
  !*    neq                    : Number of equations
  !*    nn                     : Number of nodes in the mesh
  !*    npes                   : Number of processors used in the simulations
  !*    nr                     : Number of restrained nodes
  !*    numpe                  : Processor ID number
  !*
  !*    The following scalar real argument has the INTENT(IN) attribute:
  !*
  !*    tload                  : Total applied load
  !*
  !*    The following real array argument has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*  CREATION DATE
  !*    26.02.2010
  !*  COPYRIGHT
  !*    (c) University of Manchester 2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: numpe,npes,nn,nr,neq,meshgen
  REAL(iwp), INTENT(IN)          :: timest(:),tload

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  
  CHARACTER(LEN=50)              :: fname
  INTEGER                        :: i          ! loop counter
 
  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".res"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

!------------------------------------------------------------------------------
! 3. Write basic details about the problem
!------------------------------------------------------------------------------

    WRITE(11,'(A,I12)')    "Number of processors used:                  ",npes 
    WRITE(11,'(A,I12)')    "Number of nodes in the mesh:                ",nn
    WRITE(11,'(A,I12)')    "Number of nodes that were restrained:       ",nr
    WRITE(11,'(A,I12)')    "Number of equations solved:                 ",neq
    WRITE(11,'(A,E12.4)')  "Total load applied:                         ",tload

!------------------------------------------------------------------------------
! 4. Output timing data
!------------------------------------------------------------------------------

    IF(meshgen == 1) THEN
      WRITE(11,'(/A)') "The mesh used the Smith and Griffiths numbering scheme"
      WRITE(11,'(A)') "Note: the ParaFEM-Viewer expects the Abaqus format"
    ELSE
      WRITE(11,'(/A)') "The mesh used the Abaqus numbering scheme"
    END IF
    WRITE(11,'(/3A)')   "Program section execution times                   ",  &
                        "Seconds  ", "%Total    "
    WRITE(11,'(A,F12.6,F8.2)') "Setup                                        ",&
                           timest(2)-timest(1),                                &
                          ((timest(2)-timest(1))/(timest(11)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute steering array                       ",&
                           timest(3)-timest(2),                                &
                          ((timest(3)-timest(2))/(timest(11)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute interprocessor communication tables  ",&
                           timest(4)-timest(3),                                &
                          ((timest(4)-timest(3))/(timest(11)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Allocate neq_pp arrays                       ",&
                           timest(5)-timest(4),                                &
                          ((timest(5)-timest(4))/(timest(11)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Calculate diagonal mass matrix               ",&
                           timest(6)-timest(5),                                &
                          ((timest(6)-timest(5))/(timest(11)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read applied forces + assign to equations    ",&
                           timest(7)-timest(6),                                &
                          ((timest(7)-timest(6))/(timest(11)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Element stress-strain relationship           ",&
                          timest(9),                                           &
                         ((timest(9))/(timest(11)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Output the results                           ",&
                          timest(10),                                           &
                         ((timest(10))/(timest(11)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A)')    "Total execution time                         ",&
                          timest(11)-timest(1),"  100.00"
    CLOSE(11)

  END IF 
  
  RETURN
  END SUBROUTINE WRITE_P1210
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE WRITE_XX1(ewt,fixed_freedoms,iters,job_name,loaded_nodes,neq,nn, &
                       npes,nr,numpe,timest,tload)

  !/****f* output/write_xx1
  !*  NAME
  !*    SUBROUTINE: write_xx1
  !*  SYNOPSIS
  !*    Usage:      CALL write_xx1(ewt,fixed_freedoms,iters,job_name,         &
  !*                               loaded_nodes,neq,nn,npes,nr,numpe,timest,  &
  !*                               tload)
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and 
  !*    some performance data
  !*  INPUTS
  !*    The following scalar integers have the INTENT(IN) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    iters                  : Number of PCG iterations taken to solve problem
  !*    loaded_nodes           : Number of loaded_nodes
  !*    neq                    : Total number of equations in the mesh
  !*    nn                     : Number of nodes in the mesh
  !*    npes                   : Number of processors used in the simulations
  !*    nr                     : Number of restrained nodes in the mesh
  !*    numpe                  : Processor number
  !*
  !*    The following scalar real has the INTENT(IN) attribute:
  !*
  !*    tload                  : Total applied load
  !*    ewt                    : Total strain energy
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : Job name used to name output file
  !*
  !*    The following dynamic real array has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*    Based on Smith I.M. and Griffiths D.V. "Programming the Finite Element
  !*    Method", Edition 4, Wiley, 2004.
  !*  CREATION DATE
  !*    24.07.2014
  !*  COPYRIGHT
  !*    (c) University of Manchester 2014-2015
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: numpe,npes,nn,nr,neq,iters
  INTEGER, INTENT(IN)            :: fixed_freedoms,loaded_nodes
  REAL(iwp), INTENT(IN)          :: timest(:),tload,ewt

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  
  CHARACTER(LEN=50)              :: fname
  INTEGER                        :: i          ! loop counter
 
! IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".res"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

!------------------------------------------------------------------------------
! 2. Write basic details about the problem
!------------------------------------------------------------------------------

    WRITE(11,'(/A)')   "BASIC JOB DATA                                  "     
 
    WRITE(11,'(A,I12)')    "Number of processors used                   ",npes 
    WRITE(11,'(A,I12)')    "Number of nodes in the mesh                 ",nn
    WRITE(11,'(A,I12)')    "Number of nodes that were restrained        ",nr
    WRITE(11,'(A,I12)')    "Number of equations solved                  ",neq
    WRITE(11,'(A,I12)')    "Number of PCG iterations                    ",iters
    IF(loaded_nodes > 0) THEN
      WRITE(11,'(A,I12)')    "Number of loaded nodes                      ",   &
                              loaded_nodes 
      WRITE(11,'(A,E12.4)')  "Total load applied                          ",   &
                              tload
    END IF
    IF(fixed_freedoms > 0) THEN
      WRITE(11,'(A,I12)')    "Number of fixed displacements               ",   &
                              fixed_freedoms 
    END IF
    WRITE(11,'(A,E12.4)')  "Total strain energy                         ",ewt

!------------------------------------------------------------------------------
! 3. Output timing data
!------------------------------------------------------------------------------

    WRITE(11,'(/3A)')   "PROGRAM SECTION EXECUTION TIMES                  ",  &
                        "SECONDS  ", "%TOTAL    "
    WRITE(11,'(A,F12.6,F8.2)') "Setup                                       ",&
                           timest(2)-timest(1),                               &
                           ((timest(2)-timest(1))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read element steering array                 ",&
                           timest(3)-timest(2),                               &
                           ((timest(3)-timest(2))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Convert Abaqus to S&G node ordering         ",&
                           timest(4)-timest(3),                               &
                           ((timest(4)-timest(3))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read nodal coordinates                      ",&
                           timest(5)-timest(4),                               &
                           ((timest(5)-timest(4))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read restrained nodes                       ",&
                           timest(6)-timest(5),                               &
                           ((timest(6)-timest(5))/(timest(17)-timest(1)))*100                             
    WRITE(11,'(A,F12.6,F8.2)') "Read material types                         ",&
                           timest(7)-timest(6),                               &
                           ((timest(7)-timest(6))/(timest(17)-timest(1)))*100                             
    WRITE(11,'(A,F12.6,F8.2)') "Compute steering array and neq              ",&
                           timest(8)-timest(7),                               &
                          ((timest(8)-timest(7))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute interprocessor communication tables ",&
                           timest(9)-timest(8),                               &
                          ((timest(9)-timest(8))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Allocate neq_pp arrays                      ",&
                           timest(10)-timest(9),                              &
                          ((timest(10)-timest(9))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute element stiffness matrices          ",&
                            timest(11)-timest(10),                            &
                          ((timest(11)-timest(10))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Build the preconditioner                    ",&
                           timest(12)-timest(11),                             &
                          ((timest(12)-timest(11))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Get starting r                              ",&
                           timest(13)-timest(12),                             &
                          ((timest(13)-timest(12))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Solve equations                             ",&
                           timest(14)-timest(13),                             &
                           ((timest(14)-timest(13))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Calculate nodes_pp                          ",&
                           timest(15)-timest(14),                             &
                          ((timest(15)-timest(14))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Scatter the nodes                           ",&
                           timest(16)-timest(15),                             &
                          ((timest(16)-timest(15))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Output displacements                        ",&
                           timest(17)-timest(16),                             &
                          ((timest(17)-timest(16))/(timest(17)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute strain energy                       ",&
                           timest(19)-timest(18),                             &
                          ((timest(19)-timest(18))/(timest(19)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A/)')   "Total execution time                        ",&
                          timest(19)-timest(1),"  100.00"
    CLOSE(11)
    
! END IF
  
  RETURN

  END SUBROUTINE WRITE_XX1

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE WRITE_XX4(fixed_freedoms,iters,job_name,loaded_nodes,            &
                       mult_method,neq,nn,npes,nr,numpe,timest,tload,         &
                       vox_storkm)

  !/****f* output/write_xx4
  !*  NAME
  !*    SUBROUTINE: write_xx4
  !*  SYNOPSIS
  !*    Usage:      CALL write_xx4(fixed_freedoms,iters,job_name,             &
  !*                               loaded_nodes,mult_method,neq,nn,npes,nr,   &
  !*                               numpe,timest,tload,vox_storkm)
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and 
  !*    some performance data
  !*  INPUTS
  !*    The following scalar integers have the INTENT(IN) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    iters                  : Number of PCG iterations taken to solve problem
  !*    loaded_nodes           : Number of loaded_nodes
  !*    mult_method            : Method of matrix-vector multiplication
  !*    neq                    : Total number of equations in the mesh
  !*    nn                     : Number of nodes in the mesh
  !*    npes                   : Number of processors used in the simulations
  !*    nr                     : Number of restrained nodes in the mesh
  !*    numpe                  : Processor number
  !*
  !*    The following scalar real has the INTENT(IN) attribute:
  !*
  !*    tload                  : Total applied load
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : Job name used to name output file
  !*
  !*    The following scalar logical has the INTENT(IN) attribute:
  !*
  !*    vox_storkm             : True if using a voxel mesh
  !*
  !*    The following dynamic real array has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*    Based on Smith I.M. and Griffiths D.V. "Programming the Finite Element
  !*    Method", Edition 4, Wiley, 2004.
  !*  CREATION DATE
  !*    06.02.2019
  !*  COPYRIGHT
  !*    (c) University of Manchester 2019
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: numpe,npes,nn,nr,neq,iters
  INTEGER, INTENT(IN)            :: fixed_freedoms,loaded_nodes,mult_method
  REAL(iwp), INTENT(IN)          :: timest(:),tload
  LOGICAL, INTENT(IN)            :: vox_storkm

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  
  CHARACTER(LEN=50)              :: fname
  CHARACTER(LEN=60)              :: date
  CHARACTER(LEN=24)              :: fdate      ! instrinsic procedure
  CHARACTER(LEN=80)              :: hostname
  INTEGER                        :: i          ! loop counter
  INTEGER*4                      :: status
  INTEGER*4                      :: hostnm     ! intrinsic procedure
 
  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".res"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

!------------------------------------------------------------------------------
! 2. Retrieve information about the environment
!------------------------------------------------------------------------------

!   CALL get_environment_variable("VARIABLE",variable)
    
    WRITE(11,'(A)') fdate() 

    status = hostnm(hostname)
 
    IF(status==0) THEN
      WRITE(11,'(/2A)')  "This job run on host ", hostname 
    ELSE
      WRITE(11,'(/A,A)') "This job run on host ", hostname 
      WRITE(11,'(A,I2)') "Routine HOSTNM returned error code", status 
    END IF

!------------------------------------------------------------------------------
! 3. Report which case was used
!------------------------------------------------------------------------------

    IF(mult_method==0) THEN
      IF (vox_storkm) THEN
        WRITE(11,'(A)') "Test case 0: Matmul on CPU for voxel meshes"
      ELSE
        WRITE(11,'(A)') "Test case 0: Matmul on CPU for standard meshes"
      END IF
    END IF
    IF(mult_method==1) WRITE(11,'(A)') "Test case 1: Matmul on GPU naive"
    IF(mult_method==2) WRITE(11,'(A)') "Test case 2: Matmul on GPU local memory"
    IF(mult_method==3) WRITE(11,'(A)') "Test case 3: Matmul on GPU AMD clBlas"
    IF(mult_method==4) THEN
      WRITE(11,'(A)') "Test case 4: Matmul on GPU AMD clBLAS Asynchronous"
    END IF
 
!------------------------------------------------------------------------------
! 4. Write basic details about the problem
!------------------------------------------------------------------------------

    WRITE(11,'(/A)')       "BASIC JOB DATA                                  "  
    WRITE(11,'(A,I12)')    "Number of processors used                   ",npes 
    WRITE(11,'(A,I12)')    "Number of nodes in the mesh                 ",nn
    WRITE(11,'(A,I12)')    "Number of nodes that were restrained        ",nr
    WRITE(11,'(A,I12)')    "Number of equations solved                  ",neq
    WRITE(11,'(A,I12)')    "Number of PCG iterations                    ",iters
    IF(loaded_nodes > 0) THEN
      WRITE(11,'(A,I12)')    "Number of loaded nodes                      ",   &
                              loaded_nodes 
      WRITE(11,'(A,E12.4)')  "Total load applied                          ",   &
                              tload
    END IF
    IF(fixed_freedoms > 0) THEN
      WRITE(11,'(A,I12)')    "Number of fixed displacements               ",   &
                              fixed_freedoms 
    END IF

!------------------------------------------------------------------------------
! 5. Output timing data
!------------------------------------------------------------------------------

    WRITE(11,'(/3A)')   "PROGRAM SECTION EXECUTION TIMES                  ",  &
                        "SECONDS  ", "%TOTAL    "
    WRITE(11,'(A,F12.6,F8.2)') "Setup                                       ",&
                           timest(2)-timest(1),                               &
                           ((timest(2)-timest(1))/(timest(24)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read element steering array                 ",&
                           timest(3)-timest(2),                               &
                           ((timest(3)-timest(2))/(timest(24)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Convert Abaqus to S&G node ordering         ",&
                           timest(4)-timest(3),                               &
                           ((timest(4)-timest(3))/(timest(24)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read nodal coordinates                      ",&
                           timest(5)-timest(4),                               &
                           ((timest(5)-timest(4))/(timest(24)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read restrained nodes                       ",&
                           timest(6)-timest(5),                               &
                           ((timest(6)-timest(5))/(timest(24)-timest(1)))*100                             
    WRITE(11,'(A,F12.6,F8.2)') "Compute steering array and neq              ",&
                           timest(7)-timest(6),                               &
                          ((timest(7)-timest(6))/(timest(24)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute interprocessor communication tables ",&
                           timest(8)-timest(7),                               &
                          ((timest(8)-timest(7))/(timest(24)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Allocate neq_pp arrays                      ",&
                           timest(9)-timest(8),                               &
                          ((timest(9)-timest(8))/(timest(24)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute element stiffness matrices          ",&
                            timest(10)-timest(9),                             &
                          ((timest(10)-timest(9))/(timest(24)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Build the preconditioner                    ",&
                           timest(11)-timest(10),                             &
                          ((timest(11)-timest(10))/(timest(24)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Get starting r                              ",&
                           timest(12)-timest(11),                             &
                          ((timest(12)-timest(11))/(timest(24)-timest(1)))*100  

    IF(mult_method==1 .OR. mult_method==2) THEN
      WRITE(11,'(A,F12.6,F8.2)')                                              &
        "Initiallise GPU                             ",                       &
         timest(14)-timest(13),                                               &
       ((timest(14)-timest(13))/(timest(24)-timest(1)))*100  
      WRITE(11,'(A,F12.6,F8.2)')                                              &
        "Copy KM data to GPU                         ",                       &
         timest(15)-timest(14),                                               &
       ((timest(15)-timest(14))/(timest(24)-timest(1)))*100  
      WRITE(11,'(A,F12.6,F8.2)')                                              &
        "Transfer vectors to/from GPU                ",                       &
         timest(19),                                                          &
        (timest(19))/(timest(24)-timest(1))*100  
      WRITE(11,'(A,F12.6,F8.2)')                                              &
        "Matrix-vector multiplication on GPU         ",                       &
         timest(21),                                                          &
        (timest(21))/(timest(24)-timest(1))*100  
    ELSE
      WRITE(11,'(A,F12.6,F8.2)')                                              &
        "Matrix-vector multiplication on CPU         ",                       &
         timest(17),                                                          &
        (timest(17))/(timest(24)-timest(1))*100  
    END IF
    WRITE(11,'(A,F12.6,F8.2)') "Output results                              ",&
                           timest(24)-timest(23),                             &
                          ((timest(24)-timest(23))/(timest(24)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)')                                                &
        "TOTAL: Solve equations                      ",                       &
         timest(23)-timest(12),                                               &
       ((timest(23)-timest(12))/(timest(24)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A/)')  "TOTAL: Program execution time               ", &
                          timest(24)-timest(1),"  100.00"
    CLOSE(11)
    
  END IF
  
  RETURN

  END SUBROUTINE WRITE_XX4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE WRITE_XX10(fixed_freedoms,iters,job_name,loaded_nodes,neq,nn,    &
                        npes,nr,numpe,timest,tload)

  !/****f* output/write_xx10
  !*  NAME
  !*    SUBROUTINE: write_xx10
  !*  SYNOPSIS
  !*    Usage:      CALL write_xx10(fixed_freedoms,iters,job_name,loaded_nodes,&
  !*                                neq,nn,npes,nr,numpe,timest,tload)
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and 
  !*    some performance data
  !*  INPUTS
  !*    The following scalar integers have the INTENT(IN) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    iters                  : Number of PCG iterations taken to solve problem
  !*    loaded_nodes           : Number of loaded_nodes
  !*    neq                    : Total number of equations in the mesh
  !*    nn                     : Number of nodes in the mesh
  !*    npes                   : Number of processors used in the simulations
  !*    nr                     : Number of restrained nodes in the mesh
  !*    numpe                  : Processor number
  !*
  !*    The following scalar real has the INTENT(IN) attribute:
  !*
  !*    tload                  : Total applied load
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : Job name used to name output file
  !*
  !*    The following dynamic real array has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*    Based on Smith I.M. and Griffiths D.V. "Programming the Finite Element
  !*    Method", Edition 4, Wiley, 2004.
  !*  CREATION DATE
  !*    07.02.2019
  !*  COPYRIGHT
  !*    (c) University of Manchester 2019
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: numpe,npes,nn,nr,neq,iters
  INTEGER, INTENT(IN)            :: fixed_freedoms,loaded_nodes
  REAL(iwp), INTENT(IN)          :: timest(:),tload

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  
  CHARACTER(LEN=50)              :: fname
  CHARACTER(LEN=60)              :: date
  CHARACTER(LEN=60)              :: hostname
  INTEGER                        :: i          ! loop counter
 
  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".res"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

!------------------------------------------------------------------------------
! 2. Retrieve information about the environment
!
!    Run the following command before running the job to get the date and
!    time the job was executed or submitted.
!
!    $export MYDATE=$(date)
!------------------------------------------------------------------------------

    CALL get_environment_variable("MYDATE",date)
    CALL get_environment_variable("HOSTNAME",hostname)
    
    WRITE(11,'(/A)')  TRIM(date)
    WRITE(11,'(/2A)') "This job run on host ",  TRIM(hostname)

!------------------------------------------------------------------------------
! 3. Write basic details about the problem
!------------------------------------------------------------------------------

    WRITE(11,'(/A)')   "BASIC JOB DATA                                  "     
 
    WRITE(11,'(A,I12)')    "Number of processors used                   ",npes 
    WRITE(11,'(A,I12)')    "Number of nodes in the mesh                 ",nn
    WRITE(11,'(A,I12)')    "Number of nodes that were restrained        ",nr
    WRITE(11,'(A,I12)')    "Number of equations solved                  ",neq
    WRITE(11,'(A,I12)')    "Number of PCG iterations                    ",iters
    IF(loaded_nodes > 0) THEN
      WRITE(11,'(A,I12)')    "Number of loaded nodes                      ",   &
                              loaded_nodes 
      WRITE(11,'(A,E12.4)')  "Total load applied                          ",   &
                              tload
    END IF
    IF(fixed_freedoms > 0) THEN
      WRITE(11,'(A,I12)')    "Number of fixed displacements               ",   &
                              fixed_freedoms 
    END IF

!------------------------------------------------------------------------------
! 3. Output timing data
!------------------------------------------------------------------------------

    WRITE(11,'(/3A)')   "PROGRAM SECTION EXECUTION TIMES                  ",  &
                        "SECONDS  ", "%TOTAL    "
    WRITE(11,'(A,F12.6,F8.2)') "Setup                                       ",&
                           timest(2)-timest(1),                               &
                           ((timest(2)-timest(1))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read element steering array                 ",&
                           timest(3)-timest(2),                               &
                           ((timest(3)-timest(2))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Convert Abaqus to S&G node ordering         ",&
                           timest(4)-timest(3),                               &
                           ((timest(4)-timest(3))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read nodal coordinates                      ",&
                           timest(5)-timest(4),                               &
                           ((timest(5)-timest(4))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read restrained nodes                       ",&
                           timest(6)-timest(5),                               &
                           ((timest(6)-timest(5))/(timest(14)-timest(1)))*100                             
    WRITE(11,'(A,F12.6,F8.2)') "Compute steering array and neq              ",&
                           timest(7)-timest(6),                               &
                          ((timest(7)-timest(6))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute interprocessor communication tables ",&
                           timest(8)-timest(7),                               &
                          ((timest(8)-timest(7))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Allocate neq_pp arrays                      ",&
                           timest(9)-timest(8),                               &
                          ((timest(9)-timest(8))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute element stiffness matrices          ",&
                            timest(10)-timest(9),                              &
                          ((timest(10)-timest(9))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Build the preconditioner                    ",&
                           timest(11)-timest(10),                               &
                          ((timest(11)-timest(10))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Get starting r                              ",&
                           timest(12)-timest(11),                               &
                          ((timest(12)-timest(11))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Solve equations                             ",&
                           timest(13)-timest(12),                               &
                           ((timest(13)-timest(12))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Output results                              ",&
                           timest(14)-timest(13),                              &
                          ((timest(14)-timest(13))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A/)')  "Total execution time                        ",  &
                          timest(14)-timest(1),"  100.00"
    CLOSE(11)
    
  END IF
  
  RETURN

  END SUBROUTINE WRITE_XX10

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE WRITE_XX12(fixed_freedoms,job_name,loaded_freedoms,neq,nn,npes,  &
                        nr,numpe,timest,tload,iters,iters_tot,tol,val0,ntime, &
                        timesteps_int,timesteps_real)
  
  !/****f* output/write_xx12
  !*  NAME
  !*    SUBROUTINE: write_xx12
  !*  SYNOPSIS
  !*    Usage:      CALL write_xx12(fixed_freedoms,job_name,loaded_freedoms,  &
  !*                                neq,nn,npes,nr,numpe,timest,tload,iters,  &
  !*                                iters_tot,tol,val0,ntime,timesteps_int,   &
  !*                                timesteps_real)
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem and 
  !*    some performance data
  !*  INPUTS
  !*    The following scalar integers have the INTENT(IN) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    loaded_freedoms        : Number of loaded_nodes
  !*    neq                    : Total number of equations in the mesh
  !*    nn                     : Number of nodes in the mesh
  !*    npes                   : Number of processors used in the simulations
  !*    nr                     : Number of restrained nodes in the mesh
  !*    numpe                  : Processor number
  !*    iters                  : Number of PCG iterations in last time step
  !*    iters_tot              : Total number of PCG iterations
  !*    ntime                  : Number of time sections
  !*
  !*    The following scalar real have the INTENT(IN) attribute:
  !*
  !*    tload                  : Total applied load
  !*    tol                    : Stopping criterion for PCG iterations
  !*    val0                   : Initial global temperature
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : Job name used to name output file
  !*
  !*    The following dynamic real array has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*    timesteps_real(:,:)    : Holds timestep size for time sections
  !*
  !*    The following dynamic integer array has the INTENT(IN) attribute:
  !*
  !*    iters(:)               : Holds iteration information for time sections
  !*    iters_tot(:)           : Holds iteration count for time sections
  !*    timesteps_int(:,:)     : Holds timestep count for time sections
  !*
  !*  AUTHOR
  !*    Lee Margetts, Llion Evans
  !*    Based on Smith I.M. and Griffiths D.V. "Programming the Finite Element
  !*    Method", Edition 4, Wiley, 2004.
  !*  CREATION DATE
  !*    30.07.2014
  !*  COPYRIGHT
  !*    (c) University of Manchester 2014
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  
  
  IMPLICIT NONE
  
  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: numpe,npes,nn,nr,neq,ntime
  INTEGER, INTENT(IN)            :: fixed_freedoms,loaded_freedoms
  INTEGER, INTENT(IN)            :: iters(:),iters_tot(:),timesteps_int(:,:)
  REAL(iwp), INTENT(IN)          :: tload,tol,val0
  REAL(iwp), INTENT(IN)          :: timest(:),timesteps_real(:,:)
  
!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  
  CHARACTER(LEN=50)              :: fname
  INTEGER                        :: i          ! loop counter
  INTEGER                        :: nstep,npri,npri_chk
  REAL(iwp)                      :: dtim,time_tot
  REAL(iwp),ALLOCATABLE          :: totals(:)
  
  IF(numpe==1) THEN
    
    ALLOCATE(totals(5))
    totals=0.0_iwp
    
    fname       = job_name(1:INDEX(job_name, " ")-1) // ".res"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     
    
!------------------------------------------------------------------------------
! 2. Write basic details about the problem
!------------------------------------------------------------------------------
    
    WRITE(11,'(/A)')   "BASIC JOB DATA                                  "     
    
    WRITE(11,'(A,I12)')    "Number of processors used                   ",npes
    WRITE(11,'(A,I12)')    "Number of nodes in the mesh                 ",nn
    WRITE(11,'(A,I12)')    "Number of equations solved                  ",neq
    WRITE(11,'(A,E12.4)')  "Stopping criterion for PCG iterations       ",tol
    
    WRITE(11,'(/A)')   "BOUNDARY CONDITION DATA                         "
    
    WRITE(11,'(A,I12)')    "Number of nodes that were restrained        ",nr
    IF(loaded_freedoms > 0) THEN
      WRITE(11,'(A,I12)')    "Number of loaded freedoms                   ",  &
                              loaded_freedoms 
      WRITE(11,'(A,E12.4)')  "Total energy delivered (load)               ",  &
                              tload
    END IF
    WRITE(11,'(A,E12.4)')  "Initial global temperature                  ",val0
    IF(fixed_freedoms > 0) THEN
      WRITE(11,'(A,I12)')    "Number of fixed nodal values                ",  &
                              fixed_freedoms 
    END IF
    
    WRITE(11,'(/A)')   "TRANSIENT STATE DETAILS                         "
    WRITE(11,'(A,A)')    "Section   dt               nstep        npri",      &
                         "   npri_chk   time        #iters*   Tot#iters"
    time_tot = 0
    DO i=1,ntime
      dtim      = timesteps_real(i,1)
      nstep     = timesteps_int(i,1)
      npri      = timesteps_int(i,2)
      npri_chk  = timesteps_int(i,3)
      WRITE(11,'(I7,A,E12.4,A,I9,A,I9,A,I11,A,E12.4,A,I6,A,I9)')              &
            i," ",dtim,"   ",nstep,"   ",npri,"",npri_chk," ",                &
            dtim*nstep,"   ",iters(i),"   ",iters_tot(i)
      time_tot = time_tot + dtim*nstep
      totals(1) = totals(1) + nstep
      totals(2) = totals(2) + nstep/npri
      totals(3) = totals(3) + int((nstep/npri)/npri_chk)*1.0
      totals(4) = totals(4) + dtim*nstep
      totals(5) = totals(5) + iters_tot(i)
    END DO
    WRITE(11,'(A,A)')  "----------------------------------------------------",&
                       "-------------------------------------"
    WRITE(11,'(A,I9,A,I9,A,I8,A,E12.4,A,I9)')  "Totals                 ",     &
          int(totals(1)),"   ",          int(totals(2)),"   ",int(totals(3)), &
          " ",totals(4),"            ",int(totals(5))
    WRITE(11,'(A,A)')  "----------------------------------------------------",&
                       "-------------------------------------"
    WRITE(11,'(A,A)')  "*Number of iterations to complete last timestep of",  &
                     " section"
    
!------------------------------------------------------------------------------
! 3. Output timing data
!------------------------------------------------------------------------------
    
    WRITE(11,'(/3A)')   "PROGRAM SECTION EXECUTION TIMES                  ",  &
                        "SECONDS  ", "%TOTAL    "
    WRITE(11,'(A,F12.6,F8.2)') "Setup                                       ",&
                           timest(2)-timest(1),                               &
                           ((timest(2)-timest(1))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read element steering array                 ",&
                           timest(3)-timest(2),                               &
                           ((timest(3)-timest(2))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Convert Abaqus to S&G node ordering         ",&
                           timest(4)-timest(3),                               &
                           ((timest(4)-timest(3))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read nodal coordinates                      ",&
                           timest(5)-timest(4),                               &
                           ((timest(5)-timest(4))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read restrained nodes                       ",&
                           timest(6)-timest(5),                               &
                           ((timest(6)-timest(5))/(timest(14)-timest(1)))*100                             
    WRITE(11,'(A,F12.6,F8.2)') "Compute steering array and neq              ",&
                           timest(7)-timest(6),                               &
                          ((timest(7)-timest(6))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute interprocessor communication tables ",&
                           timest(8)-timest(7),                               &
                          ((timest(8)-timest(7))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Allocate neq_pp arrays                      ",&
                           timest(9)-timest(8),                               &
                          ((timest(9)-timest(8))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute element stiffness matrices          ",&
                            timest(10)-timest(9),                             &
                          ((timest(10)-timest(9))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Build the preconditioner                    ",&
                           timest(11)-timest(10),                             &
                          ((timest(11)-timest(10))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Get starting r                              ",&
                           timest(12)-timest(11),                             &
                          ((timest(12)-timest(11))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Solve equations                             ",&
                           timest(13)-timest(12),                             &
                           ((timest(13)-timest(12))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Output results                              ",&
                           timest(14)-timest(13),                             &
                          ((timest(14)-timest(13))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A/)')  "Total execution time                        ", &
                          timest(14)-timest(1),"  100.00"
    CLOSE(11)
    
  END IF
  
  RETURN
  
  END SUBROUTINE WRITE_XX12

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE WRITE_RFEMSOLVE(fixed_freedoms,iters,job_name,loaded_nodes,      &
                             mises,neq,nn,nodecount,npes,nr,numpe,timest,tload)

  !/****f* output/write_rfemsolve
  !*  NAME
  !*    SUBROUTINE: write_rfemsolve
  !*  SYNOPSIS
  !*    Usage:      CALL write_rfemsolve(fixed_freedoms,iters,job_name,       &
  !*                                     loaded_nodes,mises,neq,nn,nodecount, &
  !*                                     npes,nr,numpe,timest,tload)
  !*  FUNCTION
  !*    Master processor writes out brief details about the problem, some 
  !*    performance data and the key results of the analysis
  !*  INPUTS
  !*    The following scalar integers have the INTENT(IN) attribute:
  !*
  !*    fixed_freedoms         : Number of fixed displacements
  !*    iters                  : Number of PCG iterations taken to solve problem
  !*    loaded_nodes           : Number of loaded_nodes
  !*    neq                    : Total number of equations in the mesh
  !*    nn                     : Number of nodes in the mesh
  !*    nodecount              : Number of nodes with value above threshold
  !*    npes                   : Number of processors used in the simulations
  !*    nr                     : Number of restrained nodes in the mesh
  !*    numpe                  : Processor number
  !*
  !*    The following scalar real has the INTENT(IN) attribute:
  !*
  !*    mises                  : Threshold mises value
  !*    tload                  : Total applied load
  !*
  !*    The following scalar character has the INTENT(IN) attribute:
  !*
  !*    job_name               : Job name used to name output file
  !*
  !*    The following dynamic real array has the INTENT(IN) attribute:
  !*
  !*    timest(:)              : Holds timing information
  !*
  !*  AUTHOR
  !*    Lee Margetts
  !*    Based on Smith I.M. and Griffiths D.V. "Programming the Finite Element
  !*    Method", Edition 4, Wiley, 2004.
  !*  CREATION DATE
  !*    12.06.2012
  !*  COPYRIGHT
  !*    (c) University of Manchester 2012
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*
  !*/  
  
  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN)  :: job_name
  INTEGER, INTENT(IN)            :: numpe,npes,nn,nodecount,nr,neq,iters
  INTEGER, INTENT(IN)            :: fixed_freedoms,loaded_nodes
  REAL(iwp), INTENT(IN)          :: timest(:),tload,mises

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  
  CHARACTER(LEN=50)              :: fname
  INTEGER                        :: i          ! loop counter
  REAL(iwp)                      :: nodecountpercent
 
  IF(numpe==1) THEN

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".res"
    OPEN(11,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

    fname       = job_name(1:INDEX(job_name, " ")-1) // ".tnc"
    OPEN(31,FILE=fname,STATUS='REPLACE',ACTION='WRITE')     

!------------------------------------------------------------------------------
! 2. Write basic details about the problem
!------------------------------------------------------------------------------

    WRITE(11,'(/A)')   "BASIC JOB DATA                                  "     
 
    WRITE(11,'(A,I12)')    "Number of processors used                   ",npes 
    WRITE(11,'(A,I12)')    "Number of nodes in the mesh                 ",nn
    WRITE(11,'(A,I12)')    "Number of nodes that were restrained        ",nr
    WRITE(11,'(A,I12)')    "Number of equations solved                  ",neq
    WRITE(11,'(A,I12)')    "Number of PCG iterations                    ",iters
    IF(loaded_nodes > 0) THEN
      WRITE(11,'(A,I12)')    "Number of loaded nodes                      ",   &
                              loaded_nodes 
      WRITE(11,'(A,E12.4)')  "Total load applied                          ",   &
                              tload
    END IF
    IF(fixed_freedoms > 0) THEN
      WRITE(11,'(A,I12)')    "Number of fixed displacements               ",   &
                              fixed_freedoms 
    END IF

    WRITE(11,'(A,F12.6)')  "Threshold values for von Mises stress       ",mises
    WRITE(11,'(A,I12)')    "Number of nodes greater than threshold      ",     &
                            nodecount
    nodecountpercent = (100.0 / nn) * nodecount
    WRITE(31,'(I12,1F8.2)')       nodecount, nodecountpercent

!-------------------------------------------------------------------------------
! 3. Output timing data
!-------------------------------------------------------------------------------

    WRITE(11,'(/3A)')   "PROGRAM SECTION EXECUTION TIMES                  ",   &
                        "SECONDS  ", "%TOTAL    "
    WRITE(11,'(A,F12.6,F8.2)') "Setup                                       ", &
                           timest(2)-timest(1),                                &
                           ((timest(2)-timest(1))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read element steering array                 ", &
                           timest(3)-timest(2),                                &
                           ((timest(3)-timest(2))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Convert Abaqus to S&G node ordering         ", &
                           timest(4)-timest(3),                                &
                           ((timest(4)-timest(3))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read nodal coordinates                      ", &
                           timest(5)-timest(4),                                &
                           ((timest(5)-timest(4))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Read restrained nodes                       ", &
                           timest(6)-timest(5),                                &
                           ((timest(6)-timest(5))/(timest(14)-timest(1)))*100                             
    WRITE(11,'(A,F12.6,F8.2)') "Compute steering array and neq              ", &
                           timest(7)-timest(6),                                &
                          ((timest(7)-timest(6))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute interprocessor communication tables ", &
                           timest(8)-timest(7),                                &
                          ((timest(8)-timest(7))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Allocate neq_pp arrays                      ", &
                           timest(9)-timest(8),                                &
                          ((timest(9)-timest(8))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute element stiffness matrices          ", &
                            timest(10)-timest(9),                              &
                          ((timest(10)-timest(9))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Build the preconditioner                    ", &
                           timest(11)-timest(10),                              &
                          ((timest(11)-timest(10))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Get starting r                              ", &
                           timest(12)-timest(11),                              &
                          ((timest(12)-timest(11))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Solve equations                             ", &
                           timest(13)-timest(12),                              &
                           ((timest(13)-timest(12))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Output results                              ", &
                           timest(14)-timest(13),                              &
                          ((timest(14)-timest(13))/(timest(14)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A/)')  "Total execution time                        ",  &
                          timest(14)-timest(1),"  100.00"
   
    WRITE(11,'(/,A)') "Please cite the following references:"
    WRITE(11,'(/,A)')   "DOI:10.1007/s11831-014-9139-3"
    WRITE(11,'(A)')     "DOI:10.1016/j.nucmat.2015.05.058"
 
    CLOSE(11)
    CLOSE(31)
    
  END IF
  
  RETURN

  END SUBROUTINE WRITE_RFEMSOLVE

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
  SUBROUTINE WRITE_NODAL_VARIABLE(text,filnum,iload,nodes_pp,npes,numpe, &
                                  numvar,stress)

  !/****f* output/write_nodal_variable
  !*  NAME
  !*    SUBROUTINE: write_nodal_variable
  !*  SYNOPSIS
  !*    Usage:      CALL write_nodal_variable(text,filnum,iload,nodes_pp, &
  !*                                          numpe,numvar,stress) 
  !*  FUNCTION
  !*    Write the values of a nodal variable to a file. This subroutine is 
  !*    parallel and requires MPI. The master processor collects all the data
  !*    from the slave processors.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    text                    : Character
  !*                            : Text indicating the variable to write
  !*
  !*    filnum                  : Integer
  !*                            : File number to write
  !*
  !*    iload                   : Integer
  !*                            : Load step number
  !*
  !*    nodes_pp                : Integer
  !*                            : Number of nodes assigned to a process
  !*
  !*    npes                    : Integer
  !*                            : Number of processes
  !*
  !*    numpe                   : Integer
  !*                            : Process number
  !*
  !*    numvar                  : Integer
  !*                            : Number of components of the variable
  !*                              (1-scalar, 3-vector, 6-tensor)
  !*
  !*    stress(nodes_pp*numvar) : Real
  !*                            : Nodal variables to print
  !*                             
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*  AUTHOR
  !*    F. Calvo
  !*    L. Margetts
  !*  CREATION DATE
  !*    04.10.2007
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2010
  !******
  !*
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN) :: text
  INTEGER, INTENT(IN)           :: filnum, iload, nodes_pp, npes, numpe, numvar
  REAL(iwp), INTENT(IN)         :: stress(:)
  INTEGER                       :: i, j, idx1, nod_r, bufsize1, bufsize2
  INTEGER                       :: ier, iproc, n, bufsize
  INTEGER                       :: statu(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: get(:)
  REAL(iwp), ALLOCATABLE        :: stress_r(:)

!------------------------------------------------------------------------------
! 1. Master processor writes out the results for its own assigned nodes
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    WRITE(filnum,'(a)')text
    WRITE(filnum,*)iload
    DO i = 1,nodes_pp
      idx1 = (i-1)*numvar + 1
      IF (numvar==1) THEN
        WRITE(filnum,2001)i,(stress(j),j=idx1,idx1+numvar-1)
      END IF
      IF (numvar==3) THEN
        WRITE(filnum,2003)i,(stress(j),j=idx1,idx1+numvar-1)
      END IF
      IF (numvar==4) THEN
        WRITE(filnum,2004)i,(stress(j),j=idx1,idx1+numvar-1)
      END IF
      IF (numvar==6) THEN
        WRITE(filnum,2006)i,(stress(j),j=idx1,idx1+numvar-1)
      END IF
    END DO
  END IF    
  
!------------------------------------------------------------------------------
! 2. Allocate arrays involved in communications
!------------------------------------------------------------------------------

  ALLOCATE(get(npes))
  ALLOCATE(stress_r(nodes_pp*numvar)) 

!------------------------------------------------------------------------------
! 3. Master processor populates the array "get" containing "nodes_pp" of 
!    every processor. Slave processors send this number to the master processor
!------------------------------------------------------------------------------

  get = 0
  get(1) = nodes_pp

  bufsize = 1

  DO i = 2,npes
    IF(numpe==i) THEN
      CALL MPI_SEND(nodes_pp,bufsize,MPI_INTEGER,0,i,MPI_COMM_WORLD,ier)
    END IF
    IF(numpe==1) THEN
      CALL MPI_RECV(nod_r,bufsize,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
      get(i) = nod_r
    END IF
  END DO

!------------------------------------------------------------------------------
! 4. Master processor receives the nodal variables from the other processors 
!    and writes them
!------------------------------------------------------------------------------

  DO iproc = 2,npes
    IF (numpe==iproc) THEN
      bufsize1 = nodes_pp*numvar
      CALL MPI_SEND(stress,bufsize1,MPI_REAL8,0,iproc,MPI_COMM_WORLD,ier)
    END IF
    IF (numpe==1) THEN
      bufsize2 = get(iproc)*numvar
      CALL MPI_RECV(stress_r,bufsize2,MPI_REAL8,iproc-1,iproc, &
                    MPI_COMM_WORLD,statu,ier)
      n = 1
      DO j = 2,iproc
        n = n + get(j-1)
      END DO
      DO i = 1,get(iproc)
        idx1 = (i-1)*numvar + 1
        IF (numvar==1) THEN
          WRITE(filnum,2001)n-1+i,(stress_r(j),j=idx1,idx1+numvar-1)
        END IF
        IF (numvar==3) THEN
          WRITE(filnum,2003)n-1+i,(stress_r(j),j=idx1,idx1+numvar-1)
        END IF
        IF (numvar==4) THEN
          WRITE(filnum,2004)n-1+i,(stress_r(j),j=idx1,idx1+numvar-1)
        END IF
        IF (numvar==6) THEN
          WRITE(filnum,2006)n-1+i,(stress_r(j),j=idx1,idx1+numvar-1)
        END IF
      END DO
    END IF
  END DO

!------------------------------------------------------------------------------
! 5. Deallocate communication arrays
!------------------------------------------------------------------------------

  DEALLOCATE(get,stress_r)

!------------------------------------------------------------------------------
! 6. Set formats used in this subroutine
!------------------------------------------------------------------------------

  2001 FORMAT(i8,1(1p,e12.4))
  2003 FORMAT(i8,3(1p,e12.4))
  2004 FORMAT(i8,4(1p,e12.4))
  2006 FORMAT(i8,6(1p,e12.4))

  END SUBROUTINE WRITE_NODAL_VARIABLE

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE WRITE_NODAL_VARIABLE2(text,filnum,iload,nodes_pp,npes,numpe, &
                                  numvar,stress)

  !/****f* output/write_nodal_variable2
  !*  NAME
  !*    SUBROUTINE: write_nodal_variable2
  !*  SYNOPSIS
  !*    Usage:      CALL write_nodal_variable2(text,filnum,iload,nodes_pp, &
  !*                                          numpe,numvar,stress) 
  !*  FUNCTION
  !*    Write the values of a nodal variable to a file. This subroutine is 
  !*    parallel and requires MPI. The master processor collects all the data
  !*    from the slave processors.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    text                    : Character
  !*                            : Text indicating the variable to write
  !*
  !*    filnum                  : Integer
  !*                            : File number to write
  !*
  !*    iload                   : Integer
  !*                            : Load step number
  !*
  !*    nodes_pp                : Integer
  !*                            : Number of nodes assigned to a process
  !*
  !*    npes                    : Integer
  !*                            : Number of processes
  !*
  !*    numpe                   : Integer
  !*                            : Process number
  !*
  !*    numvar                  : Integer
  !*                            : Number of components of the variable
  !*                              (1-scalar, 3-vector, 6-tensor)
  !*
  !*    stress(nodes_pp*numvar) : Real
  !*                            : Nodal variables to print
  !*                             
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*  AUTHOR
  !*    F. Calvo
  !*    L. Margetts
  !*  CREATION DATE
  !*    04.10.2007
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2010
  !******
  !*
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN) :: text
  INTEGER, INTENT(IN)           :: filnum, iload, nodes_pp, npes, numpe, numvar
  REAL(iwp), INTENT(IN)         :: stress(:)
  INTEGER                       :: i, j, idx1, nod_r, bufsize1, bufsize2
  INTEGER                       :: ier, iproc, n, bufsize
  INTEGER                       :: statu(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: get(:)
  REAL(iwp), ALLOCATABLE        :: stress_r(:)

!------------------------------------------------------------------------------
! 1. Master processor writes out the results for its own assigned nodes
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    WRITE(filnum,'(a)')text
    WRITE(filnum,*)iload
    DO i = 1,nodes_pp
      idx1 = (i-1)*numvar + 1
      IF (numvar==1) THEN
        WRITE(filnum,2001)i,(stress(j),j=idx1,idx1+numvar-1)
      END IF
      IF (numvar==3) THEN
        WRITE(filnum,2003)i,(stress(j),j=idx1,idx1+numvar-1)
      END IF
      IF (numvar==4) THEN
        WRITE(filnum,2004)i,(stress(j),j=idx1,idx1+numvar-1)
      END IF
      IF (numvar==6) THEN
        WRITE(filnum,2006)i,(stress(j),j=idx1,idx1+numvar-1)
      END IF
    END DO
  END IF    
  
!------------------------------------------------------------------------------
! 2. Allocate arrays involved in communications
!------------------------------------------------------------------------------

  ALLOCATE(get(npes))
  ALLOCATE(stress_r(nodes_pp*numvar)) 

!------------------------------------------------------------------------------
! 3. Master processor populates the array "get" containing "nodes_pp" of 
!    every processor. Slave processors send this number to the master processor
!------------------------------------------------------------------------------

  get = 0
  get(1) = nodes_pp

  bufsize = 1

  DO i = 2,npes
    IF(numpe==i) THEN
      CALL MPI_SEND(nodes_pp,bufsize,MPI_INTEGER,0,i,MPI_COMM_WORLD,ier)
    END IF
    IF(numpe==1) THEN
      CALL MPI_RECV(nod_r,bufsize,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
      get(i) = nod_r
    END IF
  END DO

!------------------------------------------------------------------------------
! 4. Master processor receives the nodal variables from the other processors 
!    and writes them
!------------------------------------------------------------------------------

  DO iproc = 2,npes
    IF (numpe==iproc) THEN
      bufsize1 = nodes_pp*numvar
      CALL MPI_SEND(stress,bufsize1,MPI_REAL8,0,iproc,MPI_COMM_WORLD,ier)
    END IF
    IF (numpe==1) THEN
      bufsize2 = get(iproc)*numvar
      CALL MPI_RECV(stress_r,bufsize2,MPI_REAL8,iproc-1,iproc, &
                    MPI_COMM_WORLD,statu,ier)
      n = 1
      DO j = 2,iproc
        n = n + get(j-1)
      END DO
      DO i = 1,get(iproc)
        idx1 = (i-1)*numvar + 1
        IF (numvar==1) THEN
          WRITE(filnum,2001)n-1+i,(stress_r(j),j=idx1,idx1+numvar-1)
        END IF
        IF (numvar==3) THEN
          WRITE(filnum,2003)n-1+i,(stress_r(j),j=idx1,idx1+numvar-1)
        END IF
        IF (numvar==4) THEN
          WRITE(filnum,2004)n-1+i,(stress_r(j),j=idx1,idx1+numvar-1)
        END IF
        IF (numvar==6) THEN
          WRITE(filnum,2006)n-1+i,(stress_r(j),j=idx1,idx1+numvar-1)
        END IF
      END DO
    END IF
  END DO

!------------------------------------------------------------------------------
! 5. Deallocate communication arrays
!------------------------------------------------------------------------------

  DEALLOCATE(get,stress_r)

!------------------------------------------------------------------------------
! 6. Set formats used in this subroutine
!------------------------------------------------------------------------------

  2001 FORMAT(i8,1(1p,e19.8))
  2003 FORMAT(i8,3(1p,e19.8))
  2004 FORMAT(i8,4(1p,e19.8))
  2006 FORMAT(i8,6(1p,e19.8))

  END SUBROUTINE WRITE_NODAL_VARIABLE2

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
  SUBROUTINE WRITE_NODAL_VARIABLE_BINARY(text,filnum,iload,nodes_pp,npes,     &
                                         numpe,numvar,stress)

  !/****f* output/write_nodal_variable_binary
  !*  NAME
  !*    SUBROUTINE: write_nodal_variable_binary
  !*  SYNOPSIS
  !*    Usage:      CALL write_nodal_variable_binary(text,filnum,iload,       &
  !*                                                 nodes_pp,numpe,numvar,   &
  !*                                                 stress) 
  !*  FUNCTION
  !*    Write the values of a nodal variable to a file. This subroutine is 
  !*    parallel and requires MPI. The master processor collects all the data
  !*    from the slave processors. Binary output.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    text                    : Character
  !*                            : Text indicating the variable to write
  !*
  !*    filnum                  : Integer
  !*                            : File number to write
  !*
  !*    iload                   : Integer
  !*                            : Load step number
  !*
  !*    nodes_pp                : Integer
  !*                            : Number of nodes assigned to a process
  !*
  !*    npes                    : Integer
  !*                            : Number of processes
  !*
  !*    numpe                   : Integer
  !*                            : Process number
  !*
  !*    numvar                  : Integer
  !*                            : Number of components of the variable
  !*                              (1-scalar, 3-vector, 6-tensor)
  !*
  !*    stress(nodes_pp*numvar) : Real
  !*                            : Nodal variables to print
  !*                             
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*  AUTHOR
  !*    F. Calvo
  !*    L. Margetts
  !*  CREATION DATE
  !*    27.11.2012
  !*  COPYRIGHT
  !*    (c) University of Manchester 2012-2012
  !******
  !*
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN) :: text
  INTEGER, INTENT(IN)           :: filnum, iload, nodes_pp, npes, numpe, numvar
  REAL(iwp), INTENT(IN)         :: stress(:)
  INTEGER                       :: i, j, idx1, nod_r, bufsize1, bufsize2
  INTEGER                       :: ier, iproc, n, bufsize
  INTEGER                       :: statu(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: get(:)
  REAL(iwp), ALLOCATABLE        :: stress_r(:)

!------------------------------------------------------------------------------
! 1. Master processor writes out the results for its own assigned nodes
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    WRITE(filnum)text
    WRITE(filnum)iload
!   DO i = 1,nodes_pp
!     idx1 = (i-1)*numvar + 1
!     IF (numvar==1) THEN
        WRITE(filnum)stress
!     END IF
!     IF (numvar==3) THEN
!       WRITE(filnum)i,(stress(j),j=idx1,idx1+numvar-1)
!     END IF
!     IF (numvar==4) THEN
!       WRITE(filnum)i,(stress(j),j=idx1,idx1+numvar-1)
!     END IF
!     IF (numvar==6) THEN
!       WRITE(filnum)i,(stress(j),j=idx1,idx1+numvar-1)
!     END IF
!   END DO
  END IF    
  
!------------------------------------------------------------------------------
! 2. Allocate arrays involved in communications
!------------------------------------------------------------------------------

  ALLOCATE(get(npes))
  ALLOCATE(stress_r(nodes_pp*numvar)) 

!------------------------------------------------------------------------------
! 3. Master processor populates the array "get" containing "nodes_pp" of 
!    every processor. Slave processors send this number to the master processor
!------------------------------------------------------------------------------

  get = 0
  get(1) = nodes_pp

  bufsize = 1

  DO i = 2,npes
    IF(numpe==i) THEN
      CALL MPI_SEND(nodes_pp,bufsize,MPI_INTEGER,0,i,MPI_COMM_WORLD,ier)
    END IF
    IF(numpe==1) THEN
      CALL MPI_RECV(nod_r,bufsize,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
      get(i) = nod_r
    END IF
  END DO

!------------------------------------------------------------------------------
! 4. Master processor receives the nodal variables from the other processors 
!    and writes them
!------------------------------------------------------------------------------

  DO iproc = 2,npes
    IF (numpe==iproc) THEN
      bufsize1 = nodes_pp*numvar
      CALL MPI_SEND(stress,bufsize1,MPI_REAL8,0,iproc,MPI_COMM_WORLD,ier)
    END IF
    IF (numpe==1) THEN
      bufsize2 = get(iproc)*numvar
      CALL MPI_RECV(stress_r,bufsize2,MPI_REAL8,iproc-1,iproc, &
                    MPI_COMM_WORLD,statu,ier)
!     n = 1
!     DO j = 2,iproc
!       n = n + get(j-1)
!     END DO
!     DO i = 1,get(iproc)
!       idx1 = (i-1)*numvar + 1
!       IF (numvar==1) THEN
          WRITE(filnum)stress_r !(1:(get(iproc)))
!       END IF
!       IF (numvar==3) THEN
!         WRITE(filnum)n-1+i,(stress_r(j),j=idx1,idx1+numvar-1)
!       END IF
!       IF (numvar==4) THEN
!         WRITE(filnum)n-1+i,(stress_r(j),j=idx1,idx1+numvar-1)
!       END IF
!       IF (numvar==6) THEN
!         WRITE(filnum)n-1+i,(stress_r(j),j=idx1,idx1+numvar-1)
!       END IF
!     END DO
    END IF
  END DO

  IF(numpe==1) THEN
    DO i=1,npes
      WRITE(filnum+1,*) get(i)
    END DO
  END IF


!------------------------------------------------------------------------------
! 5. Deallocate communication arrays
!------------------------------------------------------------------------------

  DEALLOCATE(get,stress_r)

!------------------------------------------------------------------------------
! 6. Set formats used in this subroutine
!------------------------------------------------------------------------------

! 2001 FORMAT(i8,1(1p,e12.4))
! 2003 FORMAT(i8,3(1p,e12.4))
! 2004 FORMAT(i8,4(1p,e12.4))
! 2006 FORMAT(i8,6(1p,e12.4))

  END SUBROUTINE WRITE_NODAL_VARIABLE_BINARY
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
 
  SUBROUTINE WRITE_X_PP(text,filnum,iload,nodes_pp,npes,     &
                                         numpe,numvar,stress)

  !/****f* output/write_x_pp
  !*  NAME
  !*    SUBROUTINE: write_x_pp
  !*  SYNOPSIS
  !*    Usage:      CALL write_x_pp(text,filnum,iload,nodes_pp,              &
  !*                                                 numpe,numvar,stress)
  !*  FUNCTION
  !*    Write the values of a nodal variable to a file. This subroutine is 
  !*    parallel and requires MPI. The master processor collects all the data
  !*    from the slave processors. Binary output.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    text                    : Character
  !*                            : Text indicating the variable to write
  !*
  !*    filnum                  : Integer
  !*                            : File number to write
  !*
  !*    iload                   : Integer
  !*                            : Load step number
  !*
  !*    nodes_pp                : Integer
  !*                            : Number of nodes assigned to a process
  !*
  !*    npes                    : Integer
  !*                            : Number of processes
  !*
  !*    numpe                   : Integer
  !*                            : Process number
  !*
  !*    numvar                  : Integer
  !*                            : Number of components of the variable
  !*                              (1-scalar, 3-vector, 6-tensor)
  !*
  !*    stress(nodes_pp*numvar) : Real
  !*                            : Nodal variables to print
  !*                             
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*  AUTHOR
  !*    F. Calvo
  !*    L. Margetts
  !*  CREATION DATE
  !*    27.11.2012
  !*  COPYRIGHT
  !*    (c) University of Manchester 2012-2012
  !******
  !*
  !*/

  IMPLICIT NONE

  CHARACTER(LEN=50), INTENT(IN) :: text
  INTEGER, INTENT(IN)           :: filnum, iload, nodes_pp, npes, numpe, numvar
  REAL(iwp), INTENT(IN)         :: stress(:)
  INTEGER                       :: i, j, idx1, nod_r, bufsize1, bufsize2
  INTEGER                       :: ier, iproc, n, bufsize
  INTEGER                       :: statu(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: get(:)
  REAL(iwp), ALLOCATABLE        :: stress_r(:)

!------------------------------------------------------------------------------
! 1. Master processor writes out the results for its own assigned nodes
!------------------------------------------------------------------------------

  IF (numpe==1) THEN
    WRITE(filnum)text
    WRITE(filnum)iload
    WRITE(filnum)stress
  END IF    
  
!------------------------------------------------------------------------------
! 2. Allocate arrays involved in communications
!------------------------------------------------------------------------------

  ALLOCATE(get(npes))
  ALLOCATE(stress_r(nodes_pp*numvar)) 

!------------------------------------------------------------------------------
! 3. Master processor populates the array "get" containing "nodes_pp" of 
!    every processor. Slave processors send this number to the master processor
!------------------------------------------------------------------------------

  get = 0
  get(1) = nodes_pp

  bufsize = 1

  DO i = 2,npes
    IF(numpe==i) THEN
      CALL MPI_SEND(nodes_pp,bufsize,MPI_INTEGER,0,i,MPI_COMM_WORLD,ier)
    END IF
    IF(numpe==1) THEN
      CALL MPI_RECV(nod_r,bufsize,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
      get(i) = nod_r
    END IF
  END DO

!------------------------------------------------------------------------------
! 4. Master processor receives the nodal variables from the other processors 
!    and writes them
!------------------------------------------------------------------------------

  DO iproc = 2,npes
    IF (numpe==iproc) THEN
      bufsize1 = nodes_pp*numvar
      CALL MPI_SEND(stress,bufsize1,MPI_REAL8,0,iproc,MPI_COMM_WORLD,ier)
    END IF
    IF (numpe==1) THEN
      bufsize2 = get(iproc)*numvar
      CALL MPI_RECV(stress_r,bufsize2,MPI_REAL8,iproc-1,iproc,                &
                    MPI_COMM_WORLD,statu,ier)
      WRITE(filnum)stress_r
    END IF
  END DO
  
!------------------------------------------------------------------------------
! 5. Deallocate communication arrays
!------------------------------------------------------------------------------
  
  DEALLOCATE(get,stress_r)
  
  END SUBROUTINE WRITE_X_PP
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE job_name_error(numpe,program_name)

  !/****f* output/job_name_error
  !*  NAME
  !*    SUBROUTINE: job_name_error
  !*  SYNOPSIS
  !*    Usage:      CALL job_name_error(numpe,program_name)
  !*  FUNCTION
  !*    Generates error message if commandline argument is missing.
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2012
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/


 IMPLICIT none

 INTEGER, INTENT(IN)           :: numpe        ! processor ID number
 CHARACTER(LEN=50), INTENT(IN) :: program_name ! program name 
 CHARACTER(LEN=50)             :: fname

 IF (numpe==1) THEN
   IF (TRIM(program_name)=='rfemsolve') THEN
     fname   = program_name(1:INDEX(program_name, " ")-1)//".res"
     OPEN (11, file=fname, status='replace', action='write')
     WRITE(11,'(/4A/)') "Fatal error: ", TRIM(program_name), " did not run",   &
                        " - one or more arguments missing."
     WRITE(11,'(4A/)')"Usage: ", TRIM(program_name), " job_in", " job_out"
     WRITE(11,'(4A)') "       ", TRIM(program_name), " will read job_in,",     &
                      " job_out, <job_in>.dat,"       
     WRITE(11,'(A)')  "       <job_in>.d, <job_in>.bnd, <job_in>.lds,"  
     WRITE(11,'(A/)') "       and write <job_out>.dis"
     WRITE(11,'(A/)') "       If analysis_type > 2, <job_in>.mat is expected."
     CLOSE(11)
   ELSE
     fname   = program_name(1:INDEX(program_name, " ")-1)//".res"
     OPEN (11, file=fname, status='replace', action='write')
     WRITE(11,'(/4A/)') "Fatal error: ", TRIM(program_name), " did not run",   &
                        " - job_name is missing."
     WRITE(11,'(3A/)')"Usage: ", TRIM(program_name), " job_name"
     WRITE(11,'(4A)') "       ", TRIM(program_name), " will read job_name,",   &
                      " <job_name>.dat,"       
     WRITE(11,'(A)')  "       <job_name>.d, <job_name>.bnd, <job_name>.lds,"  
     WRITE(11,'(A/)') "       and write <job_name>.dis"
     WRITE(11,'(A/)') "       If analysis_type > 2, <job_name>.mat is expected."
     CLOSE(11)
   END IF
 END IF

 CALL shutdown()
 STOP
 
 END SUBROUTINE job_name_error
 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------ 

  SUBROUTINE getFileNumber(fileNumber,count)
 
  !/****f* output/getFileNumber
  !*  NAME
  !*    SUBROUTINE: getFileNumber
  !*  SYNOPSIS
  !*    Usage:      CALL getFileNumber(fileNumber,count)
  !*  FUNCTION
  !*    Returns file number.
  !*  AUTHOR
  !*    L. Margetts
  !*  COPYRIGHT
  !*    (c) University of Manchester 2007-2011
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/

  IMPLICIT NONE
 
  INTEGER, INTENT(IN)         :: count
  INTEGER                     :: fileIndexLength
  INTEGER                     :: maxCount
  CHARACTER(*), INTENT(INOUT) :: fileNumber
        
  IF(count <    10 )                     fileIndexLength = 1
  IF(count >     9 .and. count <    100) fileIndexLength = 2
  IF(count >    99 .and. count <   1000) fileIndexLength = 3
  IF(count >   999 .and. count <  10000) fileIndexLength = 4
  IF(count >  9999 .and. count < 100000) fileIndexLength = 5 
  IF(count > 99999 .and. count < 999999) fileIndexLength = 6

  SELECT CASE (fileIndexLength)

    CASE(1)

    WRITE (fileNumber,"(I1)") count 
  
    CASE(2)

    WRITE (fileNumber,"(I2)") count

    CASE(3)

    WRITE (fileNumber,"(I3)") count
       
    CASE(4)

    WRITE (fileNumber,"(I4)") count

    CASE(5)

    WRITE (fileNumber,"(I5)") count
        
    CASE(6)

    WRITE (fileNumber,"(I6)") count

    CASE DEFAULT

    PRINT *, " "
    PRINT *, "[Main] Warning"
    PRINT *, "[Main] Too many files"
    PRINT *, "[Main] File number 999,999 has been overwritten with new data"
    PRINT *, " "
       
    maxCount = 999999
    WRITE (fileNumber,"(I6)") maxCount

  END SELECT

  END SUBROUTINE getFileNumber
 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------ 
 
  SUBROUTINE DISMSH_ENSI_P(filnum,iload,nodes_pp,npes,numpe, &
                                  numvar,stress)

  !/****f* output/write_nodal_variable
  !*  NAME
  !*    SUBROUTINE: write_nodal_variable
  !*  SYNOPSIS
  !*    Usage:      CALL write_nodal_variable(filnum,iload,nodes_pp, &
  !*                                          numpe,numvar,stress) 
  !*  FUNCTION
  !*    Write the values of a nodal variable to a file. This subroutine is 
  !*    parallel and requires MPI. The master processor collects all the data
  !*    from the slave processors.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    filnum                  : Integer
  !*                            : File number to write
  !*
  !*    iload                   : Integer
  !*                            : Load step number
  !*
  !*    nodes_pp                : Integer
  !*                            : Number of nodes assigned to a process
  !*
  !*    npes                    : Integer
  !*                            : Number of processes
  !*
  !*    numpe                   : Integer
  !*                            : Process number
  !*
  !*    numvar                  : Integer
  !*                            : Number of components of the variable
  !*                              (1-scalar, 3-vector, 6-tensor)
  !*
  !*    stress(nodes_pp*numvar) : Real
  !*                            : Nodal variables to print
  !*                             
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*  AUTHOR
  !*    L. Margetts
  !*  CREATION DATE
  !*    11.01.2013
  !*  COPYRIGHT
  !*    (c) University of Manchester 2013
  !******
  !*
  !*/

  IMPLICIT NONE

  INTEGER, INTENT(IN)           :: filnum, iload, nodes_pp, npes, numpe, numvar
  REAL(iwp), INTENT(IN)         :: stress(:)
  INTEGER                       :: i, j, idx1, nod_r, bufsize1, bufsize2
  INTEGER                       :: ier, iproc, n, bufsize
  INTEGER                       :: statu(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: get(:)
  REAL(iwp), ALLOCATABLE        :: stress_r(:)

!------------------------------------------------------------------------------
! 1. Master processor writes out the results for its own assigned nodes
!------------------------------------------------------------------------------

  IF(numpe==1) THEN
    DO i = 1,nodes_pp
      WRITE(filnum,'(e12.5)') stress(i)
    END DO
  END IF    
  
!------------------------------------------------------------------------------
! 2. Allocate arrays involved in communications
!------------------------------------------------------------------------------

  ALLOCATE(get(npes))
  ALLOCATE(stress_r(nodes_pp*numvar)) 

!------------------------------------------------------------------------------
! 3. Master processor populates the array "get" containing "nodes_pp" of 
!    every processor. Slave processors send this number to the master processor
!------------------------------------------------------------------------------

  get = 0
  get(1) = nodes_pp

  bufsize = 1

  DO i = 2,npes
    IF(numpe==i) THEN
      CALL MPI_SEND(nodes_pp,bufsize,MPI_INTEGER,0,i,MPI_COMM_WORLD,ier)
    END IF
    IF(numpe==1) THEN
      CALL MPI_RECV(nod_r,bufsize,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
      get(i) = nod_r
    END IF
  END DO

!------------------------------------------------------------------------------
! 4. Master processor receives the nodal variables from the other processors 
!    and writes them
!------------------------------------------------------------------------------

  DO iproc = 2,npes
    IF (numpe==iproc) THEN
      bufsize1 = nodes_pp*numvar
      CALL MPI_SEND(stress,bufsize1,MPI_REAL8,0,iproc,MPI_COMM_WORLD,ier)
    END IF
    IF (numpe==1) THEN
      bufsize2 = get(iproc)*numvar
      CALL MPI_RECV(stress_r,bufsize2,MPI_REAL8,iproc-1,iproc, &
                    MPI_COMM_WORLD,statu,ier)
      n = 1
      DO j = 2,iproc
        n = n + get(j-1)
      END DO
      DO i = 1,get(iproc)
        WRITE(filnum,'(e12.5)') stress_r(i)
      END DO
    END IF
  END DO

!------------------------------------------------------------------------------
! 5. Deallocate communication arrays
!------------------------------------------------------------------------------

  DEALLOCATE(get,stress_r)

  END SUBROUTINE DISMSH_ENSI_P







      SUBROUTINE SWAP_F8 (FLOAT8)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN OUT) :: FLOAT8

      INTEGER(KIND=1), DIMENSION(8) :: BYTE_ARR, BYTE_ARR_TMP
      INTEGER :: I


      BYTE_ARR = TRANSFER (FLOAT8, BYTE_ARR)

      BYTE_ARR_TMP = BYTE_ARR

      DO I = 1, 8
         BYTE_ARR(I) = BYTE_ARR_TMP(9-I)
      END DO

      FLOAT8 = TRANSFER (BYTE_ARR, FLOAT8)

      RETURN

      END SUBROUTINE SWAP_F8



  SUBROUTINE DISMSH_ENSI_PB(filnum,iload,nodes_pp,npes,numpe, &
                                  numvar,stress)

  !/****f* output/dismsh_ensi_pb
  !*  NAME
  !*    SUBROUTINE: dismsh_ensi_pb
  !*  SYNOPSIS
  !*    Usage:      CALL dismsh_ensi_pb(filnum,iload,nodes_pp, &
  !*                                          numpe,numvar,stress) 
  !*  FUNCTION
  !*    ! DO NOT USE - UNDER DEVELOPMENT !
  !*    Write the values of a nodal variable to a binary file. This subroutine is 
  !*    parallel and requires MPI. The master processor collects all the data
  !*    from the slave processors. Currently used in development program xx12.f90.
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    filnum                  : Integer
  !*                            : File number to write
  !*
  !*    iload                   : Integer
  !*                            : Load step number
  !*
  !*    nodes_pp                : Integer
  !*                            : Number of nodes assigned to a process
  !*
  !*    npes                    : Integer
  !*                            : Number of processes
  !*
  !*    numpe                   : Integer
  !*                            : Process number
  !*
  !*    numvar                  : Integer
  !*                            : Number of components of the variable
  !*                              (1-scalar, 3-vector, 6-tensor)
  !*
  !*    stress(nodes_pp*numvar) : Real
  !*                            : Nodal variables to print
  !*                             
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*  AUTHOR
  !*    L. Margetts
  !*  CREATION DATE
  !*    11.01.2013
  !*  COPYRIGHT
  !*    (c) University of Manchester 2013
  !******
  !*
  !*/

  IMPLICIT NONE

  INTEGER, INTENT(IN)           :: filnum, iload, nodes_pp, npes, numpe, numvar
  REAL(iwp), INTENT(IN)         :: stress(:)
  INTEGER                       :: i, j, idx1, nod_r, bufsize1, bufsize2
  INTEGER                       :: ier, iproc, n, bufsize
  INTEGER                       :: statu(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: get(:)
  REAL*4                        :: stress_single
  REAL*4                        :: stress_single_little
  REAL(iwp), ALLOCATABLE        :: stress_r(:)

!------------------------------------------------------------------------------
! 1. Master processor writes out the results for its own assigned nodes
!------------------------------------------------------------------------------

  IF(numpe==1) THEN
    DO i = 1,nodes_pp
      stress_single = stress(i)
      WRITE(filnum) stress_single
    END DO
  END IF    
  
!------------------------------------------------------------------------------
! 2. Allocate arrays involved in communications
!------------------------------------------------------------------------------

  ALLOCATE(get(npes))
  ALLOCATE(stress_r(nodes_pp*numvar)) 

!------------------------------------------------------------------------------
! 3. Master processor populates the array "get" containing "nodes_pp" of 
!    every processor. Slave processors send this number to the master processor
!------------------------------------------------------------------------------

  get = 0
  get(1) = nodes_pp

  bufsize = 1

  DO i = 2,npes
    IF(numpe==i) THEN
      CALL MPI_SEND(nodes_pp,bufsize,MPI_INTEGER,0,i,MPI_COMM_WORLD,ier)
    END IF
    IF(numpe==1) THEN
      CALL MPI_RECV(nod_r,bufsize,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
      get(i) = nod_r
    END IF
  END DO

!------------------------------------------------------------------------------
! 4. Master processor receives the nodal variables from the other processors 
!    and writes them
!------------------------------------------------------------------------------

  DO iproc = 2,npes
    IF (numpe==iproc) THEN
      bufsize1 = nodes_pp*numvar
      CALL MPI_SEND(stress,bufsize1,MPI_REAL4,0,iproc,MPI_COMM_WORLD,ier)
    END IF
    IF (numpe==1) THEN
      bufsize2 = get(iproc)*numvar
      CALL MPI_RECV(stress_r,bufsize2,MPI_REAL4,iproc-1,iproc, &
                    MPI_COMM_WORLD,statu,ier)
      n = 1
      DO j = 2,iproc
        n = n + get(j-1)
      END DO
      DO i = 1,get(iproc)
        stress_single = stress_r(i)
        WRITE(filnum) stress_single
      END DO
    END IF
  END DO

!------------------------------------------------------------------------------
! 5. Deallocate communication arrays
!------------------------------------------------------------------------------

  DEALLOCATE(get,stress_r)


  END SUBROUTINE DISMSH_ENSI_PB
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  
  SUBROUTINE DISMSH_ENSI_PB2(filnum,iload,nodes_pp,npes,numpe, &
                                  numvar,stress)

  !/****f* output/dismsh_ensi_pb2
  !*  NAME
  !*    SUBROUTINE: dismsh_ensi_pb2
  !*  SYNOPSIS
  !*    Usage:      CALL dismsh_ensi_pb2(filnum,iload,nodes_pp, &
  !*                                          numpe,numvar,stress) 
  !*  FUNCTION
  !*    ! DO NOT USE - UNDER DEVELOPMENT !
  !*    Write the values of a nodal variable to a binary file. This subroutine is 
  !*    parallel and requires MPI. The master processor collects all the data
  !*    from the slave processors.
  !*
  !*    This version writes out C binary for use in Paraview. The test program
  !*    is xx13.f90
  !*
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    filnum                  : Integer
  !*                            : File number to write
  !*
  !*    iload                   : Integer
  !*                            : Load step number
  !*
  !*    nodes_pp                : Integer
  !*                            : Number of nodes assigned to a process
  !*
  !*    npes                    : Integer
  !*                            : Number of processes
  !*
  !*    numpe                   : Integer
  !*                            : Process number
  !*
  !*    numvar                  : Integer
  !*                            : Number of components of the variable
  !*                              (1-scalar, 3-vector, 6-tensor)
  !*
  !*    stress(nodes_pp*numvar) : Real
  !*                            : Nodal variables to print
  !*                             
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*  AUTHOR
  !*    L. Margetts
  !*  CREATION DATE
  !*    23.07.2014
  !*  COPYRIGHT
  !*    (c) University of Manchester 2014
  !******
  !*
  !*/

  USE, INTRINSIC :: ISO_C_BINDING
  
  IMPLICIT NONE

  INTEGER, INTENT(IN)           :: filnum, iload, nodes_pp, npes, numpe, numvar
  REAL(iwp), INTENT(IN)         :: stress(:)
  INTEGER                       :: i, j, idx1, nod_r, bufsize1, bufsize2
  INTEGER                       :: ier, iproc, n, bufsize
  INTEGER                       :: statu(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: get(:)
  REAL(iwp), ALLOCATABLE        :: stress_r(:)

!------------------------------------------------------------------------------
! 1. Master processor writes out the results for its own assigned nodes
!------------------------------------------------------------------------------

  IF(numpe==1) THEN
!   DO i = 1,nodes_pp
      WRITE(filnum) real(stress,kind=c_float)
!   END DO
  END IF    
  
!------------------------------------------------------------------------------
! 2. Allocate arrays involved in communications
!------------------------------------------------------------------------------

  ALLOCATE(get(npes))
  ALLOCATE(stress_r(nodes_pp*numvar)) 

!------------------------------------------------------------------------------
! 3. Master processor populates the array "get" containing "nodes_pp" of 
!    every processor. Slave processors send this number to the master processor
!------------------------------------------------------------------------------

  get = 0
  get(1) = nodes_pp

  bufsize = 1

  DO i = 2,npes
    IF(numpe==i) THEN
      CALL MPI_SEND(nodes_pp,bufsize,MPI_INTEGER,0,i,MPI_COMM_WORLD,ier)
    END IF
    IF(numpe==1) THEN
      CALL MPI_RECV(nod_r,bufsize,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
      get(i) = nod_r
    END IF
  END DO

!------------------------------------------------------------------------------
! 4. Master processor receives the nodal variables from the other processors 
!    and writes them
!------------------------------------------------------------------------------

  DO iproc = 2,npes
    IF (numpe==iproc) THEN
      bufsize1 = nodes_pp*numvar
      CALL MPI_SEND(stress,bufsize1,MPI_REAL8,0,iproc,MPI_COMM_WORLD,ier)
    END IF
    IF (numpe==1) THEN
      bufsize2 = get(iproc)*numvar
      CALL MPI_RECV(stress_r,bufsize2,MPI_REAL8,iproc-1,iproc, &
                    MPI_COMM_WORLD,statu,ier)
!     n = 1
!     DO j = 2,iproc
!       n = n + get(j-1)
!     END DO

!     DO i = 1,get(iproc)
        n = get(iproc)
        WRITE(filnum) real(stress_r(1:n),kind=c_float)
!     END DO
    END IF
  END DO

!------------------------------------------------------------------------------
! 5. Deallocate communication arrays
!------------------------------------------------------------------------------

  DEALLOCATE(get,stress_r)


END SUBROUTINE DISMSH_ENSI_PB2

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  
  SUBROUTINE DISMSH_ENSI_PB3(filnum,iload,nodes_pp,npes,numpe, &
                                  numvar,stress)

  !/****f* output/dismsh_ensi_pb3
  !*  NAME
  !*    SUBROUTINE: dismsh_ensi_pb3
  !*  SYNOPSIS
  !*    Usage:      CALL dismsh_ensi_pb3(filnum,iload,nodes_pp, &
  !*                                          numpe,numvar,stress) 
  !*  FUNCTION
  !*    ! DO NOT USE - UNDER DEVELOPMENT !
  !*    Write the values of a nodal variable to a binary file. This subroutine is 
  !*    parallel and requires MPI. The master processor collects all the data
  !*    from the slave processors.
  !*
  !*    This version writes out C binary for use in Paraview. The test program
  !*    is xx13.f90
  !*
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    filnum                  : Integer
  !*                            : File number to write
  !*
  !*    iload                   : Integer
  !*                            : Load step number
  !*
  !*    nodes_pp                : Integer
  !*                            : Number of nodes assigned to a process
  !*
  !*    npes                    : Integer
  !*                            : Number of processes
  !*
  !*    numpe                   : Integer
  !*                            : Process number
  !*
  !*    numvar                  : Integer
  !*                            : Number of components of the variable
  !*                              (1-scalar, 3-vector, 6-tensor)
  !*
  !*    stress(nodes_pp*numvar) : Real
  !*                            : Nodal variables to print
  !*                             
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*  AUTHOR
  !*    L. Margetts
  !*  CREATION DATE
  !*    23.07.2014
  !*  COPYRIGHT
  !*    (c) University of Manchester 2014
  !******
  !*
  !*/

  USE, INTRINSIC :: ISO_C_BINDING
  
  IMPLICIT NONE

  INTEGER, INTENT(IN)           :: filnum, iload, nodes_pp, npes, numpe, numvar
  REAL(iwp), INTENT(IN)         :: stress(:)
  INTEGER                       :: i, j, idx1, nod_r, bufsize1, bufsize2
  INTEGER                       :: ier, iproc, n, bufsize
  INTEGER                       :: statu(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: get(:)
  REAL(iwp), ALLOCATABLE        :: stress_r(:)
  REAL(KIND=C_DOUBLE), ALLOCATABLE :: stress_r_c(:)

!------------------------------------------------------------------------------
! 1. Master processor writes out the results for its own assigned nodes
!------------------------------------------------------------------------------

  IF(numpe==1) THEN
!   DO i = 1,nodes_pp
      WRITE(filnum) real(stress,kind=c_float)
!   END DO
  END IF    
  
!------------------------------------------------------------------------------
! 2. Allocate arrays involved in communications
!------------------------------------------------------------------------------

  ALLOCATE(get(npes))
  ALLOCATE(stress_r(nodes_pp*numvar)) 
  ALLOCATE(stress_r_c(nodes_pp*numvar)) 

!------------------------------------------------------------------------------
! 3. Master processor populates the array "get" containing "nodes_pp" of 
!    every processor. Slave processors send this number to the master processor
!------------------------------------------------------------------------------

  get = 0
  get(1) = nodes_pp

  bufsize = 1

  DO i = 2,npes
    IF(numpe==i) THEN
      CALL MPI_SEND(nodes_pp,bufsize,MPI_INTEGER,0,i,MPI_COMM_WORLD,ier)
    END IF
    IF(numpe==1) THEN
      CALL MPI_RECV(nod_r,bufsize,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
      get(i) = nod_r
    END IF
  END DO

!------------------------------------------------------------------------------
! 4. Master processor receives the nodal variables from the other processors 
!    and writes them
!------------------------------------------------------------------------------

  DO iproc = 2,npes
    IF (numpe==iproc) THEN
      bufsize1 = nodes_pp*numvar
      CALL MPI_SEND(stress,bufsize1,MPI_REAL8,0,iproc,MPI_COMM_WORLD,ier)
    END IF
    IF (numpe==1) THEN
      bufsize2 = get(iproc)*numvar
      CALL MPI_RECV(stress_r,bufsize2,MPI_REAL8,iproc-1,iproc, &
                    MPI_COMM_WORLD,statu,ier)
!     n = 1
!     DO j = 2,iproc
!       n = n + get(j-1)
!     END DO

!     DO i = 1,get(iproc)
        n = get(iproc)
!        WRITE(filnum) real(stress_r(1:n),kind=c_float)
        stress_r_c = stress_r
        WRITE(filnum) real(stress_r_c(1:n),kind=c_double)
!     END DO
    END IF
  END DO

!------------------------------------------------------------------------------
! 5. Deallocate communication arrays
!------------------------------------------------------------------------------

  DEALLOCATE(get,stress_r)
  
  END SUBROUTINE DISMSH_ENSI_PB3

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  
  SUBROUTINE DISMSH_ENSI_PB2_INT(filnum,iload,nodes_pp,npes,numpe, &
                                  numvar,stress)

  !/****f* output/dismsh_ensi_pb2_int
  !*  NAME
  !*    SUBROUTINE: dismsh_ensi_pb2_int
  !*  SYNOPSIS
  !*    Usage:      CALL dismsh_ensi_pb2_int(filnum,iload,nodes_pp,npes &
  !*                                          numpe,numvar,stress) 
  !*  FUNCTION
  !*    ! DO NOT USE - UNDER DEVELOPMENT !
  !*    Write the values of a nodal variable to a binary file. This subroutine is 
  !*    parallel and requires MPI. The master processor collects all the data
  !*    from the slave processors.
  !*
  !*    This version writes out C binary for use in Paraview. The test program
  !*    is xx12.f90
  !*
  !*  INPUTS
  !*    The following arguments have the INTENT(IN) attribute:
  !*
  !*    filnum                  : Integer
  !*                            : File number to write
  !*
  !*    iload                   : Integer
  !*                            : Load step number
  !*
  !*    nodes_pp                : Integer
  !*                            : Number of nodes assigned to a process
  !*
  !*    npes                    : Integer
  !*                            : Number of processes
  !*
  !*    numpe                   : Integer
  !*                            : Process number
  !*
  !*    numvar                  : Integer
  !*                            : Number of components of the variable
  !*                              (1-scalar, 3-vector, 6-tensor)
  !*
  !*    stress(nodes_pp*numvar) : Integer
  !*                            : Nodal variables to print
  !*                             
  !*
  !*    The following arguments have the INTENT(OUT) attribute:
  !*
  !*  AUTHOR(s)
  !*    L. Margetts
  !*    Ll.M. Evans
  !*  CREATION DATE
  !*    23.07.2014
  !*  COPYRIGHT
  !*    (c) University of Manchester 2014
  !******
  !*
  !*/

  USE, INTRINSIC :: ISO_C_BINDING
  
  IMPLICIT NONE

  INTEGER, INTENT(IN)           :: filnum, iload, nodes_pp, npes, numpe, numvar, stress(:)
  INTEGER                       :: i, j, idx1, nod_r, bufsize1, bufsize2
  INTEGER                       :: ier, iproc, n, bufsize
  INTEGER                       :: statu(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE          :: get(:),stress_r(:)

!------------------------------------------------------------------------------
! 1. Master processor writes out the results for its own assigned nodes
!------------------------------------------------------------------------------

  IF(numpe==1) THEN
!   DO i = 1,nodes_pp
      WRITE(filnum) int(stress,kind=c_int)
!   END DO
  END IF    
  
!------------------------------------------------------------------------------
! 2. Allocate arrays involved in communications
!------------------------------------------------------------------------------

  ALLOCATE(get(npes))
  ALLOCATE(stress_r(nodes_pp*numvar)) 

!------------------------------------------------------------------------------
! 3. Master processor populates the array "get" containing "nodes_pp" of 
!    every processor. Slave processors send this number to the master processor
!------------------------------------------------------------------------------

  get = 0
  get(1) = nodes_pp

  bufsize = 1

  DO i = 2,npes
    IF(numpe==i) THEN
      CALL MPI_SEND(nodes_pp,bufsize,MPI_INTEGER,0,i,MPI_COMM_WORLD,ier)
    END IF
    IF(numpe==1) THEN
      CALL MPI_RECV(nod_r,bufsize,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
      get(i) = nod_r
    END IF
  END DO

!------------------------------------------------------------------------------
! 4. Master processor receives the nodal variables from the other processors 
!    and writes them
!------------------------------------------------------------------------------

  DO iproc = 2,npes
    IF (numpe==iproc) THEN
      bufsize1 = nodes_pp*numvar
      CALL MPI_SEND(stress,bufsize1,MPI_INTEGER,0,iproc,MPI_COMM_WORLD,ier)
    END IF
    IF (numpe==1) THEN
      bufsize2 = get(iproc)*numvar
      CALL MPI_RECV(stress_r,bufsize2,MPI_INTEGER,iproc-1,iproc, &
                    MPI_COMM_WORLD,statu,ier)
!     n = 1
!     DO j = 2,iproc
!       n = n + get(j-1)
!     END DO

!     DO i = 1,get(iproc)
        n = get(iproc)
        WRITE(filnum) int(stress_r(1:n),kind=c_int)
!     END DO
    END IF
  END DO

!------------------------------------------------------------------------------
! 5. Deallocate communication arrays
!------------------------------------------------------------------------------

  DEALLOCATE(get,stress_r)


END SUBROUTINE DISMSH_ENSI_PB2_INT

END MODULE OUTPUT
