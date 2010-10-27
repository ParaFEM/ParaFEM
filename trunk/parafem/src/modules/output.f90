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
  !*    Subroutine             Purpose
  !*    WRITE_P121             Writes out basic program data and timing info
  !*    WRITE_P126             Writes out basic program data and timing info  
  !*    WRITE_P129             Writes out basic program data and timing info
  !*    WRITE_NODAL_VARIABLE   Writes out results computed at the nodes
  !*    JOB_NAME_ERROR         Writes error message if job_name is missing
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
  !*    (c) University of Manchester 2010
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

    WRITE(11,'(/A)')   "PROGRAM SECTION EXECUTION TIMES                 "     
 
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
                           ((timest(2)-timest(1))/(timest(10)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute steering array and neq              ",&
                           timest(3)-timest(2),                               &
                          ((timest(3)-timest(2))/(timest(10)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute interprocessor communication tables ",&
                           timest(4)-timest(3),                               &
                          ((timest(4)-timest(3))/(timest(10)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Allocate neq_pp arrays                      ",&
                           timest(5)-timest(4),                               &
                          ((timest(5)-timest(4))/(timest(10)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Compute element stiffness matrices          ",&
                            timest(6)-timest(5),                              &
                          ((timest(6)-timest(5))/(timest(10)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Build the preconditioner                    ",&
                           timest(7)-timest(6),                               &
                          ((timest(7)-timest(6))/(timest(10)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Get starting r                              ",&
                           timest(8)-timest(7),                               &
                          ((timest(8)-timest(7))/(timest(10)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Solve equations                             ",&
                           timest(9)-timest(8),                               &
                           ((timest(9)-timest(8))/(timest(10)-timest(1)))*100  
    WRITE(11,'(A,F12.6,F8.2)') "Output results                              ",&
                           timest(10)-timest(9),                              &
                          ((timest(10)-timest(9))/(timest(10)-timest(1)))*100  
    WRITE(11,'(A,F12.6,A/)')  "Total execution time                        ",  &
                          timest(10)-timest(1),"  100.00"
    CLOSE(11)
    
  END IF
  
  RETURN
  END SUBROUTINE WRITE_P121

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
  !*    (c) University of Manchester 2007-2010
  !******
  !*  Place remarks that should not be included in the documentation here.
  !*/


 IMPLICIT none

 INTEGER, INTENT(IN)           :: numpe        ! processor ID number
 CHARACTER(LEN=50), INTENT(IN) :: program_name ! program name 
 CHARACTER(LEN=50)             :: fname

 IF (numpe==1) THEN
   fname   = program_name(1:INDEX(program_name, " ")-1)//".res"
   OPEN (11, file=fname, status='replace', action='write')
   WRITE(11,'(/4A/)') "Fatal error: ", TRIM(program_name), " did not run",     &
                      " - job_name is missing."
   WRITE(11,'(3A/)')  "Usage: ", TRIM(program_name), " job_name"
   WRITE(11,'(4A)')   "       ", TRIM(program_name), " will read job_name,",   &
                      " <job_name>.dat,"       
   WRITE(11,'(A)')    "       <job_name>.d, <job_name>.bnd, <job_name>.lds,"    
   WRITE(11,'(A/)')   "       and write <job_name>.dis"
   WRITE(11,'(A/)')   "       If analysis_type > 2, <job_name>.mat is expected."
   CLOSE(11)
 END IF

 CALL shutdown()
 STOP
 
 END SUBROUTINE job_name_error
 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------ 

END MODULE OUTPUT
