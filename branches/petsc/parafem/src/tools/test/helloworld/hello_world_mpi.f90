PROGRAM hello_world_mpi
 
! Program hello_world_mpi runs the classic "Hello World" test using
! multiple MPI processes.

include 'mpif.h'

INTEGER   :: numtasks, rank, ierr, rc, len, i
CHARACTER* (MPI_MAX_PROCESSOR_NAME) name

CALL MPI_INIT(ierr)
IF (ierr .ne. MPI_SUCCESS) THEN
  PRINT *, "Error starting MPI program. Terminating."
  CALL MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
END IF

! Find the number of cpus this job is using

CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

! Find the rank of this cpu

CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

! Find the name of this cpu

CALL MPI_GET_PROCESSOR_NAME(name, len, ierr)
IF(ierr .ne. MPI_SUCCESS) THEN
  PRINT *, "Error getting the name of the cpu. Terminating program."
  CALL MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
END IF

PRINT "('helloworld.f90: Number of tasks=',I3,' My rank=',I3,' My name=',A,& 
        '')",numtasks,rank,trim(name)

! Tell the MPI library to release all resources it is using:

CALL MPI_FINALIZE(ierr)

END PROGRAM

