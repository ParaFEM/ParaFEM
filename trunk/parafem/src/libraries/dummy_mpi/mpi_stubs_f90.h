!
!  Dummy parameters for MPI F90 stubs
!
  integer mpi_comm_world
  parameter ( mpi_comm_world = 0 )
!
!  Return values
!
  integer mpi_failure
  parameter ( mpi_failure = 1 )
  integer mpi_success
  parameter ( mpi_success = 0 )
!
!  recv message status
!
  integer mpi_status_size
  parameter ( mpi_status_size = 4 )
  integer mpi_source
  parameter ( mpi_source = 1 )
  integer mpi_tag
  parameter ( mpi_tag = 2 )
  integer mpi_count
  parameter ( mpi_count = 3 )
!
!  recv flags
!
  integer mpi_any_source
  parameter ( mpi_any_source = -1 )
  integer mpi_any_tag
  parameter ( mpi_any_tag = -1 )
!
!  data types and sizes
!
  integer MPI_INTEGER
  parameter ( MPI_INTEGER = 28 )
  integer MPI_REAL
  parameter ( MPI_REAL = 11 )
  integer mpi_double_precision
  parameter ( mpi_double_precision = 3 )
  integer mpi_logical
  parameter ( mpi_logical = 4 )
  integer mpi_character
  parameter ( mpi_character = 5 )
!
!  allreduce operations
!
  INTEGER MPI_SUM, MPI_MAX, MPI_MIN, MPI_PROD
  PARAMETER (MPI_MAX=100,MPI_MIN=101,MPI_SUM=102,MPI_PROD=103)

!
!  timer
!
  double precision mpi_wtime

