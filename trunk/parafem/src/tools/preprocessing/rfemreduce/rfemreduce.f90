PROGRAM rfemreduce

!/****h* tools/preprocessing/rfem
!*  NAME
!*    PROGRAM rfemreduce
!*  FUNCTION
!*    
!*    
!*    
!*  AUTHOR
!*    Louise Lever
!*  COPYRIGHT
!*    (c) University of Manchester 2012
!****
!*/

  IMPLICIT NONE
  
!------------------------------------------------------------------------------
! 1. Declare variables used in the main program
!------------------------------------------------------------------------------

  INTEGER                :: argc,iargc
  CHARACTER(LEN=100)     :: job_name,fname,program_name
  CHARACTER(LEN=100)     :: num_instances_arg,num_bins_arg
  INTEGER                :: num_instances, num_bins
  INTEGER                :: i, bin
  CHARACTER(LEN=15)      :: inst
  REAL                   :: bin_width
  
!------------------------------------------------------------------------------
! 2. Declare dynamic arrays
!------------------------------------------------------------------------------   

  INTEGER, ALLOCATABLE   :: inst_count(:)
  REAL, ALLOCATABLE      :: inst_count_percent(:)
  INTEGER, ALLOCATABLE   :: bin_count(:)

!------------------------------------------------------------------------------
! 3. Read model_job_name and rfem_job_name from the command line
!------------------------------------------------------------------------------

  argc = iargc()
  IF( argc /= 3 ) THEN
     PRINT*
     PRINT*, "Usage:  rfemreduce <job_name> <num-instances> <num-bins>"
     PRINT*
     PRINT*, "        program expects as input:"
     PRINT*, "          <job_name>-###.tnc"
     PRINT*
     PRINT*, "        and outputs:"
     PRINT*, "          <job-name>.hist"
     PRINT*
     STOP
  END IF
  CALL GETARG(1, job_name)
  CALL GETARG(2, num_instances_arg)
  CALL GETARG(3, num_bins_arg)

  READ(num_instances_arg,*) num_instances
  READ(num_bins_arg,*) num_bins
  bin_width = 100.0 / num_bins
  
!------------------------------------------------------------------------------
! 4. Read Threshold Node Count result data files and store in inst_count
!    and bin the data to get histogram
!------------------------------------------------------------------------------

  ALLOCATE(inst_count(num_instances))
  ALLOCATE(inst_count_percent(num_instances))
  ALLOCATE(bin_count(num_bins))

  bin_count = 0

  DO i=1,num_instances
     WRITE(inst,'(1I0)') i
     fname = job_name(1:LEN_TRIM(job_name)) // "-" // inst(1:LEN_TRIM(inst)) // ".tnc"
     OPEN(10, file=fname, status='old', action='read')
     READ(10,*) inst_count(i), inst_count_percent(i)

     bin = min(num_bins,int(inst_count_percent(i) / bin_width) + 1)
     bin_count(bin) = bin_count(bin) + 1

     CLOSE(10)
  END DO

!------------------------------------------------------------------------------
! 6. Output histogram data for num_bins bins
!------------------------------------------------------------------------------

  fname = job_name(1:INDEX(job_name," ")-1) // ".hist"
  OPEN(20,FILE=fname,STATUS='REPLACE',ACTION='WRITE')
  
  WRITE(20,'(1I12)') num_bins

  DO i = 1,num_bins
     WRITE(20,'(1I12)') bin_count(i)
  END DO

  CLOSE(20)
  
!------------------------------------------------------------------------------
! . Cleanup
!------------------------------------------------------------------------------

  DEALLOCATE(inst_count)
  DEALLOCATE(inst_count_percent)
  DEALLOCATE(bin_count)
  
END PROGRAM rfemreduce
