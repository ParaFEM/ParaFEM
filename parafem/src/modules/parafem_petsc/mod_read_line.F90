MODULE mod_read_line
  USE iso_fortran_env
  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE read_line(unit,line,iostat)      
    !/****f* input/read_line
    !*  NAME
    !*    SUBROUTINE: read_line
    !*  SYNOPSIS
    !*    Usage:      CALL read_line(unit,line,iostat)
    !*  FUNCTION
    !*    This subroutine reads a line of any length into a deferred-length
    !*    character string
    !*  ARGUMENTS
    !*    INTENT(IN)
    !*  INPUTS
    !*
    !*    unit               : Integer
    !*                         Formatted input unit
    !*
    !*    INTENT(OUT)
    !*
    !*    line               : Character string
    !*                         The line of input
    !*    iostat             : Integer
    !*                         == 0          : success
    !*                          > 0          : error occurred
    !*                         == IOSTAT_END : success but end-of-file before
    !*                                         end-of-record
    !*
    !*  AUTHOR
    !*    Mark Filipiak 
    !*  CREATION DATE
    !*    17.03.2016
    !*  MODIFICATION HISTORY
    !*    Version 1, 17.03.2016, Mark Filipiak
    !*  COPYRIGHT
    !*    (c) University of Edinburgh 2016
    !******
    !*  Place remarks that should not be included in the documentation here.
    !*
    !*  The iostat argument is not optional:  you have to get to the end of the
    !*  record to read the line.
    !*
    !*  There is no test for running out of memory.
    !*/

    INTEGER,      INTENT(IN)               :: unit
    CHARACTER(:), INTENT(OUT), ALLOCATABLE :: line
    INTEGER,      INTENT(OUT), OPTIONAL    :: iostat

    CHARACTER(1024) :: buffer
    INTEGER         :: size
    
    iostat = 0
    line = ""
    DO
      READ (unit,'(A)',ADVANCE='NO',IOSTAT=iostat,SIZE=size) buffer
      IF (iostat > 0) THEN
        RETURN
      END IF
      line = line // buffer(:size)
      ! If end-of-record or end-of-file then the line has been read
      ! successfully.
      IF (iostat < 0) THEN
        IF (iostat == IOSTAT_EOR) iostat = 0
        RETURN
      END IF
    END DO
  END SUBROUTINE read_line
    
END MODULE mod_read_line
