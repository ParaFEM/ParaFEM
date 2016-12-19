PROGRAM normdiff
  IMPLICIT NONE
  INTEGER, PARAMETER :: string_length = 1024
  CHARACTER(len=string_length) :: filename1,filename2,numvar_s,reltol_s
  CHARACTER(len=string_length) :: info1,info2,line1,line2
  INTEGER :: numvar,f1,f2,stat1,stat2,iload1,iload2,i1,i2,i,j
  DOUBLE PRECISION, DIMENSION(6) :: var1,var2,diff
  DOUBLE PRECISION :: norm1,norm2,mean_norm,absdiff,reldiff,reltol

  reltol = 1.0e-2

  IF (command_argument_count() /= 4) THEN
    WRITE(0,*) "normdiff:  Usage is"
    WRITE(0,*) "normdiff file1 file2 numvar reltol"
    STOP 2
  END IF
  CALL get_command_argument(1,filename1)
  CALL get_command_argument(2,filename2)
  CALL get_command_argument(3,numvar_s)
  READ(numvar_s,*) numvar
  CALL get_command_argument(4,reltol_s)
  READ(reltol_s,*) reltol

  IF (numvar==1) THEN
  ELSE IF (numvar==3) THEN
  ELSE IF (numvar==4) THEN
    WRITE(0,*) "normdiff:  Don't know what to do when numvar = 4"
    STOP 1
  ELSE IF (numvar==6) THEN
  END IF
  
  IF (filename1 == filename2) THEN
    WRITE(0,*) "normdiff:  file1 and file2 cannot be the same file"
    STOP 1
  END IF

  OPEN(newunit=f1,file=filename1,status='old',action='read',iostat=stat1)
  OPEN(newunit=f2,file=filename2,status='old',action='read',iostat=stat2)

  info1 = ""
  info2 = ""

  READ(f1,'(a)',iostat=stat1) info1
  READ(f2,'(a)',iostat=stat2) info2
  IF (info1 /= info2) THEN
    WRITE(0,'(a)') "normdiff:  info1 /= info2"
    STOP 1
  END IF
  
  outer: DO
    READ(f1,'(i0)',iostat=stat1) iload1
    READ(f2,'(i0)',iostat=stat2) iload2
    IF (iload1 /= iload2) THEN
      WRITE(0,'(a)') "normdiff:  iload1 /= iload2"
      STOP 1
    END IF
    
    inner: DO
      READ(f1,'(a)',iostat=stat1) line1
      READ(f2,'(a)',iostat=stat2) line2

      IF (stat1 < 0 .AND. stat2 < 0) THEN
        ! both files ended at same place
        EXIT outer
      ELSE IF (stat1 < 0 .OR. stat2 < 0) THEN
        ! one file ended before the other
        WRITE(*,'(a)') "normdiff:  one file longer than the other"
        EXIT outer
      END IF
      
      IF (line1 == info1 .AND. line2 == info2) THEN
        EXIT inner
      END IF

      IF (numvar==1) THEN
        READ(line1,2001,iostat=stat1)i1,var1(1:numvar)
        READ(line2,2001,iostat=stat2)i2,var2(1:numvar)
      ELSE IF (numvar==3) THEN
        READ(line1,2003,iostat=stat1)i1,var1(1:numvar)
        READ(line2,2003,iostat=stat2)i2,var2(1:numvar)
      ELSE IF (numvar==4) THEN
        READ(line1,2004,iostat=stat1)i1,var1(1:numvar)
        READ(line2,2004,iostat=stat2)i2,var2(1:numvar)
      ELSE IF (numvar==6) THEN
        READ(line1,2006,iostat=stat1)i1,var1(1:numvar)
        READ(line2,2006,iostat=stat2)i2,var2(1:numvar)
      END IF
2001  FORMAT(i8,1(1p,e12.4))
2003  FORMAT(i8,3(1p,e12.4))
2004  FORMAT(i8,4(1p,e12.4))
2006  FORMAT(i8,6(1p,e12.4))

      IF (stat1 > 0) THEN
        WRITE(0,'(a)') "normdiff:  read error on file1"
        STOP 1
      END IF
      IF (stat2 > 0) THEN
        WRITE(0,'(a)') "normdiff:  read error on file2"
        STOP 1
      END IF

      diff(1:numvar) = (var1(1:numvar) - var2(1:numvar))
      IF (numvar==1) THEN
        norm1 = SQRT(SUM(var1(1:1)**2))
        norm2 = SQRT(SUM(var2(1:1)**2))
        absdiff = SQRT(SUM(diff(1:1)**2))
      ELSE IF (numvar==3) THEN
        norm1 = SQRT(SUM(var1(1:3)**2))
        norm2 = SQRT(SUM(var2(1:3)**2))
        absdiff = SQRT(SUM(diff(1:3)**2))
      ELSE IF (numvar==4) THEN
      ELSE IF (numvar==6) THEN
        norm1 = SQRT(SUM(var1(1:3)**2) + 2 * SUM(var1(4:6)**2))
        norm2 = SQRT(SUM(var2(1:3)**2) + 2 * SUM(var2(4:6)**2))
        absdiff = SQRT(SUM(diff(1:3)**2) + 2 * SUM(diff(4:6)**2))
      END IF
      mean_norm = (norm1 + norm2) / 2
      IF (mean_norm /= 0.0) THEN
        reldiff = absdiff / mean_norm
      ELSE
        reldiff = 0.0
      END IF
      IF (reldiff > reltol) THEN
        WRITE(*,'(8(1p,e14.6))') absdiff,mean_norm,var1(1:numvar)
        WRITE(*,'(28x,6(1p,e14.6))')               var2(1:numvar)
      END IF
    END DO inner
  END DO outer
  
  CLOSE(f1)
  CLOSE(f2)
   
END PROGRAM normdiff
