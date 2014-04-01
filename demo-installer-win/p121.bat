ECHO -------------------------------------------------------------------------------
ECHO DEPENDING ON THE SPEED OF YOUR PC THIS MAY TAKE 5-10 MINS TO COMPLETE
ECHO IF PROMPTED TO ALLOW FIREWALL ACCESS, PLEASE DO SO
ECHO -------------------------------------------------------------------------------
SET PARAFEM_HOME=@INSTDIR@
PATH=%PARAFEM_HOME%\OpenMPI_v1.6.2-win32\bin;%PATH%
mpirun --mca btl self,sm -np 1 p121.exe ..\examples\p121\demo\p121_small
rem mpirun --mca btl self,sm -np %NUMBER_OF_PROCESSORS% p121.exe ..\examples\p121\demo\p121_small

ECHO.
ECHO -------------------------------------------------------------------------------
ECHO Execution complete.  You can now view the output in ParaView
ECHO -------------------------------------------------------------------------------
pause
"%PARAFEM_HOME%\ParaView 4.1.0\bin\paraview" --state=%PARAFEM_HOME%\examples\p121\demo\p121_small.pvsm