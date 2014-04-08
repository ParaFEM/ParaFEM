@ECHO OFF
PATH=%PARAFEM_HOME%\OpenMPI_v1.6.2-win32\bin;%PATH%
ECHO
ECHO -------------------------------------------------------------------------------
ECHO  %1 - IF PROMPTED TO ALLOW FIREWALL ACCESS, PLEASE DO SO
ECHO -------------------------------------------------------------------------------
SET PARAFEM_HOME=C:\ParaFEM_Demo
cd "%PARAFEM_HOME%\bin"

mpirun --mca btl self,sm -np 1 %1.exe ..\examples\%1\demo\%1_demo
rem mpirun --mca btl self,sm -np %NUMBER_OF_PROCESSORS% ..\examples\%1\demo\%1_demo

ECHO.
ECHO -------------------------------------------------------------------------------
ECHO Execution complete.  You can now view the output in ParaView
ECHO -------------------------------------------------------------------------------
pause
"%PARAFEM_HOME%\ParaView 4.1.0\bin\paraview" --state=%PARAFEM_HOME%\examples\%1\demo\%1_demo.pvsm

