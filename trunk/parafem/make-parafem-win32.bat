echo
rem  *** EDIT THE NEXT THREE LINES IF NECESSARY ***
rem
rem
set PARAFEM=d:\parafem\parafem
set G95=c:\g95\bin\g95
set AR=c:\g95\bin\ar
set VER=1
rem
rem 
rem *** CLEAN ***
rem
rem
cd /D %PARAFEM%\bin
del *.a *.o *.mod *.exe
cd /D %PARAFEM%\lib
del *.a
rem
rem
rem *** BUILD DUMMY_MPI LIBRARY ***
rem
rem
cd /D %PARAFEM%\src\libraries\dummy_mpi
del *.a *.o *.mod
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check mpi_stubs.f90
%AR% -r libmpi_stubs.a mpi_stubs.o
move libmpi_stubs.a %PARAFEM%\lib
move mpi_stubs.mod %PARAFEM%\include
del *.o
rem
rem
rem *** BUILD GAF77 LIBRARY ***
rem
rem
cd /D %PARAFEM%\src\libraries\gaf77
del *.a *.o *.mod
%G95% -c -r8 -i4 -ftrace=full -Wall -fbounds-check *.f
%AR% -r libgaf77.a *.o
move libgaf77.a %PARAFEM%\lib
del *.o
rem
rem *** BUILD PARAFEM LIBRARY ***
rem
cd /D %PARAFEM%\src\modules\shared
del *.a *.o *.mod
cd /D %PARAFEM%\src\modules\mpi
del *.a *.o *.mod
cd /D %PARAFEM%\src\modules\shared
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check precision.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check global_variables.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check geometry.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check elements.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check steering.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check timing.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check partition.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check plasticity.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check fluid.f90
move *.o %PARAFEM%\src\modules\mpi
move *.mod %PARAFEM%\src\modules\mpi
cd /D %PARAFEM%\src\modules\mpi
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check mpi_wrapper.f90 -I../../include
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check mp_interface.f90 -I../../libraries/dummy_mpi/include
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check maths.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check input.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check output.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check gather_scatter.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check loading.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check large_strain.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check pcg.f90
%G95% -c -r8 -i4 -Wall -ftrace=full -fbounds-check bicg.f90
%AR% -r libparafem.a *.o
move libparafem.a %PARAFEM%\lib
move *.mod %PARAFEM%\include
del *.o 
rem
rem
rem *** BUILD RFEMCUBE ***
rem
rem
cd /D %PARAFEM%\src\tools\preprocessing\rfemcube
del rfemcube.exe rfemcube.o
%G95% -r8 -i4 -Wall -ftrace=full -fbounds-check rfemcube.f90 -L%PARAFEM%\lib -lparafem -I%PARAFEM%\include -o rfemcube.exe
del *.o *.a
move rfemcube.exe %PARAFEM%\bin
rem
rem
rem *** BUILD RFEMBC ***
rem
rem
cd /D %PARAFEM%\src\tools\preprocessing\rfembc
del rfembc.exe rfembc.o
%G95% -r8 -i4 -Wall -ftrace=full -fbounds-check rfembc.f90 -L%PARAFEM%\lib -lparafem -I%PARAFEM%\include -o rfembc.exe
del *.o *.a
move rfembc.exe %PARAFEM%\bin
rem
rem 
rem *** BUILD RFEMFIELD ***
rem
rem
cd /D %PARAFEM%\src\tools\preprocessing\rfemfield
del rfemfield.exe rfemfield.o
%G95% -r8 -i4 -Wall -ftrace=full -fbounds-check rfemfield.f90 -L%PARAFEM%\lib -lgaf77 -lparafem -I%PARAFEM%\include -o rfemfield.exe
del *.o *.a
move rfemfield.exe %PARAFEM%\bin
rem
rem *** BUILD RFEMSOLVE ***
rem
cd /D %PARAFEM%\src\programs\rfem
%G95% -r8 -i4 -Wall -ftrace=full -fbounds-check rfemsolve.f90 -L%PARAFEM%\lib -lmpi_stubs -lparafem -I%PARAFEM%\include -o rfemsolve.exe
del *.o *.a
move rfemsolve.exe %PARAFEM%\bin
cd /D %PARAFEM%

:end
