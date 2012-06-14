echo
rem  *** EDIT THE NEXT THREE LINES IF NECESSARY ***
rem
rem
set PARAFEM=c:\parafem\parafem
set G95=c:\MinGW\g95\bin\g95
set AR=c:\MinGW\g95\bin\ar
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
%G95% -c -r8 mpi_stubs.f90
%AR% -r libmpi_stubs.a mpi_stubs.o 
move libmpi_stubs.a %PARAFEM%\lib
del *.o
rem
rem
rem *** BUILD GAF77 LIBRARY ***
rem
rem
cd /D %PARAFEM%\src\libraries\gaf77
del *.a *.o *.mod
%G95% -c -r8 *.f
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
%G95% -c -r8 precision.f90
%G95% -c -r8 global_variables.f90
%G95% -c -r8 geometry.f90
%G95% -c -r8 elements.f90
%G95% -c -r8 steering.f90
%G95% -c -r8 timing.f90
%G95% -c -r8 partition.f90
%G95% -c -r8 plasticity.f90
%G95% -c -r8 fluid.f90
move *.o %PARAFEM%\src\modules\mpi
move *.mod %PARAFEM%\src\modules\mpi
cd /D %PARAFEM%\src\modules\mpi
%G95% -c -r8 mp_interface.f90
%G95% -c -r8 maths.f90
%G95% -c -r8 input.f90
%G95% -c -r8 output.f90
%G95% -c -r8 gather_scatter.f90
%G95% -c -r8 loading.f90
%G95% -c -r8 large_strain.f90
%G95% -c -r8 pcg.f90
%G95% -c -r8 bicg.f90
%AR% -r libparafem.a *.o
move libparafem.a %PARAFEM%\lib
move *.mod %PARAFEM%\include
del *.o 
rem
rem
rem *** BUILD RFEMBC ***
rem
rem
cd /D %PARAFEM%\src\tools\preprocessing\rfembc
del rfembc.exe rfembc.o
%G95% -r8 rfembc.f90 -L%PARAFEM%\lib -lparafem -I%PARAFEM%\include -o rfembc.exe
del *.o *.a
move rfembc.exe %PARAFEM%\bin
rem
rem 
rem *** BUILD RFEMFIELD ***
rem
rem
cd /D %PARAFEM%\src\tools\preprocessing\rfemfield
del rfemfield.exe rfemfield.o
%G95% -r8 rfemfield.f90 -L%PARAFEM%\lib -lgaf77 -lparafem -I%PARAFEM%\include -o rfemfield.exe
del *.o *.a
move rfemfield.exe %PARAFEM%\bin
rem
rem *** BUILD RFEMSOLVE ***
rem
cd /D %PARAFEM%\src\programs\rfem
%G95% -r8 rfemsolve.f90 -L%PARAFEM%\lib -lmpi_stubs -lparafem -I%PARAFEM%\include -o rfemsolve.exe
del *.o *.a
move rfemsolve.exe %PARAFEM%\bin
cd /D %PARAFEM%

