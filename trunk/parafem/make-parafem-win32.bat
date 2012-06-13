echo
rem  *** EDIT THE NEXT THREE LINES IF NECESSARY ***
rem
rem
set PARAFEM=c:\parafem\parafem
set G95=c:\g95
set VER=1
rem
rem 
rem *** BUILD DUMMY_MPI LIBRARY ***
rem
rem
cd /D %PARAFEM%\bin
del *.a *.o *.mod rfemsolve.exe
cd /D %PARAFEM%\src\libraries\dummy_mpi
del *.a *.o *.mod
g95 -c -r8 mpi_stubs.f90
ar -r libmpi_stubs.a mpi_stubs.o 
move libmpi_stubs.a %PARAFEM%\bin
move *.o %PARAFEM%\bin
rem
rem
rem *** BUILD GAF77 LIBRARY ***
rem
rem
cd /D %PARAFEM%\src\libraries\gaf77
del *.a *.o *.mod
g95 -c -r8 *.f
ar -r gaf77.a *.o
move gaf77.a %PARAFEM%\bin
move *.o %PARAFEM%\bin
rem
rem *** BUILD PARAFEM LIBRARY ***
rem
cd /D %PARAFEM%\src\modules\shared
del *.a *.o *.mod
cd /D %PARAFEM%\src\modules\mpi
del *.a *.o *.mod
cd /D %PARAFEM%\src\modules\shared
g95 -c -r8 precision.f90
g95 -c -r8 global_variables.f90
g95 -c -r8 geometry.f90
g95 -c -r8 elements.f90
g95 -c -r8 steering.f90
g95 -c -r8 timing.f90
g95 -c -r8 partition.f90
g95 -c -r8 plasticity.f90
g95 -c -r8 fluid.f90
move *.o %PARAFEM%\src\modules\mpi
move *.mod %PARAFEM%\src\modules\mpi
cd /D %PARAFEM%\src\modules\mpi
g95 -c -r8 mp_interface.f90
g95 -c -r8 maths.f90
g95 -c -r8 input.f90
g95 -c -r8 output.f90
g95 -c -r8 gather_scatter.f90
g95 -c -r8 loading.f90
g95 -c -r8 large_strain.f90
g95 -c -r8 pcg.f90
g95 -c -r8 bicg.f90
ar -r libparafem.a *.o
move libparafem.a %PARAFEM%\bin
move *.mod %PARAFEM%\src\programs\rfem
del *.o 
rem 
rem *** BUILD RFEMFIELD ***
rem
rem *** BUILD RFEMSOLVE ***
rem
cd /D %PARAFEM%\src\programs\rfem
g95 -r8 rfemsolve.f90 -L%PARAFEM%\bin\ -lmpi_stubs -lparafem -I%PARAFEM%\src\programs\rfem -o rfemsolve.exe
del *.o *.a
move rfemsolve.exe %PARAFEM%\bin
cd /D %PARAFEM%

