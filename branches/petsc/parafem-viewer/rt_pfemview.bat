@echo off
setlocal
echo "Starting ParaFEM Viewer ..."

set XP_LICENSE_SERVER=thomas.rcs.manchester.ac.uk:33333
set XP_FEATURE=VIZ_EXPRESS
set MACHINE=pc
set PATH=runtime\bin\%MACHINE%;%PATH%
cd %FEAR_HOME%
bin\%MACHINE%\express -novcp %1 %2 %3
