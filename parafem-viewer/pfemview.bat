@echo off

echo "Starting ParaFEM Viewer ..."

cd %FEAR_HOME%
bin\%MACHINE%\express apps\pfemview.v %1 %2 %3


