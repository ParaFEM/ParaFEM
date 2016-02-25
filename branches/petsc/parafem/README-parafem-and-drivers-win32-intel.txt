===============================================================================
Building ParaFEM and 5th Edition example driver programs on Microsoft Windows
J.Robinson@software.ac.uk
===============================================================================

OS Prerequisites
----------------

Microsoft Windows XP or newer (Tested on Windows 7, 8 and 2008 Server).
	
Software Prerequisites
----------------------

1. Download and install the Microsoft Visual Studio 2010 Shell (Isolated) 
Redistributable Package [1]. This is required by to use the Intel Fortran 
compiler at the command line.

2. Download and install MS Visual Studio C++ 2010 Express [2]. This provides 
the linker and libraries required by ifort.

3. For performance reasons, the preferred Fortran compiler is ifort.  
Download and install Intel Visual Fortran Studio XE for Windows [3].  (This is 
also referred to as "Intel Visual Fortran Composer XE"). A 30 day trial version
of the compiler is available. You can accept the default options and ignore any
warnings about missing debugger extensions.  

4. Download and install Open MPI 1.6.x for Windows (Native, *not* Cygwin) [4].
When prompted, allow the installer to add OpenMPI to the system PATH.

Build Procedure
---------------

This assumes you have checked out the ParaFEM Subversion trunk to C:\parafem

1. Open a Visual Studio command prompt thus:
	Start -> Intel Parallel Studio XE 2013 -> IA-32 Visual Studio 2010 mode
2. Type:
	cd C:\parafem\parafem
	make-parafem-and-drivers-win32-intel.bat

Set the variable BUILD_GROUP_5ED to include the drivers you want built
	
If the build proceeds correctly, the driver (example) file executables will be 
placed into C:\parafem\parafem\bin

References
----------

[1] http://www.microsoft.com/en-gb/download/details.aspx?id=1366
[2] http://www.visualstudio.com/downloads/download-visual-studio-vs#DownloadFamilies_4
[3]	http://software.intel.com/en-us/intel-fortran-studio-xe-evaluation-options
[4] http://www.open-mpi.org/software/ompi/v1.6/
