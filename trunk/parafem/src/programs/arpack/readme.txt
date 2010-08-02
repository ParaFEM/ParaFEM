
  PROGRAM: p128ar.f90

  Program p128ar.f90 computes the eigenvalues and eigenvectors of a 3d elastic
  solid. The program is a modified version of p128 published in Smith I.M. and 
  Griffiths D.V., "Programming the Finite Element Method", 4th Edition, Wiley, 
  2004.
    
    Usage: p128ar <job_name>
  
  
  NOTES
  
  p128ar.f90 requires ARPACK. This distribution includes a copy of ARPACK. 
  The original copy of ARPACK was called:

    arpack96.tar.gz

  Two files have been altered, to allow ARPACK to be built with ParaFEM on a 
  CRAY XT4:

    /ARPACK/Makefile and
    /ARPACK/ARmake.inc

  To unzip, untar and build type:

    $ gunzip arpack96-pf.tar.gz
    $ tar xvf arpack96-pf.tar
    $ cd ARPACK

  In ARmake.inc, modify the $(HOME) directory to reflect where you have 
  installed ParaFEM. Then type:

    $ make lib

  The library will be placed in the /parafem/lib directory.

  
  AUTHOR: lee.margetts@manchester.ac.uk

