
  PROGRAM: xx6.f90

  Program xx6.f90 computes the eigenvalues and eigenvectors of a 3d elastic
  solid. The program is a modified version of p128 published in Smith I.M. and 
  Griffiths D.V., "Programming the Finite Element Method", 4th Edition, Wiley, 
  2004.
    
    Usage: xx6 <job_name>
  
  
  NOTES
  
  xx6.f90 requires ARPACK. This distribution includes a copy of ARPACK. 
  The original copy of ARPACK was called:

    arpack96.tar.gz

  Two files have been altered, to allow ARPACK to be built with ParaFEM on a 
  CRAY XE6:

    /ARPACK/Makefile and
    /ARPACK/ARmake.inc
  
  AUTHOR: lee.margetts@manchester.ac.uk

