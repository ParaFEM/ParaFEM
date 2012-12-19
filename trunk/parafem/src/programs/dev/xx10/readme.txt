
  PROGRAM: xx10.f90

  xx10 performs the 3D analysis of a linear elastic solid. The program is a 
  modified version of p121.f90 published in Smith I.M. and Griffiths D.V., 
  "Programming the Finite Element Method", 4th Edition, Wiley, 2004.

  This program supports GPU acceleration using OpenCL.

  Uses a special case where only one element stiffness matrix is formed
  leading to a matrix-matrix kernal in PCG which replaces a loop of element
  based matrix-vector multiplications. 
  
  *** Warning: Under development ***
  
    Usage: xx10 <job_name>

   
  AUTHOR(S)

  lee.margetts@manchester.ac.uk
  george.leaver@manchester.ac.uk
