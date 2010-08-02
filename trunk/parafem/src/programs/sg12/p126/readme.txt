
  PROGRAM: p126.f90

  p126 solves the steady state Navier Stokes equations for an incompressible
  viscous fluid using BiCGSTAB(l). The program is a modified version of the
  one published in Smith I.M. and Griffiths D.V., "Programming the Finite 
  Element Method", 4th Edition, Wiley, 2004.

    Usage: p126 <job_name>

   
  NOTES

  o The results are output into a file named <job_name>.upvw 
    To view the model in the ParaFEM-Viewer, the .dis file is created using
    the program upvw2dis.f90 in /parafem/src/programs/utilities
    The velocities are written as a first step in .dis and the pressures as 
    a second step. Pressures are only computed on the corner nodes, so do 
    not visualize properly on a mesh with 20-node bricks.


  AUTHOR

  lee.margetts@manchester.ac.uk

