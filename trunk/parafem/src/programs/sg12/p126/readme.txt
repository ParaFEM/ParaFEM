
NOTES ABOUT PROGRAM P161

Program p161 appears to be working properly. There are some issues:

  1. The "built-in" problem from Smith and Griffiths does not appear to
     apply the lid velocity correctly.

  2. The results are output into a file named <job_name>.upvw 
     To view the model in the ParaFEM-Viewer, the .dis file is created using
     the program upvw2dis.f90 in /parafem/src/programs/utilities
     The velocities are written as a first step in .dis and the pressures as 
     a second step. Pressures are only computed on the corner nodes, so do 
     not visualize properly on a mesh with 20-node bricks.


