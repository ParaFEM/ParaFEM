
The files in this directory are input files for the program sg12mg.f90

Using the input file <filename>.mg, program sg12mg.f90 will create a full 
input deck for program p121. 

  1) ed4-<size>-<element>.mg
 
     Will create an input deck for the example test problem in the 4th edition
     of Smith and Griffiths "Programming the Finite Element Method."

  2) bsq-<size>-<element>.mg
 
     Will create an input deck for the Boussinesq problem of a point load
     on an infinite elastic half space.

<element> denotes element type. The options are C3D8 (8 node hexahedron) and 
C3D20 (20 node hexahedron). The codes are those used in the commercial software
package Abaqus.

