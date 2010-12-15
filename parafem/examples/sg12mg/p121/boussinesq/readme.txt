
EXAMPLES: bsq-<size>-<element>.mg

The files in this directory are input files for the program sg12mg.f90

Using the input file <filename>.mg, program sg12mg.f90 will create a full 
input deck for program p121. 

     bsq-<size>-<element>.mg
 
     Will create an input deck for the Boussinesq problem of a point load
     applied to an infinite elastic half-space.
     
     The command:
     
       $source ./run-all-bsq.sh
     
     Will create all the input decks available in this directory.     

<element> denotes element type. The options are C3D8 (8 node hexahedron) and 
C3D20 (20 node hexahedron). The codes are those used in the commercial software
package Abaqus.


AUTHOR

Lee Margetts

