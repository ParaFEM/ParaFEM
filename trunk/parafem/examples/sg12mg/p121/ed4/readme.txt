EXAMPLES: ed4-<size>-<element>.mg

The files in this directory are input files for the program sg12mg.f90

Using the input file <filename>.mg, program sg12mg.f90 will create a full 
input deck for program p121. 

     ed4-<size>-<element>.mg
 
     Will create an input deck for the problem described in Chapter 12 of
     Smith and Griffiths "Programming the Finite Element Method".
     
     The commands:
    
       $chmod +x run-all-ed4.sh 
       $source ./run-all-ed4.sh
     
     Will create all the input decks available in this directory.

<element> denotes element type. The options are C3D8 (8 node hexahedron) and 
C3D20 (20 node hexahedron). The codes are those used in the commercial software
package Abaqus.


AUTHOR

lee.margetts@manchester.ac.uk
