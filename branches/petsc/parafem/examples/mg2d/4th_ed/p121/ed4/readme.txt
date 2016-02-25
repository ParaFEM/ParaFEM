EXAMPLES: ed4-<size>-<element>.mg

The files in this directory are input files for the program mg2d.f90

Using the input file <filename>.mg, program mg2d.f90 will create a full 
input deck for program p121. 

     ed4-<size>-<element>.mg
 
     Will create an input deck for the problem described in Chapter 12 of
     Smith and Griffiths "Programming the Finite Element Method".
     
     The commands:
    
       chmod +x run-all-ed4.sh 
       source ./run-all-ed4.sh
     
     Will create all the input decks available in this directory.

<element> denotes element type. The options are C3D8 (8 node hexahedron) and 
C3D20 (20 node hexahedron). The element codes (C3D8 and C3D20) are those used
in the commercial software package Abaqus.


SIZE "GPU"

The files ed4-gpu1-c3d20.mg and so on have been generated to try on GPUs,
maximizing usage of different on board RAM memory sizes.

Memory is estimated by multiplying elements (45x45x45 = 91,125) by
the storage required per element (60x60 + 2x60 = 3720) by the number of bytes
required for each double precision real number (8).

This gives 91,125 x 3720 x 8 = 2.7GB for ed4-gpu4-c3d20.mg


AUTHOR

lee.margetts@manchester.ac.uk
