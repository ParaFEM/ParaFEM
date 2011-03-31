
  PROGRAM: inp2d, inp2d.awk

  inp2d and inp2d.awk are a simple SHELL script and AWK script that converts simple
  Abaqus Input Decks into the required ParaFEM files ready for job submission. The
  simple shell script can be used, or the AWK script can be called via AWK itself.

    Usage: inp2d <filename.inp>
  
           e.g., inp2d cylinder.inp
 
    Alternative usage: awk -f inp2d.awk <filename.inp>
  
           e.g., awk -f inp2d.awk cylinder.inp

 NOTES 

  o The program has not yet been finished
  o By default, the element node ordering follows the Abaqus convention
  o The script currently handles C3D4, C3D8 and C3D8I element types
  o It processes 'Initial' step Displacement/Rotation Boundary Conditions (BC)
    o Though rotations are not currently supported by ParaFEM
  o It processes BCs where the DOF are given as a range
    o e.g., *Boundary
            ZMIN, 1, 1
  o And it processes BCs where the keywords ENCASTRE and PINNED are used
  o For BCs in subsequent Steps, fixed displacements are processed, using a DOF range with values
    o e.g., *Boundary
            ZMAX, 3, 3, -0.0452

  AUTHOR

  louise.lever@manchester.ac.uk
