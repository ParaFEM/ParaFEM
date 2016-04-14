
  PROGRAM: cylinder.f90

  Cylinder is a simple program that generates cylindrical meshes for 
  use with CGPACK. The version here has been adapted to output in binary 
  Ensight Gold format.

   Usage: 
 
     ./cylinder  radius  length  error  label  number_elements

   Example: 

    ./cylinder 1 1 0.01 Newfile 2
    
  NOTES 

  o The program has not yet been finished
  o By default, the element node ordering follows the Abaqus convention


  AUTHORS

  Anton Shterenlikht, University of Bristol 
  J. David Arregui, University of Manchester
  Lee Margetts, University of Manchester
