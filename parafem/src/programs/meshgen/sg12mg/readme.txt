
  PROGRAM: sg12mg.f90

  sg12mg is a simple program that generates the built-in meshes for the 10 
  programs in Smith I.M. and Griffiths D.V., "Programming the Finite Element
  Method", 4th Edition, Wiley, 2004.

    Usage: sg12mg <job_name>
  
           e.g., sg12mg ed4-small-c3d8
           e.g., sg12mg ed4-small-c3d8.mg           
  NOTES 

  o The program has not yet been finished
  o By default, the element node ordering follows the Abaqus convention
  o The program will accept a job name with or without the .mg file extension, however the
    file name itself must still include the extension.

  AUTHOR

  lee.margetts@manchester.ac.uk
  louise.lever@manchester.ac.uk
