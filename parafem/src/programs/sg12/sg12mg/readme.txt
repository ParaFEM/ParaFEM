
  PROGRAM: sg12mg.f90

  sg12mg is a simple program that generates the built-in meshes for the 10 
  programs in Smith I.M. and Griffiths D.V., "Programming the Finite Element
  Method", 4th Edition, Wiley, 2004.

    Usage: sg12mg <job_name>
  
           eg. sg12mg ed4-small-c3d8
           
  NOTES 

  o The program has not yet been finished
  o By default, the element node ordering follows the Abaqus convention


  POSSIBLE ERROR MESSAGES
  
  o The following error message is given with the PGF90 compiler if you include 
    the file extension in the <job_name>.
    
    
    " PGFIO-F-209/OPEN/unit=10/'OLD' specified for file which does not exist.
      File name = ed4-small-c3d8.mg.mg                              
      In source file sg12mg.f90, at line number 72 "
  
  
  AUTHOR

  lee.margetts@manchester.ac.uk

