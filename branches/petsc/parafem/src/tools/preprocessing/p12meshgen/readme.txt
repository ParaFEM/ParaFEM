
  PROGRAM: p12meshgen.f90

  p12meshgen is a simple program that generates the built-in meshes for 
  the 10 programs in Smith IM, Griffiths DV, and Margetts L "Programming
  the Finite Element Method", 5th Edition, Wiley, 2013.

    Usage: p12meshgen <job_name>
  
           p12meshgen p121_5

  NOTES 

  o The program has not yet been finished
  o By default, the element node ordering follows the Abaqus convention
  o The program will accept a job name with or without the .mg file 
    extension, however the file name itself must still include the 
    extension.

  AUTHOR

  lee.margetts@manchester.ac.uk
  louise.lever@manchester.ac.uk
