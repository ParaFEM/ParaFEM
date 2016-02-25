
  *** READ ONLY PROGRAM ***
  
  PROGRAM: p123.f90

  p123 performs the three dimensional analysis of Laplace's equation using
  8-node bricks. This program is identical to the one published in Smith IM, 
  Griffiths DV and Margetts L "Programming the Finite Element Method", 5th 
  Edition, Wiley, 2014. Please do not commit any changes to this program.

    Usage: p123 <job_name>

  This was originally program p123 in the 4th Edition of the book and is 
  similar to program xx11 which is still under development by Llion Evans in 
  folder /dev

  ERRATUM

  In the book, the argument NUMVAR in the subroutine SCATTER_NODES (line 176)
  is given the value NDIM(=3). This is incorrect and has been replaced with 
  NODOF(=1).
  
  AUTHOR

  lee.margetts@manchester.ac.uk

