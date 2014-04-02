
  *** READ ONLY PROGRAM ***
  
  PROGRAM: p125.f90

  p125 performs the 3D explicit analysis of the transient conduction equation. 
  This program is identical to the one published in Smith IM, Griffiths DV and 
  Margetts L "Programming the Finite Element Method", 5th Edition, Wiley, 2014. 
  Please do not commit any changes to this program.

    Usage: p125 <job_name>
  

  ERRATUM

  In the book, the argument NUMVAR in the subroutine SCATTER_NODES (line 110)
  is given the value NDIM(=3). This is incorrect and has been replaced with
  NODOF(=1).

  AUTHOR

  lee.margetts@manchester.ac.uk

