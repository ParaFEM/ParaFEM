
  *** READ ONLY PROGRAM ***
  
  PROGRAM: p124.f90

  p124 performs the 3D analysis of the transient heat conduction equation. 
  This program is identical to the one published in Smith IM, Griffiths DV and 
  Margetts L "Programming the Finite Element Method", 5th Edition, Wiley, 2014.
  Please do not commit any changes to this program.

    Usage: p124 <job_name>

  This was originally program p124 in the 4th Edition of the book and is 
  similar to program xx12 which is still under development by Llion Evans in 
  folder /dev

  ERRATUM

  On some platforms, the original program fails because loaded_freedoms_pp
  has not been set to zero. The program has been updated on line 132:

    IF(loaded_freedoms==0) loaded_freedoms_pp=0

  AUTHOR

  lee.margetts@manchester.ac.uk

