
  PROGRAM: prog83.f90

  prog83 is a general purpose BEM program for solving elasticity problems. 
  This version is parallel and uses BiCGStab(l). It must be linked to an MPI
  library on compilation. The program is a modified modified version of the 
  one published in Beer G., Smith I.M. and Duenser C., "The Boundary Element 
  Method with Programming", Springer-Verlag, 2008.

    Usage: [mpirun -np 2] prog83

  The statement in square brackets depends on the MPI implementation.
  
  
  REFERENCES
    
  If you use these programs or modify them for your research, please cite the 
  following papers:

    Smith, I.M, Margetts L. Beer G. and Duesner C. (2007) "Parallelising the 
      Boundary Element Method Using ParaFEM", Numerical Models in Geomechanics -
      NUMOG X, pp177-181.
  
    Smith, I.M. and Margetts L. (2007) "Parallel Boundary Element Analysis of 
      Tunnels, Proceedings EUROTUN 2007, Vienna, Austria.
    

  MODIFICATIONS

  There are several differences between this version and the version in the 
  book. 

    o utility.f90           -> split over several modules (maths.f90)
    o global_variables1.f90 -> global_variables.f90
    o gather_scatter6.f90   -> gather_scatter.f90
    o nels, ielpe, timest(:) no longer declared as global_variables 
    o subroutine calc_nels_pp requires argument "nels"    
     
  AUTHOR

  ian.smith@manchester.ac.uk
  lee.margetts@manchester.ac.uk

