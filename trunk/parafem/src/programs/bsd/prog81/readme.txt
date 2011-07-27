 
 PROGRAM: prog81.f90

 prog81.f90 is a general purpose Boundary Element Method program for solving
 elasticity and potential problems. It is a modified version of a program with
 the same name from Beer, Smith and Duenser, "The Boundary Element Method with
 Programming", Springer, 2008.


   Usage: prog81   (Requires an input deck called prog81.dat)


 The program is a serial program and does not use message passing.

 Known issues: The program crashes with a segmentation fault when compiled 
               using the PGI compiler on HECToR. It seems to work fine when 
               compiled using the Pathscale compiler on the same machine.

  REFERENCES
    
  If you use these programs or modify them for your research, please cite the 
  following paper:
 
    Smith, I.M, Margetts L. Beer G. and Duesner C. (2007) "Parallelising the 
      Boundary Element Method Using ParaFEM", Numerical Models in Geomechanics -
      NUMOG X, pp177-181.
    
    Smith, I.M. and Margetts L. (2007) "Parallel Boundary Element Analysis of 
      Tunnels, Proceedings EUROTUN 2007, Vienna, Austria.
      
 AUTHORS

 ian.smith@manchester.ac.uk
 lee.margetts@manchester.ac.uk
