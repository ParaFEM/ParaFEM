PROGRAM: gaitfem.f90

This program solves a 3D small strain elasticity problem, in serial, using a 
direct solver. 

  Usage: gaitfem <job_name>


ABOUT

The purpose of the program is to contain and test a series of subroutines that
will be embedded within the GaitSym program developed by Dr Bill Sellers at 
the University of Manchester.

GaitSym runs multiple realisations of a genetic algorithm across systems with
many thousands of processors. The next stage of development is to embed a 
simple FEM capability that will run independently on each processor for each
realisation.

Each input deck for GaitSym will have a single FE model that will be 
subjected to loads that change in magnitude and direction over the walking 
cycle. This means that a two stage solver with separate LU factorization and 
solution phases is an appropriate strategy.

  
AUTHOR

lee.margetts@manchester.ac.uk

