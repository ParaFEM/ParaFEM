 DESCRIPTION:	The files in this directory, xx14_<size>.mg,
		are input files for the program p12meshgen.f90.

       USAGE:	Using the input file <filename>.mg,
		program p12meshgen.f90 will create a full 
		input deck for program xx14:

		p12meshgen xx14_medium

		by default, the p12meshgen executable can be found
		under parafem/bin:

		../../../../parafem/bin/p12meshgen xx14_medium
 
		Then the executable xx14 can be run
		on the produced input deck.
		Use like this (change options as required):

		mpirun -np 16 xx14 xx14_medium

  REFERENCES:	See the problem described on pages 536-541 in
		Chapter 12 of Smith, Griffiths and Margetts
		"Programming the Finite Element Method",
		5th Edition, Wiley, 2014

LAST UPDATED:	3-OCT-2014

      AUTHOR:	Anton Shterenlikht, mexas@bris.ac.uk
