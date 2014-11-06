
  PROGRAM: xx1.f90

  xx1 performs the 3D analysis of a linear elastic solid and adds support for 
  UMAT integration. The program is a modified version of p121.f90 published in 
  Smith I.M., Griffiths D.V. and Margetts L., "Programming the Finite Element 
  Method", 5th Edition, Wiley, 2014.

    Usage: xx1 <job_name>

  Note that this program is under development and therefore subject to 
  frequent modifications and updates.
  
  INPUT DATA
  
  The program is currently being modified to read and write the model data in 
  the C Binary Ensight Gold format. The expected files are:
  

  1. <job_name>.bin.ensi.case          Ensight Gold format case file 
                                    
     An ASCII file that contains basic details about the job and a list of 
     associated files.
                                                
  3. <job_name>.bin.ensi.geo           Ensight Gold format geometry data
  
     A C binary file that contains the coordinates of the nodes G_COORD and 
     the element steering array G_NUM. 

  4. <job_name>.bin.ensi.MATID         Ensight Gold format material type      
  
     A C binary file that lists the element property type ETYPE for each 
     element. 
                           
  5. <job_name>.bin.ensi.NDBND         Ensight Gold format restraints    
  
     A C binary file that lists the restraints for all nodes. The convention 
     is given on page 20 of the 5th edition of the text book.
                                                
  6. <job_name>.bin.ensi.NDLDS         Ensight Gold format loads         
  
     A C binary file that lists nodal loads for all nodes as a vector quantity
     for each node.
                                                
  7. <job_name>.bin.ensi.NDFIX         Ensight Gold format fixed nodes   
  
     A C binary file that lists fixed values for all nodes, for example fixed 
     displacements.

  8. <job_name>.bin.ensi.dat           ParaFEM format control data
                                   
     An ASCII file that contains the basic control data required by program 
     xx1.f90. This is not an Ensight Gold file. The "bin.ensi" naming
     convention is followed for convenience.

  9. <job_name>.bin.ensi.mat           ParaFEM format materials data
                                   
     An ASCII file that contains the materials data required by program 
     xx1.f90. This is not an Ensight Gold file. The "bin.ensi" naming
     convention is followed for convenience.
     
     The file contains the following variables where NMATS is the number of
     materials and NVALS is the number of values for those materials:
     
     KEYWORD NMATS NVALS
     LABEL1 LABEL2 ... LABEL_N
     MATERIALTYPE_1 VALUE_1 VALUE_2 ... VALUE_N
     MATERIALTYPE_2 VALUE_1 VALUE_2 ... VALUE_N
     .
     .
     .
     MATERIALTYPE_N VALUE_1 VALUE_2 ... VALUE_N
     
     An example follows for a model with two elastic materials:
     
     *MATERIAL  2  2
     e v
     1 104000. 0.3
     2 370000. 0.4
     
  OUTPUT DATA
  
  10. <job_name>.bin.ensi.DISPL-000001  Ensight Gold format displacements
  
      A C binary file that contains the computed displacements at each node.
                                    
  11. <job_name>.bin.ensi.res           ParaFEM format results 
  
      An ASCII filecontaining basic job information as well as a break down of
      time spent in the main program sections.
  
  NOTES

  The build scripts do not currently cope well with mixed FORTRAN77 and
  modern Fortran. The makefile may need customising for different platforms.
  This will be addressed.
   
  AUTHOR(S)

  lee.margetts@manchester.ac.uk
  louise.lever@manchester.ac.uk
