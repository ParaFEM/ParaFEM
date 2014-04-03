
  P124 Demonstrator File Set
  
  This set of files is for use with the demonstrator application. Program p124 
  reads the "P124 Input Files" and outputs the "P124 Output Files". 
  The "P124 ParaView Files" are required to visualize the model using ParaView.
  
  To load the model into ParaView with pre-assigned settings, open the 
  p124_demo.pvsm file. Note that the "pvsm" file uses an absolute path to the
  files, so a workaround is required.
  
  To load the model into ParaView with default settings, open the 
  p124_demo.ensi.case file.
  

  P124 Input Files
  
  p124_demo.dat                 ! control data
  p124_demo.d                   ! geometry
  p124_demo.bnd                 ! boundary conditions
  

  P124 Output Files 
  
  p124_demo.res                 ! summary output
  p124_demo.ensi.NDTTR-******   ! nodal temperature
  
  
  P124 Paraview Files
  
  p124_demo.pvsm                ! Paraview settings
  p124_demo.ensi.case           ! Paraview case file
  p124_demo.ensi.geo            ! geometry
  p124_demo.ensi.MATID          ! material numbers
  p124_demo.ensi.NDBND          ! restrained nodes
  
