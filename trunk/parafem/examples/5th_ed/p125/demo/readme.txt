
  P125 Demonstrator File Set
  
  This set of files is for use with the demonstrator application. Program p125 
  reads the "P125 Input Files" and outputs the "P125 Output Files". 
  The "P125 ParaView Files" are required to visualize the model using ParaView.
  
  To load the model into ParaView with pre-assigned settings, open the 
  p125_demo.pvsm file. Note that the "pvsm" file uses an absolute path to the
  files, so a workaround is required.
  
  To load the model into ParaView with default settings, open the 
  p125_demo.ensi.case file.
  

  P125 Input Files
  
  p125_demo.dat                 ! control data
  p125_demo.d                   ! geometry
  p125_demo.bnd                 ! boundary conditions
  

  P125 Output Files 
  
  p125_demo.res                 ! summary output
  p125_demo.ensi.NDPRE-******   ! nodal pressure
  
  
  P125 Paraview Files
  
  p125_demo.pvsm                ! Paraview settings
  p125_demo.ensi.case           ! Paraview case file
  p125_demo.ensi.geo            ! geometry
  p125_demo.ensi.MATID          ! material numbers
  p125_demo.ensi.NDBND          ! restrained nodes
  
