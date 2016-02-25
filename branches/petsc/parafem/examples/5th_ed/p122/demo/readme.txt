
  P122 Demonstrator File Set
  
  This set of files is for use with the demonstrator application. Program p122 
  reads the "P122 Input Files" and outputs the "P122 Output Files". 
  The "P122 ParaView Files" are required to visualize the model using ParaView.
  
  To load the model into ParaView with pre-assigned settings, open the 
  p122_demo.pvsm file. Note that the "pvsm" file uses an absolute path to the
  files, so a workaround is required.
  
  To load the model into ParaView with default settings, open the 
  p122_demo.ensi.case file.
  

  P122 Input Files
  
  p122_demo.dat                 ! control data
  p122_demo.d                   ! geometry
  p122_demo.bnd                 ! boundary conditions
  p122_demo.fix                 ! fixed displacements
  

  P122 Output Files 
  
  p122_demo.res                 ! summary output
  p122_demo.ensi.DISPL-000001   ! displacements
  
  
  P122 Paraview Files
  
  p122_demo.pvsm                ! Paraview settings
  p122_demo.ensi.case           ! Paraview case file
  p122_demo.ensi.geo            ! geometry
  p122_demo.ensi.MATID          ! material numbers
  p122_demo.ensi.NDBND          ! restrained nodes
  p122_demo.ensi.NDFIX          ! fixed displacements
  
