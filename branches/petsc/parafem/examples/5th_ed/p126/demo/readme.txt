
  P126 Demonstrator File Set
  
  This set of files is for use with the demonstrator application. Program p126 
  reads the "P126 Input Files" and outputs the "P126 Output Files". 
  The "P126 ParaView Files" are required to visualize the model using ParaView.
  
  To load the model into ParaView with pre-assigned settings, open the 
  p126_demo.pvsm file. Note that the "pvsm" file uses an absolute path to the
  files, so a workaround is required.
  
  To load the model into ParaView with default settings, open the 
  p126_demo.ensi.case file.
  

  P126 Input Files
  
  p126_demo.dat                 ! control data
  p126_demo.d                   ! geometry
  p126_demo.bnd                 ! boundary conditions
  p126_demo.lid                 ! velocity of lid
  

  P126 Output Files 
  
  p126_demo.res                 ! summary output
  p126_demo.ensi.VEL            ! velocities
  
  
  P126 Paraview Files
  
  p126_demo.pvsm                ! Paraview settings
  p126_demo.ensi.case           ! Paraview case file
  p126_demo.ensi.geo            ! geometry
  p126_demo.ensi.MATID          ! material numbers
  p126_demo.ensi.NDBND          ! restrained nodes
  p126_demo.ensi.NDLDS          ! applied forces
  
