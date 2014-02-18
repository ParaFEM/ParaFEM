
  P126 Demonstrator File Set
  
  This set of files is for use with the demonstrator application. Program p126 
  reads the "P126 Input Files" and outputs the "P126 Output Files". 
  The "P126 ParaView Files" are required to visualize the model using ParaView.
  
  To load the model into ParaView with pre-assigned settings, open the 
  p126_tiny.pvsm file. Note that the "pvsm" file uses an absolute path to the
  files, so a workaround is required.
  
  To load the model into ParaView with default settings, open the 
  p126_tiny.ensi.case file.
  

  P126 Input Files
  
  p126_tiny.dat                 ! control data
  p126_tiny.d                   ! geometry
  p126_tiny.bnd                 ! boundary conditions
  p126_tiny.lid                 ! velocity of lid
  

  P126 Output Files 
  
  p126_tiny.res                 ! summary output
  p126_tiny.ensi.VEL            ! velocities
  
  
  P126 Paraview Files
  
  p126_tiny.pvsm                ! Paraview settings
  p126_tiny.ensi.case           ! Paraview case file
  p126_tiny.ensi.geo            ! geometry
  p126_tiny.ensi.MATID          ! material numbers
  p126_tiny.ensi.NDBND          ! restrained nodes
  p126_tiny.ensi.NDLDS          ! applied forces
  
