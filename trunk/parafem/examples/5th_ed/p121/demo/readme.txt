
  P121 Demonstrator File Set
  
  This set of files is for use with the demonstrator application. Program p121 
  reads the "P121 Input Files" and outputs the "P121 Output Files". 
  The "P121 ParaView Files" are required to visualize the model using ParaView.
  
  To load the model into ParaView with pre-assigned settings, open the 
  p121_demo.pvsm file. Note that the "pvsm" file uses an absolute path to the
  files, so a workaround is required.
  
  To load the model into ParaView with default settings, open the 
  p121_demo.ensi.case file.
  

  P121 Input Files
  
  p121_demo.dat                 ! control data
  p121_demo.d                   ! geometry
  p121_demo.bnd                 ! boundary conditions
  p121_demo.lds                 ! applied forces
  

  P121 Output Files 
  
  p121_demo.res                 ! summary output
  p121_demo.ensi.DISPL-000001   ! displacements
  
  
  P121 Paraview Files
  
  p121_demo.pvsm                ! Paraview settings
  p121_demo.ensi.case           ! Paraview case file
  p121_demo.ensi.geo            ! geometry
  p121_demo.ensi.MATID          ! material numbers
  p121_demo.ensi.NDBND          ! restrained nodes
  p121_demo.ensi.NDLDS          ! applied forces
  
