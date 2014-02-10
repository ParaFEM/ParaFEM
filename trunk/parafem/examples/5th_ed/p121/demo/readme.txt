
  P121 Demonstrator File Set
  
  This set of files is for use with the demonstrator application. Program p121 
  reads the "P121 Input Files" and outputs the "P121 Output Files". 
  The "P121 ParaView Files" are required to visualize the model using ParaView.
  
  To load the model into ParaView with pre-assigned settings, open the 
  p121_small.pvsm file. Note that the "pvsm" file uses an absolute path to the
  files, so a workaround is required.
  
  To load the model into ParaView with default settings, open the 
  p121_small.ensi.case file.
  

  P121 Input Files
  
  p121_small.dat                 ! control data
  p121_small.d                   ! geometry
  p121_small.bnd                 ! boundary conditions
  p121_small.lds                 ! applied forces
  

  P121 Output Files 
  
  p121_small.res                 ! summary output
  p121_small.ensi.DISPL-000001   ! displacements
  
  
  P121 Paraview Files
  
  p121_small.pvsm                ! Paraview settings
  p121_small.ensi.case           ! Paraview case file
  p121_small.ensi.geo            ! geometry
  p121_small.ensi.MATID          ! material numbers
  p121_small.ensi.NDBND          ! restrained nodes
  p121_small.ensi.NDLDS          ! applied forces
  
