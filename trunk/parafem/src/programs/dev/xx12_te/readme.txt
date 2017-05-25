
  PROGRAM: xx12_te.f90

  xx12_te performs the 3D analysis of a linear elastic solid due to thermal loading.
  The program is a modified version of rfemsolve_te, which in turn is a modified
  version of p121.f90 published in Smith I.M. and Griffiths D.V., "Programming the 
  Finite Element Method", 4th Edition, Wiley, 2004.

  This program supports multiple material types.
  
  *** Warning: Under development ***
  
    Usage: xx12_te <job_name_in> <job_name_out>

  AUTHORS

  lee.margetts@manchester.ac.uk
  jose.arreguimena@postgrad.manchester.ac.uk
  llion.evans@ccfe.ac.uk

  -------
  BELOW INSTRUCTIONS ARE FOR RFEMSOLVE AND MAY NOT HOLD FOR MODIFIED VERSION
  
  WINDOWS BUILD INSTRUCTIONS
  
  1. Open the command prompt
  2. Navigate to: D:\parafem\parafem
  3. SVN update
  4. Type: make-parafem-win32.bat
  5. Add the executables to the PATH: set PATH=%PATH%;D:\parafem\parafem\bin
  
  CREATE INPUT DECK
  
  1. Use program rfemcube.exe with input *.rf from /examples/rfem
  2. Use program rfembc.exe with input created from rfemcube.exe
  3. Use program rfemfield.exe with the output created by (1) and (2)
  4. Repeat (3) for as many instances as required
  
  SOLVE THE EQUATIONS
  
  1. Use Windows prompt and use rfemsolve.exe
  
  POSTPROCESSING
  
  1. Open cygwin navigate to the working directory
  2. Use the instruction: export PATH=$PATH:/cygdrive/d/parafem/parafem/bin
  3. Copy the postprocessing tool pf2ensi to /bin -> cp src/tools/postprocessing/pf2ensi/pf2ensi* bin
  4. Transform the output data of rfemsolve with the instruction -> < pf2ense> <output-rfemsolve>
  5. Use ParaView and open the file with extension .case
    
  ISSUES
  
  1. The argument <mises-threshold> needs to be removed from rfemcube.exe
  