
  PROGRAM: rfemcube.f90

  rfemcube is a simple program that generates a simple box mesh based on the RFEM input .rf file
  It allows the cube to be generated independently of the cube in the rfemfield code. It will
  output to a given model job name and will include the required Vom Mises Stress Threshold value
  in the model's .dat file
  
    Usage: rfemcube <in_rfem_name> <out_model_name> <mises-threshold>
  
  NOTES 

  o The program has not yet been finished

  AUTHOR

  louise.lever@manchester.ac.uk
