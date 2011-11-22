
  PROGRAM: pf2ensi, pf2ensi.geo.awk, pf2ensi.var.awk

  pf2ensi, pf2ensi.geo.awk and pf2ensi.var.awk are one SHELL script and two AWK
  scripts that convert ParaFEM inout and output files into the EnSight Gold ASCII
  format, which can then be read into visualization tools such as ParaView.
  The simple shell script must be used to convert all of the files and generate
  the EnSight Gold CASE file.

    Usage: pf2ensi <prefix-filename>
  
           e.g., pf2ensi cylinder
 
 NOTES 

  o Development is not yet complete but this is a working version
  o Partial data such as Bound Nodes, Loaded Nodes and Fixed Nodes is handled but
    the standard ParaView application fails to correctly load such data
    
  AUTHOR

  louise.lever@manchester.ac.uk
