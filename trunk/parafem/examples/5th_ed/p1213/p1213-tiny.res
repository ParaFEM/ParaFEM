
BASIC JOB DATA                                  
Number of processors used                              1
Number of nodes in the mesh                          756
Number of nodes that were restrained                 396
Number of equations solved                          1640
Number of PCG iterations                              59
Number of loaded nodes                                 8
Total load applied                           -0.1000E+03

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.135839    0.58
Read element steering array                     0.028122    0.12
Convert Abaqus to S&G node ordering             0.000015    0.00
Read nodal coordinates                          0.008496    0.04
Read restrained nodes                           0.060420    0.26
Compute steering array and neq                  0.028118    0.12
Compute interprocessor communication tables     0.000518    0.00
Allocate neq_pp arrays                          0.000067    0.00
Compute element stiffness matrices              0.038425    0.16
Build the preconditioner                        0.000166    0.00
Get starting r                                  0.036657    0.16
Solve equations                                 0.036750    0.16
Output results                                 23.134413   98.41
Total execution time                           23.508006  100.00

