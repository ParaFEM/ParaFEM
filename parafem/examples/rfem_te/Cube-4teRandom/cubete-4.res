
BASIC JOB DATA                                  
Number of processors used                              1
Number of nodes in the mesh                          125
Number of nodes that were restrained                   0
Number of equations solved                           375
Number of PCG iterations                              37
Threshold values for von Mises stress           0.000000
Number of nodes greater than threshold                 0

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.000000    0.00
Read element steering array                     0.000000    0.00
Convert Abaqus to S&G node ordering             0.000000    0.00
Read nodal coordinates                          0.000000    0.00
Read restrained nodes                           0.000000    0.00
Compute steering array and neq                  0.000000    0.00
Compute interprocessor communication tables     0.000000    0.00
Allocate neq_pp arrays                          0.000000    0.00
Compute element stiffness matrices              0.031200   33.33
Build the preconditioner                        0.000000    0.00
Get starting r                                  0.000000    0.00
Solve equations                                 0.046800   50.00
Output results                                  0.015600   16.67
Total execution time                            0.093601  100.00

