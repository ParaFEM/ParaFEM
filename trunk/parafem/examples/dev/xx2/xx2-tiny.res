
BASIC JOB DATA                                  
Number of processors used                              2
Number of nodes in the mesh                          756
Number of nodes that were restrained                 396
Number of equations solved                          1640
Number of PCG iterations                              59
Number of loaded nodes                                 8
Total load applied                           -0.1000E+03

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.097719   21.57
Read element steering array                     0.078987   17.43
Convert Abaqus to S&G node ordering             0.000011    0.00
Read nodal coordinates                          0.002360    0.52
Read restrained nodes                           0.055670   12.29
Compute steering array and neq                  0.051566   11.38
Compute interprocessor communication tables     0.000369    0.08
Allocate neq_pp arrays                          0.000036    0.01
Compute element stiffness matrices              0.020336    4.49
Build the preconditioner                        0.000090    0.02
Get starting r                                  0.057675   12.73
Solve equations                                 0.022268    4.91
Output results                                  0.066032   14.57
Total execution time                            0.453119  100.00

