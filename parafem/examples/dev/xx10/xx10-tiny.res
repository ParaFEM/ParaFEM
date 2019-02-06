
BASIC JOB DATA                                  
Number of processors used                              1
Number of nodes in the mesh                          756
Number of nodes that were restrained                 396
Number of equations solved                          1640
Number of PCG iterations                              79
Number of loaded nodes                                 8
Total load applied                           -0.1000E+03

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.142979   10.07
Read element steering array                     0.001000    0.07
Convert Abaqus to S&G node ordering             0.000000    0.00
Read nodal coordinates                          0.002000    0.14
Read restrained nodes                           0.000999    0.07
Compute steering array and neq                  0.004999    0.35
Compute interprocessor communication tables     0.000000    0.00
Allocate neq_pp arrays                          0.000000    0.00
Compute element stiffness matrices              0.000000    0.00
Build the preconditioner                        0.001000    0.07
Get starting r                                  0.000000    0.00
Solve equations                                 1.170823   82.46
Output results                                  0.095985    6.76
Total execution time                            1.419785  100.00

