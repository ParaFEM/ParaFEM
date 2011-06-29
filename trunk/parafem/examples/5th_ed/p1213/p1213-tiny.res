
BASIC JOB DATA                                  
Number of processors used                              1
Number of nodes in the mesh                          756
Number of nodes that were restrained                 396
Number of equations solved                          1640
Number of PCG iterations                              91
Number of loaded nodes                                 8
Total load applied                           -0.1000E+03

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.107576   21.25
Read element steering array                     0.047052    9.30
Convert Abaqus to S&G node ordering             0.000013    0.00
Read nodal coordinates                          0.023442    4.63
Read restrained nodes                           0.008142    1.61
Compute steering array and neq                  0.024305    4.80
Compute interprocessor communication tables     0.000521    0.10
Allocate neq_pp arrays                          0.000068    0.01
Compute element stiffness matrices              0.038424    7.59
Build the preconditioner                        0.000151    0.03
Get starting r                                  0.016520    3.26
Solve equations                                 0.057134   11.29
Output results                                  0.182780   36.11
Total execution time                            0.506128  100.00

