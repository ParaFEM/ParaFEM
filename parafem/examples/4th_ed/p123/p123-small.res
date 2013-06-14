
BASIC JOB DATA                                  
Number of processors used                              1
Number of nodes in the mesh                      1030301
Number of nodes that were restrained               30301
Number of equations solved                       1000000
Number of PCG iterations                             121
Number of loaded freedoms                              1
Total load applied                            0.1000E+02

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.044002    0.07
Read element steering array                     2.840178    4.84
Convert Abaqus to S&G node ordering             0.012001    0.02
Read nodal coordinates                          1.456091    2.48
Read restrained nodes                           0.016001    0.03
Compute steering array and neq                  0.636039    1.08
Compute interprocessor communication tables     0.428027    0.73
Allocate neq_pp arrays                          0.028002    0.05
Compute element stiffness matrices             13.212826   22.50
Build the preconditioner                        0.220013    0.37
Get starting r                                  0.000000    0.00
Solve equations                                39.830490   67.83
Output results                                  0.000000    0.00
Total execution time                           58.723670  100.00

