
BASIC JOB DATA                                  
Number of processors used                              1
Number of nodes in the mesh                          756
Number of nodes that were restrained                 396
Number of equations solved                          1640
Number of PCG iterations                              80
Number of loaded nodes                                 8
Total load applied                           -0.1000E+03

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.058842   27.85
Read element steering array                     0.008638    4.09
Convert Abaqus to S&G node ordering             0.000015    0.01
Read nodal coordinates                          0.002376    1.12
Read restrained nodes                           0.001881    0.89
Compute steering array and neq                  0.013867    6.56
Compute interprocessor communication tables     0.000492    0.23
Allocate neq_pp arrays                          0.000063    0.03
Compute element stiffness matrices              0.039452   18.68
Build the preconditioner                        0.000155    0.07
Get starting r                                  0.001825    0.86
Solve equations                                 0.049988   23.66
Output results                                  0.033650   15.93
Total execution time                            0.211244  100.00

