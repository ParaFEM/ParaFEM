
BASIC JOB DATA                                  
Number of processors used                              4
Number of nodes in the mesh                          756
Number of nodes that were restrained                 396
Number of equations solved                          1640
Number of PCG iterations                              79
Number of loaded nodes                                 8
Total load applied                           -0.1000E+03

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.062667   46.39
Read element steering array                     0.003946    2.92
Convert Abaqus to S&G node ordering             0.000009    0.01
Read nodal coordinates                          0.002335    1.73
Read restrained nodes                           0.001563    1.16
Compute steering array and neq                  0.003758    2.78
Compute interprocessor communication tables     0.000330    0.24
Allocate neq_pp arrays                          0.000022    0.02
Compute element stiffness matrices              0.010405    7.70
Build the preconditioner                        0.000065    0.05
Get starting r                                  0.001371    1.01
Solve equations                                 0.019961   14.78
Output results                                  0.028659   21.21
Total execution time                            0.135091  100.00

