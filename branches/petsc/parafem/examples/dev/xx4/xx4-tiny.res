
BASIC JOB DATA                                  
Number of processors used                              4
Number of nodes in the mesh                          756
Number of nodes that were restrained                 396
Number of equations solved                          1640
Number of PCG iterations                              79
Number of loaded nodes                                 8
Total load applied                           -0.1000E+03

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.063976   46.89
Read element steering array                     0.003891    2.85
Convert Abaqus to S&G node ordering             0.000009    0.01
Read nodal coordinates                          0.002318    1.70
Read restrained nodes                           0.001787    1.31
Compute steering array and neq                  0.003758    2.75
Compute interprocessor communication tables     0.000344    0.25
Allocate neq_pp arrays                          0.000028    0.02
Compute element stiffness matrices              0.010424    7.64
Build the preconditioner                        0.000064    0.05
Get starting r                                  0.001204    0.88
Solve equations                                 0.019866   14.56
Output results                                  0.028755   21.08
Total execution time                            0.136424  100.00

