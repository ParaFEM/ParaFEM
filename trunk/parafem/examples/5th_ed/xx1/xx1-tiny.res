
BASIC JOB DATA                                  
Number of processors used                              4
Number of nodes in the mesh                          756
Number of nodes that were restrained                 396
Number of equations solved                          1640
Number of PCG iterations                              80
Number of loaded nodes                                 8
Total load applied                           -0.1000E+03

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.365610   56.44
Read element steering array                     0.042986    6.64
Convert Abaqus to S&G node ordering             0.000006    0.00
Read nodal coordinates                          0.002737    0.42
Read restrained nodes                           0.002363    0.36
Compute steering array and neq                  0.000110    0.02
Compute interprocessor communication tables     0.000344    0.05
Allocate neq_pp arrays                          0.000029    0.00
Compute element stiffness matrices              0.010424    1.61
Build the preconditioner                        0.000077    0.01
Get starting r                                  0.043404    6.70
Solve equations                                 0.020720    3.20
Output results                                  0.158968   24.54
Total execution time                            0.647778  100.00

