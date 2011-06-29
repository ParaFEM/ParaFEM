
BASIC JOB DATA                                  
Number of processors used                              1
Number of nodes in the mesh                          756
Number of nodes that were restrained                 396
Number of equations solved                          1640
Number of PCG iterations                              80
Number of loaded nodes                                 8
Total load applied                           -0.1000E+03

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.058049   27.96
Read element steering array                     0.003822    1.84
Convert Abaqus to S&G node ordering             0.000013    0.01
Read nodal coordinates                          0.002375    1.14
Read restrained nodes                           0.001588    0.76
Compute steering array and neq                  0.014662    7.06
Compute interprocessor communication tables     0.000488    0.24
Allocate neq_pp arrays                          0.000066    0.03
Compute element stiffness matrices              0.038531   18.56
Build the preconditioner                        0.000164    0.08
Get starting r                                  0.004800    2.31
Solve equations                                 0.050369   24.26
Output results                                  0.032696   15.75
Total execution time                            0.207623  100.00

