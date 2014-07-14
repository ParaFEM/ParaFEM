
BASIC JOB DATA                                  
Number of processors used                              1
Number of nodes in the mesh                          125
Number of nodes that were restrained                   0
Number of equations solved                           125
Number of PCG iterations                              10
Number of fixed displacements                          1

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.030000    0.01
Read element steering array                     0.020000    0.01
Convert Abaqus to S&G node ordering             0.000000    0.00
Read nodal coordinates                          0.000000    0.00
Read restrained nodes                           0.000000    0.00
Compute steering array and neq                  0.020000    0.01
Compute interprocessor communication tables     0.000000    0.00
Allocate neq_pp arrays                          0.000000    0.00
Compute element stiffness matrices              0.010000    0.00
Build the preconditioner                        0.010000    0.00
Get starting r                                  0.040000    0.02
Solve equations                               237.930000   98.75
Output results                                  2.870000    1.19
Total execution time                          240.930000  100.00

