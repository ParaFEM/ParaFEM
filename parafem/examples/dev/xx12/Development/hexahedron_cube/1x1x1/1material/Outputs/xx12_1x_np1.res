
BASIC JOB DATA                                  
Number of processors used                              1
Number of nodes in the mesh                            8
Number of equations solved                             8
Stopping criterion for PCG iterations         0.1000E-06

BOUNDARY CONDITION DATA                         
Number of nodes that were restrained                   0
Number of loaded freedoms                              2
Total energy delivered (load)                 0.2000E+04
Initial global temperature                    0.0000E+00
Number of fixed nodal values                           2

TRANSIENT STATE DETAILS                         
Section   dt               nstep        npri   npri_chk   time        #iters*   Tot#iters
      1   0.1250E-01         800           4        400   0.1000E+02        1        1135
-----------------------------------------------------------------------------------------
Totals                       800         200          0   0.1000E+02                 1135
-----------------------------------------------------------------------------------------
*Number of iterations to complete last timestep of section

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.124800   30.77
Read element steering array                     0.000000    0.00
Convert Abaqus to S&G node ordering             0.000000    0.00
Read nodal coordinates                          0.000000    0.00
Read restrained nodes                           0.000000    0.00
Compute steering array and neq                  0.000000    0.00
Compute interprocessor communication tables     0.000000    0.00
Allocate neq_pp arrays                          0.000000    0.00
Compute element stiffness matrices              0.000000    0.00
Build the preconditioner                        0.000000    0.00
Get starting r                                  0.015600    3.85
Solve equations                                 0.031200    7.69
Output results                                  0.234002   57.69
Total execution time                            0.405602  100.00

