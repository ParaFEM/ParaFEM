
BASIC JOB DATA                                  
Number of processors used                              2
Number of nodes in the mesh                           27
Number of equations solved                            27
Stopping criterion for PCG iterations         0.1000E-06

BOUNDARY CONDITION DATA                         
Number of nodes that were restrained                   0
Number of loaded freedoms                              2
Total energy delivered (load)                 0.2000E+04
Initial global temperature                    0.0000E+00
Number of fixed nodal values                           2

TRANSIENT STATE DETAILS                         
Section   dt               nstep        npri   npri_chk   time        #iters*   Tot#iters
      1   0.1250E-01         800           4        400   0.1000E+02        1        3012
-----------------------------------------------------------------------------------------
Totals                       800         200          0   0.1000E+02                 3012
-----------------------------------------------------------------------------------------
*Number of iterations to complete last timestep of section

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.124800   27.59
Read element steering array                     0.000000    0.00
Convert Abaqus to S&G node ordering             0.000000    0.00
Read nodal coordinates                          0.000000    0.00
Read restrained nodes                           0.000000    0.00
Compute steering array and neq                  0.000000    0.00
Compute interprocessor communication tables     0.000000    0.00
Allocate neq_pp arrays                          0.000000    0.00
Compute element stiffness matrices              0.000000    0.00
Build the preconditioner                        0.000000    0.00
Get starting r                                  0.000000    0.00
Solve equations                                 0.078000   17.24
Output results                                  0.249602   55.17
Total execution time                            0.452402  100.00

