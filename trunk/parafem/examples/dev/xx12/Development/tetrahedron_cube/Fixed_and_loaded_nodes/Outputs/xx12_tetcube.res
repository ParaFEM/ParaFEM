
BASIC JOB DATA                                  
Number of processors used                              2
Number of nodes in the mesh                          155
Number of equations solved                           155
Stopping criterion for PCG iterations         0.1000E-04

BOUNDARY CONDITION DATA                         
Number of nodes that were restrained                   0
Number of loaded freedoms                             25
Total energy delivered (load)                 0.2500E+04
Initial global temperature                    0.0000E+00
Number of fixed nodal values                          25

TRANSIENT STATE DETAILS                         
Section   dt               nstep        npri   npri_chk   time        #iters*   Tot#iters
      1   0.1000E+00          10           1          1   0.1000E+01       21         208
-----------------------------------------------------------------------------------------
Totals                        10          10         10   0.1000E+01                  208
-----------------------------------------------------------------------------------------
*Number of iterations to complete last timestep of section

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.124800   61.54
Read element steering array                     0.000000    0.00
Convert Abaqus to S&G node ordering             0.000000    0.00
Read nodal coordinates                          0.000000    0.00
Read restrained nodes                           0.000000    0.00
Compute steering array and neq                  0.000000    0.00
Compute interprocessor communication tables     0.000000    0.00
Allocate neq_pp arrays                          0.000000    0.00
Compute element stiffness matrices              0.015600    7.69
Build the preconditioner                        0.000000    0.00
Get starting r                                  0.000000    0.00
Solve equations                                 0.015600    7.69
Output results                                  0.046800   23.08
Total execution time                            0.202800  100.00

