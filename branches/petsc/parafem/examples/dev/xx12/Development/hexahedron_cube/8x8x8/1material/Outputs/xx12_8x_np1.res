
BASIC JOB DATA                                  
Number of processors used                              2
Number of nodes in the mesh                          729
Number of equations solved                           729
Stopping criterion for PCG iterations         0.1000E-06

BOUNDARY CONDITION DATA                         
Number of nodes that were restrained                   0
Number of loaded freedoms                              2
Total energy delivered (load)                 0.2000E+04
Initial global temperature                    0.0000E+00
Number of fixed nodal values                           2

TRANSIENT STATE DETAILS                         
Section   dt               nstep        npri   npri_chk   time        #iters*   Tot#iters
      1   0.1250E-01         800           4        400   0.1000E+02        5        6098
-----------------------------------------------------------------------------------------
Totals                       800         200          0   0.1000E+02                 6098
-----------------------------------------------------------------------------------------
*Number of iterations to complete last timestep of section

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.171600    5.95
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
Solve equations                                 1.060809   36.76
Output results                                  1.653609   57.30
Total execution time                            2.886018  100.00

