
BASIC JOB DATA                                  
Number of processors used                              2
Number of nodes in the mesh                          125
Number of equations solved                           125
Stopping criterion for PCG iterations         0.1000E-05

BOUNDARY CONDITION DATA                         
Number of nodes that were restrained                   0
Initial global temperature                    0.2000E+03
Number of fixed nodal values                           1

TRANSIENT STATE DETAILS                         
Section   dt               nstep        npri   npri_chk   time        #iters*   Tot#iters
      1   0.1000E-02       20000          20      10000   0.2000E+02        1       48653
-----------------------------------------------------------------------------------------
Totals                     20000        1000          0   0.2000E+02                48653
-----------------------------------------------------------------------------------------
*Number of iterations to complete last timestep of section

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.218400    2.54
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
Solve equations                                 5.928039   68.97
Output results                                  2.449215   28.49
Total execution time                            8.595654  100.00

