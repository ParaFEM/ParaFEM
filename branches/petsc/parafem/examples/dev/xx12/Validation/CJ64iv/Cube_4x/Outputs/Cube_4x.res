
BASIC JOB DATA                                  
Number of processors used                              2
Number of nodes in the mesh                          125
Number of equations solved                           125
Stopping criterion for PCG iterations         0.1000E-07

BOUNDARY CONDITION DATA                         
Number of nodes that were restrained                   0
Initial global temperature                    0.1000E+01
Number of fixed nodal values                          98

TRANSIENT STATE DETAILS                         
Section   dt               nstep        npri   npri_chk   time        #iters*   Tot#iters
      1   0.5000E-07          10           1          1   0.5000E-06        5          46
      2   0.5000E-06           9           1          1   0.4500E-05        5          45
      3   0.5000E-05           9           1          1   0.4500E-04        4          44
      4   0.5000E-04           9           1          1   0.4500E-03        5          45
      5   0.5000E-03           9           1          1   0.4500E-02        5          45
      6   0.5000E-02           9           1          1   0.4500E-01        5          45
-----------------------------------------------------------------------------------------
Totals                        55          55         55   0.5000E-01                  270
-----------------------------------------------------------------------------------------
*Number of iterations to complete last timestep of section

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.140400   56.25
Read element steering array                     0.000000    0.00
Convert Abaqus to S&G node ordering             0.000000    0.00
Read nodal coordinates                          0.000000    0.00
Read restrained nodes                           0.000000    0.00
Compute steering array and neq                  0.000000    0.00
Compute interprocessor communication tables     0.000000    0.00
Allocate neq_pp arrays                          0.000000    0.00
Compute element stiffness matrices              0.046800   18.75
Build the preconditioner                        0.000000    0.00
Get starting r                                  0.000000    0.00
Solve equations                                 0.031200   12.50
Output results                                  0.031200   12.50
Total execution time                            0.249600  100.00

