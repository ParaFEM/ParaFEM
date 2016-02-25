
BASIC JOB DATA                                  
Number of processors used                              2
Number of nodes in the mesh                            5
Number of equations solved                             5
Stopping criterion for PCG iterations         0.1000E-05

BOUNDARY CONDITION DATA                         
Number of nodes that were restrained                   0
Initial global temperature                    0.2000E+03
Number of fixed nodal values                           1

TRANSIENT STATE DETAILS                         
Section   dt               nstep        npri   npri_chk   time        #iters*   Tot#iters
      1   0.1000E-02       20000          20      10000   0.2000E+02        5      100001
-----------------------------------------------------------------------------------------
Totals                     20000        1000          0   0.2000E+02               100001
-----------------------------------------------------------------------------------------
*Number of iterations to complete last timestep of section

PROGRAM SECTION EXECUTION TIMES                  SECONDS  %TOTAL    
Setup                                           0.187200    4.80
Read element steering array                     0.015600    0.40
Convert Abaqus to S&G node ordering             0.000000    0.00
Read nodal coordinates                          0.000000    0.00
Read restrained nodes                           0.000000    0.00
Compute steering array and neq                  0.000000    0.00
Compute interprocessor communication tables     0.000000    0.00
Allocate neq_pp arrays                          0.000000    0.00
Compute element stiffness matrices              0.000000    0.00
Build the preconditioner                        0.000000    0.00
Get starting r                                  0.000000    0.00
Solve equations                                 2.464821   63.20
Output results                                  1.232403   31.60
Total execution time                            3.900024  100.00

