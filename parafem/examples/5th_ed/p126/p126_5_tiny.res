
Number of processors used:                             1
Number of nodes in the mesh:                        4961
Number of nodes that were restrained:               4232
Number of equations solved:                        11138

Total BiCGSTAB(l) iterations:                        723

                  iterations        norm
               1         105  0.3423E+02
               2          70  0.1032E+02
               3          78  0.2931E+01
               4          74  0.7937E+00
               5          79  0.2532E+00
               6          74  0.1134E+00
               7          81  0.6559E-01
               8          81  0.3384E-01
               9          81  0.1901E-01

Program section execution times                   Seconds  %Total    
Setup                                            0.107580    0.21
Compute coordinates and steering array           0.002594    0.01
Compute interprocessor communication tables      0.003401    0.01
Allocate neq_pp arrays                           0.003478    0.01
Organize fixed nodes                             0.001426    0.00
Element stiffness integration                    1.802384    3.54
Build the preconditioner                         0.015210    0.03
Solve equations                                 48.897118   95.92
Output results                                   0.038831    0.08
Total execution time                            50.974704  100.00
