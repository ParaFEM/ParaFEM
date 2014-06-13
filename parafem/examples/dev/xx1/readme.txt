Program xx1 supports either load control or displacement control, but not
both at the same time.

If loaded_nodes > 0, the program will read xx1-tiny.lds

If fixed_freedoms > 0, the program will read xx1-tiny.fix

CONTROL FILE FORMAT

<job_name>.dat contains the variables named below:

ELEMENT
MESH
PARTITION
NELS NN NR NIP NOD FIXED_FREEDOMS LOADED_NODES
E V TOL LIMIT

An example follows:

xx1-tiny.dat

'hexahedron'
2
1
125  756  396  8  20  8  0
0.1000E+03  0.3000E+00  0.1000E-04  200


