Explanation of the variables used in the p124_5_size.mg files.

program
'p124'

iotype    nels   nxe  nze  nip
'parafem' 125000 50   50   8

aa   bb   cc   kx  ky  kz  rho cp
0.02 0.02 0.02 1.0 1.0 1.0 1.0 1.0

dtim  nsteps theta
0.01  150    0.5

npri  tol     limit  val0
10    0.0001  100    100.0

np_types  loaded_freedoms  fixed_freedoms
1         0                0 

aa                x-dimension of elements
bb                y-dimension of elements 
cc                z-dimension of elements
cp                heat capacity
dtim              time step
fixed_freedoms    number of fixed freedoms
iotype            file format
kx                conductivity in x-direction
ky                conductivity in y-direction
kz                conductivity in z-direction
limit             iteration ceiling
loaded_freedoms   number of loaded freedoms
nels              number of elements in mesh
nip               number of integrating points
np_types          number of property types
npri              print interval
nstep             number of time steps in analysis
nxe               elements in x-direction
nze               elements in z-direction
program           driver program type
rho               density
theta             parameter in 'theta' integrator
tol               convergence tolerance
val0              initial value

