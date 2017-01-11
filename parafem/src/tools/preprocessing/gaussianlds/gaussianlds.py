# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 09:25:35 2016

@author: evansl
"""

###############################################################################
# SECTION 0: Header info
###############################################################################

import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Check status, if = ON, additional calcs to show working (slower)
# Set to 'OFF' for quicker output of .lds
check = 'ON'

# Define constants to describe Multivariate Gaussian Distribution (MGD)
sigma = 10.892  # In order to achieve 10% of peak power at 5 mm from centre
rad = 7.5       # As specified by NETZSCH
PI = 4.0*np.arctan(1.0)
area_total = PI*rad**2
E = 5800        # Energy (mJ) deposited by laser over 15 mm2 during pulse
Fn = 1          # Scaling factor

'''
# Simon's example for error catching
try:
    with open('gegsg','wb') as f:
        f.write('3')
except:
    raise ValueError('started')
'''

fname = 'KFoam_CAD_VHR'      # Project name
pth = fname+'/'

# Set number of timesteps for laser pulse
n_t = 4000

# Set .time parameters
ntime = 6           # Total number of timesteps
nstep_t = 20000     # Number of substeps per timestep
npri = 50           # Substep output frequency
npri_chk = nstep_t  #Until checkpointing is fixed unused

# Set which plane to use i.e. xy, xz or yz
# NEED TO SET UP WAY TO AUTO SELECT CORRECT PLANE
# Find smallest bounds (xmax-xmin,ymax-ymin,zmax-zmin)
#plane = 'xy'

###############################################################################
# SECTION 1: Create array of x, y points to plot analytical solution
###############################################################################

d=20        # Width of domain
step=1      # Step size (dx)
rmin=-d/2
rmax=d/2
nstep=int(d/step)

# Initialise arrays
r = np.zeros(shape=(nstep,nstep),dtype=float)
r = np.arange(rmin, rmax, step) # Populate r with x,y coords
x = np.zeros(shape=((nstep+1)**2),dtype=float)
y = np.zeros(shape=((nstep+1)**2),dtype=float)
mgd = np.zeros(shape=((nstep+1)**2),dtype=float)

k=0
for i in range(0,nstep+1):
    for j in range(0,nstep+1):
        #print (i,j,k)
        x[k] = rmin + (j * step)
        y[k] = rmin + (i * step)
        r = np.sqrt((x[k])**2 + (y[k])**2)
        mgd[k] = (Fn/(2*PI*(sigma**2)))*np.exp((-1/2)*((r/sigma)**2))
        k=k+1

cmap=cm.coolwarm
plt.scatter(x, y, c=mgd, cmap=cmap, linewidths=0.0)
plt.axis('equal')
plt.suptitle('Analytical flux at nodal locations (mW/mm^2)')
plt.colorbar()
plt.show()

###############################################################################
# SECTION 2: Read in nodal data and mesh
###############################################################################

# Read in nodal set
#fname = 'nset.dat'
ext = '.nset'
nset = np.fromfile(pth+fname+ext,dtype=int,count=-1,sep=" ")
length_nset = int(len(nset))

# Superceeded by section below
'''
## Read in nodal coordinates
ext = '.coords'
coords = np.fromfile(pth+fname+ext,dtype=float,count=-1,sep=" ")
# Reshape coords array
ndim = 3
length = int(len(coords)/(ndim+1))
coords.shape=(length,ndim+1)
'''

# Read in nodal coordinates
ext = '.ensi.geo'
with open(pth+fname+ext) as f:
    for _ in range(9):
        next(f)
    for line in range(1):
        length = int(next(f))
        coords = np.zeros(shape=(length,4),dtype=float)
        #next(f)
    for line in range(0,length):
        coords[line,0] = line+1
        coords[line,1] = next(f)
    for line in range(0,length):
        coords[line,2] = next(f)
    for line in range(0,length):
        coords[line,3] = next(f)

plane_ = np.zeros(shape=(3),dtype=float)
plane_[0] = max(coords[:,1])-min(coords[:,1])
plane_[1] = max(coords[:,2])-min(coords[:,2])
plane_[2] = max(coords[:,3])-min(coords[:,3])
i=np.where(plane_ == plane_.min())
i=int(i[0])
if i==0:
    plane='yz'
if i==1:
    plane='xz'
if i==2:
    plane='xy'

# Cases depending on which plane is used
if plane=='xy': # do nothing
    print ('plane==xy')
if plane=='xz': # swap y and z
    print ('plane==xz')
    temp=coords
    coords[:,2]=temp[:,3]
    coords[:,3]=temp[:,2]
if plane=='yz': # swap x and z
    print ('plane==yz')
    temp=coords
    coords[:,1]=temp[:,3]
    coords[:,3]=temp[:,1]

# Place x,y coords for nset into points
points=np.zeros(shape=(length_nset,3),dtype=float)
for i,node in enumerate(nset):
    points[i,0] = coords[node-1,1]
    points[i,1] = coords[node-1,2]
    points[i,2] = coords[node-1,3]

# Calculate central point from bounds
xc = (max(coords[:,1]) + min(coords[:,1])) / 2
yc = (max(coords[:,2]) + min(coords[:,2])) / 2

# Use Delaunay to triangulate from point cloud
from scipy.spatial import Delaunay
tri = Delaunay(points[:,0:2])
plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
plt.plot(points[:,0], points[:,1], 'o')
plt.axis('equal')
plt.suptitle('Meshed nodal points')
for i in range(0,length_nset):
    plt.text(points[i,0]+0.4, points[i,1]+0.4, str(nset[i]), fontsize=12, color="red")
plt.show()

# Place calculated triangles into tri_surf array connected with steering array
steering_array = tri.simplices
tri_surf = points[tri.simplices]
ntri = int(len(steering_array))

###############################################################################
# SECTION 3: Calculate sum of normalised nodes
# Sum of loads tends to 1 as R tends to infinity
###############################################################################

# Initialise arrays
c_points = np.zeros(shape=(ntri,2),dtype=float)
c_lds = np.zeros(shape=(ntri),dtype=float)
dc = np.zeros(shape=(3),dtype=float)
lds = np.zeros(shape=(length_nset),dtype=float)
area_node = np.zeros(shape=(length_nset),dtype=float)
area_face_sum = 0.0

# Loop round each triangle in surface mesh
for i,tri_coord in enumerate(tri_surf):
    # Get nodal coords for triangle
    x1 = tri_coord[0,0]
    x2 = tri_coord[1,0]
    x3 = tri_coord[2,0]
    y1 = tri_coord[0,1]
    y2 = tri_coord[1,1]
    y3 = tri_coord[2,1]
    # Calculate absolute and fractional triangle area, keep tally of total area
    area_face = 0.5 * (abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)))
    area_frac = area_face/area_total
    area_face_sum = area_face_sum + area_face
    
    # Calculate centroid coords
    c_points[i,0] = (x1+x2+x3)/3
    c_points[i,1] = (y1+y2+y3)/3
    # Calculate distance of nodes from centroid
    dc[0] = np.sqrt((x1-c_points[i,0])**2 + (y1-c_points[i,1])**2)
    dc[1] = np.sqrt((x2-c_points[i,0])**2 + (y2-c_points[i,1])**2)
    dc[2] = np.sqrt((x3-c_points[i,0])**2 + (y3-c_points[i,1])**2)
    dsum = sum(dc)
    # Calculate distance of centroid from laser centre
    c_r = np.sqrt((c_points[i,0]-xc)**2 + (c_points[i,1]-yc)**2)
    
    # Calculate laser intensity at centroid location
    c_mgd = (Fn/(2*PI*(sigma**2)))*np.exp((-1/2)*((c_r/sigma)**2))
    # Save load for whole triangle surface area
    c_lds[i] = c_mgd*area_frac
    
    # Loop around points in triangle
    for k in range(0,3):
        # Calculate nodal load as fraction of triangle load
        # Sum into array, nodes will have contributions from multiple triangles
        lds[steering_array[i,k]] = lds[steering_array[i,k]] + c_lds[i]*dc[k]/dsum
        # Calculate area contribution for node from triangle
        area_node[steering_array[i,k]] = area_node[steering_array[i,k]] + area_face*(dc[k]/dsum)

# Plot loads at centroid locations
plt.scatter(c_points[:,0], c_points[:,1], c=c_lds[:], cmap=cmap, linewidths=0.0)
plt.axis('equal')
plt.suptitle('Triangle centroid loads (mW or mJ/s)')
plt.colorbar()
plt.show()

# Plot loads at nodal locations
plt.scatter(points[:,0], points[:,1], c=lds[:], cmap=cmap, linewidths=0.0)
plt.axis('equal')
plt.suptitle('Applied nodal loads (mW or mJ/s)')
plt.colorbar()
plt.show()

# Back calculate applied flux from nodal loads to check
check = np.zeros(shape=(length_nset),dtype=float)
for i in range(0,length_nset):
        check[i] = (lds[i] / area_node[i])*area_total

# Plot nodal flux
plt.scatter(points[:,0], points[:,1], c=check[:], cmap=cmap, linewidths=0.0)
plt.axis('equal')
plt.suptitle('Back calculated nodal flux (mW/mm^2)')
plt.colorbar()
plt.show()

# Report check results
print('')
print('sum area_face =',area_face_sum)
print('sum area_node =',sum(area_node))
# Below not needed: constant total laser area, not sample area
#print('area_total = PI*rad**2 =',area_total)
print('')
print('sum lds =',sum(lds))
print('#nds =',length_nset)
print('c_mgd*area_face_sum/area_total =',c_mgd*area_face_sum/area_total)
print('')
print('analytical flux at edge =',(Fn/(2*PI*(sigma**2)))*np.exp((-1/2)*(((12.7/2)/sigma)**2)))
print('check min = ',min(check))
print('analytical flux at centre =',(Fn/(2*PI*(sigma**2)))*np.exp((-1/2)*((0.0/sigma)**2)))
print('check max =',max(check))
print('')

###############################################################################
# SECTION 4: Calculate scaling factor to give peak flux at centre
###############################################################################

'''
# Old way of calculating intensity volume (approximation)
V_6pt4 = 1 - np.exp(-(1/2)*((6.4/sigma)**2))
V_7pt5 = 1 - np.exp(-(1/2)*((7.5/sigma)**2))
E_exp = (V_6pt4/V_7pt5)*E
# Value found iteratively
#Fn = Fn*5489751074
'''

# Correct way using error function
V_6pt4 = special.erf(np.sqrt(6.4)/(np.sqrt(2)*sigma))
V_7pt5 = special.erf(np.sqrt(7.5)/(np.sqrt(2)*sigma))
E_exp = (V_6pt4/V_7pt5)*E

# Read in laser power profile
fname_laser = 'lfa_amp.dat'
temp = np.fromfile(fname_laser,dtype=float,count=-1,sep=" ")
# Reshape amp array
amp = np.zeros(shape=(int(len(temp)/2),2),dtype=float)
k=0
for i in range(0,int(len(temp)/2)):
    for j in range(0,2):
        #print(j)
        amp[i,j]=temp[k]
        #print(k)
        k=k+1

# Interpolate laser amplitude to give desired temporal resolution
from scipy.interpolate import interp1d
amp_temp=np.zeros(shape=(len(amp[:,0])+1,2),dtype=float)
# Add row for t=0
amp_temp[1:len(amp[:,0])+1,:]=amp[0:len(amp[:,0]),:]
# Interpolate wrt t
amp_interp = interp1d(amp_temp[:,0],amp_temp[:,1])
# Set number of timesteps for laser pulse
#n_t = 2000 # Moved to top of code
t_max=0.0008
print('dt =',t_max/n_t,'(s)')
# Set total number of timesteps (zeros after pulse)
N_t = 500000
amp_out=np.zeros(shape=(N_t,2),dtype=float)
# ParaFEM bug, replace zeros with small number
#amp_out[:,1]=1.0E-14
for i in range(1,n_t):
    amp_out[i,0]=(t_max/n_t)*i
    amp_out[i,1]=amp_interp((t_max/n_t)*i)

# Output to file
ext = '.amp'
np.savetxt(pth+fname+ext, amp_out[:,1], fmt='%g')

#Calculate area under curve measured experimentally (s)    
laser = np.trapz(amp[:,1],x=amp[:,0])

# Plot laser amplitude vs time
plt.scatter(amp_out[0:n_t,0], amp_out[0:n_t,1])
plt.suptitle('Laser amplitude vs time')
plt.show()

amp         = None
amp_temp    = None
amp_interp  = None
amp_out     = None

# Above function reads in and calculates value below from text file
#laser = 0.000187172
Fn = E_exp / (laser*sum(lds))

if check=='OFF':
    Fn = 32034552174.016979

###############################################################################
# SECTION 5: Calculate nodal loads
###############################################################################

# Initialise arrays
c_points = np.zeros(shape=(ntri,2),dtype=float)
c_lds = np.zeros(shape=(ntri),dtype=float)
dc = np.zeros(shape=(3),dtype=float)
lds = np.zeros(shape=(length_nset),dtype=float)
area_node = np.zeros(shape=(length_nset),dtype=float)
area_face_sum = 0.0

# Loop round each triangle in surface mesh
for i,tri_coord in enumerate(tri_surf):
    # Get nodal coords for triangle
    x1 = tri_coord[0,0]
    x2 = tri_coord[1,0]
    x3 = tri_coord[2,0]
    y1 = tri_coord[0,1]
    y2 = tri_coord[1,1]
    y3 = tri_coord[2,1]
    # Calculate absolute and fractional triangle area, keep tally of total area
    area_face = 0.5 * (abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)))
    area_frac = area_face/area_total
    area_face_sum = area_face_sum + area_face
    
    # Calculate centroid coords
    c_points[i,0] = (x1+x2+x3)/3
    c_points[i,1] = (y1+y2+y3)/3
    # Calculate distance of nodes from centroid
    dc[0] = np.sqrt((x1-c_points[i,0])**2 + (y1-c_points[i,1])**2)
    dc[1] = np.sqrt((x2-c_points[i,0])**2 + (y2-c_points[i,1])**2)
    dc[2] = np.sqrt((x3-c_points[i,0])**2 + (y3-c_points[i,1])**2)
    dsum = sum(dc)
    # Calculate distance of centroid from laser centre
    c_r = np.sqrt((c_points[i,0]-xc)**2 + (c_points[i,1]-yc)**2)
    
    # Calculate laser intensity at centroid location
    c_mgd = (Fn/(2*PI*(sigma**2)))*np.exp((-1/2)*((c_r/sigma)**2))
    # Save load for whole triangle surface area
    c_lds[i] = c_mgd*area_frac
    
    # Loop around points in triangle
    for k in range(0,3):
        # Calculate nodal load as fraction of triangle load
        # Sum into array, nodes will have contributions from multiple triangles
        lds[steering_array[i,k]] = lds[steering_array[i,k]] + c_lds[i]*dc[k]/dsum
        # Calculate area contribution for node from triangle
        area_node[steering_array[i,k]] = area_node[steering_array[i,k]] + area_face*(dc[k]/dsum)

# Plot loads at centroid locations
plt.scatter(c_points[:,0], c_points[:,1], c=c_lds[:], cmap=cmap, linewidths=0.0)
plt.axis('equal')
plt.suptitle('Triangle centroid loads (mW or mJ/s)')
plt.colorbar()
plt.show()

# Plot loads at nodal locations
plt.scatter(points[:,0], points[:,1], c=lds[:], cmap=cmap, linewidths=0.0)
plt.axis('equal')
plt.suptitle('Applied nodal loads (mW or mJ/s)')
plt.colorbar()
plt.show()

# Back calculate applied flux from nodal loads to check
check = np.zeros(shape=(length_nset),dtype=float)
for i in range(0,length_nset):
        check[i] = (lds[i] / area_node[i])*area_total

# Plot nodal flux
plt.scatter(points[:,0], points[:,1], c=check[:], cmap=cmap, linewidths=0.0)
plt.axis('equal')
plt.suptitle('Back calculated nodal flux (mW/mm^2)')
plt.colorbar()
plt.show()

# Below not needed: constant total laser area, not sample area
#print('area_total = PI*rad**2 =',area_total)
print('')
print('sum lds =',sum(lds))
print('#nds =',length_nset)
print('c_mgd*area_face_sum/area_total =',c_mgd*area_face_sum/area_total)
print('')
print('analytical flux at edge =',(Fn/(2*PI*(sigma**2)))*np.exp((-1/2)*(((12.7/2)/sigma)**2)))
print('check min = ',min(check))
print('analytical flux at centre =',(Fn/(2*PI*(sigma**2)))*np.exp((-1/2)*((0.0/sigma)**2)))
print('check max =',max(check))
print('')

###############################################################################
# SECTION 6: Output loads to file
###############################################################################

# Output to file
# ParaFEM format
lds_print = np.zeros(shape=(length_nset,4),dtype=float)
lds_print[:,0] = nset
lds_print[:,1] = lds
ext = '.lds'
np.savetxt(pth+fname+ext, lds_print, fmt='%i %g %f %f')
# ENSI format
lds_print=np.zeros(shape=(length),dtype=float)
for i,loaded_node in enumerate(nset):
    lds_print[loaded_node-1] = lds[i]
ext = '.ensi.NDLDS'
np.savetxt(pth+fname+ext, lds_print, fmt='%g',newline='\n',comments='',
           header='Alya Ensight Gold --- Scalar per-node variable file\n'
           'part\n1\ncoordinates')
# Output calculated fluxes in ENSI (for vis purposes)
lds_print=np.zeros(shape=(length),dtype=float)
for i,loaded_node in enumerate(nset):
    lds_print[loaded_node-1] = check[i]
ext = '.ensi.NDFLX'
np.savetxt(pth+fname+ext, lds_print, fmt='%g',newline='\n',comments='',
           header='Alya Ensight Gold --- Scalar per-node variable file\n'
           'part\n1\ncoordinates')


# Output .time file, parameters set at start of program
ext = '.time'
time_print=np.zeros(shape=(ntime,5),dtype=float)
for i in range(0,ntime):
    #print(i)
    time_print[i,0]=i+1
    #print(2**(i))
    #time_print[i,1]=(t_max/n_t)*(2**(i))
    time_print[i,1]=(t_max/n_t)*(2**(i))
    time_print[i,2]=nstep_t
    time_print[i,3]=npri
    time_print[i,4]=npri_chk
np.savetxt(pth+fname+ext, time_print, fmt='%i %g %i %i %i',newline='\n',comments='',
           header='*TIMESTEPS '+str(ntime)+' 1 3\ndtim nstep npri npri_chk')
print('dt_0                    = ',t_max/n_t)
print('Number of time steps    = ',ntime)
print('Number of time substeps = ',nstep_t*ntime)
print('Total time (s)          = ',sum(time_print[:,1])*(nstep_t))
print('Number of printed steps = ',int((nstep_t/npri)*ntime))

# Add in sum lds over timeplpl
