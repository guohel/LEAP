import sys
import os
import time
import numpy as np
from leapctype import *
import math
from typing import List, Union

leapct = tomographicModels()


def rad2deg(rad):
    deg = rad / np.pi * 180
    return deg


def deg2rad(deg):
    rad = deg * np.pi / 180
    return rad


'''
This demo script shows another example of the modular-beam geometry, except here
the source and detectors are arranged around a sphere which mimics what one may have in a flash CT configuration
like the MEFCT at ARL: https://pubs.aip.org/aip/acp/article-pdf/doi/10.1063/1.5045031/14167597/160032_1_online.pdf
'''

# Specify the number of detector columns which is used below
# Scale the number of angles and the detector pixel size with N
numCols = 368
numDetector = 1
numSubAngle = 50
numAngles = numSubAngle * numDetector  # numSubAngle*numDetector interval
# pixelSize = 0.65*512/float(numCols)*2.0*11.0/14.0
pixelSize = 0.4
# Set the number of detector rows
numRows = numCols  # 76#256
sod = 200.0
sdd = 600.0

# leapct.set_conebeam(numAngles, numRows, numCols, pixelSize, pixelSize, 0.5*(numRows-1), 0.5*(numCols-1), leapct.setAngleArray(numAngles, 360.0), sod, sdd)
# Set the scanner geometry
sourcePositions = np.ascontiguousarray(np.zeros((numAngles, 3)).astype(np.float32), dtype=np.float32)
moduleCenters = np.ascontiguousarray(np.zeros((numAngles, 3)).astype(np.float32), dtype=np.float32)
colVectors = np.ascontiguousarray(np.zeros((numAngles, 3)).astype(np.float32), dtype=np.float32)
rowVectors = np.ascontiguousarray(np.zeros((numAngles, 3)).astype(np.float32), dtype=np.float32)

# geometry setting

pos_z = np.arange(-25, 25, 1)
print(pos_z)
anglesRatio_z = sod / pos_z
print(anglesRatio_z)
angles = anglesRatio_z
for n in range(len(anglesRatio_z)):
    angles[n] = np.arctan(anglesRatio_z[n])

print(angles)
ind = 0
for n in range(len(anglesRatio_z)):
    sourcePositions[ind, 2] = 0
    sourcePositions[ind, 1] = sod
    sourcePositions[ind, 0] = -25 + n
    ind += 1
'''
for n in range(len(anglesRatio_z)): 
     sourcePositions[ind,0] = 0
     sourcePositions[ind,1] = -sod
     sourcePositions[ind,2] =-25 + n
     ind+=1  
'''
ind2 = 0
for n in range(len(anglesRatio_z)):
    moduleCenters[ind2, 2] = 0
    moduleCenters[ind2, 1] = -(sdd - sod)
    moduleCenters[ind2, 0] = -(sdd - sod) / np.tan(angles[n])
    ind2 += 1
'''
for n in range(len(anglesRatio_z)): 
     moduleCenters[ind2,0] = 0
     moduleCenters[ind2,1] = (sdd-sod)
     moduleCenters[ind2,2] = (sdd-sod) / np.tan(angles[n])
     ind2+=1  
'''
print(sourcePositions)
print(moduleCenters)

for n in range(len(anglesRatio_z)):
    rowVectors[n, 0] = -1
    rowVectors[n, 1] = 0
    rowVectors[n, 2] = 0

ind3 = 0
for n in range(len(anglesRatio_z)):
    colVectors[ind3, :] = np.cross(sourcePositions[ind3, :], rowVectors[ind3, :])
    ind3 += 1
'''    
for n in range(len(anglesRatio_z)): 
    colVectors[ind3,:] = np.cross(-sourcePositions[ind3,:], rowVectors[ind3,:])
    ind3+=1
'''

# print(colVectors)
# print(rowVectors)
# print(sourcePositions)
# print(moduleCenters)
leapct.set_modularbeam(numAngles, numRows, numCols, pixelSize, pixelSize, sourcePositions, moduleCenters, rowVectors,
                       colVectors)

# Set the volume parameters
leapct.set_default_volume()
# leapct.set_volume(leapct.get_numX(), leapct.get_numY(), leapct.get_numZ(), leapct.get_voxelWidth(), leapct.get_voxelHeight(), 0, -3, 0)


# Get the source positions, detector positions, and detector orientation for all projections
sourcePositions = leapct.get_sourcePositions()
moduleCenters = leapct.get_moduleCenters()
rowVectors = leapct.get_rowVectors()
colVectors = leapct.get_colVectors()

# Trouble-Shooting Functions
leapct.print_parameters()
# leapct.sketch_system()
# quit()

# Allocate space for the projections and the volume
g = leapct.allocateProjections()
f = leapct.allocateVolume()

# Specify simplified FORBILD head phantom
# leapct.set_FORBILD(f,True)

leapct.addObject(f, 4, np.array([0.0, 0.0, 0.0]), 120.0 * np.array([1.0, 1.0, 0.1]), 0.02, None, None, 3)
z_step = 4.0

leapct.addObject(f, 1, np.array([-60.0, 0.0, -3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([-40.0, 0.0, -3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([-20.0, 0.0, -3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 0.0, -3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([20.0, 0.0, -3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([40.0, 0.0, -3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([60.0, 0.0, -3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)

leapct.addObject(f, 1, np.array([0.0, -60.0, -2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, -40.0, -2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, -20.0, -2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 0.0, -2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 20.0, -2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 40.0, -2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 60.0, -2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)

leapct.addObject(f, 1, np.array([-60.0, 0.0, -z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([-40.0, 0.0, -z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([-20.0, 0.0, -z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 0.0, -z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([20.0, 0.0, -z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([40.0, 0.0, -z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([60.0, 0.0, -z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)

leapct.addObject(f, 1, np.array([0.0, -60.0, 0.0]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, -40.0, 0.0]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, -20.0, 0.0]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 0.0, 0.0]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 20.0, 0.0]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 40.0, 0.0]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 60.0, 0.0]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)

leapct.addObject(f, 1, np.array([-60.0, 0.0, z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([-40.0, 0.0, z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([-20.0, 0.0, z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 0.0, z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([20.0, 0.0, z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([40.0, 0.0, z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([60.0, 0.0, z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)

leapct.addObject(f, 1, np.array([0.0, -60.0, 2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, -40.0, 2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, -20.0, 2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 0.0, 2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 20.0, 2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 40.0, 2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 60.0, 2 * z_step]), np.array([80.0, 1.1, 1.0]), 0.04, None, None, 3)

leapct.addObject(f, 1, np.array([-60.0, 0.0, 3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([-40.0, 0.0, 3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([-20.0, 0.0, 3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([0.0, 0.0, 3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([20.0, 0.0, 3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([40.0, 0.0, 3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)
leapct.addObject(f, 1, np.array([60.0, 0.0, 3 * z_step]), np.array([1.1, 80.0, 1.0]), 0.04, None, None, 3)

leapct.display(f)
# quit()
# "Simulate" projection data
startTime = time.time()
leapct.project(g, f)
print('Forward Projection Elapsed Time: ' + str(time.time() - startTime))
# leapct.display(g)
# quit()


# Reset the volume array to zero, otherwise iterative reconstruction algorithm will start their iterations
# with the true result which is cheating
f[:] = 0.0

# leapct.display(g)
# quit()
# Copy data to GPU
# Comment this section out to revert back to multi-GPU solution
# with CPU-GPU data transfers to see when each case is advantageous
if has_torch:
    device_name = "cuda:" + str(leapct.get_gpu())
    device = torch.device(device_name)
    g = torch.from_numpy(g).to(device)
    f = torch.from_numpy(f).to(device)

# Reconstruct the data
startTime = time.time()
# leapct.backproject(g,f)
# leapct.FBP(g,f)
# leapct.inconsistencyReconstruction(g,f)
# leapct.print_cost = True
filters = filterSequence(1.0e1)
filters.append(TV(leapct, delta=0.02 / 40.0))
# filters.append(TV(leapct, delta=0.02))
# leapct.LS(g,f,50,'SQS')
# leapct.RLS(g,f,100,filters, 'SQS')
# leapct.RWLS(g,f,20,filters,None,'SQS')
leapct.RDLS(g, f, 10, filters, 1.0, True, 1)

print('Reconstruction Elapsed Time: ' + str(time.time() - startTime))

leapct.display(f)
