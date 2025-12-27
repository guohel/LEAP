import sys
import os
import time
import numpy as np
from leapctype import *
import math
from typing import List, Union

def batch_rotate_z(vectors, degrees):
    """
    批量绕Z轴旋转多个三维向量

    参数:
        vectors: 形状为 (n, 3) 的数组，每行是一个三维向量
        degrees: 旋转角度（度）

    返回:
        旋转后的向量数组
    """
    theta = math.radians(degrees)

    # 旋转矩阵
    rotation_matrix = np.array([
        [math.cos(theta), -math.sin(theta), 0],
        [math.sin(theta), math.cos(theta), 0],
        [0, 0, 1]
    ])

    # 批量旋转
    return np.dot(vectors, rotation_matrix.T).astype(np.float32)


def batch_rotate_y(vectors, degrees):
    """批量绕Y轴旋转"""
    theta = math.radians(degrees)

    rotation_matrix = np.array([
        [math.cos(theta), 0, math.sin(theta)],
        [0, 1, 0],
        [-math.sin(theta), 0, math.cos(theta)]
    ])

    return np.dot(vectors, rotation_matrix.T).astype(np.float32)


def batch_rotate_x(vectors, degrees):
    """批量绕X轴旋转"""
    theta = math.radians(degrees)

    rotation_matrix = np.array([
        [1, 0, 0],
        [0, math.cos(theta), -math.sin(theta)],
        [0, math.sin(theta), math.cos(theta)]
    ])

    return np.dot(vectors, rotation_matrix.T).astype(np.float32)


def generate_angle_array(
        start_angle: float,
        sweep_angle: float,
        num_points: int,
        return_numpy: bool = True
) -> Union[List[float], np.ndarray]:
    """
    优化的角度数组生成函数（使用向量化操作，效率更高）

    参数:
        start_angle: 起始角度（度）
        sweep_angle: 需要转过的总角度（度）
        num_points: 生成的角度数量
        return_numpy: 是否返回numpy数组，默认为True

    返回:
        浮点角度数组
    """
    if num_points <= 0:
        raise ValueError("角度数量必须大于0")

    # 使用numpy的linspace生成线性间隔的角度
    angles = np.linspace(start_angle, start_angle + sweep_angle, num_points, dtype=np.float32)

    # 使用向量化操作归一化到0-360范围
    normalized_angles = angles % 360.0

    # 根据参数决定返回类型
    return normalized_angles if return_numpy else normalized_angles.tolist()



leapct = tomographicModels()
leapct.about()

'''
This demo script simulates and reconstructs cone-beam laminography data.

Note that in LEAP the sources and detectors move-not the object.  Thus one will need to specify the source and detector positions
from the point of view of a fixed object.

For this we will use the modular-beam geometry.  One could model a limited number of cases with the cone-beam
geometry by shifting the detecto/r (centerRow) and the volume (offsetZ), but this is pretty limited.  Modular-beam
geometries allow more custom setups because you can specify the source location, detector location, and detector orientation
for all projections anywhere you want.  For convenience, we will start with a cone-beam geometry, convert it to modular-beam
and then make the adjustments for a laminography setup.
'''

# Specify the number of detector columns which is used below
# Scale the number of angles and the detector pixel size with numCols
numCols = 1024
numAngles = int(360)
pixelSize = 0.2

# Set the number of detector rows
numRows = 768


# Set the scanner geometry
sod = 28.451 # source-to-object distance (mm)
sdd = 608 # source-to-detector distance (mm)
# leapct.set_conebeam(numAngles, numRows, numCols, pixelSize, pixelSize, 0.5*(numRows-1), 0.5*(numCols-1), leapct.setAngleArray(numAngles, numAngles), sod, sdd)
angs = generate_angle_array(0, numAngles - 1, numAngles) # 260度测试
angs = np.ascontiguousarray(angs[::-1])
leapct.set_conebeam(numAngles, numRows, numCols, pixelSize, pixelSize, 0.5*(numRows-1), 0.5*(numCols-1), angs, sod, sdd)
# leapct.sketch_system()

voxel = pixelSize / (sdd / sod)

# Set the lamonography angle which is the rotation of the axis of rotation from the z-axis
# 不知道上面的表达是否有误，实际测试laminographyAngle是与X轴的夹角，如设置10°，则射线源跟平板的位置上升或者下降很小
laminographyAngle = 26.2 # degrees

# Switch to modular-beam coordinates
leapct.convert_to_modularbeam()

# Get the source positions, detector positions, and detector orientation for all projections
# convert_to_modularbeam之后sourcePositions、moduleCenters、rowVecs、colVecs才会有值，否则全部都为0
sourcePositions = leapct.get_sourcePositions()
moduleCenters = leapct.get_moduleCenters()
rowVecs = leapct.get_rowVectors()
colVecs = leapct.get_colVectors()
# leapct.sketch_system('yz', [0, 45, 90, 135, 180, 250])

# 平面CT平板绕Z轴做圆周运动的时候，平板是不会自旋
# colVecs = np.full((numAngles, 3), [0.0, -1.0, 0.0], dtype=np.float32)
# rowVecs = np.full((numAngles, 3), [1.0, 0.0, 0.0], dtype=np.float32)

# colVecs = batch_rotate_z(colVecs, -2)
# rowVecs = batch_rotate_z(rowVecs, -2)

# Shift the source up and the detector down so that the source and detector are
# aiming down by "laminographyAngle" degrees
# one could also rotate the "rowVecs" parameter is necessary, but note that if this
# is rotated more than 5 degrees the FBP reconstruction algorithms will not work
# and one will be required to reconstruct with an iterative method
sourcePositions[:,2] = np.tan((90-laminographyAngle)*np.pi/180.0)*sod
moduleCenters[:,2] = -np.tan((90-laminographyAngle)*np.pi/180.0)*(sdd-sod)
colVecs = np.full((numAngles, 3), [0.0, -1.0, 0.0], dtype=np.float32)
rowVecs = np.full((numAngles, 3), [1.0, 0.0, 0.0], dtype=np.float32)

# leapct.set_modularbeam(numAngles, numRows, numCols, pixelSize, pixelSize, sourcePositions, moduleCenters, rowVecs, colVecs)
# leapct.sketch_system('yz', [0, 45, 90, 135, 180, 250])

# Use the following 7 lines to rotate the detector as well
from scipy.spatial.transform import Rotation as R
# sin_theta = np.sin(-0.5*laminographyAngle*np.pi/180.0)
# cos_theta = np.cos(-0.5*laminographyAngle*np.pi/180.0)
sin_theta = np.sin(-0.5*90*np.pi/180.0)
cos_theta = np.cos(-0.5*90*np.pi/180.0)
for n in range(numAngles):
    q = np.append(colVecs[n,:].copy()*sin_theta, cos_theta)
    A = R.from_quat(q).as_matrix()
    # rowVecs[n,:] = np.matmul(A, rowVecs[n,:])

# Now re-set the modular-beam geometry with the modified source and detector locations
leapct.set_modularbeam(numAngles, numRows, numCols, pixelSize, pixelSize, sourcePositions, moduleCenters, rowVecs, colVecs)
# leapct.sketch_system('yz', [0, 45, 90, 135, 180, 225, 270, 315])

# Set the volume parameters.
# It is best to do this after the CT geometry is set

# leapct.set_volume(400, 580, 2*int(np.ceil(2.15 / voxel+1)), voxel, voxel, 0, 1.5, 0)
# leapct.set_volume(900, 900, 462, voxel, voxel, 0, 0, 0)
# leapct.set_diameterFOV(5)
leapct.set_default_volume(1.5)
leapct.set_numZ(40)
leapct.set_offsetZ(-0.66)
# leapct.set_numZ(2*int(np.ceil(2.15 / leapct.get_voxelHeight()+1))) # reduce the number of z-slices to those that just occupy the object
# leapct.set_diameterFOV(1.0e16)

# If you want to specify the volume yourself, use this function:
#leapct.set_volume(numX, numY, numZ, voxelWidth=None, voxelHeight=None, offsetX=None, offsetY=None, offsetZ=None):

# Trouble-Shooting Functions
# Print the parameters to the screen.  We also plot 5 of the projections to ensure
# that the geometry was set properly
leapct.print_parameters()
leapct.sketch_system('yz', [0, 45, 90, 135, 180, 250])


# Allocate space for the projections and the volume
# You don't have to use these functions; they are provided just for convenience
# All you need is for the data to be C contiguous float32 arrays with the right dimensions

g = leapct.allocateProjections() # shape is numAngles, numRows, numCols
# data_dir =  '/mnt/d/20-51-38-366-14-40-361-3/normalized_tif/'
data_dir =  '/mnt/d/20-51-38-366-14-40-361-3/'
g = leapct.load_projections(data_dir + 'proj.tif')
# 对数化并进行归一化
# I_0 = 14000.0
I_0 = 65535.0
g[g>I_0] = I_0
# g = g / 65535.0
g = -np.log(g/I_0)


# leapct.display(g)
f = leapct.allocateVolume() # shape is numZ, numY, numX


# Add noise to the data (just for demonstration purposes)
# I_0 = 50000.0
# g[:] = -np.log(np.random.poisson(I_0*np.exp(-g))/I_0)

# Reset the volume array to zero, otherwise iterative reconstruction algorithm will start their iterations
# with the true result which is cheating
f[:] = 0.0

# Copy data to GPU
#'''
# Comment this section out to revert back to multi-GPU solution
# with CPU-GPU data transfers to see when each case is advantageous
if has_torch:
    device_name = "cuda:" + str(leapct.get_gpu())
    device = torch.device(device_name)
    g = torch.from_numpy(g).to(device)
    f = torch.from_numpy(f).to(device)
#'''

# Reconstruct the data
# We will start with an FBP reconstruction and refine it with an iterative method to remove some of the
# so-called "cone-beam" artifacts.  There are many iterative reconstruction algorithms to choose from.
startTime = time.time()
#leapct.backproject(g,f)
# leapct.FBP(g,f)    # Error: FBP only implemented for modular geometries whose rowVectors are aligned with the z-axis
# leapct.display(f)
filters = filterSequence(1.0e0)
# filters.append(TV(leapct, delta=0.02/20.0))
filters.append(TV(leapct, delta=0.0000002/20.0))
# filters.append(histogramSparsity(leapct, mus=[0.3, 0.5],weight=0.01))
# leapct.ASDPOCS(g,f,1,360,1,filters)
leapct.SART(g,f,1,360)
# leapct.OSEM(g,f,1,360)
# leapct.LS(g,f,50,'SQS')
# leapct.RWLS(g,f,10,filters,None,'SARR')
# leapct.set_numTVneighbors(6)  # 提速
# leapct.RWLS(g,f,10,filters,None,'SQS') # 同样的图片，同样的参数，RWLS算法z方向体积的大小会影响重建质量，z越大越清晰
# leapct.RDLS(g,f,50,filters,1.0,True,1)
# leapct.MLTR(g,f,10,10,filters)
print('Reconstruction Elapsed Time: ' + str(time.time()-startTime))

# leapct.display(f)
# Post Reconstruction Smoothing (optional)
# Here are some optional post reconstruction noise filters that can be applied
# Try uncommenting out these lines to test how they work
# startTime = time.time()
# leapct.diffuse(f,0.00002/20.0,8)
# leapct.MedianFilter(f)
#leapct.BlurFilter(f,2.0)
# print('Post-Processing Elapsed Time: ' + str(time.time()-startTime))

# Display the result with napari
# leapct.save_volume(os.path.join('/mnt/d/20-51-38-366-14-40-361-3/recon/', 'recon.tif'), f)
leapct.display(f)
#import matplotlib.pyplot as plt
#plt.imshow(np.squeeze(f[f.shape[0]//2,:,:]), cmap='gray')
#plt.show()
