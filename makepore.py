##################################################################################################################################
# Making Boundaries from inside out
# WORKS, so far - need to adjust the tolerance value to make sense
##################################################################################################################################
import numpy as np
import math
from scipy.interpolate import make_interp_spline
from matplotlib import colors
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.rcParams['text.usetex'] = True
#get coords of point
filename = "ionpore3.pdb"
Lattice = np.zeros([160,165,122])
NAtoms = 0
XYZFile = np.zeros([3])
f = open(filename,'r')
lines = f.readlines()
for i, line in enumerate(lines):
    vline = line.rstrip()
    words = vline.split()
    if (words[0] == 'ATOM'):
        xcoord = math.floor(float(vline[30:38]))
        ycoord = math.floor(float(vline[38:46]))
        zcoord = math.floor(float(vline[46:54]))
        Lattice[xcoord][ycoord][zcoord] = 1
        newrow = [xcoord, ycoord, zcoord]
        XYZFile = np.vstack([XYZFile, newrow])
f.close()
XYZFile = np.delete(XYZFile, 0, 0) #removing the first zeros line
XYZFile = np.unique(XYZFile, axis=0) #getting unique values only
# Constructing the input point data
x = XYZFile[:,0]
y = XYZFile[:,1]
z = XYZFile[:,2]
xmax = int(np.amax(XYZFile[:,0]))
ymax = int(np.amax(XYZFile[:,1]))
zmax = int(np.amax(XYZFile[:,2]))
xmin = int(np.amin(XYZFile[:,0]))
ymin = int(np.amin(XYZFile[:,1]))
zmin = int(np.amin(XYZFile[:,2]))
#getting the avg point in every z plane
#for every z plane it will store the x and y coord of the midpt
zMids = np.zeros([122,3], dtype=int)
for zc in range (zmin,zmax+1):
    sliceatoms = sumx = sumy = 0
    for i in range (len(z)):
        if (z[i] == zc):
            sliceatoms += 1
            sumx += x[i]
            sumy += y[i]
    zMids[zc][0] = int(sumx/sliceatoms) if sliceatoms > 0 else -1
    zMids[zc][1] = int(sumy/sliceatoms) if sliceatoms > 0 else -1
    zMids[zc][2] = int(zc) if sliceatoms > 0 else -1
#remove all zero points in zMids
zMids2 = zMids[~np.all(zMids == 0, axis=1)]
xmax2 = int(np.amax(zMids2[:,0]))
ymax2 = int(np.amax(zMids2[:,1]))
zmax2 = int(np.amax(zMids2[:,2]))
xmin2 = int(np.amin(zMids2[:,0]))
ymin2 = int(np.amin(zMids2[:,1]))
zmin2 = int(np.amin(zMids2[:,2]))
# print(zMids)
#Find the atom closest to the midpoint, and fill any sites in the lattice that are closer than that
for zc in range (len(zMids)):
    if (zmin+5 <= zc <= zmax2-17):
        zslice = zMids[zc][2]
        if (zslice != 0):
            dmin = 1000 #Large enough region
            #Find the smallest distance (the largest empty r)
            for a in range (len(XYZFile)):
                if (XYZFile[a][2] == zslice):
                    temp_d = math.sqrt(((zMids[zc][0]-XYZFile[a][0])**2)+((zMids[zc][1]-XYZFile[a][1])**2))
                    dmin = temp_d if temp_d <= dmin else dmin
            # Now fill 2's in place of the empty regions
            # for i in range (xmin,xmax+1):
            #     for j in range (ymin,ymax+1):
            # Reduce Size by one unit to eliminate the last ends
            for i in range (xmin,xmax):
                for j in range (ymin,ymax):
                    temp_r = math.sqrt(((zMids[zc][0]-i)**2)+((zMids[zc][1]-j)**2))
                    Lattice[i][j][zc] = 2 if ((temp_r <= dmin) and (Lattice[i][j][zc] == 0)) else Lattice[i][j][zc]
# Any Lattice site with 2+ pore neighbours is also a pore
tolerance = 2
ep = 0
for t in range (tolerance):
    for i in range (xmin2,xmax2):
        for j in range (ymin2,ymax2):
            for k in range (zmin2+5,zmax2-17):
                flags = np.array([Lattice[i-1][j][k], Lattice[i+1][j][k], Lattice[i-2][j][k], Lattice[i+2][j][k],
                                  Lattice[i][j-1][k], Lattice[i][j+1][k], Lattice[i][j-2][k], Lattice[i][j+2][k],
                                  Lattice[i][j][k-1], Lattice[i][j][k+1], Lattice[i][j][k-2], Lattice[i][j][k+2]], dtype = int)
                p = flags.tolist().count(2) #counts number of pore points
                q = flags.tolist().count(1) #counts number of atom points
                if (p==1):
                    ep +=1 #counts number of end points
                Lattice[i][j][k] = 2 if ((Lattice[i][j][k] == 0) and (p >= 3) and (q == 0)) else Lattice[i][j][k]

Pore = np.zeros([3])
for i in range (160):
    for j in range (165):
        for k in range (122):
            if (Lattice[i][j][k] == 2):
                nr = [i,j,k]
                Pore = np.vstack([Pore, nr])
Pore = np.delete(Pore, 0, 0) #removing the first zeros line
Pore = np.unique(Pore, axis=0) #getting unique values only
print (len(Pore)) # 5712 with tol=1 p>=2
print (ep) # 1541
o = open('PoreCoordinates1f.xyz','w')
for i in range (len(Pore)):
    o.write("H" + "   " + str(Pore[i][0]) + "   " + str(Pore[i][1]) + "   " + str(Pore[i][2]) + "\n")
o.close()

xarr = np.zeros([xmax2-xmin2-4],dtype=int)
yarr = np.zeros([ymax2-ymin2-4],dtype=int)
zarr = np.zeros([zmax2-zmin2-4],dtype=int)
count = 0 
for i in range (xmin2+2,xmax2-2):
    count = 0 
    for j in range (ymin2+2,ymax2-2):
        for k in range (zmin2+2,zmax2-2):
            if (Lattice[i][j][k] == 2):
                count = count+1
    xarr[i-xmin2-2] = count
for j in range (ymin2+2,ymax2-2):
    count = 0 
    for i in range (xmin2+2,xmax2-2):
        for k in range (zmin2+2,zmax2-2):
            if (Lattice[i][j][k] == 2):
                count = count+1
    yarr[j-ymin2-2] = count
for k in range (zmin2+2,zmax2-2):
    count = 0 
    for j in range (ymin2+2,ymax2-2):
        for i in range (xmin2+2,xmax2-2):
            if (Lattice[i][j][k] == 2):
                count = count+1
    zarr[k-zmin2-2] = count
#index of cross section to be taken
xcross = np.argmax(xarr)
ycross = np.argmax(yarr)
zcross = np.argmax(zarr)
x_coordinate = [i+xmin2+2 for i in range(len(xarr))]
y_coordinate = [i+ymin2+2 for i in range(len(yarr))]
z_coordinate = [i+zmin2+2 for i in range(len(zarr))]
x_coordinate = np.array(x_coordinate)
y_coordinate = np.array(y_coordinate)
z_coordinate = np.array(z_coordinate)
# xSpline = make_interp_spline(x_coordinate, xarr)
# ySpline = make_interp_spline(y_coordinate, yarr)
# zSpline = make_interp_spline(z_coordinate, zarr)
# xofx = np.linspace(x_coordinate.min(), x_coordinate.max(), 500)
# yofx = xSpline(xofx)
# xofy = np.linspace(y_coordinate.min(), y_coordinate.max(), 500)
# yofy = ySpline(xofy)
# xofz = np.linspace(z_coordinate.min(), z_coordinate.max(), 500)
# yofz = zSpline(xofz)
plt.plot(x_coordinate,xarr, label='normal to X', linewidth=3)
plt.plot(y_coordinate,yarr, label='normal to Y', linewidth=3)
plt.plot(z_coordinate,zarr, label='normal to Z', linewidth=3)
# plt.plot(xofx,yofx, label='normal to X', linewidth=3)
# plt.plot(xofy,yofy, label='normal to Y', linewidth=3)
# plt.plot(xofz,yofz, label='normal to Z', linewidth=3)
# plt.title(r'\textbf{Number of Lattice Points per Cross Section}', fontsize=40)
plt.xlabel(r'Normal Coordinate / \(Na^+\) diameters \newline', fontsize=40)
plt.ylabel(r'Number of Lattice Points', fontsize=40)
plt.xticks(fontsize=40)
plt.yticks(fontsize=40)
# plt.legend()
plt.legend(loc=2, prop={'size': 40})
plt.show()

# # Create 2D cross sections
# XY = np.ones([len(xarr),len(yarr)],dtype=int)
# YZ = np.ones([len(yarr),len(zarr)],dtype=int)
# ZX = np.ones([len(zarr),len(xarr)],dtype=int)
# for i in range (xmin2+2,xmax2-2):
#     for j in range (ymin2+2,ymax2-2):
#         if (Lattice[i][j][zcross] == 2):
#             XY[i-xmin2-2][j-ymin2-2] = 0
# for j in range (ymin2+2,ymax2-2):
#     for k in range (zmin2+2,zmax2-2):
#         if (Lattice[i][j][k] == 2):
#             YZ[j-ymin2-2][k-zmin2-2] = 0
# for k in range (zmin2+2,zmax2-2):
#     for i in range (xmin2+2,xmax2-2):
#         if (Lattice[i][j][k] == 2):
#             ZX[k-zmin2-2][i-xmin2-2] = 0
# # plt.imshow(XY, cmap='hot', interpolation='nearest')
# # plt.imshow(YZ, cmap='hot', interpolation='nearest')
# # plt.imshow(ZX, cmap='hot', interpolation='nearest')
# # plt.show()

# # # Create 2D projections
# XY = np.ones([160,165],dtype=int)
# YZ = np.ones([165,122],dtype=int)
# ZX = np.ones([122,160],dtype=int)
# for i in range (160):
#     for j in range (165):
#         for k in range (122):
#             if (Lattice[i][j][k] == 2):
#                 XY[i][j] = 0
#                 YZ[j][k] = 0
#                 ZX[k][i] = 0
# #Fill zeros inside the projections
# for i in range (1,len(XY)-1):
#     for j in range (1,len(XY[0])-1):
#         flags = np.array([XY[i-1][j], XY[i+1][j], XY[i][j-1], XY[i][j+1]], dtype = int)
#         p = flags.tolist().count(0) #counts number of pore points
#         XY[i][j] = 0 if ((XY[i][j] == 1) and (p >= 3)) else XY[i][j]
# for i in range (1,len(YZ)-1):
#     for j in range (1,len(YZ[0])-1):
#         flags = np.array([YZ[i-1][j], YZ[i+1][j], YZ[i][j-1], YZ[i][j+1]], dtype = int)
#         p = flags.tolist().count(0) #counts number of pore points
#         YZ[i][j] = 0 if ((YZ[i][j] == 1) and (p >= 3)) else YZ[i][j]
# for i in range (1,len(ZX)-1):
#     for j in range (1,len(ZX[0])-1):
#         flags = np.array([ZX[i-1][j], ZX[i+1][j], ZX[i][j-1], ZX[i][j+1]], dtype = int)
#         p = flags.tolist().count(0) #counts number of pore points
#         ZX[i][j] = 0 if ((ZX[i][j] == 1) and (p >= 3)) else ZX[i][j]
# XY = np.delete(XY, slice(105,-1), 0)
# XY = np.delete(XY, slice(0,65), 0)
# XY = np.delete(XY, slice(95,-1), 1)
# XY = np.delete(XY, slice(0,55), 1)
# YZ = np.delete(YZ, slice(100,-1), 0)
# YZ = np.delete(YZ, slice(0,60), 0)
# YZ = np.delete(YZ, slice(100,-1), 1)
# YZ = np.delete(YZ, slice(0,20), 1)
# ZX = np.delete(ZX, slice(100,-1), 0)
# ZX = np.delete(ZX, slice(0,20), 0)
# ZX = np.delete(ZX, slice(100,-1), 1)
# ZX = np.delete(ZX, slice(0,60), 1)
# # Make HeatMap of XY, YZ.T, ZX
# # XY
# plt.imshow((1-XY), cmap='hot', interpolation='nearest')
# ax = plt.gca()
# ax.invert_yaxis()
# # Major ticks at very 10%
# xXY = len(XY[0])
# yXY = len(XY)
# ax.set_xticks(np.arange(0, xXY, 5))
# ax.set_yticks(np.arange(0, yXY, 5))
# # Labels for major ticks
# ax.set_xticklabels(np.arange(0, xXY, 5))
# ax.set_yticklabels(np.arange(0, yXY, 5))
# # # Minor ticks at every 5%
# ax.set_xticks(np.arange(-.5, xXY, 1), minor=True)
# ax.set_yticks(np.arange(-.5, yXY, 1), minor=True)
# # Gridlines based on minor ticks
# ax.grid(which='minor', color='b', linestyle='-', linewidth=0.5)
# # plt.title(r'\(XY-\)Projection of the Lattice, in normalized coordinates', fontsize=25)
# plt.xlabel(r'Lattice Index in \(x\)', fontsize=25)
# plt.ylabel(r'Lattice Index in \(y\)', fontsize=25)
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# # plt.savefig('PoreProjXY1.png', dpi=1000)
# plt.show()
# # YZ.T
# plt.imshow((1-YZ.T), cmap='hot', interpolation='nearest')
# ax = plt.gca()
# ax.invert_yaxis()
# # Major ticks at very 10%
# xYZT = len(YZ.T[0])
# yYZT = len(YZ.T)
# ax.set_xticks(np.arange(0, xYZT, 5))
# ax.set_yticks(np.arange(0, yYZT, 5))
# # Labels for major ticks
# ax.set_xticklabels(np.arange(0, xYZT, 5))
# ax.set_yticklabels(np.arange(0, yYZT, 5))
# # # Minor ticks at every 5%
# ax.set_xticks(np.arange(-.5, xYZT, 1), minor=True)
# ax.set_yticks(np.arange(-.5, yYZT, 1), minor=True)
# # Gridlines based on minor ticks
# ax.grid(which='minor', color='b', linestyle='-', linewidth=0.5)
# # plt.title(r'\(ZY-\)Projection of the Lattice, in normalized coordinates', fontsize=25)
# plt.xlabel(r'Lattice Index in \(y\)', fontsize=25)
# plt.ylabel(r'Lattice Index in \(z\)', fontsize=25)
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# # plt.savefig('PoreProjZY1.png', dpi=1000)
# plt.show()
# # ZX
# plt.imshow((1-ZX), cmap='hot', interpolation='nearest')
# ax = plt.gca()
# ax.invert_yaxis()
# # Major ticks at very 5
# xZX = len(ZX[0])
# yZX = len(ZX)
# ax.set_xticks(np.arange(0, xZX, 5))
# ax.set_yticks(np.arange(0, yZX, 5))
# # Labels for major ticks
# ax.set_xticklabels(np.arange(0, xZX, 5))
# ax.set_yticklabels(np.arange(0, yZX, 5))
# # # Minor ticks at every 1
# ax.set_xticks(np.arange(-.5, xZX, 1), minor=True)
# ax.set_yticks(np.arange(-.5, yZX, 1), minor=True)
# # Gridlines based on minor ticks
# ax.grid(which='minor', color='b', linestyle='-', linewidth=0.5)
# # plt.title(r'\(ZX-\)Projection of the Lattice, in normalized coordinates', fontsize=25)
# plt.xlabel(r'Lattice Index in \(x\)', fontsize=25)
# plt.ylabel(r'Lattice Index in \(z\)', fontsize=25)
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# # plt.savefig('PoreProjZX1.png', dpi=1000)
# plt.show()