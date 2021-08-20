#! /usr/bin/python3
"""
Generate/plot a polygon with radius that spirals in or out

If output of image is desired, include any string as a 
command line argument.

"""
import sys

import matplotlib.pyplot as plt
import numpy as np

#### SpiraGon Parameters ####
xCenter = 4.0
yCenter = 6.0
                   # Reverse values for inward spiral
radius0 = 2.0      # Initial radius
radiusF = 5.5      # Final radius
nVert = 7          # Number of vertices
nTurns = 3         # Number of turns
rotateAngl = 5.0   # Deg
nPts = 300
#############################

if len(sys.argv) == 2:
    outFlag = True
else:
    outFlag = False

fig = plt.figure(facecolor='lightgrey')
ax = fig.add_subplot(111)
ax.patch.set_facecolor('black')

ax.set_aspect('equal')
ax.grid(True, which='both')

# Insert X and Y axes
ax.axhline(y=0, color='w')
ax.axvline(x=0, color='w')

# Remove #'s to suppress tick labels
#ax.set_xticklabels([])
#ax.set_yticklabels([])

# Make vertices with spiraling
ttlAngle = 2.0*np.pi*nTurns
# Rate of change of the radius
b = (radiusF-radius0)/(ttlAngle)
# Initial and final angle (rad)
theta0 = rotateAngl*np.pi/180.0
thetaF = ttlAngle + theta0

ttlEdges = nTurns*nVert

# Angle of each vertex
t = np.linspace(theta0, thetaF, ttlEdges+1, dtype=float)

# x-y coordinates of vertices
vx = (radius0+b*t)*np.cos(t)+xCenter
vy = (radius0+b*t)*np.sin(t)+yCenter

# Determine number of points on each edge
ptsPerEdge = int(nPts/(ttlEdges))
ttlPts = ptsPerEdge*(ttlEdges)
leftover = nPts - ttlPts
perEdge = [ptsPerEdge]*(ttlEdges)
for n in range(leftover):
    perEdge[n] = perEdge[n] + 1

# (x,y) coordinates of individual points on the polygon
xx = []
yy = []
for n in range(len(vx)-1):
    npX = np.linspace(vx[n], vx[n+1], perEdge[n]+1, dtype=float)
    npY = np.linspace(vy[n], vy[n+1], perEdge[n]+1, dtype=float)
    edgeX = npX.tolist()
    edgeY = npY.tolist()
    xx = xx + edgeX[:-1]
    yy = yy + edgeY[:-1]

print('len(xx) = ', len(xx))

print('xmin = {0:8.2f}  xmax = {1:8.2f}'.format(min(xx), max(xx)))
print('ymin = {0:8.2f}  ymax = {1:8.2f}'.format(min(yy), max(yy)))

# s is marker size
ax.scatter(xx, yy, marker='o', s=4.0, color='red')

# Conditionally output image
if outFlag == True:
    fig.savefig('genSpiraGon.png', facecolor='lightgrey', format='png')
    dataFile = 'genSpiraGon.txt'
    np.savetxt(dataFile, np.array([xx, yy]).T, fmt='%+8.3f')

plt.show()
