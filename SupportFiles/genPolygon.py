#! /usr/bin/python3
"""
Generate/plot a regular polygon.

If output of image is desired, include any string as a 
command line argument.

"""
import sys

import math
import matplotlib.pyplot as plt
import numpy as np

##### Polygon Parameters #####
xCenter = 0.0
yCenter = 0.0
nPts = 20
Radius = 5.0
nVert = 7
rotateAngl = 10.0  # Deg

"""
The following combinations will generate a regular polygon with
the bottom edge horizontal (i.e. "flat bottom").

Vertices Rotation Angle
   3	     90.000
   4	     45.000
   5	     18.000
   6  	      0.000
   7	     38.571
   8	     22.500
   9	     10.000
  10	      0.000
  11	     24.545
  12	     15.000
  13	      6.923
  14	      0.000
  15	     18.000
  16	     11.250
  17	      5.294
  18	      0.000
  19	     14.211
  20	      9.000
"""
##############################

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

# Make vertices
vertAngles = np.linspace(0.0, 2*np.pi, nVert+1, dtype=float)
cc = Radius*(np.cos(vertAngles))
ss = Radius*(np.sin(vertAngles))
c = cc.tolist()
s = ss.tolist()

# Rotate the vertices
vx = []
vy = []
angRotate = math.radians(rotateAngl)
for n in range(len(c)):
    dist = math.sqrt(c[n]*c[n] + s[n]*s[n])
    xyAng = math.atan2(s[n], c[n])
    ttlAng = angRotate + xyAng
    x = dist*math.cos(ttlAng)
    vx.append(x+xCenter)
    y = dist*math.sin(ttlAng)
    vy.append(y+yCenter)
vx.append(vx[0])
vy.append(vy[0])

#dataFile = 'genPolygon_vxvy.txt'
#np.savetxt(dataFile, np.array([vx, vy]).T, fmt='%+8.3f')

# Determine number of points on each edge
ptsPerEdge = int(nPts/nVert)
ttlPts = ptsPerEdge*nVert
leftover = nPts - ttlPts
nEdge = [ptsPerEdge]*(nVert)
for n in range(leftover):
    nEdge[n] = nEdge[n] + 1

# (x,y) coordinates of individual points on the polygon
xx = []
yy = []
for n in range(nVert):
    npX = np.linspace(vx[n], vx[n+1], nEdge[n]+1, dtype=float)
    npY = np.linspace(vy[n], vy[n+1], nEdge[n]+1, dtype=float)
    edgeX = npX.tolist()
    edgeY = npY.tolist()
    xx = xx + edgeX[:-1]
    yy = yy + edgeY[:-1]

print('len(xx) = ', len(xx))

print('xmin = {0:8.2f}  xmax = {1:8.2f}'.format(min(xx), max(xx)))
print('ymin = {0:8.2f}  ymax = {1:8.2f}'.format(min(yy), max(yy)))

# s is marker size
ax.scatter(xx, yy, s=20.0, marker='o', color='red')

# Conditionally output image
if outFlag == True:
    fig.savefig('genPolygon.png', facecolor='lightgrey', format='png')
    dataFile = 'genPolygon.txt'
    np.savetxt(dataFile, np.array([xx, yy]).T, fmt='%+8.3f')

plt.show()
