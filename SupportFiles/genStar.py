#! /usr/bin/python3
"""
Generate/plot a star with outward/inward spikes.

If output of image is desired, include any string as a 
command line argument.

"""
import sys

import math
import matplotlib.pyplot as plt
import numpy as np

##### Star Parameters #####
xCenter = 2.0
yCenter = 3.0
nPts = 300
radius0 = 3.0      # Initial spike radius
radiusF = 6.0      # Final spike radius
nSpikes = 9
rotateAngl = 20.0  # Initial rotation angle, degrees
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

# Convert initial rotation angle to radians
theta0 = rotateAngl*np.pi/180.0

# Make spike vertices from initial radius
nSpikesCirc0 = np.linspace(0.0, 2*np.pi, nSpikes+1, dtype=float)
nSpikesCirc0 = nSpikesCirc0 + theta0
xx0 = radius0*(np.cos(nSpikesCirc0))
yy0 = radius0*(np.sin(nSpikesCirc0))
x0 = xx0.tolist()
y0 = yy0.tolist()

# Half the angle between spikes vertices from initial radius
halfAngle = (nSpikesCirc0[1] - nSpikesCirc0[0])/2.0

# Make spike vertices from final radius
nSpikesCircF = nSpikesCirc0 + halfAngle
xxF = radiusF*(np.cos(nSpikesCircF))
yyF = radiusF*(np.sin(nSpikesCircF))
xF = xxF.tolist()
yF = yyF.tolist()

vx = []
vy = []
# Combine the two sets of vertices
for n in range(len(x0)):
    vx.append(x0[n] + xCenter)
    vx.append(xF[n] + xCenter)
    vy.append(y0[n] + yCenter)
    vy.append(yF[n] + yCenter)

# Determine number of points on each spike edge
nSpikes2 = nSpikes*2
ptsPerEdge = int(nPts/(nSpikes2))
ttlPts = ptsPerEdge*nSpikes2
leftover = nPts - ttlPts
nEdge = [ptsPerEdge]*(nSpikes2)
for n in range(leftover):
    nEdge[n] = nEdge[n] + 1

# (x,y) coordinates of individual points on each edge
xx = []
yy = []
for n in range(nSpikes2):
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
ax.scatter(xx, yy, s=10.0, marker='o', color='red')

# Conditionally output image and (x,y) data
if outFlag == True:
    fig.savefig('genStar.png', facecolor='lightgrey', format='png')
    dataFile = 'genStar.txt'
    np.savetxt(dataFile, np.array([xx, yy]).T, fmt='%+8.3f')

plt.show()
