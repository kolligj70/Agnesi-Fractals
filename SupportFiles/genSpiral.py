#! /usr/bin/python3
"""
Generate/plot a spiral.

If output of image is desired, include any string as a 
command line argument.

"""
import sys

import matplotlib.pyplot as plt
import numpy as np

#### Spiral Parameters ####
xCenter = 2.0
yCenter = 4.0
nPts = 300
                   # Reverse values for inward spiral
radius0 = 2.0      # Initial radius
radiusF = 5.0      # Final radius
nTurns = 3         # Number of turns
theta0_deg = 45.0  # Initial rotation angle, degrees
###########################

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

# Radius growth rate across total number of turns
b = (radiusF-radius0)/(2.0*np.pi*nTurns)

# Convert initial rotation angle to radians
theta0 = theta0_deg*np.pi/180.0
# Determine final rotation angle across total number of turns
thetaF = 2.0*np.pi*nTurns + theta0

# Array of incremental angles
t = np.arange(theta0, thetaF, (thetaF-theta0)/nPts, dtype=float)

# (x,y) coordinates of individual points on the spiral
x = (radius0+b*t)*np.cos(t)+xCenter
y = (radius0+b*t)*np.sin(t)+yCenter

print('len(x) = ', len(x))

print('xmin = {0:8.2f}  xmax = {1:8.2f}'.format(min(x), max(x)))
print('ymin = {0:8.2f}  ymax = {1:8.2f}'.format(min(y), max(y)))

# s is marker size
ax.scatter(x, y, marker='o', s=4.0, color='red')

# Conditionally output image
if outFlag == True:
    fig.savefig('genSpiral.png', facecolor='lightgrey', format='png')

plt.show()
