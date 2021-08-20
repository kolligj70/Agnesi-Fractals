#!/usr/bin/python3
"""
Generate/plot a circle with radius that varies as a sine wave.

If output of image is desired, include any string as a 
command line argument.

"""
import sys

import matplotlib.pyplot as plt
import numpy as np

#### Circle-Sine Parameters ####
xCenter = 5.0
yCenter = 2.0
Radius = 5.0
nPts = 300
freq = 11 # Number of "lumps" on the circle
amplitude = 0.75
################################

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

# Incremental angles around the circle
angle = np.linspace(0, 2*np.pi, nPts, dtype=float)

# Treat the full circumference as 1 second, i.e. basis for Hz
# Array of points in "one second" 
t = np.linspace(0, 1.0, nPts, dtype=float)

# Vary the radius according to sine wave frequency and amplitude
delRadius = amplitude*np.sin(2.0*np.pi*freq*t)

# Convert to x-y coordinates
x = (Radius+delRadius)*np.cos(angle)
x = x + xCenter
y = (Radius+delRadius)*np.sin(angle)
y = y + yCenter

print('len(x) = ', len(x))
print('xmin = {0:8.2f}  xmax = {1:8.2f}'.format(min(x), max(x)))
print('ymin = {0:8.2f}  ymax = {1:8.2f}'.format(min(y), max(y)))

# s is marker size
ax.scatter(x, y, marker='o', s=4.0, color='red')

# Conditionally output image
if outFlag == True:
    fig.savefig('genCircleSine.png', facecolor='lightgrey', format='png')

plt.show()
