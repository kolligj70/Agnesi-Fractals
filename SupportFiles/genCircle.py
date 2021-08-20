#!/usr/bin/python3
"""
Generate/plot a circle.

If output of image is desired, include any string as a 
command line argument.

"""
import sys
import matplotlib.pyplot as plt
import numpy as np

#### Circle Parameters ####
xCenter = 2.0
yCenter = -4.0
nPts = 300
Radius = 5.0
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

# Generate incremental angles around full 360 deg
angle = np.linspace(0, 2*np.pi, nPts, dtype=float)

# (x,y) coordinates of individual points on the circle
x = (Radius)*np.cos(angle)
x = x + xCenter
y = (Radius)*np.sin(angle)
y = y + yCenter

print('len(x) = ', len(x))
# Output min/max values to the screen
print('tmin = {0:8.2f}  tmax = {1:8.2f}'.format(min(x), max(x)))
print('ymin = {0:8.2f}  ymax = {1:8.2f}'.format(min(y), max(y)))

# lw is linewidth
#ax.plot(x, y, lw=4, color='red')

# s is marker size
ax.scatter(x, y, s=4.0, marker='o', color='red')

# Conditionally output image
if outFlag == True:
    fig.savefig('genCircle.png', facecolor='lightgrey', format='png')

plt.show()
