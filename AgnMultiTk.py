#!/usr/bin/python3
"""
Witch of Agnesi Images with Variations - See the full paper referenced in
the Menu->Overview for a detailed description of the original algorithm.

This program provides various modifications to the referenced algorithm,
including: choice of several Base Figures besides a circle; additional
algorithms for the iterative calculation of (x,y) points; several methods
for plotting subsets of the full data set; marker choice, along with
marker size setting, and the ability to taper the marker size; choice of
either standard or custom colormaps.

"""
import os
import sys
import string
import math
import datetime
from collections import OrderedDict

import tkinter as tk
from tkinter import filedialog
from tkinter import ttk

from PIL import ImageTk

import matplotlib as mpl
import matplotlib.colors as mc

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)

import numpy as np

import AgnMultiConfig as xfig

cmFilePath = ''

# Debug Dictionary
# debug1 - Enable plot window toolbar
# debug2 - Output marker tapering array to file
# debug3 - Output BaseFigure settings, Algorithm selection at start of "Plot"
# debug4 - Output BaseFigure x-y data to file
# debug5 - Print Color List Colormap text file used for plot
# debug6 - Print the converted value from the ckForFract function
debugDict = {"debug1": False,
             "debug2": False,
             "debug3": False,
             "debug4": False,
             "debug5": False,
             "debug6": False}


# Binding function for the Iterative Algorithm tab selection
def on_tabAlg_selected(event):
    selected_tab = event.widget.select()
    tab_text = event.widget.tab(selected_tab, "text")
    varAlgorithm.set(tab_text)


# Binding function for the Base Figure tab selection
def on_tabBaseFig_selected(event):
    selected_tab = event.widget.select()
    tab_text = event.widget.tab(selected_tab, "text")
    varBaseFig.set(tab_text)


# Iterative Algorithm calculation
def functX(a, x):
    alg = varAlgorithm.get()
    if alg == 'Version1':
        x2 = x*x
        num = 1 + a*(x + x2 + x*x2)
        den = 1.0 + x2
    elif alg == 'Version2':
        ax = a*x
        num = 1.0 + ax + ax*ax + ax*ax*ax
        den = 1.0 + x*x
    else:
        ax = a*x
        num = 1.0 - ax + ax*ax - ax*ax*ax
        den = 1.0 + x*x
    fx = num/den
    return fx


# Update the x-coordinate
def updateX(x, y, a, b):
    fx = functX(a, x)
    xnext = b*y + fx
    return xnext


# Update the y-coordinate
def updateY(xn, x, a):
    ynext = -x + functX(a, xn)
    return ynext


# Make the Base Figure: Circle
def makeCircle(xCenter, yCenter, Radius, N2make):
    xx = []
    yy = []
    for i in range(N2make):
        aRad = 2.0*np.pi*i/N2make   # Angle in radians
        x = Radius*math.cos(aRad) + xCenter
        xx.append(x)
        y = Radius*math.sin(aRad) + yCenter
        yy.append(y)
    return(xx, yy)


# Make the Base Figure: CircleSine
def makeCircleSine(xCenter, yCenter, Radius, N2make, Cycles, Amplitude):
    # Angle around the circle
    angle = np.linspace(0, 2*np.pi, N2make)
    # Generate sine wave as a function of frequency and amplitude
    # aka the radius variation
    t = np.linspace(0, 1.0, N2make)
    delRadius = Amplitude*np.sin(2.0*np.pi*Cycles*t)

    x = (Radius+delRadius) * np.cos(angle)
    x = x + xCenter
    y = (Radius+delRadius) * np.sin(angle)
    y = y + yCenter

    xx = x.tolist()
    yy = y.tolist()

    return(xx, yy)


# Make the Base Figure: Star
def makeStar(xCenter, yCenter, Radius0, RadiusF, nSpikes,
             degRotate, nPts):
    # Convert initial rotation angle to radians
    theta0 = degRotate*np.pi/180.0

    # Make Star vertices from initial radius
    nSpikesCirc0 = np.linspace(0.0, 2*np.pi, nSpikes+1, dtype=float)
    nSpikesCirc0 = nSpikesCirc0 + theta0
    xx0 = Radius0*(np.cos(nSpikesCirc0))
    yy0 = Radius0*(np.sin(nSpikesCirc0))
    x0 = xx0.tolist()
    y0 = yy0.tolist()

    # Half the angle between spikes vertices using initial radius
    halfAngle = (nSpikesCirc0[1] - nSpikesCirc0[0])/2.0

    # Make spike vertices from final radius
    nSpikesCircF = nSpikesCirc0 + halfAngle
    xxF = RadiusF*(np.cos(nSpikesCircF))
    yyF = RadiusF*(np.sin(nSpikesCircF))
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

    return(xx, yy)


# Make the Base Figure: Spiral
def makeSpiral(xCenter, yCenter, Inner, Outer, NTurns, degRotate, N2make):
    xx = []
    yy = []

    ttlAngle = 2.0*np.pi*NTurns
    # Radius incremental change
    b = (Outer-Inner)/(ttlAngle)
    theta0 = degRotate*np.pi/180.0
    thetaF = ttlAngle + theta0

    t = theta0
    tDel = (thetaF - theta0)/N2make
    for i in range(N2make):
        x = (Inner+b*t)*np.cos(t)+xCenter
        xx.append(x)
        y = (Inner+b*t)*np.sin(t)+yCenter
        yy.append(y)
        t = t + tDel

    return(xx, yy)


# Make the Base Figure: Polygon
def makePolygon(xCenter, yCenter, Radius, nVert, degRotate, N2make):
    vx = []
    vy = []

    # Make vertices
    vertAngles = np.linspace(0.0, 2*np.pi, nVert+1, dtype=float)
    cc = Radius*(np.cos(vertAngles))
    ss = Radius*(np.sin(vertAngles))
    c = cc.tolist()
    s = ss.tolist()

    if degRotate != 0.0:
        angRotate = math.radians(degRotate)
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
    else:
        for n in range(nVert):
            vx.append(c[n]+xCenter)
            vy.append(s[n]+yCenter)

        vx.append(c[0]+xCenter)
        vy.append(s[0]+yCenter)

    # Populate edges (i.e. lines between vertices)
    ptsPerEdge = int(N2make/nVert)
    ttlPts = ptsPerEdge*nVert
    leftover = N2make - ttlPts
    nEdge = [ptsPerEdge]*(nVert)
    for n in range(leftover):
        nEdge[n] = nEdge[n] + 1

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

    return(xx, yy)


# Make Base Figure: SpiraGon
def makeSpiraGon(xCenter, yCenter, Inner, Outer, nTurns, nVerts,
                 degRotate, nPts):

    # Make vertices with spiraling
    ttlAngle = 2.0*np.pi*nTurns
    # Rate of change of the radius
    b = (Outer-Inner)/(ttlAngle)
    theta0 = degRotate*np.pi/180.0
    thetaF = ttlAngle + theta0

    ttlEdges = nTurns*nVerts

    # Angle of each vertex
    t = np.linspace(theta0, thetaF, ttlEdges+1, dtype=float)
    # x-y coordinates of vertices
    vx = (Inner+b*t)*np.cos(t)+xCenter
    vy = (Inner+b*t)*np.sin(t)+yCenter

    # Determine number of points on each edge
    ptsPerEdge = int(nPts/ttlEdges)
    ttlPts = ptsPerEdge*ttlEdges
    leftover = nPts - ttlPts
    perEdge = [ptsPerEdge]*(ttlEdges)
    for n in range(leftover):
        perEdge[n] = perEdge[n] + 1

    # Populate edges (i.e. lines between vertices)
    xx = []
    yy = []
    for n in range(len(vx)-1):
        npX = np.linspace(vx[n], vx[n+1], perEdge[n]+1, dtype=float)
        npY = np.linspace(vy[n], vy[n+1], perEdge[n]+1, dtype=float)
        edgeX = npX.tolist()
        edgeY = npY.tolist()
        xx = xx + edgeX[:-1]
        yy = yy + edgeY[:-1]

    return(xx, yy)


# Calculate the Orbit for each point of the Base Figure
def makeOrbit(x0, y0, a, b, N):
    xx = [x0]
    yy = [y0]

    for n in range(N-1):
        x = xx[-1]
        y = yy[-1]
        xn = updateX(x, y, a, b)
        yn = updateY(xn, x, a)
        xx.append(xn)
        yy.append(yn)
    return(xx, yy)


# Create plot of data
def makePlot():
    global cmFilePath

    def ckNumSpikes(N):
        if N < 3:
            print('Erroneous number of spikes. Defaulting to 3')
            N = 3
        return(N)

    def ckNumVertices(N):
        if N < 3:
            print('Erroneous number of vertices. Defaulting to 3')
            N = 3
        return(N)

    def ckForFract(strVal):
        # Check for fraction or decimal entry
        cnt = strVal.count('/')
        if cnt == 0:
            fract = float(strVal)
        elif cnt == 1:
            strVal.strip()
            nd = strVal.split('/')
            num = float(nd[0])
            den = float(nd[1])
            fract = num/den
        else:
            print('Erroneous entry box value. Returning')
            return
        strStr = '{0:<8.5e}\n'.format(fract)
        if debugDict["debug6"]:
            print(strVal, strStr)
        return(strStr, fract)

    def SaveDat():
        # Plot file names = concatenation of parts of: base name
        # and datTime.
        picFile = xfig.baseName + datTime + '.png'
        FigureCanvasTkAgg.print_png(canvas, picFile)

        algVersStr = 'algVers = {0:s}\n'.format(varAlgorithm.get())
        boolVal = saveRaw.get()
        if boolVal is True:
            strVal = 'True'
        else:
            strVal = 'False'

        saveRawStr = 'saveRaw = {0:s}\n'.format(strVal)

        fullText = baseFigText + aParamStr + bParamStr + itersPerOrbitStr + \
            algVersStr + pltStartStr + pltStopStr + pltStrideStr + \
            markerStr + markerSizeStr + markerTaperStr + \
            cMapTypeStr + cMapStrOut + \
            figScaleStr + saveRawStr

        infoFile = xfig.baseName + datTime + '.agd'
        # Write plot information to file
        with open(infoFile, "w") as outfile:
            outfile.write(fullText)

        # Conditionally save x-y data
        if boolVal is True:
            dataFile = 'rawAgnDat' + datTime + '.txt'
            np.savetxt(dataFile, np.array([xPlot, yPlot]).T, fmt='%+2.12e')

    # Extract GUI data
    strVal = aParamSV.get()
    strStr, aParam = ckForFract(strVal)
    aParamStr = 'aParam = ' + strStr
    strVal = bParamSV.get()
    strStr, bParam = ckForFract(strVal)
    bParamStr = 'bParam = ' + strStr

    itersPerOrbit = int(nOrbitPtsSV.get())
    itersPerOrbitStr = 'itersPerOrbit = {0:<8d}\n'.format(itersPerOrbit)

    pltStart = int(pltStartSV.get())
    pltStartStr = 'pltStart = {0:<8d}\n'.format(pltStart)
    pltStop = int(pltStopSV.get())
    pltStopStr = 'pltStop = {0:<8d}\n'.format(pltStop)
    pltStride = int(pltStrideSV.get())
    pltStrideStr = 'pltStride = {0:<8d}\n'.format(pltStride)

    marker = MarkerSV.get()
    marker = marker.strip()
    markerStr = 'marker = {0:s}\n'.format(marker)

    markerSize = float(MarkerSizeSV.get())
    markerSizeStr = 'markerSize = {0:<8.2f}\n'.format(markerSize)

    markerTaper = MarkerTaperSV.get()
    markerTaperStr = 'markerTaper = {0:s}\n'.format(markerTaper)

    baseFig = varBaseFig.get()
    print('baseFig = ', baseFig)
    baseFigStr = baseFig + '\n'
    baseFigStr = 'baseFig = {0:s}\n'.format(baseFig)

    # Capture date/time that plot processing started
    datTime = datetime.datetime.now()
    datTime = datTime.strftime("%m%d%y_%H%M%S")
    title = 'Agn' + baseFig + '_' + datTime

    # Create Base Figure: Circle
    if baseFig == 'Circle':
        xCenter = float(xCenterCircleSV.get())
        xCenterStr = 'xCenter = {0:<8.3f}\n'.format(xCenter)
        yCenter = float(yCenterCircleSV.get())
        yCenterStr = 'yCenter = {0:<8.3f}\n'.format(yCenter)
        Radius = float(RadiusCircleSV.get())
        RadiusStr = 'Radius = {0:<8.3f}\n'.format(Radius)
        nPts = int(nPtsCircleSV.get())
        nPtsStr = 'nPts = {0:<6d}\n'.format(nPts)
        baseFigText = baseFigStr + xCenterStr + yCenterStr + \
            RadiusStr + nPtsStr

        # Make Circle: xP0, yP0 are coordinates of circle samples
        xP0, yP0 = makeCircle(xCenter, yCenter, Radius, nPts)

    # Create Base Figure: CircleSine
    elif baseFig == 'CircleSine':
        xCenter = float(xCenterCircleSineSV.get())
        xCenterStr = 'xCenter = {0:<8.3f}\n'.format(xCenter)
        yCenter = float(yCenterCircleSineSV.get())
        yCenterStr = 'yCenter = {0:<8.3f}\n'.format(yCenter)
        Radius = float(RadiusCircleSineSV.get())
        RadiusStr = 'Radius = {0:<8.3f}\n'.format(Radius)
        nPts = int(nPtsCircleSineSV.get())
        nPtsStr = 'nPts = {0:<6d}\n'.format(nPts)
        Cycles = float(CyclesCircleSineSV.get())
        CyclesStr = 'Cycles = {0:<8.3f}\n'.format(Cycles)
        Amplitude = float(AmplitudeCircleSineSV.get())
        AmplitudeStr = 'Amplitude = {0:<8.3f}\n'.format(Amplitude)
        baseFigText = baseFigStr + xCenterStr + yCenterStr + RadiusStr + \
            nPtsStr + CyclesStr + AmplitudeStr

        # Make CircleSine: xP0, yP0 are coordinates of circle samples
        xP0, yP0 = makeCircleSine(
                        xCenter, yCenter, Radius, nPts, Cycles, Amplitude)

    # Create Base Figure: Star
    elif baseFig == 'Star':
        xCenter = float(xCenterStarSV.get())
        xCenterStr = 'xCenter = {0:<8.3f}\n'.format(xCenter)
        yCenter = float(yCenterStarSV.get())
        yCenterStr = 'yCenter = {0:<8.3f}\n'.format(yCenter)
        StartRadius = float(StartRadiusStarSV.get())
        StartRadiusStr = 'StartRadius = {0:<8.3f}\n'.format(StartRadius)
        EndRadius = float(EndRadiusStarSV.get())
        EndRadiusStr = 'EndRadius = {0:<8.3f}\n'.format(EndRadius)
        nSpikes = int(nSpikesStarSV.get())
        nSpikes = ckNumSpikes(nSpikes)
        nSpikesStr = 'nSpikes = {0:<3d}\n'.format(nSpikes)
        initialRotation = float(RotateStarSV.get())
        initialRotationStr = \
            'initialRotation = {0:<8.3f}\n'.format(initialRotation)
        nPts = int(nPtsStarSV.get())
        nPtsStr = 'nPts = {0:<6d}\n'.format(nPts)
        baseFigText = baseFigStr + xCenterStr + yCenterStr + \
                      StartRadiusStr + EndRadiusStr + \
                      nSpikesStr + initialRotationStr + nPtsStr

        # Make Star: xP0, yP0 are coordinates of circle samples
        xP0, yP0 = makeStar(xCenter, yCenter, StartRadius, EndRadius,
                            nSpikes, initialRotation, nPts)

    # Create Base Figure: Spiral
    elif baseFig == 'Spiral':
        xCenter = float(xCenterSpiralSV.get())
        xCenterStr = 'xCenter = {0:<8.3f}\n'.format(xCenter)
        yCenter = float(yCenterSpiralSV.get())
        yCenterStr = 'yCenter = {0:<8.3f}\n'.format(yCenter)
        StartRadius = float(StartRadiusSpiralSV.get())
        StartRadiusStr = 'StartRadius = {0:<8.3f}\n'.format(StartRadius)
        EndRadius = float(EndRadiusSpiralSV.get())
        EndRadiusStr = 'EndRadius = {0:<8.3f}\n'.format(EndRadius)
        nTurns = int(nTurnsSpiralSV.get())
        nTurnsStr = 'nTurns = {0:<3d}\n'.format(nTurns)
        initialRotation = float(RotateSpiralSV.get())
        initialRotationStr = \
            'initialRotation = {0:<8.3f}\n'.format(initialRotation)
        nPts = int(nPtsSpiralSV.get())
        nPtsStr = 'nPts = {0:<6d}\n'.format(nPts)

        baseFigText = baseFigStr + xCenterStr + yCenterStr + StartRadiusStr + \
            EndRadiusStr + nTurnsStr + initialRotationStr + nPtsStr

        # Make Spiral: xP0, yP0 are coordinates of spiral samples
        xP0, yP0 = makeSpiral(xCenter, yCenter, StartRadius, EndRadius,
                              nTurns, initialRotation, nPts)

    # Create Base Figure: Polygon
    elif baseFig == 'Polygon':
        xCenter = float(xCenterPolygonSV.get())
        xCenterStr = 'xCenter = {0:<8.3f}\n'.format(xCenter)
        yCenter = float(yCenterPolygonSV.get())
        yCenterStr = 'yCenter = {0:<8.3f}\n'.format(yCenter)
        Radius = float(RadiusPolygonSV.get())
        RadiusStr = 'Radius = {0:<8.3f}\n'.format(Radius)
        nVert = int(nVertPolygonSV.get())
        nVert = ckNumVertices(nVert)
        nVertStr = 'nVert = {0:<3d}\n'.format(nVert)
        initialRotation = float(RotatePolygonSV.get())
        initialRotationStr = \
            'initialRotation = {0:<8.3f}\n'.format(initialRotation)
        nPts = int(nPtsPolygonSV.get())
        nPtsStr = 'nPts = {0:<6d}\n'.format(nPts)

        baseFigText = baseFigStr + xCenterStr + yCenterStr + RadiusStr + \
            nVertStr + initialRotationStr + nPtsStr

        # Make Polygon: xP0, yP0 are coordinates of polygon samples
        xP0, yP0 = makePolygon(xCenter, yCenter, Radius, nVert,
                               initialRotation, nPts)

    # Create Base Figure: SpiraGon
    elif baseFig == 'SpiraGon':
        xCenter = float(xCenterSpiraGonSV.get())
        xCenterStr = 'xCenter = {0:<8.3f}\n'.format(xCenter)
        yCenter = float(yCenterSpiraGonSV.get())
        yCenterStr = 'yCenter = {0:<8.3f}\n'.format(yCenter)
        StartRadius = float(StartRadiusSpiraGonSV.get())
        StartRadiusStr = 'StartRadius = {0:<8.3f}\n'.format(StartRadius)
        EndRadius = float(EndRadiusSpiraGonSV.get())
        EndRadiusStr = 'EndRadius = {0:<8.3f}\n'.format(EndRadius)
        nTurns = int(nTurnsSpiraGonSV.get())
        nTurnsStr = 'nTurns = {0:<3d}\n'.format(nTurns)
        nVert = int(nVertSpiraGonSV.get())
        nVert = ckNumVertices(nVert)
        nVertStr = 'nVert = {0:<3d}\n'.format(nVert)
        initialRotation = float(RotateSpiraGonSV.get())
        initialRotationStr = \
            'initialRotation = {0:<8.3f}\n'.format(initialRotation)
        nPts = int(nPtsSpiraGonSV.get())
        nPtsStr = 'nPts = {0:<6d}\n'.format(nPts)

        baseFigText = baseFigStr + xCenterStr + yCenterStr + StartRadiusStr + \
            EndRadiusStr + nTurnsStr + nVertStr + initialRotationStr + \
            nPtsStr

        # Make SpiraGon: xP0, yP0 are coordinates of SpiraGon samples
        xP0, yP0 = makeSpiraGon(xCenter, yCenter, StartRadius, EndRadius,
                                nTurns, nVert, initialRotation, nPts)

    # Else default to Circle on error
    else:
        print('Unknown Base Figure. Defaulting to Circle')
        baseFig = 'Circle'
        xCenter = 0.0
        xCenterStr = 'xCenter = {0:<8.3f}\n'.format(xCenter)
        yCenter = 0.0
        yCenterStr = 'xCenter = {0:<8.3f}\n'.format(yCenter)
        Radius = 5.0
        RadiusStr = 'Radius = {0:<8.3f}\n'.format(Radius)
        nPts = 360
        nPtsStr = 'nPts = {0:<6d}\n'.format(nPts)

        # Make circle: xP0, yP0 are coordinates of circle samples
        xP0, yP0 = makeCircle(xCenter, yCenter, Radius, nPts)

    if debugDict["debug3"]:
        alg = varAlgorithm.get()
        print('\n' + alg)
        print(baseFigText)

    if debugDict["debug4"]:
        np.savetxt('baseFig.txt', np.array([xP0, yP0]).T, fmt='%+2.12e')

    # Make orbit for each point on the Base Figure
    xCalc = []
    yCalc = []
    for N in range(len(xP0)):
        xxx, yyy = makeOrbit(xP0[N], yP0[N], aParam, bParam, itersPerOrbit)
        for i in range(len(xxx)):
            xCalc.append(xxx[i])
            yCalc.append(yyy[i])

    # Apply start/stop restrictions to plot data
    xPlot = xCalc[pltStart:pltStop+1]
    yPlot = yCalc[pltStart:pltStop+1]

    # Apply stride to data
    if pltStride < 1:
        print('Erroneous Plot Stride value. Defaulting to 1')
        pltStride = 1
    xPlot = xPlot[0::pltStride]
    yPlot = yPlot[0::pltStride]

    # Find minimum and maximum x,y data
    xmin = min(xPlot)
    xmax = max(xPlot)
    ymin = min(yPlot)
    ymax = max(yPlot)

    # Create plot window
    window = tk.Toplevel()
    window.geometry('+200+30')

    tk.Label(window, text=title).grid(row=0, column=0, sticky=tk.EW)

    # Set plot window extents according to data ranges and
    # extra based on configuration data
    xExtent = xmax-xmin
    yExtent = ymax-ymin
    xWindowXtra = (xfig.plotExpander)*xExtent
    yWindowXtra = (xfig.plotExpander)*yExtent

    figScale = varFigureScaling.get()
    figScaleStr = 'figScale = {0:s}\n'.format(figScale)

    f = Figure(figsize=(xfig.windowWidth, xfig.windowHeight), dpi=xfig.dpi)
    figPlot = f.add_subplot(111)
    if (figScale == 'EQUAL'):
        # Set equal x and y axes
        figPlot.set_aspect('equal', 'box')
    else:
        # Set x-axis and y-axis according to x&y min/max data
        figPlot.set_xlim([xmin-xWindowXtra, xmax+xWindowXtra])
        figPlot.set_ylim([ymin-yWindowXtra, ymax+yWindowXtra])

    # Set plot facecolor
    figPlot.patch.set_facecolor(xfig.faceColor)

    # Minimize the border surrounding the plot area
    f.set_tight_layout(True)

    # Turn off grid labels
    figPlot.xaxis.set_visible(False)
    figPlot.yaxis.set_visible(False)

    # Output min/max x and y values to screen
    totalCalculations = nPts*itersPerOrbit
    print('Total Calculations = {0:<8d}'.format(totalCalculations))
    print('xmin = {0:8.1f}  xmax = {1:8.1f}'.format(xmin, xmax))
    print('ymin = {0:8.1f}  ymax = {1:8.1f}'.format(ymin, ymax))

    # Create data array for taper scaling
    xLen = len(xPlot)
    xxt = np.arange(xLen)

    # Set array of marker tapering scaling
    if markerTaper == "None":
        scaler = markerSize
    elif markerTaper == "Zero2Full":
        scaler = markerSize*(xxt/xLen)
    elif markerTaper == "Full2Zero":
        scaler = np.flipud(markerSize*(xxt/xLen))
    else:
        print("Erroneous markerScale. Defaulting to no scaling")
        scaler = markerSize

    if debugDict["debug2"]:
        np.savetxt('taper.txt', np.array([scaler]).T, fmt='%+10.5f')

    # Set up specified CMAP
    cMapType = varcMapType.get()
    cMapTypeStr = 'cMapType = {0:s}\n'.format(cMapType)
    if (cMapType == "STD"):
        # Here to apply Standard colormap to plot
        cMap = cMapStd.get()
        cMap = cMap.strip()
        cMapStrOut = 'cMapStd = {0:s}\n'.format(cMap)
        figPlot.scatter(xPlot, yPlot, s=scaler, c=xxt, cmap=cMap,
                        marker=marker, edgecolor='black')
    elif (cMapType == "CUST"):
        cMapListStr = cMapList.get().strip()
        if (debugDict["debug5"]):
            print("LIST colorFile ", cMapListStr)
        if cmFilePath == '':
            cmFilePath = './' + cMapListStr
        # Here to generate Color List colormap and apply to plot
        # Get Color List custom colormap file
        with open(cmFilePath, "r") as inFile:
            lines = []
            for line in inFile:
                if (debugDict["debug5"]):
                    # Suppress the extra newlines
                    print(line.rstrip())
                # Discard lines starting with #, or zero-length lines
                if (line.startswith('#')) or (len(line.split()) == 0):
                    continue
                else:
                    # Strip trailing newlines
                    line = line.rstrip("\n")
                    lines.append(line)

        for line in lines:
            lineSpl = line.split('=')
            # Parse lines with edge data
            if (lineSpl[0].strip() == 'edgeList'):
                edgeListStr = lineSpl[1].strip()
                edgeListStr = \
                    edgeListStr.strip('][').replace(",", "").split()
                edgeList = [float(x) for x in edgeListStr]
            # Parse lines with colors (names of hex codes)
            elif (lineSpl[0].strip() == 'rgbList'):
                rgbListStr = lineSpl[1].strip()
                rgbListStr = rgbListStr.strip('][').replace(",", "").split()
            else:
                continue

        # Check Color List for proper hex format, and named colors
        rgbOK = ckColorList(rgbListStr)
        if rgbOK is False:
            print("makePlot: Color List error")
            return
        # Check Edge List initial/final values and increasing values
        edgeOK = ckEdgeVals(edgeList)
        if edgeOK is False:
            print("makePlot: Edge List error")
            return

        # Adapt Edge List to data array
        bounds = [i*(len(xPlot)) for i in edgeList]
        # Generate colormap and plot
        cmap = mpl.colors.ListedColormap(rgbListStr)
        cMapStrOut = 'cMapList = {0:s}\n'.format(cMapListStr)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N, clip=True)
        figPlot.scatter(xPlot, yPlot, s=scaler, c=xxt, cmap=cmap,
                        marker=marker, edgecolor='black', norm=norm)
    else:
        print("makePlot: unknown colormap type")
        return

    # Create canvas to hold plot
    canvas = FigureCanvasTkAgg(f, master=window)  # A tk.DrawingArea.
    canvas._tkcanvas.grid(row=1, column=0, sticky=tk.N)
    canvas.draw()

    # Define plot window Save button
    btnSave = tk.Button(master=window, text='Save', command=SaveDat)
    btnSave.grid(row=2, column=1, sticky=tk.E)

    if debugDict["debug1"]:
        # Enable plot window toolbar for unit testing. Toolbar
        # data is displayed in lower-left corner of main GUI window.
        toolbarFrame = tk.Frame(master=root)
        toolbarFrame.grid(row=2, column=0)
        toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)


# Check Color List for proper hex, and recognized named colors
def ckColorList(lst):
    # Check for digit is allowable hex
    def isHex(N):
        hexDigits = string.hexdigits
        OK = True
        for d in N:
            if d not in hexDigits:
                OK = False
                break
        return (OK)

    # Get known colors dictionary
    namedColorsDict = mc.cnames
    namedColorsList = []
    # Extract known color names from dictionary
    for k in namedColorsDict.keys():
        namedColorsList.append(k)
        # Sort the color list
        namedColorsList.sort(key=str.lower)

    # Check each digit is hex and correct number of digits
    OK = True
    for N in lst:
        # Check hex value
        if N[0] == "#":
            N = N.lstrip("#")
            if ((len(N)) != 6) or (isHex(N) is False):
                msg = 'Erroneous hex number'
                print(msg)
                OK = False
        else:
            # Check named color against known colors list
            if N not in namedColorsList:
                msg = "Named color not in referenced list"
                print(msg)
                OK = False
    return(OK)


# Color List colormap - check edge value list for mandatory initial/final
# values and increasing values across whole list
def ckEdgeVals(lst):
    # Check list for monotonically increasing values
    def ckIncVals(lst):
        fltList = [float(x) for x in lst]
        OK = np.all(np.diff(fltList) > 0.0)
        return (OK)

    # Check for "equals"
    def ckClose(a, b):
        if (abs(a - b) < 1e-9):
            OK = True
        else:
            OK = False
        return (OK)

    # Check first/last edge values
    if (ckClose(float(lst[0]), 0.0) is False) or \
       (ckClose(float(lst[-1]), 1.0) is False) or \
       (ckIncVals(lst) is False):
        msg = 'First/last edge value error OR values not increasing'
        print(msg)
        OK = False
    else:
        OK = True
    return (OK)


# Make initial GUI form
def makeForm(root):
    global CustomCmapEntry
    global tabCtrl, tabAlg1, tabAlg2, tabAlg3
    global baseFigTabCtrl
    global BaseFigCircleFrame, BaseFigCircleSineFrame, BaseFigSpiralFrame
    global BaseFigPolygonFrame, BaseFigSpiraGonFrame, BaseFigStarFrame
    # Base Figure - Circle: Entries
    global xCenterCircleEntry, yCenterCircleEntry
    global RadiusCircleEntry, nPtsCircleEntry
    # Base Figure - CircleSine: Entries
    global xCenterCircleSineEntry, yCenterCircleSineEntry
    global RadiusCircleSineEntry, nPtsCircleSineEntry
    global CyclesCircleSineEntry, AmplitudeCircleSineEntry
    # Base Figure - Star: Entries
    global xCenterStarEntry, yCenterStarEntry
    global StartRadiusStarEntry, EndRadiusStarEntry
    global nSpikesStarEntry, RotateStarEntry, nPtsStarEntry
    # Base Figure - Spiral: Entries
    global xCenterSpiralEntry, yCenterSpiralEntry
    global StartRadiusSpiralEntry, EndRadiusSpiralEntry
    global nTurnsSpiralEntry, RotateSpiralEntry, nPtsSpiralEntry
    # Base Figure - Polygon: Entries
    global xCenterPolygonEntry, yCenterPolygonEntry, RadiusPolygonEntry
    global nVertPolygonEntry, RotatePolygonEntry, nPtsPolygonEntry
    # Base Figure - SpiraGon: Entries
    global xCenterSpiraGonEntry, yCenterSpiraGonEntry
    global StartRadiusSpiraGonEntry, EndRadiusSpiraGonEntry
    global nVertSpiraGonEntry, RotateSpiraGonEntry
    global nTurnsSpiraGonEntry, nPtsSpiraGonEntry
    # Iterative Calculation: Entries
    global aParamEntry, bParamEntry, nOrbitPtsEntry
    # Screen Data: Entries
    global pltStartEntry, pltStopEntry, pltStrideEntry
    global MarkerEntry, MarkerSizeEntry
    global varcMapType, StdCmapEntry, CustomCmapEntry
    global algTabCtrl, varAlgorithm

    # Function to create GUI Label/Value for specified entry
    def makeLine(dictList, txtvar, frame):
        lbl = dictList[0]
        val = dictList[1]
        row = dictList[2]
        lineLabel = tk.Label(frame,
                             width=xfig.widthRowLbl,
                             text=lbl,
                             font=xfig.fontLbl,
                             anchor='w')
        lineEntry = tk.Entry(frame,
                             textvariable=txtvar,
                             font=xfig.fontEnt,
                             width=xfig.widthRowEnt,
                             justif='left')
        return(lineLabel, lineEntry, val, row)

    ipadx = 2
    ipady = 2
    padx = 1
    pady = 2
    hgt = 660

    # Create frame for data entry
    dataFrame = tk.Frame(root, bg=xfig.backgndColor)
    dataFrame.grid(row=1, column=0,
                   sticky=tk.N, columnspan=2)

    # Section Label
    txt = '*Base Figure*'
    lbl = tk.Label(dataFrame, width=xfig.widthSectionTitle, text=txt)
    lbl.config(font=xfig.fontSectionTitle)
    lbl.grid(row=1, column=0, ipadx=ipadx, ipady=ipady,
             sticky=tk.N, columnspan=2)

    # Create tabbed window manager
    baseFigTabCtrl = ttk.Notebook(dataFrame)
    baseFigTabCtrl.grid(row=2, column=0, ipadx=ipadx, ipady=ipady,
                        sticky=tk.EW, columnspan=2)
    baseFigTabCtrl.bind("<<NotebookTabChanged>>", on_tabBaseFig_selected)

    # Circle
    # Tabbed Window: Base Fig Circle
    BaseFigCircleFrame = ttk.Frame(baseFigTabCtrl)
    baseFigTabCtrl.add(BaseFigCircleFrame, text='Circle')

    # X Coordinate - Circle Center
    dictList = xfig.datCircle['xCenterCircle']
    xCenterCircleLabel, xCenterCircleEntry, val, xCenterCircleRow = \
        makeLine(dictList, xCenterCircleSV, BaseFigCircleFrame)
    xCenterCircleSV.set(val)

    # Y Coordinate - Circle Center
    dictList = xfig.datCircle['yCenterCircle']
    yCenterCircleLabel, yCenterCircleEntry, val, yCenterCircleRow = \
        makeLine(dictList, yCenterCircleSV, BaseFigCircleFrame)
    yCenterCircleSV.set(val)

    # Circle Radius
    dictList = xfig.datCircle['Radius']
    RadiusCircleLabel, RadiusCircleEntry, val, radiusCircleRow = \
        makeLine(dictList, RadiusCircleSV, BaseFigCircleFrame)
    RadiusCircleSV.set(val)

    # Number of Points on Circle
    dictList = xfig.datCircle['nPts']
    nPtsCircleLabel, nPtsCircleEntry, val, nptsCircleRow = \
        makeLine(dictList, nPtsCircleSV, BaseFigCircleFrame)
    nPtsCircleSV.set(val)

    xCenterCircleLabel.grid(row=xCenterCircleRow, column=0, sticky=tk.W)
    xCenterCircleEntry.grid(row=xCenterCircleRow, column=1)
    yCenterCircleLabel.grid(row=yCenterCircleRow, column=0, sticky=tk.W)
    yCenterCircleEntry.grid(row=yCenterCircleRow, column=1)
    RadiusCircleLabel.grid(row=radiusCircleRow, column=0, sticky=tk.W)
    RadiusCircleEntry.grid(row=radiusCircleRow, column=1)
    nPtsCircleLabel.grid(row=nptsCircleRow, column=0, sticky=tk.W)
    nPtsCircleEntry.grid(row=nptsCircleRow, column=1)

    # CircleSine
    # Tabbed Window: Base Fig CircleSine
    BaseFigCircleSineFrame = ttk.Frame(baseFigTabCtrl)
    baseFigTabCtrl.add(BaseFigCircleSineFrame, text='CircleSine')

    # X Coordinate - CircleSine Center
    dictList = xfig.datCircleSine['xCenterCircleSine']
    xCenterCircleSineLabel, xCenterCircleSineEntry, val, xCenterCircSinRow = \
        makeLine(dictList, xCenterCircleSineSV, BaseFigCircleSineFrame)
    xCenterCircleSineSV.set(val)

    # Y Coordinate - CircleSine Center
    dictList = xfig.datCircleSine['yCenterCircleSine']
    yCenterCircleSineLabel, yCenterCircleSineEntry, val, yCenterCircSinRow = \
        makeLine(dictList, yCenterCircleSineSV, BaseFigCircleSineFrame)
    yCenterCircleSineSV.set(val)

    # CircleSine Radius
    dictList = xfig.datCircleSine['Radius']
    RadiusCircleSineLabel, RadiusCircleSineEntry, val, radiusCircSinRow = \
        makeLine(dictList, RadiusCircleSineSV, BaseFigCircleSineFrame)
    RadiusCircleSineSV.set(val)

    # Number of Points on CircleSine
    dictList = xfig.datCircleSine['nPts']
    nPtsCircleSineLabel, nPtsCircleSineEntry, val, nptsCircSinRow = \
        makeLine(dictList, nPtsCircleSineSV, BaseFigCircleSineFrame)
    nPtsCircleSineSV.set(val)

    dictList = xfig.datCircleSine['Cycles']
    CyclesCircleSineLabel, CyclesCircleSineEntry, val, cyclesCircSinRow = \
        makeLine(dictList, CyclesCircleSineSV, BaseFigCircleSineFrame)
    CyclesCircleSineSV.set(val)

    dictList = xfig.datCircleSine['Amplitude']
    AmplitudeCircleSineLabel, AmplitudeCircleSineEntry, val, ampCircSinRow = \
        makeLine(dictList, AmplitudeCircleSineSV, BaseFigCircleSineFrame)
    AmplitudeCircleSineSV.set(val)

    xCenterCircleSineLabel.grid(row=xCenterCircSinRow, column=0, sticky=tk.W)
    xCenterCircleSineEntry.grid(row=xCenterCircSinRow, column=1)

    yCenterCircleSineLabel.grid(row=yCenterCircSinRow, column=0, sticky=tk.W)
    yCenterCircleSineEntry.grid(row=yCenterCircSinRow, column=1)

    RadiusCircleSineLabel.grid(row=radiusCircSinRow, column=0, sticky=tk.W)
    RadiusCircleSineEntry.grid(row=radiusCircSinRow, column=1)

    nPtsCircleSineLabel.grid(row=nptsCircSinRow, column=0, sticky=tk.W)
    nPtsCircleSineEntry.grid(row=nptsCircSinRow, column=1)

    CyclesCircleSineLabel.grid(row=cyclesCircSinRow, column=0, sticky=tk.W)
    CyclesCircleSineEntry.grid(row=cyclesCircSinRow, column=1)

    AmplitudeCircleSineLabel.grid(row=ampCircSinRow, column=0, sticky=tk.W)
    AmplitudeCircleSineEntry.grid(row=ampCircSinRow, column=1)

    # Star
    # Tabbed Window: Base Fig Star
    BaseFigStarFrame = ttk.Frame(baseFigTabCtrl)
    baseFigTabCtrl.add(BaseFigStarFrame, text='Star')

    # X Coordinate - Star Center
    dictList = xfig.datStar['xCenterStar']
    xCenterStarLabel, xCenterStarEntry, val, xCenterStarRow = \
        makeLine(dictList, xCenterStarSV, BaseFigStarFrame)
    xCenterStarSV.set(val)

    # Y Coordinate - Star Center
    dictList = xfig.datStar['yCenterStar']
    yCenterStarLabel, yCenterStarEntry, val, yCenterStarRow = \
        makeLine(dictList, yCenterStarSV, BaseFigStarFrame)
    yCenterStarSV.set(val)

    # Star - Spike Start Radius
    dictList = xfig.datStar['StartRadius']
    StartRadiusStarLabel, StartRadiusStarEntry, val, StartRadiusStarRow = \
        makeLine(dictList, StartRadiusStarSV, BaseFigStarFrame)
    StartRadiusStarSV.set(val)

    # Star - Spike End Radius
    dictList = xfig.datStar['EndRadius']
    EndRadiusStarLabel, EndRadiusStarEntry, val, EndRadiusStarRow = \
        makeLine(dictList, EndRadiusStarSV, BaseFigStarFrame)
    EndRadiusStarSV.set(val)

    # Star - Number of Spikes
    dictList = xfig.datStar['nSpikes']
    nSpikesStarLabel, nSpikesStarEntry, val, nSpikesStarRow = \
        makeLine(dictList, nSpikesStarSV, BaseFigStarFrame)
    nSpikesStarSV.set(val)

    # Star - Initial Rotation
    dictList = xfig.datStar['degRotation']
    RotateStarLabel, RotateStarEntry, val, rotateStarRow = \
        makeLine(dictList, RotateStarSV, BaseFigStarFrame)
    RotateStarSV.set(val)

    # Number of Points on Star
    dictList = xfig.datStar['nPts']
    nPtsStarLabel, nPtsStarEntry, val, nptsStarRow = \
        makeLine(dictList, nPtsStarSV, BaseFigStarFrame)
    nPtsStarSV.set(val)

    xCenterStarLabel.grid(row=xCenterStarRow, column=0, sticky=tk.W)
    xCenterStarEntry.grid(row=xCenterStarRow, column=1)

    yCenterStarLabel.grid(row=yCenterStarRow, column=0, sticky=tk.W)
    yCenterStarEntry.grid(row=yCenterStarRow, column=1)

    StartRadiusStarLabel.grid(row=StartRadiusStarRow, column=0,
                              sticky=tk.W)
    StartRadiusStarEntry.grid(row=StartRadiusStarRow, column=1)

    EndRadiusStarLabel.grid(row=EndRadiusStarRow, column=0, sticky=tk.W)
    EndRadiusStarEntry.grid(row=EndRadiusStarRow, column=1)

    nSpikesStarLabel.grid(row=nSpikesStarRow, column=0, sticky=tk.W)
    nSpikesStarEntry.grid(row=nSpikesStarRow, column=1)

    RotateStarLabel.grid(row=rotateStarRow, column=0, sticky=tk.W)
    RotateStarEntry.grid(row=rotateStarRow, column=1)

    nPtsStarLabel.grid(row=nptsStarRow, column=0, sticky=tk.W)
    nPtsStarEntry.grid(row=nptsStarRow, column=1)

    # Spiral
    # Tabbed Window: Base Fig: Spiral
    BaseFigSpiralFrame = ttk.Frame(baseFigTabCtrl)
    baseFigTabCtrl.add(BaseFigSpiralFrame, text='Spiral')

    # X Coordinate - Spiral Center
    dictList = xfig.datSpiral['xCenterSpiral']
    xCenterSpiralLabel, xCenterSpiralEntry, val, xCenterSpiralRow = \
        makeLine(dictList, xCenterSpiralSV, BaseFigSpiralFrame)
    xCenterSpiralSV.set(val)

    # Y Coordinate - Spiral Center
    dictList = xfig.datSpiral['yCenterSpiral']
    yCenterSpiralLabel, yCenterSpiralEntry, val, yCenterSpiralRow = \
        makeLine(dictList, yCenterSpiralSV, BaseFigSpiralFrame)
    yCenterSpiralSV.set(val)

    # Spiral Starting Radius
    dictList = xfig.datSpiral['StartRadius']
    StartRadiusSpiralLabel, StartRadiusSpiralEntry, val, \
        startRadiusSpiralRow = \
        makeLine(dictList, StartRadiusSpiralSV, BaseFigSpiralFrame)
    StartRadiusSpiralSV.set(val)

    # Spiral Ending Radius
    dictList = xfig.datSpiral['EndRadius']
    EndRadiusSpiralLabel, EndRadiusSpiralEntry, val, endRadiusSpiralRow = \
        makeLine(dictList, EndRadiusSpiralSV, BaseFigSpiralFrame)
    EndRadiusSpiralSV.set(val)

    # Spiral - Number of Turns
    dictList = xfig.datSpiral['nTurns']
    nTurnsSpiralLabel, nTurnsSpiralEntry, val, nTurnsSpiralRow = \
        makeLine(dictList, nTurnsSpiralSV, BaseFigSpiralFrame)
    nTurnsSpiralSV.set(val)

    # Spiral - Initial Rotation
    dictList = xfig.datSpiral['degRotation']
    RotateSpiralLabel, RotateSpiralEntry, val, rotateSpiralRow = \
        makeLine(dictList, RotateSpiralSV, BaseFigSpiralFrame)
    RotateSpiralSV.set(val)

    # Number of Points on Spiral
    dictList = xfig.datSpiral['nPts']
    nPtsSpiralLabel, nPtsSpiralEntry, val, nptsSpiralRow = \
        makeLine(dictList, nPtsSpiralSV, BaseFigSpiralFrame)
    nPtsSpiralSV.set(val)

    xCenterSpiralLabel.grid(row=xCenterSpiralRow, column=0, sticky=tk.W)
    xCenterSpiralEntry.grid(row=xCenterSpiralRow, column=1)

    yCenterSpiralLabel.grid(row=yCenterSpiralRow, column=0, sticky=tk.W)
    yCenterSpiralEntry.grid(row=yCenterSpiralRow, column=1)

    StartRadiusSpiralLabel.grid(row=startRadiusSpiralRow, column=0,
                                sticky=tk.W)
    StartRadiusSpiralEntry.grid(row=startRadiusSpiralRow, column=1)

    EndRadiusSpiralLabel.grid(row=endRadiusSpiralRow, column=0, sticky=tk.W)
    EndRadiusSpiralEntry.grid(row=endRadiusSpiralRow, column=1)

    nTurnsSpiralLabel.grid(row=nTurnsSpiralRow, column=0, sticky=tk.W)
    nTurnsSpiralEntry.grid(row=nTurnsSpiralRow, column=1)

    RotateSpiralLabel.grid(row=rotateSpiralRow, column=0, sticky=tk.W)
    RotateSpiralEntry.grid(row=rotateSpiralRow, column=1)

    nPtsSpiralLabel.grid(row=nptsSpiralRow, column=0, sticky=tk.W)
    nPtsSpiralEntry.grid(row=nptsSpiralRow, column=1)

    # Polygon
    # Tabbed Window: Base Fig: Polygon
    BaseFigPolygonFrame = ttk.Frame(baseFigTabCtrl)
    baseFigTabCtrl.add(BaseFigPolygonFrame, text='Polygon')

    # X Coordinate - Polygon Center
    dictList = xfig.datPolygon['xCenterPolygon']
    xCenterPolygonLabel, xCenterPolygonEntry, val, xCenterPolygonRow = \
        makeLine(dictList, xCenterPolygonSV, BaseFigPolygonFrame)
    xCenterPolygonSV.set(val)

    # Y Coordinate - Polygon Center
    dictList = xfig.datPolygon['yCenterPolygon']
    yCenterPolygonLabel, yCenterPolygonEntry, val, yCenterPolygonRow = \
        makeLine(dictList, yCenterPolygonSV, BaseFigPolygonFrame)
    yCenterPolygonSV.set(val)

    # Polygon Radius
    dictList = xfig.datPolygon['Radius']
    RadiusPolygonLabel, RadiusPolygonEntry, val, radiusPolygonRow = \
        makeLine(dictList, RadiusPolygonSV, BaseFigPolygonFrame)
    RadiusPolygonSV.set(val)

    # Polygon - Number of Vertices
    dictList = xfig.datPolygon['nVert']
    nVertPolygonLabel, nVertPolygonEntry, val, nvertPolygonRow = \
        makeLine(dictList, nVertPolygonSV, BaseFigPolygonFrame)
    nVertPolygonSV.set(val)

    # Polygon - Initial Rotation
    dictList = xfig.datPolygon['degRotation']
    RotatePolygonLabel, RotatePolygonEntry, val, rotatePolygonRow = \
        makeLine(dictList, RotatePolygonSV, BaseFigPolygonFrame)
    RotatePolygonSV.set(val)

    # Number of Points on Polygon
    dictList = xfig.datPolygon['nPts']
    nPtsPolygonLabel, nPtsPolygonEntry, val, nptsPolygonRow = \
        makeLine(dictList, nPtsPolygonSV, BaseFigPolygonFrame)
    nPtsPolygonSV.set(val)

    xCenterPolygonLabel.grid(row=xCenterPolygonRow, column=0, sticky=tk.W)
    xCenterPolygonEntry.grid(row=xCenterPolygonRow, column=1)

    yCenterPolygonLabel.grid(row=yCenterPolygonRow, column=0, sticky=tk.W)
    yCenterPolygonEntry.grid(row=yCenterPolygonRow, column=1)

    RadiusPolygonLabel.grid(row=radiusPolygonRow, column=0, sticky=tk.W)
    RadiusPolygonEntry.grid(row=radiusPolygonRow, column=1)

    nVertPolygonLabel.grid(row=nvertPolygonRow, column=0, sticky=tk.W)
    nVertPolygonEntry.grid(row=nvertPolygonRow, column=1)

    RotatePolygonLabel.grid(row=rotatePolygonRow, column=0, sticky=tk.W)
    RotatePolygonEntry.grid(row=rotatePolygonRow, column=1)

    nPtsPolygonLabel.grid(row=8, column=0, sticky=tk.W)
    nPtsPolygonEntry.grid(row=8, column=1)

    # SpiraGon
    # Tabbed Window: Base Fig: SpiraGon
    BaseFigSpiraGonFrame = ttk.Frame(baseFigTabCtrl)
    baseFigTabCtrl.add(BaseFigSpiraGonFrame, text='SpiraGon')

    # X Coordinate - SpiraGon Center
    dictList = xfig.datSpiraGon['xCenterSpiraGon']
    xCenterSpiraGonLabel, xCenterSpiraGonEntry, val, xCenterSpiraGonRow = \
        makeLine(dictList, xCenterSpiraGonSV, BaseFigSpiraGonFrame)
    xCenterSpiraGonSV.set(val)

    # Y Coordinate - SpiraGon Center
    dictList = xfig.datSpiraGon['yCenterSpiraGon']
    yCenterSpiraGonLabel, yCenterSpiraGonEntry, val, yCenterSpiraGonRow = \
        makeLine(dictList, yCenterSpiraGonSV, BaseFigSpiraGonFrame)
    yCenterSpiraGonSV.set(val)

    # SpiraGon Start Radius
    dictList = xfig.datSpiraGon['StartRadius']
    StartRadiusSpiraGonLabel, StartRadiusSpiraGonEntry, val, \
        StartRadiusSpiraGonRow = \
        makeLine(dictList, StartRadiusSpiraGonSV, BaseFigSpiraGonFrame)
    StartRadiusSpiraGonSV.set(val)

    # SpiraGon End Radius
    dictList = xfig.datSpiraGon['EndRadius']
    EndRadiusSpiraGonLabel, EndRadiusSpiraGonEntry, val, \
        EndRadiusSpiraGonRow = \
        makeLine(dictList, EndRadiusSpiraGonSV, BaseFigSpiraGonFrame)
    EndRadiusSpiraGonSV.set(val)

    # SpiraGon Number of Vertices
    dictList = xfig.datSpiraGon['nVert']
    nVertSpiraGonLabel, nVertSpiraGonEntry, val, nvertSpiraGonRow = \
        makeLine(dictList, nVertSpiraGonSV, BaseFigSpiraGonFrame)
    nVertSpiraGonSV.set(val)

    # SpiraGon Number of Turns
    dictList = xfig.datSpiraGon['nTurns']
    nTurnsSpiraGonLabel, nTurnsSpiraGonEntry, val, nTurnsSpiraGonRow = \
        makeLine(dictList, nTurnsSpiraGonSV, BaseFigSpiraGonFrame)
    nTurnsSpiraGonSV.set(val)

    # SpiraGon Initial Rotation
    dictList = xfig.datSpiraGon['degRotation']
    RotateSpiraGonLabel, RotateSpiraGonEntry, val, rotateSpiraGonRow = \
        makeLine(dictList, RotateSpiraGonSV, BaseFigSpiraGonFrame)
    RotateSpiraGonSV.set(val)

    # Number of Points on SpiraGon
    dictList = xfig.datSpiraGon['nPts']
    nPtsSpiraGonLabel, nPtsSpiraGonEntry, val, nptsSpiraGonRow = \
        makeLine(dictList, nPtsSpiraGonSV, BaseFigSpiraGonFrame)
    nPtsSpiraGonSV.set(val)

    xCenterSpiraGonLabel.grid(row=xCenterSpiraGonRow, column=0, sticky=tk.W)
    xCenterSpiraGonEntry.grid(row=xCenterSpiraGonRow, column=1)

    yCenterSpiraGonLabel.grid(row=yCenterSpiraGonRow, column=0, sticky=tk.W)
    yCenterSpiraGonEntry.grid(row=yCenterSpiraGonRow, column=1)

    StartRadiusSpiraGonLabel.grid(row=StartRadiusSpiraGonRow, column=0,
                                  sticky=tk.W)
    StartRadiusSpiraGonEntry.grid(row=StartRadiusSpiraGonRow, column=1)

    EndRadiusSpiraGonLabel.grid(row=EndRadiusSpiraGonRow, column=0,
                                sticky=tk.W)
    EndRadiusSpiraGonEntry.grid(row=EndRadiusSpiraGonRow, column=1)

    nVertSpiraGonLabel.grid(row=nvertSpiraGonRow, column=0, sticky=tk.W)
    nVertSpiraGonEntry.grid(row=nvertSpiraGonRow, column=1)

    nTurnsSpiraGonLabel.grid(row=nTurnsSpiraGonRow, column=0, sticky=tk.W)
    nTurnsSpiraGonEntry.grid(row=nTurnsSpiraGonRow, column=1)

    RotateSpiraGonLabel.grid(row=rotateSpiraGonRow, column=0, sticky=tk.W)
    RotateSpiraGonEntry.grid(row=rotateSpiraGonRow, column=1)

    nPtsSpiraGonLabel.grid(row=nptsSpiraGonRow, column=0, sticky=tk.W)
    nPtsSpiraGonEntry.grid(row=nptsSpiraGonRow, column=1)

    # Separator line
    iterativeSeparatorRow = xfig.datIterative['Separator'][2]
    sep1 = ttk.Separator(dataFrame, orient='horizontal')
    sep1.grid(row=iterativeSeparatorRow, column=0, pady=1,
              stick=tk.EW, columnspan=2)

    # Section Label
    txt = '*Iterative Calculation Data*'
    iterativeSectionLblRow = xfig.datIterative['SectionLbl'][2]
    lbl = tk.Label(dataFrame, width=xfig.widthSectionTitle, text=txt)
    lbl.config(font=xfig.fontSectionTitle)
    lbl.grid(row=iterativeSectionLblRow, column=0, ipadx=ipadx, ipady=ipady,
             sticky=tk.N, columnspan=2)

    # Iterative Calculations - A Parameter
    dictList = xfig.datIterative['aParam']
    aParamLabel, aParamEntry, val, aParamRow = \
        makeLine(dictList, aParamSV, dataFrame)
    aParamSV.set(val)

    # Iterative Calculations - B Parameter
    dictList = xfig.datIterative['bParam']
    bParamLabel, bParamEntry, val, bParamRow = \
        makeLine(dictList, bParamSV, dataFrame)
    bParamSV.set(val)

    # Iterative Calculations - Points Per Orbit
    dictList = xfig.datIterative['orbitPts']
    nOrbitPtsLabel, nOrbitPtsEntry, val, orbitPtsRow = \
        makeLine(dictList, nOrbitPtsSV, dataFrame)
    nOrbitPtsSV.set(val)

    aParamLabel.grid(row=aParamRow, column=0, sticky=tk.W)
    aParamEntry.grid(row=aParamRow, column=1)
    bParamLabel.grid(row=bParamRow, column=0, sticky=tk.W)
    bParamEntry.grid(row=bParamRow, column=1)
    nOrbitPtsLabel.grid(row=orbitPtsRow, column=0, sticky=tk.W)
    nOrbitPtsEntry.grid(row=orbitPtsRow, column=1)

    # Separator line
    datScreenSeparatorRow = xfig.datScreen['Separator'][2]
    sep2 = ttk.Separator(dataFrame, orient='horizontal')
    sep2.grid(row=datScreenSeparatorRow, column=0, pady=1,
              stick=tk.EW, columnspan=2)

    # Section Label
    txt = '*Screen Data*'
    datScreenSectionLblRow = xfig.datScreen['SectionLbl'][2]
    lbl = tk.Label(dataFrame, width=xfig.widthSectionTitle, text=txt)
    lbl.config(font=xfig.fontSectionTitle)
    lbl.grid(row=datScreenSectionLblRow, column=0, ipadx=ipadx, ipady=ipady,
             sticky=tk.N, columnspan=2)

    # Screen Data - Plot Start Iteration
    dictList = xfig.datScreen['pltStart']
    pltStartLabel, pltStartEntry, val, pltStartRow = \
        makeLine(dictList, pltStartSV, dataFrame)
    pltStartSV.set(val)

    # Screen Data - Plot Stop Iteration
    dictList = xfig.datScreen['pltStop']
    pltStopLabel, pltStopEntry, val, pltStopRow = \
        makeLine(dictList, pltStopSV, dataFrame)
    pltStopSV.set(val)

    # Screen Data - Plot Stride Count
    dictList = xfig.datScreen['pltStride']
    pltStrideLabel, pltStrideEntry, val, pltStrideRow = \
        makeLine(dictList, pltStrideSV, dataFrame)
    pltStrideSV.set(val)

    # Screen Data - Marker Definition
    dictList = xfig.datScreen['Marker']
    MarkerLabel, MarkerEntry, val, markerRow = \
        makeLine(dictList, MarkerSV, dataFrame)
    MarkerSV.set(val)

    # Screen Data - Marker Size
    dictList = xfig.datScreen['MarkerSize']
    MarkerSizeLabel, MarkerSizeEntry, val, markerSizeRow = \
        makeLine(dictList, MarkerSizeSV, dataFrame)
    MarkerSizeSV.set(val)

    dictList = xfig.datScreen['MarkerTaper']
    lbl = dictList[0]
    val = dictList[1].strip()
    MarkerTaperLabel = tk.Label(dataFrame,
                                width=xfig.widthRowLbl,
                                text=lbl,
                                font=xfig.fontLbl,
                                anchor='w')
    datScreenMarkerTaperRow = dictList[2]

    markerTaperBtn1 = tk.Radiobutton(dataFrame, text='None',
                                     variable=MarkerTaperSV,
                                     value='None', font=xfig.fontRadio,
                                     width=int(xfig.widthRadio/2))
    markerTaperBtn2 = tk.Radiobutton(dataFrame, text='Zero2Full',
                                     variable=MarkerTaperSV,
                                     value='Zero2Full', font=xfig.fontRadio,
                                     width=int(xfig.widthRadio/2))
    markerTaperBtn3 = tk.Radiobutton(dataFrame, text='Full2Zero',
                                     variable=MarkerTaperSV,
                                     value='Full2Zero', font=xfig.fontRadio)
    MarkerTaperSV.set(val)
    TaperBtnRow = xfig.datScreen['TaperBtn'][2]

    pltStartLabel.grid(row=pltStartRow, column=0, sticky=tk.W)
    pltStartEntry.grid(row=pltStartRow, column=1)
    pltStopLabel.grid(row=pltStopRow, column=0, sticky=tk.W)
    pltStopEntry.grid(row=pltStopRow, column=1)
    pltStrideLabel.grid(row=pltStrideRow, column=0, sticky=tk.W)
    pltStrideEntry.grid(row=pltStrideRow, column=1)
    MarkerLabel.grid(row=markerRow, column=0, sticky=tk.W)
    MarkerEntry.grid(row=markerRow, column=1)
    MarkerSizeLabel.grid(row=markerSizeRow, column=0, sticky=tk.W)
    MarkerSizeEntry.grid(row=markerSizeRow, column=1)
    MarkerTaperLabel.grid(row=datScreenMarkerTaperRow, column=0, sticky=tk.W)
    markerTaperBtn1.grid(row=TaperBtnRow, column=0, columnspan=1, sticky=tk.W)
    markerTaperBtn2.grid(row=TaperBtnRow, column=0, columnspan=1, sticky=tk.E)
    markerTaperBtn3.grid(row=TaperBtnRow, column=1, columnspan=1, sticky=tk.W)

    # Separator line
    cmapSepRow = xfig.cmapSection['Separator'][2]
    sep3 = ttk.Separator(dataFrame, orient='horizontal')
    sep3.grid(row=cmapSepRow, column=0, pady=1, stick=tk.EW, columnspan=2)

    # Default CMAP type radio button
    varcMapType.set(xfig.cMapTypeInit)

    # Radio button for Standard Matplotlib colormap
    radioStdRow = xfig.cmapSection['radioSTD'][2]
    radioSTD = tk.Radiobutton(dataFrame,
                              text=xfig.lblStdCMAP,
                              bg='white',
                              variable=varcMapType,
                              value='STD',
                              font=xfig.fontRadio,
                              anchor=tk.W,
                              width=xfig.widthRadio)
    radioSTD.grid(row=radioStdRow, column=0, ipadx=ipadx, ipady=ipady,
                  padx=padx, pady=pady, sticky=tk.W)

    # Name of Standard Matplotlib colormap
    StdCmapEntry = tk.Entry(dataFrame,
                            textvariable=cMapStd,
                            font=xfig.fontEnt,
                            width=xfig.widthRowEnt)
    StdCmapEntry.grid(row=radioStdRow, column=1, ipadx=ipadx, ipady=ipady,
                      padx=padx, pady=pady, sticky=tk.N)
    StdCmapEntry.insert(0, ' '+xfig.defaultStdCMAP)

    # Radio button for Color List colormap
    radioCustomRow = xfig.cmapSection['radioCUST'][2]
    radioCUSTOMLIST = tk.Radiobutton(dataFrame,
                                     text=xfig.lblCustomCMAP,
                                     bg='white',
                                     variable=varcMapType,
                                     value='CUST',
                                     font=xfig.fontRadio,
                                     anchor=tk.W,
                                     width=xfig.widthRadio)
    radioCUSTOMLIST.grid(row=radioCustomRow, column=0, ipadx=ipadx,
                         ipady=ipady, padx=padx, pady=pady, sticky=tk.W)

    # Name of Custom:LIST colormap
    CustomCmapEntry = tk.Entry(dataFrame,
                               textvariable=cMapList,
                               width=xfig.widthRowEnt,
                               state=tk.NORMAL,
                               font=xfig.fontEnt)
    CustomCmapEntry.grid(row=radioCustomRow, column=1, ipadx=ipadx,
                         ipady=ipady, padx=padx, pady=pady, sticky=tk.N)
    CustomCmapEntry.insert(0, ' '+xfig.defaultCustomCMAP)

    # Separator line
    axesSepRow = xfig.axesSection['Separator'][2]
    sep4 = ttk.Separator(dataFrame, orient='horizontal')
    sep4.grid(row=axesSepRow, column=0, pady=1, stick=tk.EW, columnspan=2)

    # Default Figure axes scaling
    varFigureScaling.set(xfig.AxesInit)

    # Radio button for equal sized X&Y axes
    radioEqualXYRow = xfig.axesSection['axesEQUAL'][2]
    radioEqualXY = tk.Radiobutton(dataFrame,
                                  text=xfig.lblEqualAxes,
                                  bg='white',
                                  variable=varFigureScaling,
                                  value='EQUAL',
                                  font=xfig.fontRadio,
                                  anchor=tk.W,
                                  width=xfig.widthRadio)
    radioEqualXY.grid(row=radioEqualXYRow, column=0, ipadx=ipadx,
                      ipady=ipady, padx=padx, pady=pady, sticky=tk.W)

    # Radio button for axes sizing from X&Y MinMax
    radioMINMAXRow = xfig.axesSection['axesMINMAX'][2]
    radioMinMax = tk.Radiobutton(dataFrame,
                                 text=xfig.lblAxesMinMax,
                                 bg='white',
                                 variable=varFigureScaling,
                                 value='MINMAX',
                                 font=xfig.fontRadio,
                                 anchor=tk.W,
                                 width=xfig.widthRadio)
    radioMinMax.grid(row=radioMINMAXRow, column=0, ipadx=ipadx, ipady=ipady,
                     padx=padx, pady=pady, sticky=tk.W)

    # Separator line
    saveRawSepRow = xfig.saveRawSection['Separator'][2]
    sep5 = ttk.Separator(dataFrame, orient='horizontal')
    sep5.grid(row=saveRawSepRow, column=0, pady=1, stick=tk.EW, columnspan=1)

    # Save raw data check button
    btnSaveRawRow = xfig.saveRawSection['rawBtn'][2]
    rawCkBtn = tk.Checkbutton(
        dataFrame, text='SaveRaw', bg='white',
        variable=saveRaw, onvalue=True, offvalue=False,
        font=xfig.fontLbl, anchor=tk.W, width=xfig.widthCkBtn)
    rawCkBtn.grid(row=btnSaveRawRow, column=0, ipadx=ipadx,
                  padx=padx, sticky=tk.E)

    # Create frame for Iterative Algorithm Version
    algorithmFrame = tk.Frame(root, height=hgt, bg=xfig.backgndColor)
    algorithmFrame.grid(row=1, column=2, ipadx=ipadx, ipady=ipady,
                        sticky=tk.N, columnspan=2)

    txt = '*Iterative Algorithm*'
    lbl = tk.Label(algorithmFrame, width=xfig.widthSectionTitle, text=txt)
    lbl.config(font=xfig.fontSectionTitle)
    lbl.grid(row=2, column=0, ipadx=ipadx, ipady=ipady,
             sticky=tk.N, columnspan=2)

    # Create tabbed window manager
    algTabCtrl = ttk.Notebook(algorithmFrame)
    algTabCtrl.grid(row=3, column=0, ipadx=ipadx, ipady=ipady,
                    sticky=tk.W, columnspan=2)
    algTabCtrl.bind('<<NotebookTabChanged>>', on_tabAlg_selected)

    # Tabbed Window: Algorithm Version 1
    tabAlg1 = ttk.Frame(algTabCtrl)
    algTabCtrl.add(tabAlg1, text='Version1')
    img1 = ImageTk.PhotoImage(file=xfig.tabAlgsRefs['Version1'])
    canvas1 = tk.Canvas(tabAlg1, bd=0, width=img1.width(),
                        height=img1.height())
    canvas1.grid(row=0, column=0, ipadx=ipadx, ipady=ipady, sticky=tk.N)
    canvas1.create_image(0, 0, image=img1, anchor=tk.NW)
    # Following prevents image garbage collection after return
    # which allows display of the image.
    canvas1.img = img1

    # Tabbed Window: Algorithm Version 2
    tabAlg2 = ttk.Frame(algTabCtrl)
    algTabCtrl.add(tabAlg2, text='Version2')
    img2 = ImageTk.PhotoImage(file=xfig.tabAlgsRefs['Version2'])
    canvas2 = tk.Canvas(tabAlg2, bd=0, width=img2.width(),
                        height=img2.height())
    canvas2.grid(row=0, column=0, ipadx=ipadx, ipady=ipady, sticky=tk.N)
    canvas2.create_image(0, 0, image=img2, anchor=tk.NW)
    # Following prevents image garbage collection after return
    # which allows display of the image.
    canvas2.img = img2

    # Tabbed Window: Algorithm Version 3
    tabAlg3 = ttk.Frame(algTabCtrl)
    algTabCtrl.add(tabAlg3, text='Version3')
    img3 = ImageTk.PhotoImage(file=xfig.tabAlgsRefs['Version3'])
    canvas3 = tk.Canvas(tabAlg3, bd=0, width=img3.width(),
                        height=img3.height())
    canvas3.grid(row=0, column=0, ipadx=ipadx, ipady=ipady, sticky=tk.N)
    canvas3.create_image(0, 0, image=img3, anchor=tk.NW)
    # Following prevents image garbage collection after return
    # which allows display of the image.
    canvas3.img = img3


# Load and parse image parameters
def getImageData():
    global baseFigTabCtrl
    # Base Figure - Circle: Entries
    global xCenterCircleEntry, yCenterCircleEntry
    global RadiusCircleEntry, nPtsCircleEntry
    # Base Figure - CircleSine: Entries
    global xCenterCircleSineEntry, yCenterCircleSineEntry
    global RadiusCircleSineEntry, nPtsCircleSineEntry
    global CyclesCircleSineEntry, AmplitudeCircleSineEntry
    # Base figure - Star
    global xCenterStarEntry, yCenterStarEntry
    global StartRadiusStarEntry, EndRadiusStarEntry
    global nSpikesStarEntry, RotateStarEntry, nPtsStarEntry
    # Base Figure - Spiral: Entries
    global xCenterSpiralEntry, yCenterSpiralEntry
    global StartRadiusSpiralEntry, EndRadiusSpiralEntry
    global nTurnsSpiralEntry, RotateSpiralEntry, nPtsSpiralEntry
    # Base Figure - Polygon: Entries
    global xCenterPolygonEntry, yCenterPolygonEntry, RadiusPolygonEntry
    global nVertPolygonEntry, RotatePolygonEntry, nPtsPolygonEntry
    # Base Figure - SpiraGon: Entries
    global xCenterSpiraGonEntry, yCenterSpiraGonEntry
    global EndRadiusSpiraGonEntry, EndRadiusSpiraGonEntry
    global nVertSpiraGonEntry, RotateSpiraGonEntry
    global nTurnsSpiraGonEntry, nPtsSpiraGonEntry
    # Iterative Calculation: Entries
    global aParamEntry, bParamEntry, nOrbitPtsEntry
    # Screen Data: Entries
    global pltStartEntry, pltStopEntry, pltStrideEntry
    global MarkerEntry, MarkerSizeEntry
    global varcMapType, StdCmapEntry, CustomCmapEntry
    global algTabCtrl, varAlgorithm

    global varBaseFig

    # Load/read image parameters file
    def dataFileReader():
        dataDict = OrderedDict()
        name = filedialog.askopenfilename(
                initialdir=".",
                title="Select file",
                filetypes=(("Image Parameters", "*.agd"),
                           ("all files", "*.*")))
        if not name:
            return(None)
        else:
            with open(name, 'r') as inFile:
                for m in inFile:
                    # Discard comment lines or blank lines
                    if (m.lstrip().startswith('#')) or (len(m.split()) == 0):
                        continue
                    else:
                        # Discard trailing newline
                        m = m.rstrip('\n')
                        # Discard trailing comment
                        if (m.find("#") != -1):
                            # Found trailing comment
                            temp = m.split('#')
                            m = temp[0]
                        mSpl = m.split('=')
                        key = mSpl[0].replace(' ', '')
                        val = mSpl[1].replace(' ', '')
                        dataDict[key] = val

        return(dataDict)

    # Create dictionary of data loaded from file
    dataDict = dataFileReader()
    if dataDict is None:
        return

    def updateEntry(ident, ent):
        newVal = ' ' + dataDict[ident]
        ent.delete(0, tk.END)
        ent.insert(0, newVal)
        return(newVal)

    # Process Base Figure type and parameters
    BASEFIG = dataDict['baseFig']
    if BASEFIG == 'Circle':
        # Set tab focus to Circle
        varBaseFig.set('Circle')
        baseFigTabCtrl.select([BaseFigCircleFrame])

        newVal = updateEntry('xCenter', xCenterCircleEntry)
        xCenterCircleSV.set(newVal)

        newVal = updateEntry('yCenter', yCenterCircleEntry)
        yCenterCircleSV.set(newVal)

        newVal = updateEntry('Radius', RadiusCircleEntry)
        RadiusCircleSV.set(newVal)

        newVal = updateEntry('nPts', nPtsCircleEntry)
        nPtsCircleSV.set(newVal)

    elif BASEFIG == 'CircleSine':
        # Set tab focus to CircleSine
        varBaseFig.set('CircleSine')
        baseFigTabCtrl.select([BaseFigCircleSineFrame])

        newVal = updateEntry('xCenter', xCenterCircleSineEntry)
        xCenterCircleSineSV.set(newVal)

        newVal = updateEntry('yCenter', yCenterCircleSineEntry)
        yCenterCircleSineSV.set(newVal)

        newVal = updateEntry('Radius', RadiusCircleSineEntry)
        RadiusCircleSineSV.set(newVal)

        newVal = updateEntry('nPts', nPtsCircleSineEntry)
        nPtsCircleSineSV.set(newVal)

        newVal = updateEntry('Cycles', CyclesCircleSineEntry)
        CyclesCircleSineSV.set(newVal)

        newVal = updateEntry('Amplitude', AmplitudeCircleSineEntry)
        AmplitudeCircleSineSV.set(newVal)

    elif BASEFIG == 'Star':
        # Set tab focus to Star
        varBaseFig.set('Star')
        baseFigTabCtrl.select([BaseFigStarFrame])

        newVal = updateEntry('xCenter', xCenterStarEntry)
        xCenterStarSV.set(newVal)

        newVal = updateEntry('yCenter', yCenterStarEntry)
        yCenterStarSV.set(newVal)

        newVal = updateEntry('StartRadius', StartRadiusStarEntry)
        StartRadiusStarSV.set(newVal)

        newVal = updateEntry('EndRadius', EndRadiusStarEntry)
        EndRadiusStarSV.set(newVal)

        newVal = updateEntry('nSpikes', nSpikesStarEntry)
        nSpikesStarSV.set(newVal)

        newVal = updateEntry('initialRotation', RotateStarEntry)
        RotateStarSV.set(newVal)

        newVal = updateEntry('nPts', nPtsStarEntry)
        nPtsStarSV.set(newVal)
    elif BASEFIG == 'Spiral':
        # Set tab focus to Spiral
        varBaseFig.set('Spiral')
        baseFigTabCtrl.select([BaseFigSpiralFrame])

        newVal = updateEntry('xCenter', xCenterSpiralEntry)
        xCenterSpiralSV.set(newVal)

        newVal = updateEntry('yCenter', yCenterSpiralEntry)
        yCenterSpiralSV.set(newVal)

        newVal = updateEntry('StartRadius', StartRadiusSpiralEntry)
        StartRadiusSpiralSV.set(newVal)

        newVal = updateEntry('EndRadius', EndRadiusSpiralEntry)
        EndRadiusSpiralSV.set(newVal)

        newVal = updateEntry('nTurns', nTurnsSpiralEntry)
        nTurnsSpiralSV.set(newVal)

        newVal = updateEntry('initialRotation', RotateSpiralEntry)
        RotateSpiralSV.set(newVal)

        newVal = updateEntry('nPts', nPtsSpiralEntry)
        nPtsSpiralSV.set(newVal)

    elif BASEFIG == 'Polygon':
        # Set tab focus to a regular polygon
        varBaseFig.set('Polygon')
        baseFigTabCtrl.select([BaseFigPolygonFrame])

        newVal = updateEntry('xCenter', xCenterPolygonEntry)
        xCenterPolygonSV.set(newVal)

        newVal = updateEntry('yCenter', yCenterPolygonEntry)
        yCenterPolygonSV.set(newVal)

        newVal = updateEntry('Radius', RadiusPolygonEntry)
        RadiusPolygonSV.set(newVal)

        newVal = updateEntry('nVert', nVertPolygonEntry)
        nVertPolygonSV.set(newVal)

        newVal = updateEntry('initialRotation', RotatePolygonEntry)
        RotatePolygonSV.set(newVal)

        newVal = updateEntry('nPts', nPtsPolygonEntry)
        nPtsPolygonSV.set(newVal)

    elif BASEFIG == 'SpiraGon':
        # Set tab focus to a regular polygon with varying radius
        varBaseFig.set('SpiraGon')
        baseFigTabCtrl.select([BaseFigSpiraGonFrame])

        newVal = updateEntry('xCenter', xCenterSpiraGonEntry)
        xCenterSpiraGonSV.set(newVal)

        newVal = updateEntry('yCenter', yCenterSpiraGonEntry)
        yCenterSpiraGonSV.set(newVal)

        newVal = updateEntry('StartRadius', StartRadiusSpiraGonEntry)
        StartRadiusSpiraGonSV.set(newVal)

        newVal = updateEntry('EndRadius', EndRadiusSpiraGonEntry)
        EndRadiusSpiraGonSV.set(newVal)

        newVal = updateEntry('nVert', nVertSpiraGonEntry)
        nVertSpiraGonSV.set(newVal)

        newVal = updateEntry('nTurns', nTurnsSpiraGonEntry)
        nTurnsSpiraGonSV.set(newVal)

        newVal = updateEntry('initialRotation', RotateSpiraGonEntry)
        RotateSpiraGonSV.set(newVal)

        newVal = updateEntry('nPts', nPtsSpiraGonEntry)
        nPtsSpiraGonSV.set(newVal)
    else:
        print('Unknown Base Figure in input file. Returning')
        return

    # Process Iterative Calculations parameters
    newVal = updateEntry('aParam', aParamEntry)
    aParamSV.set(newVal)

    newVal = updateEntry('bParam', bParamEntry)
    bParamSV.set(newVal)

    newVal = updateEntry('itersPerOrbit', nOrbitPtsEntry)
    nOrbitPtsSV.set(newVal)

    # Process numerical Screen Data
    newVal = updateEntry('pltStart', pltStartEntry)
    pltStartSV.set(newVal)

    newVal = updateEntry('pltStop', pltStopEntry)
    pltStopSV.set(newVal)

    newVal = updateEntry('pltStride', pltStrideEntry)
    pltStrideSV.set(newVal)

    newVal = updateEntry('marker', MarkerEntry)
    MarkerSV.set(newVal)

    newVal = updateEntry('markerSize', MarkerSizeEntry)
    MarkerSizeSV.set(newVal)

    newVal = dataDict['markerTaper']
    MarkerTaperSV.set(newVal)

    # Process CMAP data
    newVal = dataDict['cMapType']
    varcMapType.set(newVal)
    if newVal == 'STD':
        newCmap = updateEntry('cMapStd', StdCmapEntry)
        cMapStd.set(newCmap)
    else:
        newCmap = updateEntry('cMapList', CustomCmapEntry)
        cMapList.set(newCmap)

    # Process image axes scaling
    newVal = dataDict['figScale']
    varFigureScaling.set(newVal)

    # Process Save Raw Data checkbox setting
    newVal = dataDict['saveRaw']
    if newVal.upper() == 'TRUE':
        saveRaw.set(True)
    elif newVal.upper() == 'FALSE':
        saveRaw.set(False)
    else:
        print('Unknown definition for saveRaw. Defaulting False')
        saveRaw.set(False)

    # Process the Iterative Algorithm setting
    newVal = dataDict['algVers']
    if newVal == 'Version1':
        algTabCtrl.select(tabAlg1)
    elif newVal == 'Version2':
        algTabCtrl.select(tabAlg2)
    elif newVal == 'Version3':
        algTabCtrl.select(tabAlg3)
    else:
        print('Erroneous algVers from file load. Defaulting to Version1')
        algTabCtrl.select(tabAlg1)
        varAlgorithm.set('Version1')
        return

    varAlgorithm.set(newVal)


# Load Color List colormap definition file name
def getColorListCMAP():
    global cmFileName, cmFilePath, CustomCmapEntry

    cmFilePath = filedialog.askopenfilename(
        initialdir=".",
        title="Select file",
        filetypes=(("Color List files", "*.cmL"), ("all files", "*.*")))
    if (cmFilePath):
        cmFileName = os.path.basename(cmFilePath)
        CustomCmapEntry.config(state=tk.NORMAL)
        CustomCmapEntry.delete(0, tk.END)
        CustomCmapEntry.insert(0, ' '+cmFileName)
        varcMapType.set("CUST")
    else:
        return


# Make general image for GUI display, or Help menu items
def makeImage(img):
    img = ImageTk.PhotoImage(file=img)
    top = tk.Toplevel()
    wth = img.width()
    hgt = img.height()
    x = 500
    y = 50
    top.geometry('%dx%d+%d+%d' % (wth, hgt, x, y))
    canvas = tk.Canvas(top, bd=0, width=wth, height=hgt)
    canvas.grid()

    canvas.create_image(0, 0, image=img, anchor=tk.NW)
    # Following prevents image garbage collection after return
    # which allows display of the image.
    canvas.img = img


class helpPg:
    def __init__(self, tabName, pgImg, tabMgr):
        ipadx = 2
        ipady = 2

        self.tabMgr = tabMgr
        self.helpPg = ttk.Frame(self.tabMgr)
        self.tabMgr.add(self.helpPg, text=tabName)
        self.imgHelp = ImageTk.PhotoImage(file=pgImg)
        width = self.imgHelp.width()
        height = self.imgHelp.height()
        self.helpCanvas = \
            tk.Canvas(self.helpPg, width=width,
                      height=height, borderwidth=0,
                      highlightthickness=0, background='black')
        self.helpCanvas.config(scrollregion=[0, 0, 1.2*width, 1.2*height])
        self.hbar = tk.Scrollbar(self.helpPg, orient='horizontal')
        self.hbar.grid(row=0, column=0, sticky=tk.EW, columnspan=100)
        self.hbar.config(command=self.helpCanvas.xview)
        self.vbar = tk.Scrollbar(self.helpPg, orient='vertical')
        self.vbar.grid(row=1, column=5, sticky=tk.NS, rowspan=100)
        self.vbar.config(command=self.helpCanvas.yview)

        self.helpCanvas.config(xscrollcommand=self.hbar.set,
                               yscrollcommand=self.vbar.set)
        self.helpCanvas.grid(row=1, column=0,
                             ipadx=ipadx, ipady=ipady, sticky=tk.N)
        self.helpCanvas.create_image(0, 0, image=self.imgHelp, anchor=tk.NW)
        # Following prevents image garbage collection after return
        # which allows display of the image.
        self.helpCanvas.img = self.imgHelp


# Display Help pages
def showHelpPages():
    ipadx = 2
    ipady = 2

    # Create Help window
    helpWindow = tk.Toplevel()
    helpWindow.geometry('+250+5')

    title = 'AgnMultiTk Help'
    lblHelp = tk.Label(helpWindow, text=title)
    lblHelp.grid(row=0, column=0, ipadx=ipadx, ipady=ipady,
                 sticky=tk.N, columnspan=100)

    # Create tabbed window manager
    helpTabCtrl = ttk.Notebook(helpWindow)
    helpTabCtrl.grid(row=2, column=0, ipadx=ipadx, ipady=ipady,
                     sticky=tk.EW, columnspan=5)

    # Generate Detailed Help pages
    tabHelpPg1 = helpPg('Page 1',
                        xfig.helpImgRefs['Page 1'],
                        helpTabCtrl)

    tabHelpPg2 = helpPg('Page 2',
                        xfig.helpImgRefs['Page 2'],
                        helpTabCtrl)

    tabHelpPg3 = helpPg('Page 3',
                        xfig.helpImgRefs['Page 3'],
                        helpTabCtrl)

    tabHelpPg4 = helpPg('Page 4',
                        xfig.helpImgRefs['Page 4'],
                        helpTabCtrl)

    tabHelpPg5 = helpPg('Page 5',
                        xfig.helpImgRefs['Page 5'],
                        helpTabCtrl)

    tabHelpPg6 = helpPg('Page 6',
                        xfig.helpImgRefs['Page 6'],
                        helpTabCtrl)

    tabHelpPg7 = helpPg('Page 7',
                        xfig.helpImgRefs['Page 7'],
                        helpTabCtrl)

    tabHelpPg8 = helpPg('Page 8',
                        xfig.helpImgRefs['Page 8'],
                        helpTabCtrl)

    tabHelpPg9 = helpPg('Page 9',
                        xfig.helpImgRefs['Page 9'],
                        helpTabCtrl)

    tabHelpPg10 = helpPg('Page 10',
                         xfig.helpImgRefs['Page 10'],
                         helpTabCtrl)


def showOverview():
    makeImage(xfig.miscRefs['Overview'])


def showGallery_Circle1():
    makeImage(xfig.galleryRefs['Circle1'])


def showGallery_Circle2():
    makeImage(xfig.galleryRefs['Circle2'])


def showGallery_CircleSine1():
    makeImage(xfig.galleryRefs['CircleSine1'])


def showGallery_CircleSine2():
    makeImage(xfig.galleryRefs['CircleSine2'])


def showGallery_Star1():
    makeImage(xfig.galleryRefs['Star1'])


def showGallery_Star2():
    makeImage(xfig.galleryRefs['Star2'])


def showGallery_Spiral1():
    makeImage(xfig.galleryRefs['Spiral1'])


def showGallery_Spiral2():
    makeImage(xfig.galleryRefs['Spiral2'])


def showGallery_Polygon1():
    makeImage(xfig.galleryRefs['Polygon1'])


def showGallery_Polygon2():
    makeImage(xfig.galleryRefs['Polygon2'])


def showGallery_SpiraGon1():
    makeImage(xfig.galleryRefs['SpiraGon1'])


def showGallery_SpiraGon2():
    makeImage(xfig.galleryRefs['SpiraGon2'])


# Make/display program About window
def aboutTxt():
    aboutWindow = tk.Toplevel()
    winWidth = 700
    winHeight = 400

    posRight = int(aboutWindow.winfo_screenwidth()/2 - winWidth/2)
    posDown = int(aboutWindow.winfo_screenheight()/3 - winHeight/2)
    aboutWindow.geometry('+{}+{}'.format(posRight, posDown))

    # Add thumbnail to About window
    canvas = tk.Canvas(aboutWindow, width=200, height=350)
    canvas.grid(row=0, column=0, sticky=tk.NW, padx=2, pady=2)
    photo = tk.PhotoImage(file=xfig.miscRefs['About'])
    canvas.create_image(102, 100, image=photo)
    # Attach photo to canvas since photo is a local object
    # that gets garbage collected
    canvas.image = photo

    # Insert version and history text.
    titleTxt = '    ** AgnMulti History **\n'
    txtBox2 = tk.Text(aboutWindow, height=25, width=75, padx=10)
    txtBox2.tag_configure('titleFont', font=('Liberation Serif', 20, 'bold'))
    scroll = tk.Scrollbar(aboutWindow, command=txtBox2.yview)
    txtBox2.grid(row=0, column=1)
    scroll.grid(row=0, column=2, sticky=tk.NS)
    txtBox2.insert(tk.END, titleTxt, 'titleFont')

    fileHandle = open('./SupportFiles/aboutText.txt', 'r')
    aboutText = fileHandle.read()  # Read the file to a variable
    fileHandle.close()

    txtBox2.tag_configure('textFont', font=('Liberation Serif', 12))
    txtBox2.insert(tk.END, aboutText, 'textFont')
    txtBox2.config(state=tk.DISABLED)


if __name__ == '__main__':
    # Get any command line arguments
    narg = len(sys.argv)

    if narg > 1:
        for k in range(1, narg):
            if sys.argv[k] in debugDict:
                debugDict[sys.argv[k]] = True
            else:
                print('Unrecognized command line arg: ', sys.argv[k])
                print('Exiting....')
                sys.exit()

    root = tk.Tk()
    root.title('')
    root.geometry("+80+120")
    root.configure(bg=xfig.backgndColor)

    style = ttk.Style()

    # Set Tab styling
    # "active" in settings means the mouse cursor is over the selection
    settings = {"TNotebook.Tab":
                {"configure": {"padding": [5, 1], "background": "lightgrey"},
                 "map": {"background": [("selected", "black"),
                                        ("active", "darkgrey")],
                         "foreground": [("selected", "white"),
                                        ("active", "black")]}
                }
               }

    style.theme_create("NBtabs", parent="alt", settings=settings)
    style.theme_use("NBtabs")

    # Base Figure: Circle - variables
    xCenterCircleSV = tk.StringVar()
    yCenterCircleSV = tk.StringVar()
    RadiusCircleSV = tk.StringVar()
    nPtsCircleSV = tk.StringVar()

    # Base Figure: CircleSine - variables
    xCenterCircleSineSV = tk.StringVar()
    yCenterCircleSineSV = tk.StringVar()
    RadiusCircleSineSV = tk.StringVar()
    nPtsCircleSineSV = tk.StringVar()
    CyclesCircleSineSV = tk.StringVar()
    AmplitudeCircleSineSV = tk.StringVar()

    # Base Figure: Star - variables
    xCenterStarSV = tk.StringVar()
    yCenterStarSV = tk.StringVar()
    StartRadiusStarSV = tk.StringVar()
    EndRadiusStarSV = tk.StringVar()
    nSpikesStarSV = tk.StringVar()
    RotateStarSV = tk.StringVar()
    nPtsStarSV = tk.StringVar()

    # Base Figure: Spiral - variables
    xCenterSpiralSV = tk.StringVar()
    yCenterSpiralSV = tk.StringVar()
    StartRadiusSpiralSV = tk.StringVar()
    EndRadiusSpiralSV = tk.StringVar()
    nTurnsSpiralSV = tk.StringVar()
    RotateSpiralSV = tk.StringVar()
    nPtsSpiralSV = tk.StringVar()

    # Base Figure: Polygon - variables
    xCenterPolygonSV = tk.StringVar()
    yCenterPolygonSV = tk.StringVar()
    RadiusPolygonSV = tk.StringVar()
    nVertPolygonSV = tk.StringVar()
    RotatePolygonSV = tk.StringVar()
    nPtsPolygonSV = tk.StringVar()

    # Base Figure: SpiraGon - variables
    xCenterSpiraGonSV = tk.StringVar()
    yCenterSpiraGonSV = tk.StringVar()
    StartRadiusSpiraGonSV = tk.StringVar()
    EndRadiusSpiraGonSV = tk.StringVar()
    nVertSpiraGonSV = tk.StringVar()
    RotateSpiraGonSV = tk.StringVar()
    nTurnsSpiraGonSV = tk.StringVar()
    nPtsSpiraGonSV = tk.StringVar()

    # Iterative Calculations - variables
    aParamSV = tk.StringVar()
    bParamSV = tk.StringVar()
    nOrbitPtsSV = tk.StringVar()

    # Screen Data - variables
    pltStartSV = tk.StringVar()
    pltStopSV = tk.StringVar()
    pltStrideSV = tk.StringVar()
    MarkerSV = tk.StringVar()
    MarkerSizeSV = tk.StringVar()
    MarkerTaperSV = tk.StringVar()
    cMapStd = tk.StringVar()
    cMapList = tk.StringVar()

    # Declare CMAP variable
    varcMapType = tk.StringVar()

    # Declare Algorithm variable
    varAlgorithm = tk.StringVar()

    # Declare Figure axes scaling
    varFigureScaling = tk.StringVar()

    # Declare Base Figure variable
    varBaseFig = tk.StringVar()

    # Initialize BooleanVar for raw data output when saving image
    saveRaw = tk.BooleanVar()
    flag = xfig.saveRawInit
    if flag.upper() == 'TRUE':
        saveRaw.set(True)
    elif flag.upper() == 'FALSE':
        saveRaw.set(False)
    else:
        print('Unknown definition for saveRaw. Defaulting False')
        saveRaw.set(False)

    # Make initial GUI form
    makeForm(root)

    # Make separate frame for Main window Plot and Quit buttons
    actionFrame1 = tk.Frame(root, bg=xfig.backgndColor)
    actionFrame1.grid(row=2, column=3, sticky=tk.EW)

    # Main window "Quit" button
    btnQuit = tk.Button(actionFrame1, text='   Quit   ',
                        font=xfig.fontButton,
                        command=root.destroy)
    btnQuit.grid(row=0, column=1, sticky=tk.EW)

    # Main window button to make the plot
    btnPlot = tk.Button(actionFrame1, text='   Plot   ',
                        font=xfig.fontButton,
                        command=(makePlot))
    btnPlot.grid(row=0, column=0, padx=2, pady=2)
    actionFrame1.grid(sticky=tk.EW)

    # Define menubar
    menubar = tk.Menu(root)

    # Load Menu
    loadMenu = tk.Menu(menubar, tearoff=False)
    menubar.add_cascade(
        label="Load",
        font=xfig.fontMenu,
        menu=loadMenu)
    # Select image parameters file name
    loadMenu.add_command(
        label="Image Parameters",
        font=xfig.fontMenu,
        command=getImageData)
    # Select custom CMAP:List file name
    loadMenu.add_command(
        label="CMAP:List",
        font=xfig.fontMenu,
        command=getColorListCMAP)

    # Help menu
    helpmenu = tk.Menu(menubar, tearoff=False)
    # Extra space in label prevents Mac attaching Search box to menu choices
    menubar.add_cascade(label=" Help", font=xfig.fontMenu, menu=helpmenu)

    # Overview
    helpmenu.add_command(
        label="Overview",
        font=xfig.fontMenu,
        command=showOverview)

    # Display Help pages
    helpmenu.add_command(
        label="Detailed Help",
        font=xfig.fontMenu,
        command=showHelpPages)

    # Image gallery
    galleryMenu = tk.Menu(tearoff=False)
    helpmenu.add_cascade(
        label="Gallery",
        font=xfig.fontMenu,
        menu=galleryMenu)

    galleryCircleMenu = tk.Menu(tearoff=False)
    galleryMenu.add_cascade(
        label="Circle",
        font=xfig.fontMenu,
        menu=galleryCircleMenu)
    galleryCircleMenu.add_command(
        label="Image 1",
        font=xfig.fontMenu,
        command=showGallery_Circle1)
    galleryCircleMenu.add_command(
        label="Image 2",
        font=xfig.fontMenu,
        command=showGallery_Circle2)

    galleryCircleSineMenu = tk.Menu(tearoff=False)
    galleryMenu.add_cascade(
        label="CircleSine",
        font=xfig.fontMenu,
        menu=galleryCircleSineMenu)
    galleryCircleSineMenu.add_command(
        label="Image 1",
        font=xfig.fontMenu,
        command=showGallery_CircleSine1)
    galleryCircleSineMenu.add_command(
        label="Image 2",
        font=xfig.fontMenu,
        command=showGallery_CircleSine2)

    galleryStarMenu = tk.Menu(tearoff=False)
    galleryMenu.add_cascade(
        label="Star",
        font=xfig.fontMenu,
        menu=galleryStarMenu)
    galleryStarMenu.add_command(
        label="Image 1",
        font=xfig.fontMenu,
        command=showGallery_Star1)
    galleryStarMenu.add_command(
        label="Image 2",
        font=xfig.fontMenu,
        command=showGallery_Star2)

    gallerySpiralMenu = tk.Menu(tearoff=False)
    galleryMenu.add_cascade(
        label="Spiral",
        font=xfig.fontMenu,
        menu=gallerySpiralMenu)
    gallerySpiralMenu.add_command(
        label="Image 1",
        font=xfig.fontMenu,
        command=showGallery_Spiral1)
    gallerySpiralMenu.add_command(
        label="Image 2",
        font=xfig.fontMenu,
        command=showGallery_Spiral2)

    galleryPolygonMenu = tk.Menu(tearoff=False)
    galleryMenu.add_cascade(
        label="Polygon",
        font=xfig.fontMenu,
        menu=galleryPolygonMenu)
    galleryPolygonMenu.add_command(
        label="Image 1",
        font=xfig.fontMenu,
        command=showGallery_Polygon1)
    galleryPolygonMenu.add_command(
        label="Image 2",
        font=xfig.fontMenu,
        command=showGallery_Polygon2)

    gallerySpiraGonMenu = tk.Menu(tearoff=False)
    galleryMenu.add_cascade(
        label="SpiraGon",
        font=xfig.fontMenu,
        menu=gallerySpiraGonMenu)
    gallerySpiraGonMenu.add_command(
        label="Image 1",
        font=xfig.fontMenu,
        command=showGallery_SpiraGon1)
    gallerySpiraGonMenu.add_command(
        label="Image 2",
        font=xfig.fontMenu,
        command=showGallery_SpiraGon2)

    helpmenu.add_command(
        label="About",
        font=xfig.fontMenu,
        command=aboutTxt)

    # Main loop
    root.config(menu=menubar)
    root.mainloop()
