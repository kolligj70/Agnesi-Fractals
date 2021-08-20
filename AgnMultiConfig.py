#!/usr/bin/python3

# Miscellaneous configuration data
windowWidth = 14.0
windowHeight = 9.5
dpi = 100
faceColor = 'black'
backgndColor = 'gainsboro'
baseName = 'AgnMulti_'
plotExpander = 0.025

fontSectionTitle = ("sans-serif 12 bold")
fontLbl = ("sans-serif 12")
fontEnt = ("sans-serif 12")
fontMenu = ("sans-serif 13")
fontRadio = ("sans-serif 12")
fontButton = ("sans-serif 14 bold")

widthRadio = 19
widthRowLbl = 20
widthRowEnt = 20
widthCkBtn = 19
widthSectionTitle = 25

# Label, Default, gridRow
datCircle = {"xCenterCircle" : [" X CircleCenter", " 6.0", 3],
             "yCenterCircle" : [" Y CircleCenter", " -10.0", 4],
             "Radius" : [" Radius", " 5.0", 5],
             "nPts" : [" N CirclePts", " 360", 6]}
datCircleSine = {"xCenterCircleSine" : [" X CircleSineCenter", " 6.0", 3],
                 "yCenterCircleSine" : [" Y CircleSineCenter", " -10.0", 4],
                 "Radius" : [" Radius", " 5.0", 5],
                 "nPts" : [" N CirclePts", " 360", 6],
                 "Cycles" : [" Cycles", " 13", 7],
                 "Amplitude" : [" Amplitude", " 0.5", 8]}
datSpiral = {"xCenterSpiral" : [" X SpiralCenter", " 6.0", 3], 
             "yCenterSpiral" : [" Y SpiralCenter", " -10.0", 4],
             "StartRadius" : [" Start Radius", " 4.5", 5],
             "EndRadius" : [" End Radius", " 5.5", 6],
             "nTurns" : [" Number of Turns", " 2", 7],
             "degRotation" : [" Initial Rotation", " 0.0", 8],
             "nPts" : [" N SpiralPts", " 360", 9]}
datPolygon = {"xCenterPolygon" : [" X PolygonCenter", " 6.0", 3],
              "yCenterPolygon" : [" Y PolygonCenter", " -10.0", 4],
              "Radius" : [" Radius", " 5.0", 5],
              "nVert" : [" Number of Vertices", " 3", 6],
              "degRotation" : [" Initial Rotation", " 0.0", 7],
              "nPts" : [" N PolygonPts", " 300", 8]}
datSpiraGon = {"xCenterSpiraGon" : [" X SpiraGonCenter", " 6.0", 3],
               "yCenterSpiraGon" : [" Y SpiraGonCenter", " -10.0", 4],
               "StartRadius" : [" Start Radius", " 4.0", 5],
               "EndRadius" : [" End Radius", " 5.0", 6],
               "nVert" : [" Number of Vertices", " 3", 7],
               "nTurns" : [" Number of Turns", " 2", 8],
               "degRotation" : [" Initial Rotation", " 0.0", 9],
               "nPts" : [" N SpiraGonPts", " 300", 10]}
datStar = {"xCenterStar" : [" X StarCenter", " 3.0", 3],
           "yCenterStar" : [" Y StarCenter", " -1.0", 4],
           "StartRadius" : [" Start Radius", " 4.0", 5],
           "EndRadius" : [" End Radius", " 5.0", 6],
           "nSpikes" : [" Number of Spikes", " 9", 7],
           "degRotation" : [" Initial Rotation", " 0.0", 8],
           "nPts" : [" N StarPts", " 400", 9]}
datIterative = {"Separator" : ["", "", 11],
                "SectionLbl" : ["", "", 12],
                "aParam" : [" A Parameter", " -0.995", 13],
                "bParam" : [" B Parameter", " 0.999", 14],
                "orbitPts" : [" Orbit Points", " 250", 15]}
datScreen = {"Separator" : ["", "", 16],
             "SectionLbl" : ["", "", 17],
             "pltStart" : [" Plot Start", " 0", 18],
             "pltStop" : [" Plot Stop", " 500000", 19],
             "pltStride" : [" Plot Stride", " 1", 20],
             "Marker" : [" Marker", " o", 21],
             "MarkerSize" : [" Marker Size", " 49.0", 22],
             "MarkerTaper" : [" Marker Taper", "None", 23],
             "TaperBtn" : ["", "", 24]}

cmapSection = {"Separator" : ["", "", 25],
               "radioSTD" : ["", "", 26],
               "radioCUST" : ["", "", 27]}

axesSection = {"Separator" : ["", "", 28],
               "axesEQUAL" : ["", "", 29],
               "axesMINMAX" : ["", "", 30]}

saveRawSection = {"Separator" : ["", "", 31],
                  "rawBtn" : ["", "", 32]}

# CMAP Radio Button Data
lblStdCMAP = "Standard CMAP"
defaultStdCMAP = " cool"
lblCustomCMAP = "Custom CMAP: List"
defaultCustomCMAP = " aqPu.cmL"
cMapTypeInit = "STD"

# Plot Axes Button Data
lblEqualAxes = " Equal XY Axes"
lblAxesMinMax = " Axes: XY MinMax"
AxesInit = "MINMAX"

# Save Raw Data Check Button
saveRawInit = 'FALSE'

helpImgRefs = \
    {'Page 1'  :  './SupportFiles/HelpImages/Pg1_1.png',
     'Page 2'  :  './SupportFiles/HelpImages/Pg2_1.png',
     'Page 3'  :  './SupportFiles/HelpImages/Pg3_1.png',
     'Page 4'  :  './SupportFiles/HelpImages/Pg4_1.png',
     'Page 5'  :  './SupportFiles/HelpImages/Pg5_1.png',
     'Page 6'  :  './SupportFiles/HelpImages/Pg6_1.png',
     'Page 7'  :  './SupportFiles/HelpImages/Pg7_1.png',
     'Page 8'  :  './SupportFiles/HelpImages/Pg8_1.png',
     'Page 9'  :  './SupportFiles/HelpImages/Pg9_1.png',
     'Page 10' :  './SupportFiles/HelpImages/Pg10_1.png'}
     
tabAlgsRefs = \
    {'Version1'    : './SupportFiles/Algorithm1.png',
     'Version2'    : './SupportFiles/Algorithm2.png',
     'Version3'    : './SupportFiles/Algorithm3.png'}

# Gallery References
galleryRefs = \
    {'Circle1'     : './Gallery/Circle/Gallery_Circle1.png',
     'Circle2'     : './Gallery/Circle/Gallery_Circle2.png',
     'CircleSine1' : './Gallery/CircleSine/Gallery_CircleSine1.png',
     'CircleSine2' : './Gallery/CircleSine/Gallery_CircleSine2.png',
     'Star1'       : './Gallery/Star/Gallery_Star1.png',
     'Star2'       : './Gallery/Star/Gallery_Star2.png',
     'Spiral1'     : './Gallery/Spiral/Gallery_Spiral1.png',
     'Spiral2'     : './Gallery/Spiral/Gallery_Spiral2.png',
     'Polygon1'    : './Gallery/Polygon/Gallery_Polygon1.png',
     'Polygon2'    : './Gallery/Polygon/Gallery_Polygon2.png',
     'SpiraGon1'   : './Gallery/SpiraGon/Gallery_SpiraGon1.png',
     'SpiraGon2'   : './Gallery/SpiraGon/Gallery_SpiraGon2.png'}

miscRefs = \
    {'Overview'    : './SupportFiles/Overview.png',
     'About'       : './SupportFiles/aboutImg.png'}


