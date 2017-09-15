#!/usr/bin/python

import os, sys, math

####################################################
# Modify coordinates
x = [0.0, 10.0]
z = [0.0, 0.50]
yMax = 1.0
yMin = -1.0

####################################################
# Check
if len(x) != len(z):
    print "Error in x and y.\n"
    sys.exit(0)

# Previous calculations
nPoi = len(x) # Points
nSec = nPoi - 1 # Sections
nFac = 2 * nSec # Facets

# Write the file
fileW = open('ramp.stl', 'w')
fileW.write('solid ramp\n')

for i in range(nSec):
    dz = (z[i+1]-z[i])
    dx = (x[i+1]-x[i])
    ang = math.atan2(-dz,dx)
    
    j = 1
    fileW.write('facet normal %f 0.0 %f\n' % (math.sin(ang), math.cos(ang)) )
    fileW.write('   outer loop\n')

    fileW.write('      vertex %f %f %f \n' % (x[i], yMax, z[i]) )
    fileW.write('      vertex %f %f %f \n' % (x[i], yMin, z[i]) )
    fileW.write('      vertex %f %f %f \n' % (x[i+1], yMin, z[i+1]) )
            
    fileW.write('   endloop\n')  
    fileW.write('endfacet\n')
    
    j = 2
    fileW.write('facet normal %f 0.0 %f\n' % (math.sin(ang), math.cos(ang)) )
    fileW.write('   outer loop\n')

    fileW.write('      vertex %f %f %f \n' % (x[i], yMax, z[i]) )
    fileW.write('      vertex %f %f %f \n' % (x[i+1], yMax, z[i+1]) )
    fileW.write('      vertex %f %f %f \n' % (x[i+1], yMin, z[i+1]) )

    fileW.write('   endloop\n') 
    fileW.write('endfacet\n')
    
fileW.write('endsolid ramp\n')
fileW.close()
