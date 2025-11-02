# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                                                 #
# UnitCellDistortionMapping.py                                                                    #
# Script for mapping trigonometric distortions to unit cell                                       #
# Units: Kilograms-Meters-Seconds                                                                 #
# Initiated by Dennis Hagnas on Mon Dec 30 2024                                                   #
#                                                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import csv
import math
import re
import numpy as np
#--------------------------------------------------------------------------------------------------
# Import parametric variables

with open('parameters.txt', 'r') as f:
    exec(f.read())
#--------------------------------------------------------------------------------------------------
# Import mesh data

with open('IdealMesh.inp', 'r') as f:
    idealLines = f.readlines()   
idealNodes = []
readNodes = False
for line in idealLines:
    if line.startswith('*Node'):
        readNodes = True
    elif line.startswith('*') and readNodes:
        break
    elif readNodes:
        parts = re.split(r',\s*', line.strip())
        label = int(parts[0])
        x = float(parts[1])
        y = float(parts[2])
        z = float(parts[3])
        idealNodes.append((label, x, y, z))
#--------------------------------------------------------------------------------------------------
# Trigonometric series defining the distortions

H = [
    [0.006, 0.000, 0.000],
    [0.000, 0.000, 0.000],
    [0.000, 0.000, 0.000]
]
def formula(x, z, H, M, N, l, s):
    v = 0.0
    for m in range(1, M+1):
        for n in range(1, N+1):
            H_mn = H[m-1][n-1]
            v += H_mn * math.sin((m * math.pi * x) / l) * math.sin((n * math.pi * z) / s)
    return v
#--------------------------------------------------------------------------------------------------
# Apply distortion

ey = hw + 0.5 * (tp + tf)
distortedNodes = []
for node in idealNodes:
    label, x, y, z = node
    newY = y + ((ey + y)/ey)*formula(x, z, H, M, N, l, s)
    distortedNodes.append((label, x, newY, z))
#--------------------------------------------------------------------------------------------------
# Write to input file
    
with open('DistortedMesh.inp', 'w') as f:
    readNodes = False
    for line in idealLines:
        if line.lstrip().startswith('*Preprint'):
            continue
        if line.startswith('*Node'):
            readNodes = True
            f.write(line)
            for node in distortedNodes:
                f.write('{0}, {1}, {2}, {3}\n'.format(node[0], node[1], node[2], node[3]))
        elif line.startswith('*') and readNodes:
            readNodes = False
            f.write(line)
        elif not readNodes:
            f.write(line)
# ---------------------------------------- END OF SCRIPT ---------------------------------------- #