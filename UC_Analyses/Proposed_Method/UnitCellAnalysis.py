# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                                                 #
# UnitCellAnalysis.py                                                                             #
# Parametric unit cell model                                                                      #
# Units: Kilograms-Meters-Seconds                                                                 #
# Initiated by Dennis Hagnas on Mon Dec 16 2024                                                   #
#                                                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from abaqusConstants import *
from odbAccess import *
import os
import csv
import numpy as np
import visualization
import math
import time
import subprocess

startTime = time.time()
# -------------------------------------------------------------------------------------------------
# Create new model database

mdb.Model(name = 'Model-0')
# -------------------------------------------------------------------------------------------------
# Set path to working directory

workDir = os.getcwd()
os.chdir(workDir)
# -------------------------------------------------------------------------------------------------
# Import parametric variables

with open('parameters.txt', 'r') as f:
    exec(f.read())
# -------------------------------------------------------------------------------------------------
# Create the sketch to be extruded

ez = hw + 0.5 * (tp + tf)
s1 = mdb.models['Model-0'].ConstrainedSketch(name = '__profile__', sheetSize = 1.0)
s1.Line(point1 = (-s/2, 0.0), point2 = (-s, 0.0))
s1.Line(point1 = (0.0,  0.0), point2 = (-s/2,0.0))
s1.Line(point1 = (-s/2, 0.0), point2 = (-s/2, -ez))
if bf > 0:
    s1.Line(point1 = (-((s/2)-(bf/2+ef)), -ez), point2 = (-s/2, -ez))
    if bf/2 > ef:
        s1.Line(point1 = (-((s/2)+(bf/2-ef)), -ez), point2 = (-s/2, -ez))
# -------------------------------------------------------------------------------------------------
# Extrude sketch

p1 = mdb.models['Model-0'].Part(name = 'Part-UnitCell', dimensionality = THREE_D, 
    type = DEFORMABLE_BODY)
p1.BaseShellExtrude(sketch = s1, depth = l)
# -------------------------------------------------------------------------------------------------
# Create linear elastic material

mdb.models['Model-0'].Material(name = 'Material-1')
mdb.models['Model-0'].materials['Material-1'].Elastic(table = ((E, nu), ))
# -------------------------------------------------------------------------------------------------
# Create sections

mdb.models['Model-0'].HomogeneousShellSection(name = 'Section-Plate', material = 'Material-1', 
    thickness = tp, poissonDefinition = DEFAULT)	
mdb.models['Model-0'].HomogeneousShellSection(name = 'Section-Web', material = 'Material-1', 
    thickness = tw, poissonDefinition = DEFAULT)
if tf > 0:	
    mdb.models['Model-0'].HomogeneousShellSection(name = 'Section-Flange', material = 'Material-1', 
        thickness = tf, poissonDefinition = DEFAULT)
# -------------------------------------------------------------------------------------------------
# Assign sections
	
faces = p1.faces.findAt(((0, 0, 0),), ((-s, 0, 0),) )
region = p1.Set(faces = faces, name = 'Set-Plate')
p1.SectionAssignment(region = region, sectionName = 'Section-Plate', offset = 0.0, 
    offsetType = MIDDLE_SURFACE, offsetField = '', thicknessAssignment = FROM_SECTION)
faces = p1.faces.findAt(((-s/2, -ez/2, 0),))
region = p1.Set(faces = faces, name = 'Set-Web')
p1.SectionAssignment(region = region, sectionName = 'Section-Web', offset = 0.0, 
    offsetType = MIDDLE_SURFACE, offsetField = '', thicknessAssignment = FROM_SECTION)
if bf > 0:
    if bf/2 > ef:
        faces = p1.faces.findAt(((-((s/2)+(bf/2-ef)), -ez, 0),), ((-s/2+(bf/2-ef), -ez, 0),) ) 
        region = p1.Set(faces = faces, name = 'Set-Flange')
        p1.SectionAssignment(region = region, sectionName = 'Section-Flange', offset = 0.0, 
            offsetType = MIDDLE_SURFACE, offsetField = '', thicknessAssignment = FROM_SECTION)
    else:
        faces = p1.faces.findAt(((-s/2+(bf/2), -ez, 0),) ) 
        region = p1.Set(faces = faces, name = 'Set-Flange')
        p1.SectionAssignment(region = region, sectionName = 'Section-Flange', offset = 0.0, 
            offsetType = MIDDLE_SURFACE, offsetField = '', thicknessAssignment = FROM_SECTION)       
# -------------------------------------------------------------------------------------------------
# Create assembly and orient part in global coordinate system

a1 = mdb.models['Model-0'].rootAssembly
a1.DatumCsysByDefault(CARTESIAN)
a1.Instance(name = 'Part-UnitCell-1', part = p1, dependent = OFF)
a1.rotate(instanceList = ('Part-UnitCell-1', ), axisPoint = (0.0, 0.0, 0.0), 
    axisDirection = (0.0, 1.0, 0.0), angle = 90.0)
a1.rotate(instanceList = ('Part-UnitCell-1', ), axisPoint = (0.0, 0.0, 0.0), 
    axisDirection = (1.0, 0.0, 0.0), angle = -90.0)
# -------------------------------------------------------------------------------------------------
# Generate ideal mesh

e1 = a1.instances['Part-UnitCell-1'].edges
xEdges = e1.findAt(((l/2, 0, 0),), ((l/2, s, 0),), ((l/2, s/2+(bf/2-ef),ez),))
yEdges = e1.findAt(((0, s/3, 0),), ((0, s/1.5, 0),), ((l, s/3, 0),), ((l, s/1.5, 0),))  
a1.seedEdgeByNumber(edges = xEdges, number = nx, constraint = FIXED)
a1.seedEdgeByNumber(edges = yEdges, number = ny/2, constraint = FIXED)
partInstances = (a1.instances['Part-UnitCell-1'], )
a1.generateMesh(regions = partInstances)
# -------------------------------------------------------------------------------------------------
# Create set for center node for u3 = 0 BC

n1 = a1.instances['Part-UnitCell-1'].nodes
W0 = n1.getClosest((l/2, s/2, 0), )
a1.SetFromNodeLabels(name = 'Set-W0', nodeLabels = (('Part-UnitCell-1', 
    (W0.label, )), ), unsorted = False)
# -------------------------------------------------------------------------------------------------
# Create set for reference node in the stiffener   
    
refNode = n1.getClosest((l, s/2, ez), )
a1.SetFromNodeLabels(name = 'Set-ReferenceNode', nodeLabels = (('Part-UnitCell-1', 
    (refNode.label, )), ), unsorted = False)
# -------------------------------------------------------------------------------------------------
# Define major edge node sets (corner nodes included)

tol = 1e-4

a1.Set(nodes = n1.getByBoundingBox(xMin = -tol, xMax = l+tol, yMin = s-tol, yMax = s+tol, 
    zMin = -tol, zMax = ez+tol), name = 'Set-LeftMax')
a1.Set(nodes = n1.getByBoundingBox(xMin = -tol, xMax = l+tol, yMin = -tol, yMax = tol, 
    zMin = -tol, zMax = ez+tol), name = 'Set-RightMax')
a1.Set(nodes = n1.getByBoundingBox(xMin = l-tol, xMax = l+tol, yMin = -tol, yMax = s+tol, 
    zMin = -tol, zMax = ez+tol), name = 'Set-FrontMax')
a1.Set(nodes = n1.getByBoundingBox(xMin = -tol, xMax = tol, yMin = -tol, yMax = s+tol, 
    zMin = -tol, zMax = ez+tol), name = 'Set-BackMax')    
a1.Set(nodes = n1.getByBoundingBox(xMin = -tol, xMax = tol, yMin = -tol, yMax = s+tol, 
    zMin = -tol, zMax = tol), name = 'Set-BackPlate')
a1.Set(nodes = n1.getByBoundingBox(xMin = l-tol, xMax = l+tol, yMin = -tol, yMax = s+tol, 
    zMin = -tol, zMax = tol), name = 'Set-FrontPlate')
# -------------------------------------------------------------------------------------------------
# Define minor edge node sets (corner nodes excluded)

a1.Set(nodes = n1.getByBoundingBox(xMin = tol, xMax = l-tol, yMin = s-tol, yMax = s+tol, 
    zMin = -tol, zMax = ez+tol), name = 'Set-LeftMin')
a1.Set(nodes = n1.getByBoundingBox(xMin = tol, xMax = l-tol, yMin = -tol, yMax = tol, 
    zMin = -tol, zMax = ez+tol), name = 'Set-RightMin')
a1.Set(nodes = n1.getByBoundingBox(xMin = l-tol, xMax = l+tol, yMin = tol, yMax = s-tol, 
    zMin = -tol, zMax = ez+tol), name = 'Set-FrontMin')
a1.Set(nodes = n1.getByBoundingBox(xMin = -tol, xMax = tol, yMin = tol, yMax = s-tol, 
    zMin = -tol, zMax = ez+tol), name = 'Set-BackMin')
# -------------------------------------------------------------------------------------------------
# Create a set for every edge node

for curNode in a1.sets['Set-LeftMax'].nodes:
    a1.SetFromNodeLabels(name = 'Set-LeftNode-'+str(curNode.label), 
        nodeLabels = (('Part-UnitCell-1', (curNode.label, )), ), unsorted = False)
for curNode in a1.sets['Set-RightMax'].nodes:
    a1.SetFromNodeLabels(name = 'Set-RightNode-'+str(curNode.label), 
        nodeLabels = (('Part-UnitCell-1', (curNode.label, )), ), unsorted = False)
for curNode in a1.sets['Set-FrontMax'].nodes:
    a1.SetFromNodeLabels(name = 'Set-FrontNode-'+str(curNode.label), 
        nodeLabels = (('Part-UnitCell-1', (curNode.label, )), ), unsorted = False)
for curNode in a1.sets['Set-BackMax'].nodes:
    a1.SetFromNodeLabels(name = 'Set-BackNode-'+str(curNode.label), 
        nodeLabels = (('Part-UnitCell-1', (curNode.label, )), ), unsorted = False)
# -------------------------------------------------------------------------------------------------
# Define major and minor lists containing Left-Right and Front-Back node pairs

def lrNodePairs(leftNodes, rightNodes, tol):
    matches = []
    for lNode in leftNodes:
        lx = lNode.coordinates[0]
        lz = lNode.coordinates[2]
        for rNode in rightNodes:
            rx = rNode.coordinates[0]
            rz = rNode.coordinates[2]
            if np.isclose(lx, rx, atol = tol) and np.isclose(lz, rz, atol = tol):
                matches.append((lNode.label, rNode.label))
    return matches  
def fbNodePairs(frontNodes, backNodes, tol):
    matches = []
    for fNode in frontNodes:
        fy = fNode.coordinates[1]
        fz = fNode.coordinates[2]
        for bNode in backNodes:
            by = bNode.coordinates[1]
            bz = bNode.coordinates[2]
            if np.isclose(fy, by, atol = tol) and np.isclose(fz, bz, atol = tol):
                matches.append((fNode.label, bNode.label))
    return matches
LR = lrNodePairs(a1.sets['Set-LeftMax'].nodes, a1.sets['Set-RightMax'].nodes, tol)
lr = lrNodePairs(a1.sets['Set-LeftMin'].nodes, a1.sets['Set-RightMin'].nodes, tol)
FB = fbNodePairs(a1.sets['Set-FrontMax'].nodes, a1.sets['Set-BackMax'].nodes, tol)
fb = fbNodePairs(a1.sets['Set-FrontMin'].nodes, a1.sets['Set-BackMin'].nodes, tol)
# -------------------------------------------------------------------------------------------------
# Write input file for distortion

mdb.Job(name = 'IdealMesh', model = 'Model-0', description = '', type = ANALYSIS, 
    nodalOutputPrecision = SINGLE)
mdb.jobs['IdealMesh'].writeInput(consistencyChecking = OFF)
mdb.jobs['IdealMesh'].waitForCompletion()
# -------------------------------------------------------------------------------------------------
# Execute subprocess for mapping distortions to mesh

subprocess.call(['python', 'UnitCellDistortionMapping.py'])
# -------------------------------------------------------------------------------------------------
# Import distorted mesh

mdb.ModelFromInputFile(name = 'Model-1-Neg', inputFileName = 'DistortedMesh.inp')
a2 = mdb.models['Model-1-Neg'].rootAssembly
p2 = mdb.models['Model-1-Neg'].parts.values()[0]
session.viewports['Viewport: 1'].setValues(displayedObject = p2)
session.viewports['Viewport: 1'].view.rotate(xAngle = 90, yAngle = 0, zAngle = 0, mode = MODEL)
session.viewports['Viewport: 1'].view.fitView()
# -------------------------------------------------------------------------------------------------
# Create step

mdb.models['Model-1-Neg'].StaticRiksStep(name = 'Step-1', previous = 'Initial', maxNumInc = 10000,
    initialArcInc = 0.0001, minArcInc = 1e-10, maxArcInc = 0.01, nlgeom = ON, maxLPF = 1)
# -------------------------------------------------------------------------------------------------
# Create u3 = 0 BC

mdb.models['Model-1-Neg'].DisplacementBC(name = 'BC-W0', createStepName = 'Initial', 
    region = a2.sets['SET-W0'], u3 = SET)   
# -------------------------------------------------------------------------------------------------
# Copy model

mdb.Model(name = 'Model-2-Neg', objectToCopy = mdb.models['Model-1-Neg'])
mdb.Model(name = 'Model-3-Neg', objectToCopy = mdb.models['Model-1-Neg'])
mdb.Model(name = 'Model-4-Neg', objectToCopy = mdb.models['Model-1-Neg'])
mdb.Model(name = 'Model-5-Neg', objectToCopy = mdb.models['Model-1-Neg'])
mdb.Model(name = 'Model-6-Neg', objectToCopy = mdb.models['Model-1-Neg'])
# ----------------------------------------------------------------------------------------------- #
# # # # # # # # # # # # # # # LOAD CASE 1 - EXTENSIONAL STRAIN IN X # # # # # # # # # # # # # # # #
# ----------------------------------------------------------------------------------------------- #
# Write boundary conditions for negative strain job

mdb.models['Model-1-Neg'].DisplacementBC(name = 'BC-Left', createStepName = 'Initial', 
    region = a2.sets['SET-LEFTMAX'], u2 = SET, ur1 = SET)
mdb.models['Model-1-Neg'].DisplacementBC(name = 'BC-Right', createStepName = 'Initial', 
    region = a2.sets['SET-RIGHTMAX'], u2 = SET, ur1 = SET)
mdb.models['Model-1-Neg'].DisplacementBC(name = 'BC-Front', createStepName = 'Step-1', 
    region = a2.sets['SET-FRONTMAX'], u1 = -0.005*l, ur2 = SET)
mdb.models['Model-1-Neg'].DisplacementBC(name = 'BC-Back', createStepName = 'Step-1', 
    region = a2.sets['SET-BACKMAX'], u1 = 0.005*l, ur2 = SET)    
# -------------------------------------------------------------------------------------------------
# Define equation constraints

for idx, (nodeA, nodeB) in enumerate(lr):
    mdb.models['Model-1-Neg'].Equation(name='Constraint-U1-' + str(nodeA) + '-' + str(nodeB), 
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 1), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 1)))   
    mdb.models['Model-1-Neg'].Equation(name='Constraint-U3-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 3), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 3)))    
    mdb.models['Model-1-Neg'].Equation(name='Constraint-UR2-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 5), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 5)))
for idy, (nodeA, nodeB) in enumerate(FB):   
    mdb.models['Model-1-Neg'].Equation(name='Constraint-U3-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 3), (-1.0, 'SET-BACKNODE-' + str(nodeB), 3)))   
for idy, (nodeA, nodeB) in enumerate(fb):
    mdb.models['Model-1-Neg'].Equation(name='Constraint-UR1-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 4), (-1.0, 'SET-BACKNODE-' + str(nodeB), 4)))
    mdb.models['Model-1-Neg'].Equation(name='Constraint-U2-' + str(nodeA) + '-' + str(nodeB), 
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 2), (-1.0, 'SET-BACKNODE-' + str(nodeB), 2)))    
# -------------------------------------------------------------------------------------------------
# Copy Model

mdb.Model(name = 'Model-1-Pos', objectToCopy = mdb.models['Model-1-Neg'])
# -------------------------------------------------------------------------------------------------
# Redefine boundary conditions for positive strain job

mdb.models['Model-1-Pos'].boundaryConditions['BC-Front'].setValues(u1 = 0.005*l)
mdb.models['Model-1-Pos'].boundaryConditions['BC-Back'].setValues(u1 = -0.005*l)
# -------------------------------------------------------------------------------------------------
# Create and submit Load Case 1 jobs

def submitJob(modelNo):
    job1 = 'Job-'+str(modelNo) +'-Neg'
    model1 = 'Model-'+str(modelNo) +'-Neg'
    job2 = 'Job-'+str(modelNo) +'-Pos'
    model2 = 'Model-'+str(modelNo) +'-Pos'
    mdb.Job(name = job1, model = model1, description = '', type = ANALYSIS,  
        explicitPrecision = SINGLE, nodalOutputPrecision = SINGLE, echoPrint = OFF, 
        modelPrint = OFF, contactPrint = OFF, historyPrint = OFF, numCpus = 1, numGPUs = 0)
    mdb.Job(name = job2, model = model2, description = '', type = ANALYSIS, 
        explicitPrecision = SINGLE, nodalOutputPrecision = SINGLE, echoPrint = OFF, 
        modelPrint = OFF, contactPrint = OFF, historyPrint = OFF, numCpus = 1, numGPUs = 0)
    mdb.jobs[job1].submit(consistencyChecking = OFF)
    mdb.jobs[job2].submit(consistencyChecking = OFF)
    mdb.jobs[job1].waitForCompletion()
    mdb.jobs[job2].waitForCompletion()
mdb.models['Model-1-Pos'].steps['Step-1'].setValues(maxArcInc = 0.002)
submitJob(1)
# ----------------------------------------------------------------------------------------------- #
# # # # # # # # # # # # # # # LOAD CASE 2 - EXTENSIONAL STRAIN IN Y # # # # # # # # # # # # # # # # 
# ----------------------------------------------------------------------------------------------- #
# Write boundary conditions for negative strain job

mdb.models['Model-2-Neg'].DisplacementBC(name = 'BC-Left', createStepName = 'Step-1', 
    region = a2.sets['SET-LEFTMAX'], u2 = -0.005*s, ur1 = SET)
mdb.models['Model-2-Neg'].DisplacementBC(name = 'BC-Right', createStepName = 'Step-1', 
    region = a2.sets['SET-RIGHTMAX'], u2 = 0.005*s, ur1 = SET)
mdb.models['Model-2-Neg'].DisplacementBC(name = 'BC-Front', createStepName = 'Initial', 
    region = a2.sets['SET-FRONTMAX'], u1 = SET, ur2 = SET)
mdb.models['Model-2-Neg'].DisplacementBC(name = 'BC-Back', createStepName = 'Initial', 
    region = a2.sets['SET-BACKMAX'], u1 = SET, ur2 = SET)
mdb.models['Model-2-Neg'].DisplacementBC(name = 'BC-BackPlate', createStepName = 'Initial', 
    region = a2.sets['SET-BACKPLATE'], u3 = SET)
mdb.models['Model-2-Neg'].DisplacementBC(name = 'BC-FrontPlate', createStepName = 'Initial', 
    region = a2.sets['SET-FRONTPLATE'], u3 = SET)
mdb.models['Model-2-Neg'].boundaryConditions['BC-W0'].setValues(u3 = UNSET)
# -------------------------------------------------------------------------------------------------
# Define equation constraints

for idx, (nodeA, nodeB) in enumerate(lr):
    mdb.models['Model-2-Neg'].Equation(name='Constraint-U1-' + str(nodeA) + '-' + str(nodeB), 
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 1), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 1)))
    mdb.models['Model-2-Neg'].Equation(name='Constraint-U3-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 3), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 3)))
    mdb.models['Model-2-Neg'].Equation(name='Constraint-UR2-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 5), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 5)))
for idx, (nodeA, nodeB) in enumerate(FB):
    mdb.models['Model-2-Neg'].Equation(name='Constraint-U3-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 3), (-1.0, 'SET-BACKNODE-' + str(nodeB), 3)))
for idy, (nodeA, nodeB) in enumerate(fb):
    mdb.models['Model-2-Neg'].Equation(name='Constraint-U2-' + str(nodeA) + '-' + str(nodeB), 
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 2), (-1.0, 'SET-BACKNODE-' + str(nodeB), 2)))    
    mdb.models['Model-2-Neg'].Equation(name='Constraint-UR1-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 4), (-1.0, 'SET-BACKNODE-' + str(nodeB), 4)))  
# -------------------------------------------------------------------------------------------------
# Copy Model

mdb.Model(name = 'Model-2-Pos', objectToCopy = mdb.models['Model-2-Neg'])
# -------------------------------------------------------------------------------------------------
# Redefine boundary conditions for positive strain job

mdb.models['Model-2-Pos'].boundaryConditions['BC-Left'].setValues(u2 = 0.005*s)
mdb.models['Model-2-Pos'].boundaryConditions['BC-Right'].setValues(u2 = -0.005*s)
# -------------------------------------------------------------------------------------------------
# Create and submit jobs

mdb.models['Model-2-Pos'].steps['Step-1'].setValues(maxArcInc = 0.002)
submitJob(2)
# ----------------------------------------------------------------------------------------------- #
# # # # # # # # # # # # # # # # # # LOAD CASE 3 - SHEAR STRAIN  # # # # # # # # # # # # # # # # # # 
# ----------------------------------------------------------------------------------------------- #
# Write boundary conditions for negative strain job

mdb.models['Model-3-Neg'].DisplacementBC(name = 'BC-Left', createStepName = 'Initial', 
    region = a2.sets['SET-LEFTMAX'], u1 = SET, ur1 = SET)
mdb.models['Model-3-Neg'].DisplacementBC(name = 'BC-Right', createStepName = 'Initial', 
    region = a2.sets['SET-RIGHTMAX'], u1 = SET, ur1 = SET)
mdb.models['Model-3-Neg'].DisplacementBC(name = 'BC-Front', createStepName = 'Step-1', 
    region = a2.sets['SET-FRONTMAX'], u1 = SET, u2 = -0.005*s, ur2 = SET)
mdb.models['Model-3-Neg'].DisplacementBC(name = 'BC-Back', createStepName = 'Step-1', 
    region = a2.sets['SET-BACKMAX'], u1 = SET, u2 = 0.005*s, ur2 = SET)
# -------------------------------------------------------------------------------------------------
# Define equation constraints

for idx, (nodeA, nodeB) in enumerate(lr):
    mdb.models['Model-3-Neg'].Equation(name='Constraint-U3-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 3), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 3)))  
    mdb.models['Model-3-Neg'].Equation(name='Constraint-UR2-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 5), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 5)))
for idy, (nodeA, nodeB) in enumerate(FB):
    mdb.models['Model-3-Neg'].Equation(name='Constraint-U3-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 3), (-1.0, 'SET-BACKNODE-' + str(nodeB), 3)))
for idy, (nodeA, nodeB) in enumerate(fb):        
    mdb.models['Model-3-Neg'].Equation(name='Constraint-UR1-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 4), (-1.0, 'SET-BACKNODE-' + str(nodeB), 4)))
u2 = a2.sets['SET-BACKMAX'].nodes[0] 
for curNode in a2.sets['SET-LEFTMIN'].nodes:
    x1 = curNode.coordinates[0]        
    mdb.models['Model-3-Neg'].Equation(name = 'Constraint-U2-' + str(curNode.label), terms = (
        (-1.0, 'SET-LEFTNODE-' + str (curNode.label), 2), 
        (1.0, 'SET-BACKNODE-' + str (u2.label), 2), 
        ((-2 * x1) / l, 'SET-BACKNODE-' + str (u2.label), 2)))
for curNode in a2.sets['SET-RIGHTMIN'].nodes:
    x1 = curNode.coordinates[0]        
    mdb.models['Model-3-Neg'].Equation(name = 'Constraint-U2-' + str(curNode.label), terms = (
        (-1.0, 'SET-RIGHTNODE-' + str (curNode.label), 2), 
        (1.0, 'SET-BACKNODE-' + str (u2.label), 2), 
        ((-2 * x1) / l, 'SET-BACKNODE-' + str (u2.label), 2)))
# -------------------------------------------------------------------------------------------------
# Copy Model

mdb.Model(name = 'Model-3-Pos', objectToCopy = mdb.models['Model-3-Neg'])
# -------------------------------------------------------------------------------------------------
# Redefine boundary conditions for positive strain job

mdb.models['Model-3-Pos'].boundaryConditions['BC-Front'].setValues(u2 = 0.005*s)
mdb.models['Model-3-Pos'].boundaryConditions['BC-Back'].setValues(u2 = -0.005*s)
#--------------------------------------------------------------------------------------------------
# Create and submit jobs

submitJob(3)
# ----------------------------------------------------------------------------------------------- #
# # # # # # # # # # # # # LOAD CASE 4 - UNIDIRECTIONAL BENDING STRAIN IN X  # # # # # # # # # # # #  
# ----------------------------------------------------------------------------------------------- #
# Write boundary conditions for negative strain job

thetax = 0.25 * l
mdb.models['Model-4-Neg'].DisplacementBC(name = 'BC-Right', createStepName = 'Initial', 
    region = a2.sets['SET-RIGHTMAX'], u2 = SET, ur1 = SET)
mdb.models['Model-4-Neg'].DisplacementBC(name = 'BC-Left', createStepName = 'Initial', 
    region = a2.sets['SET-LEFTMAX'], u2 = SET, ur1 = SET)
for curNode in a2.sets['SET-FRONTMAX'].nodes:
    mdb.models['Model-4-Neg'].DisplacementBC(name = 'BC-FrontNode-' + str(curNode.label), 
            createStepName = 'Step-1', region = a2.sets['SET-FRONTNODE-' + str(curNode.label)], 
            u1 = curNode.coordinates[2] * math.tan(-thetax), ur2 = -thetax)
for curNode in a2.sets['SET-BACKMAX'].nodes: 
    mdb.models['Model-4-Neg'].DisplacementBC(name = 'BC-BackNode-' + str(curNode.label), 
            createStepName = 'Step-1', region = a2.sets['SET-BACKNODE-' + str(curNode.label)], 
            u1 = curNode.coordinates[2] * math.tan(thetax), ur2 = thetax)
# -------------------------------------------------------------------------------------------------
# Define equation constraints

for idx, (nodeA, nodeB) in enumerate(lr):
    mdb.models['Model-4-Neg'].Equation(name='Constraint-U1-' + str(nodeA) + '-' + str(nodeB), 
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 1), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 1)))
    mdb.models['Model-4-Neg'].Equation(name='Constraint-U3-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 3), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 3)))    
    mdb.models['Model-4-Neg'].Equation(name='Constraint-UR2-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 5), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 5)))
for idy, (nodeA, nodeB) in enumerate(FB):
    mdb.models['Model-4-Neg'].Equation(name='Constraint-U3-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 3), (-1.0, 'SET-BACKNODE-' + str(nodeB), 3)))   
for idy, (nodeA, nodeB) in enumerate(fb):        
    mdb.models['Model-4-Neg'].Equation(name='Constraint-U2-' + str(nodeA) + '-' + str(nodeB), 
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 2), (-1.0, 'SET-BACKNODE-' + str(nodeB), 2)))
    mdb.models['Model-4-Neg'].Equation(name='Constraint-UR1-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 4), (-1.0, 'SET-BACKNODE-' + str(nodeB), 4)))
# -------------------------------------------------------------------------------------------------
# Copy Model

mdb.Model(name = 'Model-4-Pos', objectToCopy = mdb.models['Model-4-Neg'])
# -------------------------------------------------------------------------------------------------
# Redefine boundary conditions for positive strain job

for curNode in a2.sets['SET-FRONTMAX'].nodes:
    mdb.models['Model-4-Pos'].boundaryConditions['BC-FrontNode-' + str(curNode.label)].setValues(
        u1 = curNode.coordinates[2] * math.tan(thetax), ur2 = thetax) 
for curNode in a2.sets['SET-BACKMAX'].nodes:
    mdb.models['Model-4-Pos'].boundaryConditions['BC-BackNode-' + str(curNode.label)].setValues(
        u1 = curNode.coordinates[2] * math.tan(-thetax), ur2 = -thetax)
# -------------------------------------------------------------------------------------------------
# Create and submit jobs

submitJob(4)
# ----------------------------------------------------------------------------------------------- #
# # # # # # # # # # # # # LOAD CASE 5 - UNIDIRECTIONAL BENDING STRAIN IN Y  # # # # # # # # # # # #
# ----------------------------------------------------------------------------------------------- #
# Write boundary conditions for negative strain job

thetay = 0.25 * s
for curNode in a2.sets['SET-LEFTMAX'].nodes:
    mdb.models['Model-5-Neg'].DisplacementBC(name = 'BC-LeftNode-' + str(curNode.label), 
        createStepName = 'Step-1', region = a2.sets['SET-LEFTNODE-' + str(curNode.label)], 
        u2 = curNode.coordinates[2] * math.tan(thetay), ur1 = thetay) 
for curNode in a2.sets['SET-RIGHTMAX'].nodes:
    mdb.models['Model-5-Neg'].DisplacementBC(name = 'BC-RightNode-' + str(curNode.label), 
        createStepName = 'Step-1', region = a2.sets['SET-RIGHTNODE-' + str(curNode.label)], 
         u2 = curNode.coordinates[2] * math.tan(-thetay), ur1 = -thetay)
mdb.models['Model-5-Neg'].DisplacementBC(name = 'BC-Front', createStepName = 'Initial', 
    region = a2.sets['SET-FRONTMAX'], u1 = SET, ur2 = SET)
mdb.models['Model-5-Neg'].DisplacementBC(name = 'BC-Back', createStepName = 'Initial', 
    region = a2.sets['SET-BACKMAX'], u1 = SET, ur2 = SET)
# -------------------------------------------------------------------------------------------------
# Define equation constraints

for idx, (nodeA, nodeB) in enumerate(lr):
    mdb.models['Model-5-Neg'].Equation(name='Constraint-U1-' + str(nodeA) + '-' + str(nodeB), 
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 1), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 1)))
    mdb.models['Model-5-Neg'].Equation(name='Constraint-U3-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 3), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 3)))
    mdb.models['Model-5-Neg'].Equation(name='Constraint-UR2-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 5), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 5)))
for idy, (nodeA, nodeB) in enumerate(FB):
    mdb.models['Model-5-Neg'].Equation(name='Constraint-U3-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 3), (-1.0, 'SET-BACKNODE-' + str(nodeB), 3)))
for idy, (nodeA, nodeB) in enumerate(fb):
    mdb.models['Model-5-Neg'].Equation(name='Constraint-UR1-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 4), (-1.0, 'SET-BACKNODE-' + str(nodeB), 4)))    
    mdb.models['Model-5-Neg'].Equation(name='Constraint-U2-' + str(nodeA) + '-' + str(nodeB), 
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 2), (-1.0, 'SET-BACKNODE-' + str(nodeB), 2)))
# -------------------------------------------------------------------------------------------------
# Copy Model

mdb.Model(name = 'Model-5-Pos', objectToCopy = mdb.models['Model-5-Neg'])
# -------------------------------------------------------------------------------------------------
# Redefine boundary conditions for positive strain job

for curNode in a2.sets['SET-LEFTMAX'].nodes:
    mdb.models['Model-5-Pos'].boundaryConditions['BC-LeftNode-' + str(curNode.label)].setValues(
        u2 = curNode.coordinates[2] * math.tan(-thetay), ur1 = -thetay)  
for curNode in a2.sets['SET-RIGHTMAX'].nodes:
    mdb.models['Model-5-Pos'].boundaryConditions['BC-RightNode-' + str(curNode.label)].setValues(
        u2 = curNode.coordinates[2]* math.tan(thetay), ur1 = thetay)
#--------------------------------------------------------------------------------------------------
# Create and submit jobs

submitJob(5)
# ----------------------------------------------------------------------------------------------- #
# # # # # # # # # # # # # # # # # # # # LOAD CASE 6 - TWIST # # # # # # # # # # # # # # # # # # # # 
# ----------------------------------------------------------------------------------------------- #
# Write boundary conditions for negative strain job

phi = 0.25 * s
mdb.models['Model-6-Neg'].DisplacementBC(name = 'BC-Left', createStepName = 'Initial', 
    region = a2.sets['SET-LEFTMAX'], u1 = SET, ur1 = SET)
mdb.models['Model-6-Neg'].DisplacementBC(name = 'BC-Right', createStepName = 'Initial', 
    region = a2.sets['SET-RIGHTMAX'], u1 = SET, ur1 = SET)
mdb.models['Model-6-Neg'].DisplacementBC(name = 'BC-FrontPlate', createStepName = 'Initial', 
    region = a2.sets['SET-FRONTPLATE'], u1 = SET)  
for curNode in a2.sets['SET-FRONTMAX'].nodes:
    mdb.models['Model-6-Neg'].DisplacementBC(name = 'BC-FrontNode-' + str(curNode.label), 
        createStepName = 'Step-1', region = a2.sets['SET-FRONTNODE-' + str(curNode.label)], 
        u2 = curNode.coordinates[2] * math.tan(-phi), ur2 = SET) 
mdb.models['Model-6-Neg'].DisplacementBC(name = 'BC-BackPlate', createStepName = 'Initial', 
    region = a2.sets['SET-BACKPLATE'], u1 = SET)  
for curNode in a2.sets['SET-BACKMAX'].nodes:
    mdb.models['Model-6-Neg'].DisplacementBC(name = 'BC-BackNode-' + str(curNode.label), 
        createStepName = 'Step-1', region = a2.sets['SET-BACKNODE-' + str(curNode.label)], 
        u2 = curNode.coordinates[2] * math.tan(phi), ur2 = SET)        
# -------------------------------------------------------------------------------------------------
# Define equation constraints

for idx, (nodeA, nodeB) in enumerate(lr):   
    mdb.models['Model-6-Neg'].Equation(name='Constraint-U3-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 3), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 3)))
    mdb.models['Model-6-Neg'].Equation(name='Constraint-UR2-' + str(nodeA) + '-' + str(nodeB), 
        terms = ((1.0, 'SET-LEFTNODE-' + str(nodeA), 5), (-1.0, 'SET-RIGHTNODE-' + str(nodeB), 5)))
for idy, (nodeA, nodeB) in enumerate(FB):
    mdb.models['Model-6-Neg'].Equation(name='Constraint-U3-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 3), (-1.0, 'SET-BACKNODE-' + str(nodeB), 3)))
    mdb.models['Model-6-Neg'].Equation(name='Constraint-U1-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 1), (-1.0, 'SET-BACKNODE-' + str(nodeB), 1)))
        
for idy, (nodeA, nodeB) in enumerate(fb):
    mdb.models['Model-6-Neg'].Equation(name='Constraint-UR1-' + str(nodeA) + '-' + str(nodeB),
        terms = ((1.0, 'SET-FRONTNODE-' + str(nodeA), 4), (-1.0, 'SET-BACKNODE-' + str(nodeB), 4)))
# -------------------------------------------------------------------------------------------------
# Copy Model

mdb.Model(name = 'Model-6-Pos', objectToCopy = mdb.models['Model-6-Neg'])
# -------------------------------------------------------------------------------------------------
# Redefine boundary conditions for positive strain job
        
for curNode in a2.sets['SET-FRONTMAX'].nodes:
    mdb.models['Model-6-Pos'].boundaryConditions['BC-FrontNode-' + str(curNode.label)].setValues(
        u2 = curNode.coordinates[2] * math.tan(phi))
for curNode in a2.sets['SET-BACKMAX'].nodes:
    mdb.models['Model-6-Pos'].boundaryConditions['BC-BackNode-' + str(curNode.label)].setValues(
        u2 = curNode.coordinates[2] * math.tan(-phi))
# -------------------------------------------------------------------------------------------------
# Create and submit jobs

submitJob(6)
# -------------------------------------------------------------------------------------------------
# Extract generalised forces and strains

def Extraction(jobNo):
    mdb.jobs['Job-'+str(jobNo)+'-Neg'].waitForCompletion()
    odb = visualization.openOdb(path = 'Job-'+str(jobNo)+'-Neg.odb')
    with open('loadCase-'+str(jobNo)+'.csv', mode = 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        for i, frame in enumerate(reversed(odb.steps['Step-1'].frames)):
            rf1Tot, rf2Tot, rm1Tot, rm2Tot, rf12Tot, rm12Tot = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            rfField = frame.fieldOutputs['RF']
            rmField = frame.fieldOutputs['RM']
            uField  = frame.fieldOutputs['U']
            urField = frame.fieldOutputs['UR']
            for curNode in odb.rootAssembly.nodeSets['SET-FRONTMAX'].nodes[0]:
                rf1 = rfField.getSubset(region = curNode).values[0].data[0]
                rm2 = rmField.getSubset(region = curNode).values[0].data[1]
                rf12 = rfField.getSubset(region = curNode).values[0].data[1]
                u3 = uField.getSubset(region = curNode).values[0].data[2]
                x3 = curNode.coordinates[2]
                rf1Tot  += rf1
                rf12Tot += rf12
                rm2Tot  += rm2 + rf1 * x3
                rm12Tot += rf12 * x3
            for curNode in odb.rootAssembly.nodeSets['SET-LEFTMAX'].nodes[0]:
                rf2 = rfField.getSubset(region = curNode).values[0].data[1]
                rm1 = rmField.getSubset(region = curNode).values[0].data[0]
                u3 = uField.getSubset(region = curNode).values[0].data[2]
                x3 = curNode.coordinates[2]
                rf2Tot += rf2
                rm1Tot += -rm1 + rf2 * x3 
            u1  = 2 * uField.getSubset(
                region = odb.rootAssembly.nodeSets['SET-REFERENCENODE']).values[0].data[0]
            u12 = 2 * uField.getSubset(
                region = odb.rootAssembly.nodeSets['SET-REFERENCENODE']).values[0].data[1]
            u2  = 2 * uField.getSubset(
                region = odb.rootAssembly.nodeSets['SET-LEFTMAX']).values[0].data[1]
            ur1 = 2 * urField.getSubset(
                region = odb.rootAssembly.nodeSets['SET-LEFTMAX']).values[0].data[0]
            refx3 = refNode.coordinates[2]
            writer.writerow([i, rf1Tot/s, rf2Tot/l, rf12Tot/s, 
                            rm2Tot/s, rm1Tot/l, rm12Tot/s, 
                            u1/l, u2/s, u12/s, 
                            u1/(l*refx3), -ur1/s, u12/(s*refx3)])
    nFrames = len(odb.steps['Step-1'].frames)
    odb.close()
    mdb.jobs['Job-'+str(jobNo)+'-Pos'].waitForCompletion()    
    odb = visualization.openOdb(path = 'Job-'+str(jobNo)+'-Pos.odb')
    with open('loadCase-'+str(jobNo)+'.csv', mode = 'a') as file:
        writer = csv.writer(file, lineterminator = '\n')
        for i, frame in enumerate(odb.steps['Step-1'].frames):
            if i == 0:
                continue
            rf1Tot, rf2Tot, rm1Tot, rm2Tot, rf12Tot, rm12Tot = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            rfField = frame.fieldOutputs['RF']
            rmField = frame.fieldOutputs['RM']
            uField  = frame.fieldOutputs['U']
            urField = frame.fieldOutputs['UR']
            for curNode in odb.rootAssembly.nodeSets['SET-FRONTMAX'].nodes[0]:
                rf1 = rfField.getSubset(region = curNode).values[0].data[0]
                rm2 = rmField.getSubset(region = curNode).values[0].data[1]
                rf12 = rfField.getSubset(region = curNode).values[0].data[1]
                u3 = uField.getSubset(region = curNode).values[0].data[2]
                x3 = curNode.coordinates[2]
                rf1Tot += rf1
                rf12Tot += rf12
                rm2Tot += rm2 + rf1 * x3
                rm12Tot += rf12 * x3         
            for curNode in odb.rootAssembly.nodeSets['SET-LEFTMAX'].nodes[0]:
                rf2 = rfField.getSubset(region = curNode).values[0].data[1]
                rm1 = rmField.getSubset(region = curNode).values[0].data[0]
                u3 = uField.getSubset(region = curNode).values[0].data[2]
                x3 = curNode.coordinates[2]
                rf2Tot += rf2
                rm1Tot += -rm1 + rf2 * x3    
            u1  = 2 * uField.getSubset(
                region = odb.rootAssembly.nodeSets['SET-REFERENCENODE']).values[0].data[0]
            u12 = 2 * uField.getSubset(
                region = odb.rootAssembly.nodeSets['SET-REFERENCENODE']).values[0].data[1]
            u2  = 2 * uField.getSubset(
                region = odb.rootAssembly.nodeSets['SET-LEFTMAX']).values[0].data[1]
            ur1 = 2 * urField.getSubset(
                region = odb.rootAssembly.nodeSets['SET-LEFTMAX']).values[0].data[0]
            refx3 = refNode.coordinates[2]
            writer.writerow([i+nFrames-1, rf1Tot/s, rf2Tot/l, rf12Tot/s, 
                            rm2Tot/s, rm1Tot/l, rm12Tot/s, 
                            u1/l, u2/s, u12/s, 
                            u1/(l*refx3), -ur1/s, u12/(s*refx3)])
    odb.close()
for i in range(1, 7):
    Extraction(i)
# -------------------------------------------------------------------------------------------------
# Remove excess files and determine analysis time

jobs = ['Job-1-Neg', 'Job-1-Pos', 'Job-2-Neg', 'Job-2-Pos', 'Job-3-Neg', 'Job-3-Pos', 'Job-4-Neg',
        'Job-4-Pos', 'Job-5-Neg', 'Job-5-Pos', 'Job-6-Neg', 'Job-6-Pos']
extensions = ['.prt', '.com', '.msg', '.ipm', '.log']
for job in jobs:
    for extension in extensions:
        file = job + extension
        if os.path.isfile(file):
            os.remove(file) 
os.remove('IdealMesh.inp')
os.remove('DistortedMesh.inp')
endTime = time.time() 
duration = endTime - startTime
minutes = int(duration // 60)
seconds = int(duration % 60)
print("Unit cell analysis completed in {} minutes and {} seconds.".format(minutes, seconds))
# ---------------------------------------- END OF SCRIPT ---------------------------------------- #