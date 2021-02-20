# Case Study Structure 4 cols, 4 beams, real dimensions
# SI units 
# Kaatsiz - FEB20


# import OpenSeesPy
from openseespy.opensees import *

import numpy as np
import math

#import openseespy.postprocessing.ops_vis as opsv
#import matplotlib.pyplot as plt

# wipe model
wipe()

# create model
model('basic', '-ndm', 3, '-ndf', 6)

storyHeight = float(3)

# MaterialProps
fc=35 #in MPa
Ec=5000*(math.sqrt(fc))*1000 #in KPa

# print(Ec)

# adasds

# nodes

node(11,0,0,0)
node(12,8,0,0)
node(21,0,8,0)
node(22,8,8,0)
node(111,0,0,storyHeight)
node(112,8,0,storyHeight)
node(121,0,8,storyHeight)
node(122,8,8,storyHeight)
node(100,4,4,storyHeight)


# support conditions,

fixValues=[1,1,1,1,1,1] #all dofs fixed

fix(11,*fixValues)
fix(12,*fixValues)
fix(21,*fixValues)
fix(22,*fixValues)

fixCMNodes=[0,0,1,1,1,0]  #fix cm nodes where constrained

fix(100,*fixCMNodes)

# rigid diaphragm constraint for the slab
perpDirn=3 #vertical is dof 3
rNodeTag=100 #firstStoryCM
cNodes=[111,112,121,122] #slaveNodes


rigidDiaphragm(perpDirn, rNodeTag, *cNodes)



#Define floor masses
massX=39.057
massY=massX
massTheta=416.61

cmMass=[massX,massY,0,0,0,massTheta]
#assign mass to CM
mass(100,*cmMass)


### DEFINE MEMBERS

# Define geometric transformations. ---- 
# Beams: primary bending in local Z dir.
# Columns: local z is in Global X dir.

colVector=[1,0,0]
xBeamVector=[0,-1,0]
yBeamVector=[1,0,0]




colTransf=1
xBeamTransf=2
yBeamTransf=3

geomTransf('Linear', colTransf, *colVector)
geomTransf('Linear', xBeamTransf, *xBeamVector)
geomTransf('Linear', yBeamTransf, *yBeamVector)


#start with connectivity
ele1010Nodes=[111,112]
ele1020Nodes=[121,122]
ele1100Nodes=[111,121]
ele1200Nodes=[112,122]

ele1111Nodes=[11,111]
ele1121Nodes=[12,112]
ele1211Nodes=[21,121]
ele1221Nodes=[22,122]

# BEAMS
hBeam=0.60
bBeam=0.35
A_beam=hBeam*bBeam
crackConstantBeam=0.35
E_beam=Ec*crackConstantBeam
G_beam=E_beam/(2*(1+0.3)) #poisson is 0.3
J_beam=(0.229*hBeam*math.pow(bBeam,3)) #for a rectangular x section beta is approx taken as 0.229
Iy_beam=1/12*hBeam*math.pow(bBeam,3)
Iz_beam=1/12*bBeam*math.pow(hBeam,3)
# print(A_beam)


# print(E_beam)
# print(A_beam)
# print(G_beam)
# print(J_beam)
# print(Iy_beam)
# print(Iz_beam)

# asds

# XBEAMS
element('elasticBeamColumn', 1010, *ele1010Nodes, A_beam, E_beam, G_beam, J_beam, Iy_beam, Iz_beam, xBeamTransf)
element('elasticBeamColumn', 1020, *ele1020Nodes, A_beam, E_beam, G_beam, J_beam, Iy_beam, Iz_beam, xBeamTransf)
# YBEAMS
element('elasticBeamColumn', 1100, *ele1100Nodes, A_beam, E_beam, G_beam, J_beam, Iy_beam, Iz_beam, yBeamTransf)
element('elasticBeamColumn', 1200, *ele1200Nodes, A_beam, E_beam, G_beam, J_beam, Iy_beam, Iz_beam, yBeamTransf)


# COLUMNS
hColumn=0.55
bColumn=0.30
A_column=hColumn*bColumn
crackConstantColumn=0.50
E_column=Ec*crackConstantColumn
G_column=E_column/(2*(1+0.3)) #poisson is 0.3
J_column=(0.229*hColumn*math.pow(bColumn,3)) #for a rectangular x section beta is approx taken as 0.229
Iy_column=1/12*hColumn*math.pow(bColumn,3)
Iz_column=1/12*bColumn*math.pow(hColumn,3)


# Columns
element('elasticBeamColumn', 1111, *ele1111Nodes, A_column, E_column, G_column, J_column, Iy_column, Iz_column, colTransf)
element('elasticBeamColumn', 1121, *ele1121Nodes, A_column, E_column, G_column, J_column, Iy_column, Iz_column, colTransf)
element('elasticBeamColumn', 1211, *ele1211Nodes, A_column, E_column, G_column, J_column, Iy_column, Iz_column, colTransf)
element('elasticBeamColumn', 1221, *ele1221Nodes, A_column, E_column, G_column, J_column, Iy_column, Iz_column, colTransf)


#  LOADING
# Create patterin
timeSeries('Linear', 1)
pattern('Plain', 1, 1)

# ASSIGN POINT LOADS
nodalLoad=[0,0,-12.375,0,0,0] #column weights per column
load(111,*nodalLoad)
load(121,*nodalLoad)
load(112,*nodalLoad)
load(122,*nodalLoad)

# Assign beam uniform loads
Wy=-19.25 #computed total weight/4 distributed per beam

eleLoad('-ele', 1010,'-type','-beamUniform', Wy, 0.0, 0.0)
eleLoad('-ele', 1020,'-type','-beamUniform', Wy, 0.0, 0.0)
eleLoad('-ele', 1100,'-type','-beamUniform', Wy, 0.0, 0.0)
eleLoad('-ele', 1200,'-type','-beamUniform', Wy, 0.0, 0.0)



# opsv.plot_model(node_labels=1, element_labels=1, offset_nd_label=False, axis_off=0, az_el=(-60.0, 30.0), fig_wi_he=(16.0, 10.0), fig_lbrt=(0.04, 0.04, 0.96, 0.96))



# Delete the old analysis and all it's component objects
wipeAnalysis()

#ANALYSIS PARAMETERS

# create constraint handler
constraints("Transformation")

# create DOF number
numberer("RCM")

# create SOE
system("BandGeneral")

# create test object
tolerance=1*math.pow(10,-8)
test('EnergyIncr',tolerance,10)


# create integrator
integrator("LoadControl", 0.1)

# create algorithm
algorithm("Linear")

# create analysis object
analysis("Static")
# adsads
analyze(10)


#  Support Reactions
reactions()
support11=nodeReaction(11)
support12=nodeReaction(12)
support21=nodeReaction(21)
support22=nodeReaction(22)

support11=np.array(support11)
support12=np.array(support12)
support21=np.array(support21)
support22=np.array(support22)

totalSupportReaction=support11+support12+support21+support22
totalSupportReaction= np.ndarray.tolist(totalSupportReaction)

# vele=(totalSupportReaction[2])
# vele=str(vele)

print ('Total Weight is '+format(totalSupportReaction[2],"3.2f")+' kN')


#  Node Displacements - TEST OUTPUTS
node111=nodeDisp(111)
node112=nodeDisp(112)
node121=nodeDisp(121)
node122=nodeDisp(122)

cmNodeDisps=nodeDisp(100)

#beam forces - TEST OUTPUTS
beam1111Forces=eleForce(1010) #THIS GIVES GLOBAL FORCE - 12 DOFS

beam1111Data= eleResponse(1010,'force') #THE STRING here invokes setresponse() method depending on element. For an elastic beam col element only valid query is force.
# force beam column elements could give valid chordDeformation and plasticDeformation outputs using this command.



# EIGEN VALUE ANALYSIS
LambdaValues=eigen('-fullGenLapack',3)


Omega = np.sqrt(LambdaValues)
Periods=2*math.pi/Omega

Omega=np.ndarray.tolist(Omega)
Periods=np.ndarray.tolist(Periods)

# print(Omega)
# print(Periods)


# opsv.plot_model()

print ('Frequencies: '+", ".join(format(x, "2.2f") for x in Omega)+' rad/s')
print ('Periods: '+", ".join(format(x, "2.2f") for x in Periods)+' seconds')


print('done')







