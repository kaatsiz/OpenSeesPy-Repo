# Single column analysis 

# This wil evolve into a method that takes the axial load and 6 coordinates for two nodes. 
# single column spanning from x1 y1 z1 to x2 y2 z2 with an axial load at x2 y2 z2 having a value of P

# KAATSIZ MAR21


#function ARGUMENTs
nodalCoordinates = [0,0,0,0,0,3]
axialLoad=[0,0,-100,0,0,0] #apply -100kN in the direction of gravity

# PREREQUISITES

# import OpenSeesPy
from openseespy.opensees import *

import numpy as np
import math



# wipe model == this may be a required argument to be passed to the method, if we need a wipe before calculation
wipe()

# create model = these will be parametric, it could be a method as well in the future.
model('basic', '-ndm', 3, '-ndf', 6)



# MaterialProps == They will come from the material databases, will be passed to analysis
fc=35 #in MPa
Ec=5000*(math.sqrt(fc))*1000 #in KPa



# NODE DEFINITIONS 
# First node: x1 y1 z1, second node x2 y2 z2

# Nodal and other element tags are not parametic here. they can be sent from the database into these methods. Right now they are random integers.

node(11,nodalCoordinates[0],nodalCoordinates[1],nodalCoordinates[2])
node(12,nodalCoordinates[3],nodalCoordinates[4],nodalCoordinates[5])


# support conditions, This may also be a method. We will decide which nodes will be fixed and send them to related methods accordingly.

fixValues=[1,1,1,1,1,1] #all dofs fixed

fix(11,*fixValues)


## THESE ARE 3D Model Related Stuff. They will be parametretized and methodized later.
# fixCMNodes=[0,0,1,1,1,0]  #fix cm nodes where constrained

# fix(100,*fixCMNodes)

# # rigid diaphragm constraint for the slab
# perpDirn=3 #vertical is dof 3
# rNodeTag=100 #firstStoryCM
# cNodes=[111,112,121,122] #slaveNodes


# rigidDiaphragm(perpDirn, rNodeTag, *cNodes)


## THESE ARE MASS Related Stuff for 3D structures and dynamic analysis. They will be parametretized and methodized later. We will compute mass and pass the values to the methods
# #Define floor masses
# massX=39.057
# massY=massX
# massTheta=416.61

# cmMass=[massX,massY,0,0,0,massTheta]
# #assign mass to CM
# mass(100,*cmMass)


# DEFINE MEMBERS

# Define geometric transformations. ---- 
# Beams: primary bending in local Z dir.
# Columns: local z is in Global X dir.

### TRANSFORMATION STUFF = Basically it shows the orientation of the element section. We need to compute this according to geometry defined in the interface,
# these will be separate methods and basically send 3 element vectors showing the orientation. Right now it is hand calculated. colVector is the one here.

colVector=[1,0,0]
# xBeamVector=[0,-1,0]
# yBeamVector=[1,0,0]



# column transformation tag. We need to parametrize this as well.
colTransf=1
# xBeamTransf=2
# yBeamTransf=3

# and here is the transformation command. OpenSees understands the transformation as per the syntax below.
geomTransf('Linear', colTransf, *colVector)
# geomTransf('Linear', xBeamTransf, *xBeamVector)
# geomTransf('Linear', yBeamTransf, *yBeamVector)




## ASSIGNING THE COLUMN MEMBER FOR THIS METHOD.


#start with connectivity = that shows the A and B nodes of the element that we are defining. Syntax for OpenSees Coommand. While defining every member, we may need this sort of a variable list.
ele1010Nodes=[11,12]

# Right Now we are using linear elastic elements. The section properties below are defined for this member. Each member wil lhave its own methods, submethods etc.

# Some of these section properties will be computed from the section and material interfaces and will be passed to the related methods.
hBeam=0.60 #beam depth
bBeam=0.35 #beam width
A_beam=hBeam*bBeam #beam section area
crackConstantBeam=0.35 #beam section property
E_beam=Ec*crackConstantBeam #material property
G_beam=E_beam/(2*(1+0.3)) #poisson is 0.3 maeterial property
J_beam=(0.229*hBeam*math.pow(bBeam,3)) #for a rectangular x section beta is approx taken as 0.229 #section property
Iy_beam=1/12*hBeam*math.pow(bBeam,3) #section property
Iz_beam=1/12*bBeam*math.pow(hBeam,3) #section property
# print(A_beam)


# print(E_beam)
# print(A_beam)
# print(G_beam)
# print(J_beam)
# print(Iy_beam)
# print(Iz_beam)

# asds

# XBEAMS
element('elasticBeamColumn', 1010, *ele1010Nodes, A_beam, E_beam, G_beam, J_beam, Iy_beam, Iz_beam, colTransf)
# element('elasticBeamColumn', 1020, *ele1020Nodes, A_beam, E_beam, G_beam, J_beam, Iy_beam, Iz_beam, xBeamTransf)
# # YBEAMS
# element('elasticBeamColumn', 1100, *ele1100Nodes, A_beam, E_beam, G_beam, J_beam, Iy_beam, Iz_beam, yBeamTransf)
# element('elasticBeamColumn', 1200, *ele1200Nodes, A_beam, E_beam, G_beam, J_beam, Iy_beam, Iz_beam, yBeamTransf)



#  LOADING == Now this is a very separate method. We need to define a new timeSeries object and pattern for `nearly` every kind of loading. 
# Here for the nodal load we are defininf a linear timeSeries (to apply the load in a linear way) and a plain pattern (basic load acting mechanism)
# Create patterin
timeSeries('Linear', 1)
pattern('Plain', 1, 1)

# ASSIGN POINT LOADS
load(12,*axialLoad) #for this method

### Assign beam uniform loads This will be a completely different method. We will compute the distributed loads on beams and apply them accordingly
# Wy=-19.25 #computed total weight/4 distributed per beam

# eleLoad('-ele', 1010,'-type','-beamUniform', Wy, 0.0, 0.0)
# eleLoad('-ele', 1020,'-type','-beamUniform', Wy, 0.0, 0.0)
# eleLoad('-ele', 1100,'-type','-beamUniform', Wy, 0.0, 0.0)
# eleLoad('-ele', 1200,'-type','-beamUniform', Wy, 0.0, 0.0)



# Delete the old analysis and all it's component objects # may be required to be passed as an argument.
wipeAnalysis()

#ANALYSIS PARAMETERS -- now this part depends on entirely the type of analysis the user wants to do. Of course we will keep some of the options to ourselves.
# This part is definitely going to be a method such as analysisParameters etc/

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
integrator("LoadControl", 1)

# create algorithm
algorithm("Linear")

# create analysis object
analysis("Static")
# 

# This is the actual analysis command that does everyting. Here it runs once, applying the full load and writes the results to the environment. Since this is linear elastic, it is OK.
analyze(1)

# We will definitely create methods for gettin different responses from the analysis environment. Below we are getting some forces and displacements only in a non parametric manner.
#  REST is post process. Here we are getting the support reactions for the fixed node below. == if something is fixed to not to move, than forces occur.

#  Support Reactions
reactions()
support11=nodeReaction(11)


support11=np.array(support11)

# vele=(totalSupportReaction[2])
# vele=str(vele)

print ('Total Weight is '+format(support11[2],"3.2f")+' kN') # this is just a self check. We are seeing 100 kN as the vertical reaction which is equal to the applied load.

#The output for supportreatcion is beliw?
supportReaction=support11 #6 values for 3 forces in x y z directions and 3 moments.


#  Node Displacements - Lets get free node displacements under loading. This will be a diferent method as well.
topNodeDisplacement=nodeDisp(12) #3 displacements and 3 rotations for x y z directions. we will see only top node negative displacement.




#columnForces - Lets get element forces at the ends of the column. There are many ways to get internal forces for an element. We will implement most of them. This is the simplest one for visualization.
columnForce=eleForce(1010) #THIS GIVES GLOBAL FORCE - 12 DOFS = 6 for bottom end 6 for top end. only Nodal force 100 is applied, therefore 100 is shown in vertical element forces at both ends.

# beam1111Data= eleResponse(1010,'force') #THE STRING here invokes setresponse() method depending on element. For an elastic beam col element only valid query is force.
# force beam column elements could give valid chordDeformation and plasticDeformation outputs using this command.



# EIGEN VALUE ANALYSIS = This is an independent method. If mass is also defined, we can determine dynamic properties of a structure by doing the things below
# LambdaValues=eigen('-fullGenLapack',3)


# Omega = np.sqrt(LambdaValues)
# Periods=2*math.pi/Omega

# Omega=np.ndarray.tolist(Omega)
# Periods=np.ndarray.tolist(Periods)

# # print(Omega)
# # print(Periods)


# # opsv.plot_model()

# print ('Frequencies: '+", ".join(format(x, "2.2f") for x in Omega)+' rad/s')
# print ('Periods: '+", ".join(format(x, "2.2f") for x in Periods)+' seconds')


# print('done')







